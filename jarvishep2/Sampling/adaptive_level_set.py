#!/usr/bin/env python3
"""Adaptive level-set sampler (D10): gen-0 → graph crossing → barrier refine.

One algorithm, two neighbor-graph builders (Delaunay d≤3 / kNN d=4–5).
Generation-0: Bridson for d≤4, Sobol QMC for d=5. Design:
``Jarvis-Books/Jarvis-HEP V2/components/adaptive_voronoi_contour.md``.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass, field
from typing import Any, Protocol

import numpy as np

from jarvishep2.Sampling.bridson import Bridson_sampling, hypersphere_surface_sample
from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection, map_u_to_physical
from jarvishep2.Sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.Sampling.variables import load_variables
from jarvishep2.expression import CompiledExpression
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample

_EPS = float(np.finfo(np.float64).eps)
_CLIP_LO = _EPS
_CLIP_HI = 1.0 - _EPS


@dataclass
class LevelSetPoint:
    u: list[float]
    x: dict[str, float]
    f: float | None
    uuid: str = ""
    generation: int = 0


class NeighborGraph(Protocol):
    def build(self, points_u: np.ndarray) -> np.ndarray:
        """Return undirected edges as shape (m, 2) int array."""
        ...


class DelaunayGraph:
    """Exact site adjacency via Voronoi ridge_points (d ≤ 3)."""

    def build(self, points_u: np.ndarray) -> np.ndarray:
        from scipy.spatial import Voronoi

        if points_u.shape[0] < 2:
            return np.zeros((0, 2), dtype=np.int64)
        if points_u.shape[0] < points_u.shape[1] + 1:
            # Degenerate: fall back to complete graph of available pairs
            return np.array(list(itertools.combinations(range(points_u.shape[0]), 2)), dtype=np.int64)
        vor = Voronoi(points_u, qhull_options="QJ")
        ridges = np.asarray(vor.ridge_points, dtype=np.int64)
        if ridges.size == 0:
            return np.zeros((0, 2), dtype=np.int64)
        return ridges


class KNNGraph:
    """Symmetric kNN graph via cKDTree (d = 4–5 proximity mode)."""

    def __init__(self, knn_k: int) -> None:
        self.knn_k = max(1, int(knn_k))

    def build(self, points_u: np.ndarray) -> np.ndarray:
        from scipy.spatial import cKDTree

        n = int(points_u.shape[0])
        if n < 2:
            return np.zeros((0, 2), dtype=np.int64)
        k = min(self.knn_k + 1, n)
        tree = cKDTree(points_u)
        _dist, idx = tree.query(points_u, k=k)
        edges: set[tuple[int, int]] = set()
        for i in range(n):
            neighbors = np.atleast_1d(idx[i])
            for j in neighbors:
                j = int(j)
                if j == i:
                    continue
                a, b = (i, j) if i < j else (j, i)
                edges.add((a, b))
        if not edges:
            return np.zeros((0, 2), dtype=np.int64)
        return np.array(sorted(edges), dtype=np.int64)


def sample_ball(center: np.ndarray, radius: float, rng: np.random.Generator, d: int) -> np.ndarray:
    """Uniform sample in the d-ball of given radius around center."""
    if radius <= 0:
        return np.asarray(center, dtype=np.float64).copy()
    direction = rng.normal(size=d)
    norm = float(np.linalg.norm(direction))
    if norm < 1e-15:
        direction = np.ones(d)
        norm = float(np.linalg.norm(direction))
    direction /= norm
    # radius scale: r * U^{1/d}
    scale = radius * (float(rng.random()) ** (1.0 / max(d, 1)))
    return np.asarray(center, dtype=np.float64) + direction * scale


def clip_unit_cube(u: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(u, dtype=np.float64), _CLIP_LO, _CLIP_HI)


def interpolate_crossing(
    u_i: np.ndarray,
    u_j: np.ndarray,
    f_i: float,
    f_j: float,
    target: float,
) -> np.ndarray:
    denom = f_j - f_i
    if abs(denom) < 1e-15:
        return 0.5 * (u_i + u_j)
    t = (target - f_i) / denom
    t = float(np.clip(t, 0.0, 1.0))
    return u_i + t * (u_j - u_i)


class AdaptiveLevelSetSampler(CheckpointedSampler):
    """Feedback-driven level-set tracer (Method: AdaptiveLevelSet)."""

    method = "AdaptiveLevelSet"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.adaptive_level_set")
        self.vars: list = []
        self._dim = 0
        self._target_expression = ""
        self._target_value = 0.0
        self._contour_precision = 0.01
        self._function_tolerance = 0.05
        self._initial_radius = 0.08
        self._refinement_factor = 0.5
        self._max_generations = 25
        self._max_points = 5000
        self._knn_k = 8
        self._neighbor_graph_mode = "auto"
        self._max_new_per_generation = 500
        self._k_ref = 4
        self._slice_pairs: list[tuple[str, str]] | None = None
        self._simplify_tolerance: float | None = None
        self._selectionexp: str | None = None
        self._seed = 0
        self._batch_size = 16
        self._graph: NeighborGraph | None = None
        self._points: list[LevelSetPoint] = []
        self._generation = 0
        self._pending_uuids: set[str] = set()
        self._accepted_index = 0
        self._seed_seq: np.random.SeedSequence | None = None
        self._converged = False
        self._levelset: dict[str, Any] | None = None
        self._target_fn: CompiledExpression | None = None
        self._failed_regions: list[dict[str, Any]] = []
        self._uuid_to_index: dict[str, int] = {}
        self._stop_reason = ""

    # ------------------------------------------------------------------ config
    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dim = len(self.vars)
        if self._dim < 2 or self._dim > 5:
            raise ValueError(
                f"AdaptiveLevelSet requires 2 ≤ dim ≤ 5, got {self._dim}"
            )
        block = sampling.get("AdaptiveLevelSet") or sampling.get("adaptive_level_set") or {}
        if not isinstance(block, Mapping):
            raise ValueError("Sampling.AdaptiveLevelSet must be a mapping")
        self._target_expression = str(block.get("target_expression", "")).strip()
        if not self._target_expression:
            raise ValueError("AdaptiveLevelSet.target_expression is required")
        if "target_value" not in block:
            raise ValueError("AdaptiveLevelSet.target_value is required")
        self._target_value = float(block["target_value"])
        self._contour_precision = float(block.get("contour_precision", 0.01) or 0.01)
        self._function_tolerance = float(block.get("function_tolerance", 0.05) or 0.05)
        self._initial_radius = float(block.get("initial_radius", 0.08) or 0.08)
        default_factor = 0.5 if self._dim <= 3 else 0.65
        self._refinement_factor = float(
            block.get("refinement_factor", default_factor) or default_factor
        )
        self._max_generations = max(1, int(block.get("max_generations", 25) or 25))
        default_max_points = 5000 if self._dim <= 3 else 20000
        self._max_points = max(10, int(block.get("max_points", default_max_points) or default_max_points))
        self._knn_k = max(1, int(block.get("knn_k", 4 * self._dim) or (4 * self._dim)))
        self._neighbor_graph_mode = str(block.get("neighbor_graph", "auto")).strip().lower()
        self._max_new_per_generation = max(
            1,
            int(block.get("max_new_per_generation", self._max_points // 10) or max(1, self._max_points // 10)),
        )
        self._k_ref = max(1, int(block.get("k_ref", 4) or 4))
        self._selectionexp = sampling.get("selection")
        self._seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._seed_seq = np.random.SeedSequence(self._seed)
        self._compile_target()
        self._graph = self._make_graph()
        # Optional slice pairs for d≥4
        raw_pairs = block.get("slice_pairs")
        if isinstance(raw_pairs, list) and raw_pairs:
            pairs: list[tuple[str, str]] = []
            for item in raw_pairs:
                if isinstance(item, (list, tuple)) and len(item) == 2:
                    pairs.append((str(item[0]), str(item[1])))
            self._slice_pairs = pairs or None
        else:
            self._slice_pairs = None
        simp = block.get("simplify_tolerance")
        self._simplify_tolerance = float(simp) if simp is not None else None
        if self._neighbor_graph_mode == "delaunay" and self._dim >= 4:
            self._logger.warning(
                "neighbor_graph=delaunay with d=%d may be expensive", self._dim
            )

    def _compile_target(self) -> None:
        self._target_fn = self._expression_context.compile(
            self._target_expression,
            symbols=(
                "x",
                "y",
                "z",
                "shift",
                "calc_z",
                "LogL",
                "LogL_Z",
                *[str(v.name) for v in self.vars],
            ),
        )

    def _make_graph(self) -> NeighborGraph:
        mode = self._neighbor_graph_mode
        if mode == "auto":
            mode = "delaunay" if self._dim <= 3 else "knn"
        if mode == "delaunay":
            return DelaunayGraph()
        if mode == "knn":
            return KNNGraph(self._knn_k)
        raise ValueError(f"unknown neighbor_graph mode '{self._neighbor_graph_mode}'")

    def _eval_target(self, observables: Mapping[str, Any]) -> float | None:
        if self._target_fn is None:
            self._compile_target()
        assert self._target_fn is not None
        try:
            value = self._target_fn.evaluate(observables)
            if isinstance(value, np.generic):
                value = value.item()
            return float(value)
        except Exception:
            return None

    # --------------------------------------------------------------- generation-0
    def _generation_rng(self, generation: int) -> np.random.Generator:
        assert self._seed_seq is not None
        child = self._seed_seq.spawn(generation + 1)[generation]
        return np.random.default_rng(child)

    def _generation_0_points(self) -> list[np.ndarray]:
        d = self._dim
        r = self._initial_radius
        if d <= 4:
            # Bridson on unit cube
            rng = self._generation_rng(0)
            # Bridson_sampling uses global np.random; seed from child
            seed_int = int(rng.integers(0, 2**31 - 1))
            np.random.seed(seed_int)
            dims = np.ones(d, dtype=np.float64)
            points = Bridson_sampling(
                dims=dims,
                radius=r,
                k=30,
                hypersphere_sample=hypersphere_surface_sample,
            )
            # Bridson returns points in [0, dims); map to unit cube u
            us = [clip_unit_cube(np.asarray(row, dtype=np.float64) / dims) for row in points]
        else:
            # d = 5: Sobol (C10)
            from scipy.stats import qmc

            n0 = int(math.ceil((1.0 / max(r, 1e-6)) ** d))
            n0 = min(n0, max(16, self._max_points // 4))
            # next power of two for random_base2 balance
            power = max(4, int(math.ceil(math.log2(max(n0, 2)))))
            n0 = 2**power
            rng = self._generation_rng(0)
            seed_int = int(rng.integers(0, 2**31 - 1))
            engine = qmc.Sobol(d=d, scramble=True, seed=seed_int)
            raw = engine.random_base2(m=power)
            us = [clip_unit_cube(np.asarray(row, dtype=np.float64)) for row in raw]
        # selection filter
        accepted: list[np.ndarray] = []
        for u in us:
            physical = map_u_to_physical(u, self.vars)
            if self._selectionexp and not evaluate_selection(
                self._selectionexp,
                physical,
                context=self._expression_context,
            ):
                continue
            accepted.append(u)
        return accepted

    # ---------------------------------------------------------------- refinement
    def _priority_key(self, edge: tuple[int, int], f: list[float | None], points_u: np.ndarray) -> tuple:
        i, j = edge
        length = float(np.linalg.norm(points_u[i] - points_u[j]))
        fi, fj = f[i], f[j]
        df = 0.0
        if fi is not None and fj is not None:
            df = abs(float(fi) - float(fj))
        # longest first, then |Δf|, then index
        return (-length, -df, i, j)

    def _refine(
        self,
        crossing: list[tuple[int, int]],
        *,
        generation: int,
        radius: float,
    ) -> list[np.ndarray]:
        if not crossing or self._graph is None:
            return []
        points_u = np.asarray([p.u for p in self._points], dtype=np.float64)
        f_vals: list[float | None] = [p.f for p in self._points]
        ordered = sorted(
            crossing,
            key=lambda e: self._priority_key(e, f_vals, points_u),
        )
        rng = self._generation_rng(generation)
        existing = list(points_u)
        tree_pts = np.asarray(existing, dtype=np.float64) if existing else None
        from scipy.spatial import cKDTree

        tree = cKDTree(tree_pts) if tree_pts is not None and len(tree_pts) else None
        min_sep = max(radius / 2.0, 1e-9)
        new_us: list[np.ndarray] = []
        budget = self._max_new_per_generation
        for edge in ordered:
            if len(new_us) >= budget:
                break
            if len(self._points) + len(new_us) >= self._max_points:
                break
            i, j = edge
            u_i = points_u[i]
            u_j = points_u[j]
            mid = 0.5 * (u_i + u_j)
            edge_half = 0.5 * float(np.linalg.norm(u_i - u_j))
            ball_r = max(radius, edge_half)
            for _ in range(self._k_ref):
                if len(new_us) >= budget:
                    break
                cand = clip_unit_cube(sample_ball(mid, ball_r, rng, self._dim))
                physical = map_u_to_physical(cand, self.vars)
                if self._selectionexp and not evaluate_selection(
                    self._selectionexp,
                    physical,
                    context=self._expression_context,
                ):
                    continue
                # spacing reject
                if tree is not None:
                    dist, _ = tree.query(cand, k=1)
                    if float(dist) < min_sep:
                        continue
                if new_us:
                    dists = np.linalg.norm(np.asarray(new_us) - cand, axis=1)
                    if float(np.min(dists)) < min_sep:
                        continue
                new_us.append(cand)
                # grow tree lazily every few points
                if tree is None or (len(new_us) % 8 == 0):
                    all_pts = np.vstack([points_u] + new_us) if new_us else points_u
                    tree = cKDTree(all_pts)
        return new_us

    # ------------------------------------------------------------- graph / cross
    def _crossing_edges(self) -> list[tuple[int, int]]:
        assert self._graph is not None
        points_u = np.asarray([p.u for p in self._points], dtype=np.float64)
        if points_u.shape[0] < 2:
            return []
        edges = self._graph.build(points_u)
        t = self._target_value
        crossing: list[tuple[int, int]] = []
        for row in edges:
            i, j = int(row[0]), int(row[1])
            fi, fj = self._points[i].f, self._points[j].f
            # Failed/unknown f treated as crossing (conservative) but excluded from convergence
            if fi is None or fj is None:
                crossing.append((i, j))
                continue
            if (fi - t) * (fj - t) <= 0:
                crossing.append((i, j))
        return crossing

    def _check_converged(self, crossing: list[tuple[int, int]]) -> bool:
        if not crossing:
            return False
        max_len = 0.0
        max_df = 0.0
        has_known = False
        for i, j in crossing:
            fi, fj = self._points[i].f, self._points[j].f
            if fi is None or fj is None:
                continue
            has_known = True
            ui = np.asarray(self._points[i].u, dtype=np.float64)
            uj = np.asarray(self._points[j].u, dtype=np.float64)
            max_len = max(max_len, float(np.linalg.norm(ui - uj)))
            max_df = max(max_df, abs(float(fi) - float(fj)))
        if not has_known:
            return False
        return max_len < self._contour_precision and max_df < self._function_tolerance

    # --------------------------------------------------------------- submit/wait
    def _next_uuid(self) -> str:
        idx = self._accepted_index
        self._accepted_index += 1
        return deterministic_sampler_uuid(
            prefix="alevelset",
            seed=self._seed,
            sample_index=idx,
        )

    def _submit_generation(self, us: list[np.ndarray], *, generation: int) -> int:
        if not us:
            return 0
        if self.redis is None:
            raise RuntimeError("AdaptiveLevelSet requires redis")
        batch: list[Sample] = []
        for u in us:
            if len(self._points) >= self._max_points:
                break
            physical = map_u_to_physical(u, self.vars)
            uuid = self._next_uuid()
            point = LevelSetPoint(
                u=[float(x) for x in u],
                x={k: float(v) for k, v in physical.items()},
                f=None,
                uuid=uuid,
                generation=generation,
            )
            self._uuid_to_index[uuid] = len(self._points)
            self._points.append(point)
            self._pending_uuids.add(uuid)
            sample = self._build_sample(np.asarray(u, dtype=np.float64))
            sample.uuid = uuid
            batch.append(sample)
            if len(batch) >= self._batch_size:
                self._submit_group(batch)
                self._submitted_uuids.extend(s.uuid for s in batch)
                batch = []
        if batch:
            self._submit_group(batch)
            self._submitted_uuids.extend(s.uuid for s in batch)
        return len(us)

    def _wait_generation(self, *, timeout: float) -> None:
        if self.redis is None:
            raise RuntimeError("AdaptiveLevelSet requires redis")
        deadline = time.monotonic() + max(1.0, float(timeout))
        while self._pending_uuids:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(
                    f"AdaptiveLevelSet generation {self._generation} timed out "
                    f"with {len(self._pending_uuids)} pending sample(s)"
                )
            wait = max(1, min(5, int(remaining)))
            record = self.redis.pull_feedback(timeout=wait)
            if record is None:
                continue
            uuid = str(record.get("uuid", ""))
            if uuid not in self._pending_uuids:
                continue
            self._pending_uuids.discard(uuid)
            self._completed_uuids.add(uuid)
            index = self._uuid_to_index.get(uuid)
            if index is None:
                continue
            status = str(record.get("status", "Completed"))
            observables = dict(record.get("observables") or {})
            if status == "Failed":
                self._points[index].f = None
                self._failed_regions.append({"uuid": uuid, "u": list(self._points[index].u)})
                continue
            # merge physical params into observables for expression
            payload = dict(self._points[index].x)
            payload.update(observables)
            self._points[index].f = self._eval_target(payload)

    # ----------------------------------------------------------------- finalize
    def _graph_component_count(self, edges: np.ndarray, n: int) -> int:
        if n == 0:
            return 0
        parent = list(range(n))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def union(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        for row in edges:
            union(int(row[0]), int(row[1]))
        return len({find(i) for i in range(n)})

    def _chain_polylines_2d(self, crossing: list[tuple[int, int]]) -> list[list[list[float]]]:
        """Chain crossing vertices into polylines via shared endpoints (d=2)."""
        if not crossing:
            return []
        # Build crossing vertices on edges
        edge_vertices: dict[tuple[int, int], list[float]] = {}
        for i, j in crossing:
            fi, fj = self._points[i].f, self._points[j].f
            if fi is None or fj is None:
                continue
            ui = np.asarray(self._points[i].u, dtype=np.float64)
            uj = np.asarray(self._points[j].u, dtype=np.float64)
            u_star = interpolate_crossing(ui, uj, float(fi), float(fj), self._target_value)
            a, b = (i, j) if i < j else (j, i)
            edge_vertices[(a, b)] = [float(x) for x in u_star]

        # adjacency of sites via crossing edges
        adj: dict[int, list[int]] = {}
        for i, j in edge_vertices:
            adj.setdefault(i, []).append(j)
            adj.setdefault(j, []).append(i)

        # Walk components: start at degree-1 nodes, else any
        visited_edges: set[tuple[int, int]] = set()
        polylines: list[list[list[float]]] = []

        def edge_key(a: int, b: int) -> tuple[int, int]:
            return (a, b) if a < b else (b, a)

        def walk(start: int) -> list[list[float]]:
            chain: list[list[float]] = []
            prev = -1
            node = start
            while True:
                neighbors = [n for n in adj.get(node, []) if edge_key(node, n) not in visited_edges]
                if not neighbors:
                    break
                # prefer continuing, else first
                nxt = neighbors[0]
                for cand in neighbors:
                    if cand != prev:
                        nxt = cand
                        break
                key = edge_key(node, nxt)
                visited_edges.add(key)
                chain.append(edge_vertices[key])
                prev, node = node, nxt
                # stop at degree 1 after first step
                if prev != start and len(adj.get(node, [])) == 1:
                    break
                # cycle detection
                if len(chain) > len(edge_vertices) + 2:
                    break
            return chain

        starts = [n for n, ns in adj.items() if len(ns) == 1] or list(adj.keys())
        for start in starts:
            unused = any(edge_key(start, n) not in visited_edges for n in adj.get(start, []))
            if not unused:
                continue
            poly = walk(start)
            if len(poly) >= 1:
                polylines.append(poly)

        # leftover edges
        for key, vertex in edge_vertices.items():
            if key not in visited_edges:
                polylines.append([vertex])
                visited_edges.add(key)
        return polylines

    def _crossing_point_cloud(self, crossing: list[tuple[int, int]]) -> list[dict[str, Any]]:
        cloud: list[dict[str, Any]] = []
        for i, j in crossing:
            fi, fj = self._points[i].f, self._points[j].f
            if fi is None or fj is None:
                continue
            ui = np.asarray(self._points[i].u, dtype=np.float64)
            uj = np.asarray(self._points[j].u, dtype=np.float64)
            u_star = interpolate_crossing(ui, uj, float(fi), float(fj), self._target_value)
            physical = map_u_to_physical(u_star, self.vars)
            cloud.append(
                {
                    "u": [float(x) for x in u_star],
                    "x": {k: float(v) for k, v in physical.items()},
                }
            )
        return cloud

    def _slice_projections(self, cloud: list[dict[str, Any]]) -> list[dict[str, Any]]:
        names = [str(v.name) for v in self.vars]
        if self._slice_pairs:
            pairs = self._slice_pairs
        else:
            pairs = list(itertools.combinations(names, 2))
        name_to_idx = {name: idx for idx, name in enumerate(names)}
        projections: list[dict[str, Any]] = []
        for a, b in pairs:
            if a not in name_to_idx or b not in name_to_idx:
                continue
            ia, ib = name_to_idx[a], name_to_idx[b]
            pts = []
            for item in cloud:
                u = item["u"]
                pts.append(
                    {
                        a: float(u[ia]),
                        b: float(u[ib]),
                        "u_full": list(u),
                        "x": dict(item["x"]),
                    }
                )
            projections.append({"axes": [a, b], "points": pts})
        return projections

    def finalize(self, crossing: list[tuple[int, int]]) -> dict[str, Any]:
        cloud = self._crossing_point_cloud(crossing)
        points_u = np.asarray([p.u for p in self._points], dtype=np.float64)
        edges = (
            self._graph.build(points_u)
            if self._graph is not None and points_u.shape[0] >= 2
            else np.zeros((0, 2), dtype=np.int64)
        )
        components = self._graph_component_count(edges, points_u.shape[0])
        mode = "delaunay" if isinstance(self._graph, DelaunayGraph) else "knn"
        fidelity = "high" if self._dim <= 3 else "proximity-approximate"
        payload: dict[str, Any] = {
            "dim": self._dim,
            "mode": mode,
            "target_expression": self._target_expression,
            "target_value": self._target_value,
            "precision_achieved": self._contour_precision if self._converged else None,
            "function_tolerance": self._function_tolerance,
            "n_points_total": len(self._points),
            "n_generations": self._generation + 1,
            "converged": self._converged,
            "stop_reason": self._stop_reason,
            "failed_regions": list(self._failed_regions),
            "fidelity": fidelity,
            "graph_components": components,
            "crossing_points_u": [item["u"] for item in cloud],
            "crossing_points_x": [item["x"] for item in cloud],
        }
        if self._dim == 2:
            payload["polylines_u"] = self._chain_polylines_2d(crossing)
            payload["polylines_x"] = [
                [map_u_to_physical(np.asarray(pt, dtype=np.float64), self.vars) for pt in poly]
                for poly in payload["polylines_u"]
            ]
        if self._dim >= 4:
            payload["slice_projections"] = self._slice_projections(cloud)
        self._levelset = payload
        # write levelset.json
        task_root = str(
            self.config.get("task_result_dir")
            or (self.config.get("Runtime") or {}).get("task_result_dir")
            or os.getcwd()
        )
        # also try info-style keys from sampler config
        out_dir = str(self.config.get("task_result_dir") or os.getcwd())
        try:
            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, "levelset.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            self._logger.info("wrote level-set payload → %s", path)
        except OSError as exc:
            self._logger.warning("failed to write levelset.json: %s", exc)
        return payload

    # ------------------------------------------------------------------- driver
    def at_safe_barrier(self) -> bool:
        return not self._pending_uuids

    def run_adaptive(self, *, timeout: float = 3600.0) -> int:
        """Control-process entry: generation loop with feedback barriers.

        Resume-safe (D10.4 §7): if ``_points`` is already populated (checkpoint
        import), do **not** re-run generation-0. Finish any in-flight
        ``_pending_uuids`` first, then continue refine from the restored barrier.
        """
        if self.redis is None:
            raise RuntimeError("AdaptiveLevelSet.run_adaptive requires redis")
        if self._seed_seq is None:
            self._seed_seq = np.random.SeedSequence(self._seed)
        if self._graph is None:
            self._graph = self._make_graph()
        if self._target_fn is None:
            self._compile_target()

        # Already finished in a prior process (checkpointed after converge).
        if self._converged and self._levelset is not None and not self._pending_uuids:
            return 0

        total_submitted = 0
        if not self._points:
            gen0 = self._generation_0_points()
            self._generation = 0
            total_submitted += self._submit_generation(gen0, generation=0)
            self._wait_generation(timeout=timeout)
        elif self._pending_uuids:
            # Mid-generation resume: repropose_unfinished already re-queued tasks.
            self._wait_generation(timeout=timeout)

        while True:
            crossing = self._crossing_edges()
            if not crossing:
                self._stop_reason = "level-set not present in domain"
                self._converged = False
                self.finalize(crossing)
                return total_submitted
            if self._check_converged(crossing):
                self._converged = True
                self._stop_reason = "converged"
                self.finalize(crossing)
                return total_submitted
            if len(self._points) >= self._max_points:
                self._stop_reason = "max_points"
                self.finalize(crossing)
                return total_submitted
            next_gen = self._generation + 1
            if next_gen > self._max_generations:
                self._stop_reason = "max_generations"
                self._converged = self._check_converged(crossing)
                self.finalize(crossing)
                return total_submitted
            radius = self._initial_radius * (self._refinement_factor ** next_gen)
            new_us = self._refine(crossing, generation=next_gen, radius=radius)
            if not new_us:
                self._stop_reason = "no_refinement_candidates"
                self.finalize(crossing)
                return total_submitted
            self._generation = next_gen
            total_submitted += self._submit_generation(new_us, generation=next_gen)
            self._wait_generation(timeout=timeout)

    def run_distributed(self) -> int:
        """Alias so generic drivers can call run_distributed when registered."""
        return self.run_adaptive()

    def repropose_unfinished(self) -> list[str]:
        if not self._pending_uuids:
            return []
        requeued: list[str] = []
        for uuid in sorted(self._pending_uuids):
            index = self._uuid_to_index.get(uuid)
            if index is None:
                continue
            u = np.asarray(self._points[index].u, dtype=np.float64)
            sample = self._build_sample(u)
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def export_runtime_state(self) -> dict[str, Any]:
        return {
            "points": [asdict(p) for p in self._points],
            "dim": self._dim,
            "generation": self._generation,
            "pending_uuids": sorted(self._pending_uuids),
            "accepted_index": self._accepted_index,
            "seed": self._seed,
            "converged": self._converged,
            "levelset": self._levelset,
            "failed_regions": list(self._failed_regions),
            "stop_reason": self._stop_reason,
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "control_state": self._checkpoint_control_state(),
            "uuid_to_index": dict(self._uuid_to_index),
        }

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        raw_points = state.get("points") or []
        self._points = [
            LevelSetPoint(
                u=list(item.get("u") or []),
                x=dict(item.get("x") or {}),
                f=item.get("f"),
                uuid=str(item.get("uuid", "")),
                generation=int(item.get("generation", 0) or 0),
            )
            for item in raw_points
            if isinstance(item, Mapping)
        ]
        self._dim = int(state.get("dim", self._dim) or self._dim)
        self._generation = int(state.get("generation", 0) or 0)
        self._pending_uuids = set(str(u) for u in (state.get("pending_uuids") or []))
        self._accepted_index = int(state.get("accepted_index", 0) or 0)
        self._seed = int(state.get("seed", self._seed) or self._seed)
        self._seed_seq = np.random.SeedSequence(self._seed)
        self._converged = bool(state.get("converged", False))
        self._levelset = state.get("levelset")
        self._failed_regions = list(state.get("failed_regions") or [])
        self._stop_reason = str(state.get("stop_reason") or "")
        self._uuid_to_index = {
            str(k): int(v) for k, v in dict(state.get("uuid_to_index") or {}).items()
        }
        if not self._uuid_to_index:
            self._uuid_to_index = {p.uuid: i for i, p in enumerate(self._points)}
        self._import_checkpoint_control_state(state)
        if self._target_expression:
            self._compile_target()
        self._graph = self._make_graph()


__all__ = [
    "AdaptiveLevelSetSampler",
    "DelaunayGraph",
    "KNNGraph",
    "LevelSetPoint",
    "NeighborGraph",
    "clip_unit_cube",
    "interpolate_crossing",
    "sample_ball",
]
