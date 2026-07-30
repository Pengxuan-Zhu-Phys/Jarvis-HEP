#!/usr/bin/env python3
"""AdaptiveBridson sampler: outer/core densify + t_min/t_max gen gate.

Design (binding):
``Jarvis-Books/Jarvis-HEP V2/DESIGN_ADAPTIVE_BRIDSON_LIVE_BAND.md``

Control flow (intended):
  gen-0  — global Bridson with ``initial_radius``
  at fixed r_g (= one generation scale):
    1. while actionable fill is making progress (new core / closer best),
       bridge / root-corr / densify at the same scale
       (fill_pass++; generation does NOT advance)
       Persistent straddles alone do not block scale quiescence.
    2. when fill is done → inspect the 2*r_g scale band (t_min, t_max)
    3. if (t_max − t_min) < threshold → EXIT converged
       else shrink r_g, generation += 1, propose next scale
  densify centers still use absolute outer/core |f−T| bands;
  t_min/t_max (all evaluated sites within 2*r_g of best) gate next-gen and exit.

Public tuning surface:
  ``outer_half_width`` controls discovery support in f-space;
  ``min_radius`` controls the final u-space resolution.  The inner core/stop
  width is derived as ``outer_half_width / 8``.  Historical tuning keys remain
  readable for checkpoint/config compatibility, but are no longer required.
"""

from __future__ import annotations

import itertools
import json
import math
import os
import time
from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
from enum import Enum
from typing import Any, Protocol

import numpy as np

from jarvishep2.Sampling.bridson import Bridson_sampling, hypersphere_surface_sample
from jarvishep2.Sampling.feedback_sampler import FeedbackSampler
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
    """One evaluated (or pending) site in unit-cube coordinates."""

    u: list[float]
    x: dict[str, float]
    f: float | None
    uuid: str = ""
    generation: int = 0


class LoopAction(str, Enum):
    """One explicit control-flow outcome for an adaptive loop iteration."""

    CONTINUE = "continue"
    ADVANCE = "advance"
    STOP = "stop"


@dataclass(frozen=True)
class LoopDecision:
    """A named loop transition; ``STOP`` decisions always carry a reason."""

    action: LoopAction
    reason: str | None = None

    @classmethod
    def continue_fill(cls) -> "LoopDecision":
        return cls(LoopAction.CONTINUE)

    @classmethod
    def advance(cls) -> "LoopDecision":
        return cls(LoopAction.ADVANCE)

    @classmethod
    def stop(cls, reason: str) -> "LoopDecision":
        return cls(LoopAction.STOP, reason)


@dataclass(frozen=True)
class ScaleState:
    """Band state prepared once at the start of a radius-scale iteration."""

    cores: list[int]
    frontier: list[int]
    outer: list[int]
    best: int | None
    t_min: float | None
    t_max: float | None


class NeighborGraph(Protocol):
    def build(self, points_u: np.ndarray) -> np.ndarray:
        """Undirected edges as shape (m, 2) int array of site indices."""
        ...


class DelaunayGraph:
    """Site adjacency via Voronoi ridge_points (exact for d ≤ 3)."""

    def build(self, points_u: np.ndarray) -> np.ndarray:
        from scipy.spatial import Voronoi

        n = int(points_u.shape[0])
        d = int(points_u.shape[1]) if points_u.ndim == 2 else 0
        if n < 2:
            return np.zeros((0, 2), dtype=np.int64)
        # Qhull needs ≥ d+1 points for a simplex; with n == d+1 Voronoi can
        # still fail (degenerate). Fall back to complete graph for tiny clouds.
        if n <= d + 1:
            return np.array(
                list(itertools.combinations(range(n), 2)),
                dtype=np.int64,
            )
        try:
            vor = Voronoi(points_u, qhull_options="QJ")
        except Exception:
            return np.array(
                list(itertools.combinations(range(n), 2)),
                dtype=np.int64,
            )
        ridges = np.asarray(vor.ridge_points, dtype=np.int64)
        if ridges.size == 0:
            return np.zeros((0, 2), dtype=np.int64)
        return ridges


class KNNGraph:
    """Symmetric kNN graph (d = 4–5 proximity mode)."""

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
            for j in np.atleast_1d(idx[i]):
                j = int(j)
                if j == i:
                    continue
                a, b = (i, j) if i < j else (j, i)
                edges.add((a, b))
        if not edges:
            return np.zeros((0, 2), dtype=np.int64)
        return np.array(sorted(edges), dtype=np.int64)


def sample_ball(
    center: np.ndarray, radius: float, rng: np.random.Generator, d: int
) -> np.ndarray:
    """Uniform sample in the d-ball of given radius around center."""
    if radius <= 0:
        return np.asarray(center, dtype=np.float64).copy()
    direction = rng.normal(size=d)
    norm = float(np.linalg.norm(direction))
    if norm < 1e-15:
        direction = np.ones(d)
        norm = float(np.linalg.norm(direction))
    direction /= norm
    scale = radius * (float(rng.random()) ** (1.0 / max(d, 1)))
    return np.asarray(center, dtype=np.float64) + direction * scale


def clip_unit_cube(u: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(u, dtype=np.float64), _CLIP_LO, _CLIP_HI)


def _finite_f(value: Any) -> float | None:
    if value is None:
        return None
    try:
        x = float(value)
    except (TypeError, ValueError):
        return None
    if x != x or x in (float("inf"), float("-inf")):
        return None
    return x


class _SepIndex:
    """Uniform grid for O(1) Bridson exclusion queries (d ≤ 5).

    Replaces the O(N) scan in ``_try_accept_candidate`` that dominated late
    densify generations (thousands of live cores × k_ref trials).
    """

    __slots__ = ("r_sep", "r2", "cell", "dim", "buckets", "points")

    def __init__(self, r_sep: float, dim: int) -> None:
        self.r_sep = float(r_sep)
        self.r2 = self.r_sep * self.r_sep
        # Cell size = r_sep: only the 3^d neighborhood can violate exclusion.
        self.cell = max(self.r_sep, 1e-15)
        self.dim = int(dim)
        self.buckets: dict[tuple[int, ...], list[int]] = {}
        self.points: list[np.ndarray] = []

    def _key(self, u: np.ndarray) -> tuple[int, ...]:
        c = self.cell
        return tuple(int(math.floor(float(x) / c)) for x in u)

    def add(self, u: np.ndarray) -> None:
        p = np.asarray(u, dtype=np.float64)
        idx = len(self.points)
        self.points.append(p)
        key = self._key(p)
        bucket = self.buckets.get(key)
        if bucket is None:
            self.buckets[key] = [idx]
        else:
            bucket.append(idx)

    def far_enough(self, cand: np.ndarray) -> bool:
        if self.r_sep <= 0:
            return True
        key = self._key(cand)
        r2 = self.r2
        dim = self.dim
        # Iterate the Moore neighborhood without building a full product list.
        # d=2 → 9 cells; d=5 → 243 cells (still cheap vs scanning thousands).
        ranges = [(-1, 0, 1)] * dim
        from itertools import product

        for off in product(*ranges):
            nb = tuple(key[i] + off[i] for i in range(dim))
            bucket = self.buckets.get(nb)
            if not bucket:
                continue
            for j in bucket:
                dlt = cand - self.points[j]
                if float(np.dot(dlt, dlt)) < r2:
                    return False
        return True


def _thin_centers(
    positions: np.ndarray,
    *,
    min_spacing: float,
    order: Sequence[int] | None = None,
) -> list[int]:
    """Greedy spatial thinning: keep indices at least ``min_spacing`` apart.

    Used to pick densify *proposal centers* from an over-dense live-core set.
    All cores still participate in exclusion; only centers open r_g windows.
    """
    n = int(positions.shape[0])
    if n == 0:
        return []
    if min_spacing <= 0 or n == 1:
        return list(range(n)) if order is None else list(order)
    seq = list(order) if order is not None else list(range(n))
    index = _SepIndex(min_spacing, int(positions.shape[1]))
    kept: list[int] = []
    for i in seq:
        p = positions[int(i)]
        if index.far_enough(p):
            index.add(p)
            kept.append(int(i))
    return kept


class AdaptiveBridsonSampler(FeedbackSampler):
    """Live-band AdaptiveBridson (Method: AdaptiveBridson).

    See DESIGN_ADAPTIVE_BRIDSON_LIVE_BAND.md. Config block:
    ``Sampling.AdaptiveBridson``.
    """

    method = "AdaptiveBridson"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.adaptive_bridson")
        self.vars: list = []
        self._dim = 0
        self._target_expression = ""
        self._target_value = 0.0
        self._initial_radius = 0.10
        self._min_radius = 1.0 / 200.0
        self._refinement_factor = 0.5
        # Absolute half-widths in f-space (not Voronoi neighbor span).
        self._core_half_width = 0.0025  # inner band |f-T| ≤ w_core → densify centers
        self._outer_half_width = 0.02  # discovery band (anneals → w_core)
        # Optional final contour half-width (|f-T| ≤ final). None → use core_half_width.
        self._final_half_width: float | None = None
        self._final_half_width_configured: bool = False
        self._outer_shrink_factor = 0.7  # w_outer := max(w_core, factor * w_outer)
        self._core_spacing_factor = 2.0  # coverage: max core NN gap ≤ factor * r_g
        self._min_cores_for_coverage = 4  # need a few cores before scale is "covered"
        # r_g shrink only after *coverage* at this scale (append cores first).
        self._radius_shrink_mode = "on_coverage"  # on_coverage | on_fill | every_generation
        # Compat alias: threshold ≈ core stop width (also default w_core if unset).
        self._threshold: float | None = 0.0025
        self._max_generations = 25
        self._max_points = 5000
        self._max_new_per_generation = 500
        self._k_ref = 30
        self._knn_k = 8
        self._neighbor_graph_mode = "auto"
        # Gap bridging: reconnect *core* islands along the band (u-space).
        self._bridge_gaps = True
        self._bridge_span_factor = 2.5  # max core–core link = factor * r_window
        self._selectionexp: str | None = None
        self._graph: NeighborGraph | None = None
        self._points: list[LevelSetPoint] = []
        self._accepted_index = 0
        self._converged = False
        self._levelset: dict[str, Any] | None = None
        self._target_fn: CompiledExpression | None = None
        self._failed_regions: list[dict[str, Any]] = []
        self._uuid_to_index: dict[str, int] = {}
        self._stop_reason = ""
        self._radius = 0.10  # current generation window radius r_g
        self._t_min: float | None = None
        self._t_max: float | None = None
        self._live_core_indices: list[int] = []
        self._outer_indices: list[int] = []
        self._frontier_indices: list[int] = []  # outer \ core (brackets)
        self._open_brackets: int = 0
        self._bracket_support_indices: list[int] = []
        self._virtual_core_anchors: list[list[float]] = []
        self._active_virtual_cores: dict[str, dict[str, Any]] = {}
        self._retired_virtual_core_ids: set[str] = set()
        self._virtual_core_core_signature: tuple[str, ...] | None = None
        # Cores accumulate at fixed r_g; solidify only on r_g shrink.
        # Once a sample is a core it stays in the permanent core set (expand only).
        self._solidified_core_uuids: set[str] = set()
        self._accumulated_core_uuids: set[str] = set()
        self._n_radius_refines: int = 0
        # generation == number of r_g shrinks completed (scale index).
        # fill_pass counts same-r_g densify/fill rounds and does NOT advance generation.
        self._fill_pass: int = 0
        self._max_fill_passes: int = 64  # safety against infinite same-r_g fill
        self._last_densify_new: int = 0
        self._last_densify_attempted: bool = False
        # The control cloud may be pruned, but submission history must not be.
        # Otherwise a same-scale fill pass can resubmit an identical u with a
        # fresh UUID after the original point falls outside the live band.
        self._submitted_u_keys: set[tuple[float, ...]] = set()
        # A/B straddles are persistent level-set geometry, not a termination
        # condition.  Stop filling a scale after several passes without a new
        # core or a meaningfully closer best point, then let t_min/t_max decide
        # convergence versus the next (smaller-r_g) generation.
        self._quiet_fill_passes: int = 0
        self._quiet_fill_limit: int = 3
        # Persistent same-r_g geometric fronts.  A front is not discarded
        # after one unlucky/off-level proposal: it is retired only after it is
        # superseded by a core farther along its direction, reaches the cube
        # boundary, or exhausts several distinct directional attempts.
        self._active_endpoints: dict[str, dict[str, Any]] = {}
        self._retired_endpoint_ids: set[str] = set()
        self._endpoint_core_signature: tuple[str, ...] | None = None
        self._endpoint_stall_limit: int = 12
        self._endpoint_omni_probes: int = 16
        self._slice_pairs: list[tuple[str, str]] | None = None

    # ------------------------------------------------------------------ config
    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._dim = len(self.vars)
        if self._dim < 2 or self._dim > 5:
            raise ValueError(
                f"AdaptiveBridson requires 2 ≤ dim ≤ 5, got {self._dim}"
            )
        block = sampling.get("AdaptiveBridson") or sampling.get("adaptive_bridson") or {}
        if not isinstance(block, Mapping):
            raise ValueError("Sampling.AdaptiveBridson must be a mapping")
        self._target_expression = str(block.get("target_expression", "")).strip()
        if not self._target_expression:
            raise ValueError("AdaptiveBridson.target_expression is required")
        if "target_value" not in block:
            raise ValueError("AdaptiveBridson.target_value is required")
        self._target_value = float(block["target_value"])
        # Public geometry controls are intentionally limited to
        # outer_half_width and min_radius.  Keep historical keys readable so
        # old scans/checkpoints resume safely, but make the no-legacy path
        # reproduce the canonical iDM settings without YAML boilerplate.
        self._initial_radius = float(block.get("initial_radius", 0.10) or 0.10)
        if self._initial_radius <= 0:
            raise ValueError("AdaptiveBridson.initial_radius must be > 0")
        self._min_radius = float(
            block.get("min_radius", 1.0 / 200.0) or 1.0 / 200.0
        )
        if self._min_radius <= 0:
            raise ValueError("AdaptiveBridson.min_radius must be > 0")
        # Starting below the requested floor means there is simply no radius
        # refinement to perform; use the initial scale as the effective floor.
        self._min_radius = min(self._min_radius, self._initial_radius)
        self._refinement_factor = float(block.get("refinement_factor", 0.5) or 0.5)
        if not (0.0 < self._refinement_factor < 1.0):
            raise ValueError(
                "AdaptiveBridson.refinement_factor must be in (0, 1)"
            )
        # --- absolute outer / core half-widths in f ---
        # outer_half_width is public.  The canonical inner band and stop width
        # are derived from it; threshold/core_half_width remain compatibility
        # overrides for configurations written before the two-knob interface.
        outer_half_width_configured = (
            "outer_half_width" in block
            and block.get("outer_half_width") is not None
        )
        if outer_half_width_configured:
            self._outer_half_width = float(block["outer_half_width"])
        else:
            self._outer_half_width = 0.02
        if self._outer_half_width <= 0:
            raise ValueError("AdaptiveBridson.outer_half_width must be > 0")

        derived_core_half_width = self._outer_half_width / 8.0

        # threshold / function_tolerance → legacy core stop-width override.
        if "threshold" in block and block.get("threshold") is not None:
            self._threshold = float(block["threshold"])
            if self._threshold <= 0:
                raise ValueError("AdaptiveBridson.threshold must be > 0")
        elif "function_tolerance" in block and block.get("function_tolerance") is not None:
            self._threshold = float(block["function_tolerance"])
        else:
            self._threshold = None

        if "core_half_width" in block and block.get("core_half_width") is not None:
            self._core_half_width = float(block["core_half_width"])
        elif self._threshold is not None:
            self._core_half_width = float(self._threshold)
        else:
            self._core_half_width = derived_core_half_width
        if self._core_half_width <= 0:
            raise ValueError("AdaptiveBridson.core_half_width must be > 0")

        # Historical cards often supplied threshold/core_half_width alone.
        # Preserve their old implicit outer band while the public two-knob
        # path derives the core in the opposite direction (outer / 8).
        if not outer_half_width_configured and (
            self._threshold is not None
            or ("core_half_width" in block and block.get("core_half_width") is not None)
        ):
            self._outer_half_width = max(0.02, 8.0 * self._core_half_width)

        if self._outer_half_width < self._core_half_width:
            raise ValueError(
                "AdaptiveBridson.outer_half_width must be >= core_half_width"
            )

        # final_half_width: optional tighter export/stop band (null = use core).
        if "final_half_width" in block and block.get("final_half_width") is not None:
            self._final_half_width = float(block["final_half_width"])
            self._final_half_width_configured = True
            if self._final_half_width <= 0:
                raise ValueError("AdaptiveBridson.final_half_width must be > 0")
            if self._final_half_width > self._core_half_width + 1e-15:
                raise ValueError(
                    "AdaptiveBridson.final_half_width must be <= core_half_width "
                    f"(got {self._final_half_width:g} > core {self._core_half_width:g})"
                )
        else:
            self._final_half_width = None
            self._final_half_width_configured = False

        self._outer_shrink_factor = float(
            block.get("outer_shrink_factor", 0.7) or 0.7
        )
        if not (0.0 < self._outer_shrink_factor <= 1.0):
            raise ValueError(
                "AdaptiveBridson.outer_shrink_factor must be in (0, 1]"
            )
        self._core_spacing_factor = float(
            block.get("core_spacing_factor", 2.0) or 2.0
        )
        if self._core_spacing_factor <= 0:
            raise ValueError("AdaptiveBridson.core_spacing_factor must be > 0")

        mode = str(
            block.get("radius_shrink_mode", "on_coverage") or "on_coverage"
        ).strip().lower()
        if mode in {
            "on_coverage",
            "on-coverage",
            "coverage",
            "on_fill",
            "on-fill",
            "fill",
            "gated",
        }:
            # on_fill kept as alias → coverage-gated (append cores, then shrink).
            self._radius_shrink_mode = "on_coverage"
        elif mode in {"every_generation", "every-generation", "always", "legacy"}:
            self._radius_shrink_mode = "every_generation"
        else:
            raise ValueError(
                "AdaptiveBridson.radius_shrink_mode must be "
                "'on_coverage' (default), 'on_fill' (alias), or 'every_generation'"
            )
        self._min_cores_for_coverage = max(
            2,
            int(block.get("min_cores_for_coverage", self._min_cores_for_coverage) or 4),
        )

        if "max_generations" in block and block.get("max_generations") is not None:
            # Max r_g scale index (generation := number of r_g shrinks done).
            self._max_generations = max(0, int(block["max_generations"]))
        else:
            self._max_generations = 16
        if self._threshold is None:
            self._threshold = float(self._core_half_width)
        default_max_points = 50000
        self._max_points = max(
            10, int(block.get("max_points", default_max_points) or default_max_points)
        )
        self._max_new_per_generation = max(
            1,
            int(
                block.get(
                    "max_new_per_generation",
                    4000,
                )
                or 4000
            ),
        )
        # Same-r_g densify/fill rounds before giving up (not counted as generations).
        self._max_fill_passes = max(
            8,
            int(
                block.get(
                    "max_fill_passes",
                    max(64, 4 * max(1, self._max_generations)),
                )
                or max(64, 4 * max(1, self._max_generations))
            ),
        )
        self._quiet_fill_limit = max(
            1, int(block.get("quiet_fill_passes", 3) or 3)
        )
        self._endpoint_stall_limit = max(
            3, int(block.get("endpoint_stall_passes", 12) or 12)
        )
        self._endpoint_omni_probes = max(
            2 * self._dim,
            int(block.get("endpoint_omni_probes", max(16, 4 * self._dim)) or 16),
        )
        # Bridge disconnected live-core islands along the contour.
        # Default on: without this, gen-0 hits densify into separate blobs.
        if "bridge_gaps" in block:
            self._bridge_gaps = bool(block.get("bridge_gaps"))
        else:
            self._bridge_gaps = True
        self._bridge_span_factor = float(
            block.get("bridge_span_factor", self._bridge_span_factor) or 2.5
        )
        if self._bridge_span_factor < 1.0:
            raise ValueError(
                "AdaptiveBridson.bridge_span_factor must be >= 1 "
                "(max core–core bridge length = factor * r_window)"
            )
        self._k_ref = max(1, int(block.get("k_ref", 30) or 30))
        self._knn_k = max(1, int(block.get("knn_k", 4 * self._dim) or (4 * self._dim)))
        self._neighbor_graph_mode = str(
            block.get("neighbor_graph", "auto")
        ).strip().lower()
        self._selectionexp = sampling.get("selection")
        seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._init_seed_sequence(seed)
        self._compile_target()
        self._graph = self._make_graph()
        self._radius = self._initial_radius
        raw_pairs = block.get("slice_pairs")
        if isinstance(raw_pairs, list) and raw_pairs:
            pairs: list[tuple[str, str]] = []
            for item in raw_pairs:
                if isinstance(item, (list, tuple)) and len(item) == 2:
                    pairs.append((str(item[0]), str(item[1])))
            self._slice_pairs = pairs or None
        else:
            self._slice_pairs = None
        var_names = ", ".join(str(v.name) for v in self.vars) or "(none)"
        graph_mode = self._neighbor_graph_mode
        if graph_mode == "auto":
            graph_mode = f"auto→{'delaunay' if self._dim <= 3 else 'knn'}"
        self._logger.debug(
            "AdaptiveBridson configured (outer/core + root-correction, u-space)\n"
            "    dim / vars          -> %d  [%s]\n"
            "    target              -> %s = %g\n"
            "    neighbor graph      -> %s  (knn_k=%d)\n"
            "    initial_radius r0   -> %g\n"
            "    minimum radius     -> %g  (u-space Euclidean floor)\n"
            "    refinement_factor   -> %g  (r_g shrink multiplier; blue-noise sep=r_g; mode=%s)\n"
            "    outer / core half-w -> %g / %g  (outer shrink×%g)\n"
            "    bridge_gaps         -> %s  (span_factor=%g · r_g)\n"
            "    coverage            -> gap≤%g·r_g, min_cores=%d\n"
            "    max_generations     -> %d  (= max r_g shrinks / scale index)\n"
            "    max_fill_passes     -> %d  (same-r_g densify rounds, not gens)\n"
            "    quiet_fill_passes   -> %d  (no-progress passes before r_g shrink)\n"
            "    endpoint stall      -> %d  (failed same-r_g fronts before retirement)\n"
            "    endpoint omni probes-> %d  (full hypersphere per active front/pass)\n"
            "    max_points          -> %d  (max_new/pass=%d, k_ref=%d)\n"
            "    selection           -> %s",
            self._dim,
            var_names,
            self._target_expression,
            self._target_value,
            graph_mode,
            self._knn_k,
            self._initial_radius,
            self._min_radius,
            self._refinement_factor,
            self._radius_shrink_mode,
            self._outer_half_width,
            self._core_half_width,
            self._outer_shrink_factor,
            "on" if self._bridge_gaps else "off",
            self._bridge_span_factor,
            self._core_spacing_factor,
            self._min_cores_for_coverage,
            self._max_generations,
            self._max_fill_passes,
            self._quiet_fill_limit,
            self._endpoint_stall_limit,
            self._endpoint_omni_probes,
            self._max_points,
            self._max_new_per_generation,
            self._k_ref,
            self._selectionexp or "(none)",
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

    def feedback_return_spec(self) -> dict[str, Any]:
        from jarvishep2.feedback_return import resolve_feedback_return

        return resolve_feedback_return(self.config, sampler=None, method=self.method)

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
            return _finite_f(value)
        except Exception:
            return None

    # --------------------------------------------------------------- gen-0
    def _generation_0_points(self) -> list[np.ndarray]:
        d = self._dim
        r = self._initial_radius
        if d <= 4:
            fill = "Bridson"
            rng = self._generation_rng(0)
            seed_int = int(rng.integers(0, 2**31 - 1))
            np.random.seed(seed_int)
            dims = np.ones(d, dtype=np.float64)
            points = Bridson_sampling(
                dims=dims,
                radius=r,
                k=30,
                hypersphere_sample=hypersphere_surface_sample,
            )
            us = [
                clip_unit_cube(np.asarray(row, dtype=np.float64) / dims)
                for row in points
            ]
        else:
            fill = "Sobol"
            from scipy.stats import qmc

            n0 = int(math.ceil((1.0 / max(r, 1e-6)) ** d))
            n0 = min(n0, max(16, self._max_points // 4))
            power = max(4, int(math.ceil(math.log2(max(n0, 2)))))
            rng = self._generation_rng(0)
            seed_int = int(rng.integers(0, 2**31 - 1))
            engine = qmc.Sobol(d=d, scramble=True, seed=seed_int)
            raw = engine.random_base2(m=power)
            us = [clip_unit_cube(np.asarray(row, dtype=np.float64)) for row in raw]
        accepted: list[np.ndarray] = []
        n_reject = 0
        for u in us:
            physical = map_u_to_physical(u, self.vars)
            if self._selectionexp and not evaluate_selection(
                self._selectionexp,
                physical,
                context=self._expression_context,
            ):
                n_reject += 1
                continue
            accepted.append(u)
        self._logger.debug(
            "AdaptiveBridson gen-0 fill (u-space)\n"
            "    method               -> %s  (r0=%g, d=%d)\n"
            "    raw / accepted       -> %d / %d\n"
            "    selection rejected   -> %d",
            fill,
            r,
            d,
            len(us),
            len(accepted),
            n_reject,
        )
        return accepted

    # ----------------------------------------------------------- neighbors / band
    def _adjacency(self) -> dict[int, list[int]]:
        assert self._graph is not None
        n = len(self._points)
        if n < 2:
            return {i: [] for i in range(n)}
        points_u = np.asarray([p.u for p in self._points], dtype=np.float64)
        edges = self._graph.build(points_u)
        adj: dict[int, list[int]] = {i: [] for i in range(n)}
        for row in edges:
            i, j = int(row[0]), int(row[1])
            if 0 <= i < n and 0 <= j < n and i != j:
                adj[i].append(j)
                adj[j].append(i)
        for i in adj:
            adj[i] = sorted(set(adj[i]))
        return adj

    def _best_index(self) -> int | None:
        best_i: int | None = None
        best_score = float("inf")
        t = self._target_value
        for i, p in enumerate(self._points):
            fi = _finite_f(p.f)
            if fi is None:
                continue
            score = abs(fi - t)
            if score < best_score - 1e-15 or (
                abs(score - best_score) <= 1e-15 and (best_i is None or i < best_i)
            ):
                best_score = score
                best_i = i
        return best_i

    def _abs_dev(self, index: int) -> float | None:
        fi = _finite_f(self._points[index].f)
        if fi is None:
            return None
        return abs(fi - self._target_value)

    def _classify_bands(
        self,
    ) -> tuple[list[int], list[int], list[int], int | None]:
        """Split points by absolute distance to target.

        Returns
        -------
        cores : |f−T| ≤ w_core  (densify / bridge centers)
        frontier : w_core < |f−T| ≤ w_outer  (brackets / support only)
        outer : cores ∪ frontier
        best : index of min |f−T| among finite f (may be outside outer)
        """
        t = self._target_value
        w_core = float(self._core_half_width)
        w_outer = float(self._outer_half_width)
        cores: list[int] = []
        frontier: list[int] = []
        best = self._best_index()
        for i, p in enumerate(self._points):
            fi = _finite_f(p.f)
            if fi is None:
                continue
            d = abs(fi - t)
            if d <= w_core + 1e-15:
                cores.append(i)
            elif d <= w_outer + 1e-15:
                frontier.append(i)
        outer = sorted(set(cores) | set(frontier))
        return cores, frontier, outer, best

    def _radius_t_band(
        self,
    ) -> tuple[int | None, float | None, float | None]:
        """t_min/t_max from evaluated sites in the best point's 2*r_g ball.

        This is the **generation / exit gate**, not the densify-center mask.
        densify centers still use absolute outer/core |f−T| bands.

        The factor two is required by the strict r_g blue-noise invariant:
        an r_g ball would otherwise contain only the best point (apart from a
        measure-zero boundary), making its target width spuriously zero.
        """
        best = self._best_index()
        if best is None:
            return None, None, None
        best_u = np.asarray(self._points[best].u, dtype=np.float64)
        radius = 2.0 * max(0.0, float(self._radius))
        vals: list[float] = []
        for point in self._points:
            fj = _finite_f(point.f)
            if fj is not None:
                u = np.asarray(point.u, dtype=np.float64)
                if float(np.linalg.norm(u - best_u)) <= radius + 1e-15:
                    vals.append(float(fj))
        # A singleton has zero numerical width but carries no local variation
        # information, so it must never satisfy the convergence gate.
        if len(vals) < 2:
            return best, None, None
        return best, float(min(vals)), float(max(vals))

    def _neighbor_t_band(
        self,
    ) -> tuple[int | None, float | None, float | None]:
        """Compatibility alias; the band uses the current 2*r_g scale ball."""
        return self._radius_t_band()

    def _band_width_ok(
        self, t_min: float | None, t_max: float | None
    ) -> bool:
        """True when (t_max − t_min) < threshold (primary exit / no next gen)."""
        if t_min is None or t_max is None:
            return False
        tau = (
            float(self._threshold)
            if self._threshold is not None
            else float(self._core_half_width)
        )
        if tau <= 0:
            return False
        return (float(t_max) - float(t_min)) < tau

    def _is_core_index(self, index: int, core_set: set[int] | None = None) -> bool:
        if core_set is not None:
            return int(index) in core_set
        d = self._abs_dev(int(index))
        return d is not None and d <= float(self._core_half_width) + 1e-15

    def _classify_bridge_edges(
        self,
        bracket_edges: Sequence[tuple[int, int, float, float]],
        cores: Sequence[int],
        *,
        r_window: float | None = None,
    ) -> dict[str, list[tuple[int, int, float, float]]]:
        """Classify local nearest straddle edges by core endpoints.

        A: neither endpoint is core  → must root-corr / expand toward T
        B: exactly one endpoint core → expand from core along the edge
        C: both cores but u-gap large → spatial bridge (core–core densify)
        D: both cores and close      → no fill on this edge
        """
        core_set = {int(i) for i in cores}
        rw = float(r_window) if r_window is not None else float(self._radius)
        # Spatial gap threshold for core–core bridge (align with densify window).
        gap_lim = max(rw * float(self._bridge_span_factor), 2.0 * rw)
        out: dict[str, list[tuple[int, int, float, float]]] = {
            "A": [],
            "B": [],
            "C": [],
            "D": [],
        }
        for edge in bracket_edges:
            i, j, fi, fj = int(edge[0]), int(edge[1]), float(edge[2]), float(edge[3])
            ci = i in core_set
            cj = j in core_set
            if not ci and not cj:
                out["A"].append((i, j, fi, fj))
            elif ci != cj:
                out["B"].append((i, j, fi, fj))
            else:
                ui = np.asarray(self._points[i].u, dtype=np.float64)
                uj = np.asarray(self._points[j].u, dtype=np.float64)
                dist = float(np.linalg.norm(uj - ui))
                if dist > gap_lim + 1e-15:
                    out["C"].append((i, j, fi, fj))
                else:
                    out["D"].append((i, j, fi, fj))
        return out

    def _register_cores(self, cores: Sequence[int]) -> None:
        """Append-only: any current core uuid is remembered forever."""
        for i in cores:
            uuid = str(self._points[int(i)].uuid or "")
            if uuid:
                self._accumulated_core_uuids.add(uuid)

    def _accumulated_core_indices(
        self,
        *,
        t_min: float | None = None,
        t_max: float | None = None,
        require_in_band: bool = False,
    ) -> list[int]:
        """Indices of accumulated cores still in the control cloud.

        When ``require_in_band`` and (t_min,t_max) are set, only return cores
        whose f lies in [t_min, t_max] — densify must not be driven by cores
        that already left the current 2*r_g target band (prevents sideways
        replacement).
        """
        out: list[int] = []
        w = float(self._core_half_width)
        for i, p in enumerate(self._points):
            fi = _finite_f(p.f)
            is_core_now = fi is not None and abs(fi - self._target_value) <= w + 1e-15
            was_core = bool(p.uuid and p.uuid in self._accumulated_core_uuids)
            if is_core_now and p.uuid:
                self._accumulated_core_uuids.add(str(p.uuid))
            if not (is_core_now or was_core):
                continue
            if require_in_band and t_min is not None and t_max is not None:
                if fi is None or fi < float(t_min) - 1e-15 or fi > float(t_max) + 1e-15:
                    continue
            # Active densify cores must still satisfy |f−T|≤w_core (not stale uuids).
            if require_in_band and not is_core_now:
                continue
            out.append(i)
        return sorted(set(out))

    def _bridge_endpoint_indices(
        self,
        classified: Mapping[str, Sequence[tuple[int, int, float, float]]],
    ) -> list[int]:
        """A/B straddle endpoints only (for edge fill, NOT isotropic densify)."""
        ends: set[int] = set()
        for key in ("A", "B"):
            for i, j, _fi, _fj in classified.get(key, []):
                ends.add(int(i))
                ends.add(int(j))
        return sorted(ends)

    def _f_pred_on_edge(self, f_a: float, f_b: float, alpha: float) -> float:
        """Linear f prediction along edge (for near-T gating of fill points)."""
        return float(f_a + alpha * (f_b - f_a))

    def _alpha_moves_toward_t(self, f_from: float, f_to: float, alpha: float, t: float) -> bool:
        """True if stepping from f_from toward f_to by alpha reduces |f−T| (converge to level)."""
        f_pred = self._f_pred_on_edge(f_from, f_to, alpha)
        return abs(f_pred - t) < abs(f_from - t) - 1e-15

    def _expand_from_core_along_edge(
        self,
        edges_b: Sequence[tuple[int, int, float, float]],
        core_set: set[int],
        *,
        sep_index: _SepIndex,
        r_sep: float,
        budget: int,
        new_us: list[np.ndarray],
    ) -> int:
        """Case B: one end is core — place points *toward T*, not outward off-band.

        Only secant alphas that reduce |f−T| relative to the core (converge to
        the level set). Fixed large steps toward the far endpoint are forbidden
        — those expand both sides of the ridge instead of densifying the middle.
        """
        t = self._target_value
        # Predicted |f−T| must stay inside a soft outer shell around T.
        f_cap = max(float(self._outer_half_width), 3.0 * float(self._core_half_width))
        n_before = len(new_us)
        for i, j, fi, fj in edges_b:
            if len(new_us) >= budget:
                break
            i, j = int(i), int(j)
            if i in core_set and j not in core_set:
                i_core, i_other, f_core, f_other = i, j, fi, fj
            elif j in core_set and i not in core_set:
                i_core, i_other, f_core, f_other = j, i, fj, fi
            else:
                continue
            u_c = np.asarray(self._points[i_core].u, dtype=np.float64)
            u_o = np.asarray(self._points[i_other].u, dtype=np.float64)
            edge_len = float(np.linalg.norm(u_o - u_c))
            if edge_len < 1e-15:
                continue
            # Only toward-T secants (no fixed outward 0.3/0.45 steps).
            alphas = list(self._secant_alphas(f_core, f_other, t))
            seen: set[float] = set()
            edge_hits = 0
            for alpha in alphas:
                if len(new_us) >= budget or edge_hits >= 2:
                    break
                alpha = float(min(0.85, max(0.08, alpha)))
                key = round(alpha, 4)
                if key in seen:
                    continue
                seen.add(key)
                if not self._alpha_moves_toward_t(f_core, f_other, alpha, t):
                    continue
                f_pred = self._f_pred_on_edge(f_core, f_other, alpha)
                if abs(f_pred - t) > f_cap + 1e-15:
                    continue
                cand = clip_unit_cube(u_c + alpha * (u_o - u_c))
                if float(np.linalg.norm(cand - u_c)) < 0.35 * max(r_sep, 1e-4):
                    continue
                if not sep_index.far_enough(cand):
                    continue
                physical = map_u_to_physical(cand, self.vars)
                if self._selectionexp and not evaluate_selection(
                    self._selectionexp,
                    physical,
                    context=self._expression_context,
                ):
                    continue
                if self._unique_point_count() + len(new_us) >= self._max_points:
                    break
                sep_index.add(cand)
                new_us.append(cand)
                edge_hits += 1
        return len(new_us) - n_before

    def _fill_needed(
        self,
        *,
        cores: Sequence[int],
        n_local_brackets: int,
        coverage_ok: bool,
        densify_new: int | None = None,
        densify_attempted: bool = False,
        n_bridge_A: int = 0,
        n_bridge_B: int = 0,
        n_bridge_C: int = 0,
        n_endpoint_edges: int = 0,
        n_active_endpoints: int = 0,
        n_active_virtual_cores: int = 0,
        scale_quiet: bool = False,
        force_scale_close: bool = False,
    ) -> bool:
        """Whether same-r_g fill (root-corr / expand / densify / bridge) is required.

        A/B/C edges request fill only while this scale is making progress.
        Straddles normally persist along a resolved level set, so their mere
        existence must not prevent r_g from shrinking forever.
        """
        if force_scale_close:
            return False
        # Endpoint fronts are recomputed after every feedback barrier. As long
        # as a straddle root remains spatially uncovered at this r_g, advance
        # that front before considering the scale quiet.
        if int(n_endpoint_edges) > 0:
            return True
        if int(n_active_endpoints) > 0:
            return True
        if int(n_active_virtual_cores) > 0:
            return True
        # Missing/disconnected coverage is structural, so a few quiet random
        # passes must not close the scale before the remaining islands have
        # had a chance to connect.
        if not cores:
            return True
        if not coverage_ok:
            return True
        if scale_quiet:
            return False
        if int(n_bridge_A) > 0 or int(n_bridge_B) > 0:
            return True
        if int(n_bridge_C) > 0:
            return True
        # Fallback when classification not provided: any local nearest bracket.
        if int(n_local_brackets) > 0 and (
            n_bridge_A + n_bridge_B + n_bridge_C == 0
        ):
            # unclassified path — treat as need fill
            return True
        if densify_attempted and densify_new is not None:
            trickle = max(3, min(12, len(cores) // 2))
            if densify_new > trickle:
                return True
            return False
        return False

    def _edge_root_anchor(
        self, edge: tuple[int, int, float, float]
    ) -> np.ndarray:
        """Best secant estimate of the level-set root on one graph edge."""
        i, j, fi, fj = edge
        alphas = self._secant_alphas(fi, fj, self._target_value)
        alpha = min(
            alphas,
            key=lambda a: abs(
                self._f_pred_on_edge(fi, fj, float(a)) - self._target_value
            ),
        )
        alpha = float(min(0.85, max(0.08, alpha)))
        ui = np.asarray(self._points[i].u, dtype=np.float64)
        uj = np.asarray(self._points[j].u, dtype=np.float64)
        return clip_unit_cube(ui + alpha * (uj - ui))

    def _local_bracket_edges(
        self,
        edges: Sequence[tuple[int, int, float, float]],
        *,
        r_g: float | None = None,
    ) -> list[tuple[int, int, float, float]]:
        """Adjacent target-straddling edges used to define working cores."""
        limit = 2.5 * float(self._radius if r_g is None else r_g)
        local: list[tuple[int, int, float, float]] = []
        for edge in edges:
            i, j, _fi, _fj = edge
            ui = np.asarray(self._points[int(i)].u, dtype=np.float64)
            uj = np.asarray(self._points[int(j)].u, dtype=np.float64)
            if float(np.linalg.norm(uj - ui)) <= limit + 1e-15:
                local.append(edge)
        return local

    def _bracket_root_anchors(
        self,
        edges: Sequence[tuple[int, int, float, float]],
    ) -> list[np.ndarray]:
        """Secant roots of local brackets, spatially thinned at r_g."""
        if not edges:
            return []
        index = _SepIndex(max(float(self._radius), 1e-15), self._dim)
        anchors: list[np.ndarray] = []
        ranked = sorted(
            edges,
            key=lambda edge: (
                max(
                    abs(float(edge[2]) - self._target_value),
                    abs(float(edge[3]) - self._target_value),
                ),
                int(edge[0]),
                int(edge[1]),
            ),
        )
        for edge in ranked:
            anchor = self._edge_root_anchor(edge)
            if index.far_enough(anchor):
                index.add(anchor)
                anchors.append(anchor)
        return anchors

    def _virtual_core_is_covered(
        self, anchor: np.ndarray, cores: Sequence[int]
    ) -> bool:
        """Whether an evaluated core already covers a bracket root."""
        if not cores:
            return False
        core_pos = np.asarray(
            [self._points[int(index)].u for index in cores], dtype=np.float64
        )
        return bool(
            np.min(np.linalg.norm(core_pos - anchor[None, :], axis=1))
            <= 2.0 * float(self._radius) + 1e-15
        )

    def _virtual_core_id(self, anchor: np.ndarray) -> str:
        # Quantisation only supplies a deterministic identity for newly seen
        # roots. Existing roots are matched geometrically below, because the
        # secant estimate naturally moves after every feedback batch.
        scale = max(0.5 * float(self._radius), 1e-12)
        cell = tuple(int(np.floor(float(value) / scale)) for value in anchor)
        return "bracket:%d:%s" % (
            int(self._generation),
            ",".join(str(value) for value in cell),
        )

    def _refresh_virtual_cores(
        self, anchors: Sequence[np.ndarray], cores: Sequence[int]
    ) -> int:
        """Persist uncovered bracket roots across same-r_g feedback passes."""
        core_signature = tuple(
            sorted(
                str(self._points[int(index)].uuid or f"index:{int(index)}")
                for index in cores
            )
        )
        if (
            self._virtual_core_core_signature is not None
            and core_signature != self._virtual_core_core_signature
        ):
            self._retired_virtual_core_ids.clear()
        self._virtual_core_core_signature = core_signature
        remaining = [np.asarray(anchor, dtype=np.float64) for anchor in anchors]
        refreshed: dict[str, dict[str, Any]] = {}
        point_by_u = {self._u_key(point.u): point for point in self._points}

        for root_id, old in sorted(self._active_virtual_cores.items()):
            state = dict(old)
            old_anchor = np.asarray(state.get("anchor_u") or [], dtype=np.float64)
            if old_anchor.size != self._dim:
                continue
            if self._virtual_core_is_covered(old_anchor, cores):
                self._retired_virtual_core_ids.add(root_id)
                continue
            improved = False
            if bool(state.get("pending_attempt")):
                best_probe: LevelSetPoint | None = None
                best_dev = float("inf")
                for probe_u in state.get("pending_probe_us") or []:
                    point = point_by_u.get(self._u_key(probe_u))
                    if point is None:
                        continue
                    value = _finite_f(point.f)
                    if value is None:
                        continue
                    deviation = abs(float(value) - self._target_value)
                    if deviation < best_dev:
                        best_dev = deviation
                        best_probe = point
                previous_raw = state.get("best_probe_dev")
                previous = (
                    float(previous_raw)
                    if previous_raw is not None
                    else float("inf")
                )
                if best_probe is not None and best_dev < previous - 1e-15:
                    # Keep the secant root as the identity/coverage anchor, but
                    # walk the sampling centre through progressively better
                    # evaluated points. This lets a root escape a locally
                    # saturated r_g shell without promoting an off-band point
                    # to an actual core.
                    state["search_u"] = list(map(float, best_probe.u))
                    state["best_probe_dev"] = float(best_dev)
                    state["misses"] = 0
                    improved = True
                if not improved:
                    state["misses"] = int(state.get("misses", 0)) + 1
            state["pending_attempt"] = False
            state["pending_probe_us"] = []
            if (
                int(state.get("generation", self._generation))
                != int(self._generation)
                or int(state.get("misses", 0)) >= self._endpoint_stall_limit
            ):
                self._retired_virtual_core_ids.add(root_id)
                continue

            if remaining:
                distances = np.asarray(
                    [np.linalg.norm(candidate - old_anchor) for candidate in remaining]
                )
                nearest = int(np.argmin(distances))
                if float(distances[nearest]) <= float(self._radius) + 1e-15:
                    moved = remaining.pop(nearest)
                    if not self._virtual_core_is_covered(moved, cores):
                        state["anchor_u"] = moved.tolist()
            refreshed[root_id] = state

        for anchor in remaining:
            if self._virtual_core_is_covered(anchor, cores):
                continue
            root_id = self._virtual_core_id(anchor)
            if root_id in refreshed or root_id in self._retired_virtual_core_ids:
                continue
            refreshed[root_id] = {
                "id": root_id,
                "anchor_u": anchor.tolist(),
                "search_u": anchor.tolist(),
                "generation": int(self._generation),
                "attempts": 0,
                "misses": 0,
                "pending_attempt": False,
                "pending_probe_us": [],
                "best_probe_dev": None,
            }
        self._active_virtual_cores = refreshed
        return len(refreshed)

    def _virtual_core_direct_proposals(
        self,
        *,
        sep_index: _SepIndex,
        new_us: list[np.ndarray],
        budget: int,
    ) -> int:
        """Submit the secant root itself before any random exploration."""
        n_before = len(new_us)
        ordered = sorted(
            self._active_virtual_cores.items(),
            key=lambda item: (
                int(item[1].get("misses", 0)),
                (
                    float(item[1]["best_probe_dev"])
                    if item[1].get("best_probe_dev") is not None
                    else float("inf")
                ),
                item[0],
            ),
        )
        for _root_id, state in ordered:
            if len(new_us) >= budget:
                break
            anchor = np.asarray(state["anchor_u"], dtype=np.float64)
            if not self._try_accept_candidate(
                anchor,
                sep_index=sep_index,
                new_us=new_us,
                budget=budget,
            ):
                continue
            state["attempts"] = int(state.get("attempts", 0)) + 1
            state["pending_attempt"] = True
            state["pending_probe_us"] = [anchor.tolist()]
        return len(new_us) - n_before

    def _virtual_core_proposals(
        self,
        *,
        sep_index: _SepIndex,
        new_us: list[np.ndarray],
        budget: int,
        rng: np.random.Generator,
    ) -> int:
        """Give each persistent uncovered root at most two seed proposals."""
        n_before = len(new_us)
        per_anchor = 2
        d = max(1, int(self._dim))
        inner = float(self._radius)
        outer = 2.5 * inner
        for root_id in sorted(self._active_virtual_cores):
            state = self._active_virtual_cores[root_id]
            # A direct secant proposal is already the best use of this root in
            # the current feedback batch. Random fallback waits for its result.
            if bool(state.get("pending_attempt")):
                continue
            anchor = np.asarray(
                state.get("search_u") or state["anchor_u"], dtype=np.float64
            )
            # A completed trial counts even when strict blue noise rejects all
            # candidates. Otherwise a saturated root remains active forever.
            state["attempts"] = int(state.get("attempts", 0)) + 1
            state["pending_attempt"] = True
            state["pending_probe_us"] = []
            accepted = 0
            for _trial in range(max(24, 10 * per_anchor)):
                if accepted >= per_anchor or len(new_us) >= budget:
                    break
                direction = rng.normal(size=self._dim)
                norm = float(np.linalg.norm(direction))
                if norm <= 1e-15:
                    continue
                direction /= norm
                radial = (
                    inner**d
                    + float(rng.random()) * (outer**d - inner**d)
                ) ** (1.0 / d)
                if self._try_accept_candidate(
                    anchor + radial * direction,
                    sep_index=sep_index,
                    new_us=new_us,
                    budget=budget,
                ):
                    accepted += 1
                    state["pending_probe_us"].append(
                        clip_unit_cube(anchor + radial * direction).tolist()
                    )
            if len(new_us) >= budget:
                break
        return len(new_us) - n_before

    def _actionable_endpoint_edges(
        self,
        edges: Sequence[tuple[int, int, float, float]],
        cores: Sequence[int],
        *,
        r_sep: float,
    ) -> list[tuple[int, int, float, float]]:
        """Uncovered straddle roots that must advance before r_g shrinks."""
        if not edges:
            return []
        history_index = self._submitted_sep_index(r_sep)
        core_pos = (
            np.asarray([self._points[int(i)].u for i in cores], dtype=np.float64)
            if cores
            else np.empty((0, self._dim), dtype=np.float64)
        )
        min_novelty = max(0.75 * float(r_sep), 1e-4)
        active: list[tuple[int, int, float, float]] = []
        for edge in edges:
            i, j, fi, fj = edge
            ui = np.asarray(self._points[i].u, dtype=np.float64)
            uj = np.asarray(self._points[j].u, dtype=np.float64)
            for alpha in self._secant_alphas(fi, fj, self._target_value):
                alpha = float(min(0.85, max(0.08, alpha)))
                anchor = clip_unit_cube(ui + alpha * (uj - ui))
                if self._u_key(anchor) in self._submitted_u_keys:
                    continue
                if not history_index.far_enough(anchor):
                    continue
                if core_pos.size:
                    novelty = float(
                        np.min(np.linalg.norm(core_pos - anchor[None, :], axis=1))
                    )
                    if novelty <= min_novelty + 1e-15:
                        continue
                active.append(edge)
                break
        return active

    def _core_components(
        self,
        cores: Sequence[int],
        *,
        link_distance: float | None = None,
    ) -> list[list[int]]:
        """Connected core components at the current spatial resolution."""
        ids = sorted({int(i) for i in cores})
        if not ids:
            return []
        if len(ids) == 1:
            return [ids]
        radius = (
            float(self._core_spacing_factor) * float(self._radius)
            if link_distance is None
            else float(link_distance)
        )
        pts = np.asarray([self._points[i].u for i in ids], dtype=np.float64)
        parent = list(range(len(ids)))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def unite(a: int, b: int) -> None:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra

        try:
            from scipy.spatial import cKDTree

            pairs = cKDTree(pts).query_pairs(
                r=max(radius, 1e-15), output_type="ndarray"
            )
            for a, b in pairs:
                unite(int(a), int(b))
        except Exception:
            for a in range(len(ids)):
                for b in range(a + 1, len(ids)):
                    if float(np.linalg.norm(pts[a] - pts[b])) <= radius + 1e-15:
                        unite(a, b)

        groups: dict[int, list[int]] = {}
        for local, point_i in enumerate(ids):
            groups.setdefault(find(local), []).append(point_i)
        return sorted(
            (sorted(group) for group in groups.values()),
            key=lambda group: (self._u_key(self._points[group[0]].u), len(group)),
        )

    def _endpoint_id(
        self, *, kind: str, anchor: int, target: int | None = None, slot: str = ""
    ) -> str:
        anchor_id = str(self._points[int(anchor)].uuid or f"point-{int(anchor)}")
        if target is None:
            return f"{kind}:{anchor_id}:{slot}"
        target_id = str(self._points[int(target)].uuid or f"point-{int(target)}")
        return f"{kind}:{anchor_id}:{target_id}"

    def _geometric_endpoint_candidates(
        self, cores: Sequence[int]
    ) -> list[dict[str, Any]]:
        """Detect open contour ends and gap-facing component fronts.

        Local open ends have a strongly one-sided neighborhood.  For separate
        components, the nearest cross-component pair is also exposed as two
        directed fronts. Probe shells may widen after misses, but the core
        anchor advances only after feedback confirms a new on-level core.
        """
        components = self._core_components(cores)
        if not components:
            return []
        link = float(self._core_spacing_factor) * float(self._radius)
        candidates: dict[str, dict[str, Any]] = {}

        def add(kind: str, anchor: int, direction: np.ndarray, *, target: int | None = None, slot: str = "") -> None:
            norm = float(np.linalg.norm(direction))
            if norm <= 1e-15:
                return
            endpoint_id = self._endpoint_id(
                kind=kind, anchor=anchor, target=target, slot=slot
            )
            candidates[endpoint_id] = {
                "id": endpoint_id,
                "kind": kind,
                "anchor_uuid": str(
                    self._points[int(anchor)].uuid or f"point-{int(anchor)}"
                ),
                "anchor_u": list(np.asarray(self._points[int(anchor)].u, dtype=float)),
                "direction": list(np.asarray(direction, dtype=float) / norm),
                "target_uuid": (
                    None
                    if target is None
                    else str(self._points[int(target)].uuid or f"point-{int(target)}")
                ),
                "target_distance": (
                    None
                    if target is None
                    else float(
                        np.linalg.norm(
                            np.asarray(self._points[int(target)].u, dtype=float)
                            - np.asarray(self._points[int(anchor)].u, dtype=float)
                        )
                    )
                ),
                "attempts": 0,
                "misses": 0,
                "pending_attempt": False,
                "max_probe_distance": 0.0,
                "generation": int(self._generation),
            }

        # Genuine open ends: a normalized vector sum near one means all local
        # neighbors lie on one side. Closed loops have balanced neighborhoods.
        for component in components:
            if len(component) == 1:
                if len(components) == 1:
                    axis = int(self._fill_pass) % max(self._dim, 1)
                    unit = np.zeros(self._dim, dtype=np.float64)
                    unit[axis] = 1.0
                    add("terminal", component[0], unit, slot="plus")
                    add("terminal", component[0], -unit, slot="minus")
                continue
            pts = np.asarray([self._points[i].u for i in component], dtype=np.float64)
            scored: list[tuple[float, int, np.ndarray]] = []
            for local, point_i in enumerate(component):
                delta = pts - pts[local]
                dist = np.linalg.norm(delta, axis=1)
                mask = (dist > 1e-15) & (dist <= link + 1e-15)
                if not np.any(mask):
                    continue
                unit = delta[mask] / dist[mask, None]
                mean = np.mean(unit, axis=0)
                score = float(np.linalg.norm(mean))
                if score >= 0.45:
                    scored.append((score, point_i, -mean))
            scored.sort(key=lambda item: (-item[0], item[1]))
            chosen: list[int] = []
            for _score, point_i, direction in scored:
                u = np.asarray(self._points[point_i].u, dtype=np.float64)
                if any(
                    float(np.linalg.norm(u - np.asarray(self._points[j].u)))
                    < max(link, 1e-15)
                    for j in chosen
                ):
                    continue
                add("terminal", point_i, direction, slot="open")
                chosen.append(point_i)
                if len(chosen) >= 2:
                    break

        # Each component grows toward its nearest other component.  Symmetric
        # duplicates collapse by endpoint id; recomputation after feedback
        # moves the front forward rather than jumping across the entire gap.
        if len(components) > 1:
            component_pts = [
                np.asarray([self._points[i].u for i in component], dtype=np.float64)
                for component in components
            ]
            nearest: list[tuple[float, int, int] | None] = [None] * len(components)
            try:
                from scipy.spatial import cKDTree

                trees = [cKDTree(pts) for pts in component_pts]
                n_components = len(components)
                if n_components <= 64:
                    candidate_pairs = list(
                        itertools.combinations(range(n_components), 2)
                    )
                else:
                    # Avoid an O(C²) explosion for a very fragmented bootstrap
                    # cloud. Nearby component centroids provide a bounded set
                    # of exact point-to-component distance checks.
                    centroids = np.asarray(
                        [np.mean(pts, axis=0) for pts in component_pts]
                    )
                    _dist, neighbors = cKDTree(centroids).query(
                        centroids, k=min(9, n_components)
                    )
                    pair_set: set[tuple[int, int]] = set()
                    for left, row in enumerate(np.atleast_2d(neighbors)):
                        for right_raw in np.atleast_1d(row):
                            right = int(right_raw)
                            if right == left:
                                continue
                            pair_set.add(
                                (left, right) if left < right else (right, left)
                            )
                    candidate_pairs = sorted(pair_set)
                for left, right in candidate_pairs:
                    distances, local_right = trees[right].query(
                        component_pts[left], k=1
                    )
                    local_left = int(np.argmin(distances))
                    target_local = int(np.asarray(local_right)[local_left])
                    pair_lr = (
                        float(np.asarray(distances)[local_left]),
                        components[left][local_left],
                        components[right][target_local],
                    )
                    pair_rl = (pair_lr[0], pair_lr[2], pair_lr[1])
                    if nearest[left] is None or pair_lr < nearest[left]:
                        nearest[left] = pair_lr
                    if nearest[right] is None or pair_rl < nearest[right]:
                        nearest[right] = pair_rl
            except Exception:
                for left in range(len(components)):
                    for right in range(left + 1, len(components)):
                        for li, anchor in enumerate(components[left]):
                            distances = np.linalg.norm(
                                component_pts[right] - component_pts[left][li], axis=1
                            )
                            rj = int(np.argmin(distances))
                            pair_lr = (
                                float(distances[rj]), anchor, components[right][rj]
                            )
                            pair_rl = (pair_lr[0], pair_lr[2], pair_lr[1])
                            if nearest[left] is None or pair_lr < nearest[left]:
                                nearest[left] = pair_lr
                            if nearest[right] is None or pair_rl < nearest[right]:
                                nearest[right] = pair_rl
            for best_pair in nearest:
                if best_pair is None:
                    continue
                _distance, anchor, target = best_pair
                direction = np.asarray(self._points[target].u) - np.asarray(
                    self._points[anchor].u
                )
                add("gap", anchor, direction, target=target)
        return list(candidates.values())

    def _endpoint_boundary_distance(self, endpoint: Mapping[str, Any]) -> float:
        u = np.asarray(endpoint.get("anchor_u") or [], dtype=np.float64)
        direction = np.asarray(endpoint.get("direction") or [], dtype=np.float64)
        if u.size != self._dim or direction.size != self._dim:
            return 0.0
        distances: list[float] = []
        for value, step in zip(u, direction):
            if step > 1e-15:
                distances.append((1.0 - float(value)) / float(step))
            elif step < -1e-15:
                distances.append((0.0 - float(value)) / float(step))
        return min((x for x in distances if x >= 0.0), default=float("inf"))

    def _endpoint_has_advanced(
        self, endpoint: Mapping[str, Any], cores: Sequence[int]
    ) -> bool:
        anchor = np.asarray(endpoint.get("anchor_u") or [], dtype=np.float64)
        direction = np.asarray(endpoint.get("direction") or [], dtype=np.float64)
        if anchor.size != self._dim or direction.size != self._dim:
            return True
        for point_i in cores:
            delta = np.asarray(self._points[int(point_i)].u) - anchor
            forward = float(np.dot(delta, direction))
            if forward < 0.50 * float(self._radius):
                continue
            # A far component is the destination of a gap front, not evidence
            # that this anchor already advanced. Only a newly grown local step
            # may supersede the endpoint.
            max_probe = max(
                2.50 * float(self._radius),
                float(endpoint.get("max_probe_distance", 0.0) or 0.0)
                + 0.75 * float(self._radius),
            )
            target_distance = endpoint.get("target_distance")
            if target_distance is not None:
                max_probe = min(max_probe, 0.75 * float(target_distance))
            if forward > max_probe + 1e-15:
                continue
            lateral = float(np.linalg.norm(delta - forward * direction))
            if lateral <= 1.25 * float(self._radius) + 1e-15:
                return True
        return False

    def _refresh_active_endpoints(self, cores: Sequence[int]) -> int:
        """Refresh persistent fronts after a feedback barrier.

        A pending attempt is counted as a miss only if no new core advanced
        beyond its anchor. Unmatched fronts survive transient geometric graph
        changes and are retired only by explicit, conservative conditions.
        """
        core_signature = tuple(
            sorted(
                str(self._points[int(index)].uuid or f"index:{int(index)}")
                for index in cores
            )
        )
        # A retired front was exhausted only for the old core geometry. Once
        # feedback adds an actual core, rebuild all geometric fronts and allow
        # an old gap/terminal id to become active again from its new anchor.
        if (
            self._endpoint_core_signature is not None
            and core_signature != self._endpoint_core_signature
        ):
            self._retired_endpoint_ids.clear()
        self._endpoint_core_signature = core_signature
        detected = {
            str(item["id"]): item
            for item in self._geometric_endpoint_candidates(cores)
        }
        refreshed: dict[str, dict[str, Any]] = {}
        for endpoint_id, old in self._active_endpoints.items():
            state = dict(old)
            advanced = self._endpoint_has_advanced(state, cores)
            if bool(state.get("pending_attempt")) and not advanced:
                state["misses"] = int(state.get("misses", 0)) + 1
            state["pending_attempt"] = False
            retire = (
                advanced
                or self._endpoint_boundary_distance(state)
                < 0.75 * float(self._radius)
                or int(state.get("misses", 0)) >= self._endpoint_stall_limit
                or int(state.get("generation", self._generation))
                != int(self._generation)
            )
            if retire:
                self._retired_endpoint_ids.add(endpoint_id)
                continue
            if endpoint_id in detected:
                new_state = detected.pop(endpoint_id)
                new_state["attempts"] = int(state.get("attempts", 0))
                new_state["misses"] = int(state.get("misses", 0))
                new_state["max_probe_distance"] = float(
                    state.get("max_probe_distance", 0.0) or 0.0
                )
                state = new_state
            refreshed[endpoint_id] = state
        for endpoint_id, state in detected.items():
            if endpoint_id not in self._retired_endpoint_ids:
                refreshed[endpoint_id] = state
        self._active_endpoints = refreshed
        return len(refreshed)

    def _active_endpoint_proposals(
        self,
        *,
        sep_index: _SepIndex,
        new_us: list[np.ndarray],
        budget: int,
        rng: np.random.Generator,
    ) -> int:
        """Search outward from every active front with strict blue noise.

        Draw uniformly by volume from the hollow d-ball with radii
        ``[r_g, 2.5*r_g]``. Feedback moves the core anchor only when a genuine
        core is found; the endpoint set is rebuilt after the feedback barrier.
        """
        n_before = len(new_us)
        for endpoint_id in sorted(self._active_endpoints):
            if len(new_us) >= budget:
                break
            endpoint = self._active_endpoints[endpoint_id]
            anchor = np.asarray(endpoint["anchor_u"], dtype=np.float64)
            endpoint["model_guided"] = False
            endpoint["attempts"] = int(endpoint.get("attempts", 0)) + 1
            endpoint["pending_attempt"] = True
            n_front_accepted = 0
            front_cap = int(self._endpoint_omni_probes)
            d = max(1, int(self._dim))
            inner = float(self._radius)
            outer = 2.5 * inner
            max_trials = max(32, 12 * front_cap)
            for _trial in range(max_trials):
                if n_front_accepted >= front_cap or len(new_us) >= budget:
                    break
                direction = rng.normal(size=self._dim)
                norm = float(np.linalg.norm(direction))
                if norm <= 1e-15:
                    continue
                direction /= norm
                radial = (
                    inner**d
                    + float(rng.random()) * (outer**d - inner**d)
                ) ** (1.0 / d)
                cand = anchor + radial * direction
                if self._try_accept_candidate(
                    cand,
                    sep_index=sep_index,
                    new_us=new_us,
                    budget=budget,
                ):
                    endpoint["max_probe_distance"] = max(
                        float(endpoint.get("max_probe_distance", 0.0) or 0.0),
                        radial,
                    )
                    n_front_accepted += 1
        return len(new_us) - n_before

    def _reset_endpoint_state(self) -> None:
        self._active_endpoints = {}
        self._retired_endpoint_ids = set()
        self._endpoint_core_signature = None
        self._active_virtual_cores = {}
        self._retired_virtual_core_ids = set()
        self._virtual_core_core_signature = None

    @staticmethod
    def _u_key(u: Sequence[float] | np.ndarray) -> tuple[float, ...]:
        """Stable exact-ish key for permanent submitted-coordinate deduplication."""
        arr = np.asarray(u, dtype=np.float64).reshape(-1)
        return tuple(round(float(x), 14) for x in arr)

    def _unique_point_count(self) -> int:
        """Number of unique coordinates submitted over the whole run."""
        return len(self._submitted_u_keys)

    def _filter_unsubmitted(
        self, candidates: Sequence[np.ndarray]
    ) -> list[np.ndarray]:
        """Drop submitted coordinates, including points pruned from the graph."""
        unseen: list[np.ndarray] = []
        local: set[tuple[float, ...]] = set()
        for candidate in candidates:
            u = clip_unit_cube(np.asarray(candidate, dtype=np.float64))
            key = self._u_key(u)
            if key in self._submitted_u_keys or key in local:
                continue
            local.add(key)
            unseen.append(u)
        return unseen

    def _submitted_sep_index(self, min_distance: float) -> _SepIndex:
        """Blue-noise index over the complete submitted-coordinate history.

        The live control cloud is prunable, while ``_submitted_u_keys`` is
        permanent.  Seed from both so legacy/in-memory test states without a
        populated history still enforce separation from their current sites.
        """
        index = _SepIndex(max(float(min_distance), 1e-15), self._dim)
        keys = set(self._submitted_u_keys)
        keys.update(self._u_key(point.u) for point in self._points)
        for key in sorted(keys):
            index.add(np.asarray(key, dtype=np.float64))
        return index

    def _filter_blue_noise(
        self,
        candidates: Sequence[np.ndarray],
        *,
        min_distance: float | None = None,
    ) -> list[np.ndarray]:
        """Keep candidates at least ``min_distance`` from history and peers."""
        distance = (
            float(self._radius)
            if min_distance is None
            else float(min_distance)
        )
        index = self._submitted_sep_index(distance)
        accepted: list[np.ndarray] = []
        for candidate in candidates:
            u = clip_unit_cube(np.asarray(candidate, dtype=np.float64))
            if self._u_key(u) in self._submitted_u_keys:
                continue
            if not index.far_enough(u):
                continue
            index.add(u)
            accepted.append(u)
        return accepted

    def _fill_rng(self, generation: int) -> np.random.Generator:
        """Deterministic RNG stream unique to a (generation, fill_pass) pair."""
        seed = np.random.SeedSequence(
            [
                int(self._seed),
                max(0, int(generation)),
                max(0, int(self._fill_pass)),
                0xAB12,
            ]
        )
        return np.random.default_rng(seed)

    def _compute_live_band(
        self,
    ) -> tuple[int | None, float | None, float | None, list[int]]:
        """Returns (best, t_min, t_max, core_indices).

        ``t_min/t_max`` = 2*r_g scale band around best (exit/next-gen gate).
        ``core_indices`` = absolute |f−T|≤w_core densify centers.
        """
        cores, _frontier, _outer, _best_cls = self._classify_bands()
        best, t_min, t_max = self._radius_t_band()
        return best, t_min, t_max, list(cores)

    def _straddle_brackets(
        self, *, pool: Sequence[int] | None = None
    ) -> list[tuple[int, int, float, float]]:
        """Neighbor edges that straddle the target T (sign change in f−T).

        Used for secant root-correction. ``pool`` defaults to all finite points
        (needed when no cores exist yet — discovery from gen-0).
        """
        t = self._target_value
        adj = self._adjacency()
        if pool is None:
            pool_set = {
                i
                for i, p in enumerate(self._points)
                if _finite_f(p.f) is not None
            }
        else:
            pool_set = {int(i) for i in pool}
        brackets: list[tuple[int, int, float, float]] = []
        seen: set[tuple[int, int]] = set()
        for i in sorted(pool_set):
            fi = _finite_f(self._points[i].f)
            if fi is None:
                continue
            di = fi - t
            for j in adj.get(i, []):
                if j not in pool_set or j <= i:
                    continue
                fj = _finite_f(self._points[j].f)
                if fj is None:
                    continue
                dj = fj - t
                if di * dj < 0.0:  # strict straddle
                    key = (i, j)
                    if key in seen:
                        continue
                    seen.add(key)
                    brackets.append((i, j, float(fi), float(fj)))
        return brackets

    def _nearest_straddle_brackets(
        self,
        *,
        pool: Sequence[int] | None = None,
        max_distance: float | None = None,
    ) -> list[tuple[int, int, float, float]]:
        """Nearest u-space neighbour on the opposite side of target.

        Delaunay/Voronoi adjacency is useful for the general control graph,
        but pruning and the capped bracket ranking can hide an otherwise
        obvious local sign change. Query both directions (below→above and
        above→below), then retain only trustworthy, genuinely local pairs.
        This path deliberately does not use the global bracket cap.
        """
        from scipy.spatial import cKDTree

        indices = (
            sorted({int(index) for index in pool})
            if pool is not None
            else list(range(len(self._points)))
        )
        t = float(self._target_value)
        below: list[int] = []
        above: list[int] = []
        values: dict[int, float] = {}
        for index in indices:
            if not 0 <= index < len(self._points):
                continue
            value = _finite_f(self._points[index].f)
            if value is None:
                continue
            values[index] = float(value)
            if value < t:
                below.append(index)
            elif value > t:
                above.append(index)
        if not below or not above:
            return []

        distance_limit = (
            2.5 * float(self._radius)
            if max_distance is None
            else float(max_distance)
        )
        near = max(3.0 * float(self._outer_half_width), 0.05)
        far = max(15.0 * float(self._outer_half_width), 0.25)
        pairs: set[tuple[int, int]] = set()
        for sources, targets in ((below, above), (above, below)):
            target_u = np.asarray(
                [self._points[index].u for index in targets], dtype=np.float64
            )
            tree = cKDTree(target_u)
            source_u = np.asarray(
                [self._points[index].u for index in sources], dtype=np.float64
            )
            distances, neighbours = tree.query(source_u, k=1)
            for source, local_target, distance in zip(
                sources,
                np.atleast_1d(neighbours),
                np.atleast_1d(distances),
            ):
                if float(distance) > distance_limit + 1e-15:
                    continue
                target = targets[int(local_target)]
                fi, fj = values[source], values[target]
                if min(abs(fi - t), abs(fj - t)) > near + 1e-15:
                    continue
                if max(abs(fi - t), abs(fj - t)) > far + 1e-15:
                    continue
                pairs.add(
                    (source, target) if source < target else (target, source)
                )
        return [
            (i, j, values[i], values[j]) for i, j in sorted(pairs)
        ]

    @staticmethod
    def _merge_bracket_edges(
        *groups: Sequence[tuple[int, int, float, float]],
    ) -> list[tuple[int, int, float, float]]:
        merged: dict[tuple[int, int], tuple[int, int, float, float]] = {}
        for group in groups:
            for edge in group:
                key = (
                    min(int(edge[0]), int(edge[1])),
                    max(int(edge[0]), int(edge[1])),
                )
                merged[key] = edge
        return [merged[key] for key in sorted(merged)]

    @staticmethod
    def _secant_alphas(fi: float, fj: float, t: float) -> list[float]:
        """Few alphas along edge i→j aiming at f≈T (log + linear + mid).

        Avoid flooding every Delaunay chord with 7+ samples (hairball mode).
        """
        alphas: list[float] = []
        if fi > 0.0 and fj > 0.0 and t > 0.0 and fi != fj:
            try:
                alphas.append(
                    (math.log(t) - math.log(fi)) / (math.log(fj) - math.log(fi))
                )
            except (ValueError, ZeroDivisionError):
                pass
        denom = fj - fi
        if abs(denom) > 1e-30:
            alphas.append(float((t - fi) / denom))
        alphas.append(0.5)
        out: list[float] = []
        seen: set[float] = set()
        for a in alphas:
            a = float(min(0.92, max(0.08, a)))
            key = round(a, 5)
            if key in seen:
                continue
            seen.add(key)
            out.append(a)
        return out

    def _root_correct_proposals(
        self,
        brackets: Sequence[tuple[int, int, float, float]],
        *,
        sep_index: _SepIndex | None,
        r_sep: float,
        budget: int,
        r_window: float | None = None,
    ) -> list[np.ndarray]:
        """Root-correction on local nearest opposite-side pairs."""
        t = self._target_value
        new_us: list[np.ndarray] = []
        if budget <= 0 or not brackets:
            return new_us

        ranked = sorted(
            brackets,
            key=lambda edge: (
                max(
                    abs(float(edge[2]) - self._target_value),
                    abs(float(edge[3]) - self._target_value),
                ),
                int(edge[0]),
                int(edge[1]),
            ),
        )
        n_raw = len(brackets)

        # Root correction obeys the same current-r_g blue-noise exclusion as
        # ordinary Bridson proposals.  The caller normally supplies an index
        # containing the complete submission history.
        r_floor = max(float(r_sep), 1e-15)
        if sep_index is None:
            sep_index = self._submitted_sep_index(r_floor)

        f_cap = max(float(self._outer_half_width), 3.0 * float(self._core_half_width))
        n_try = 0
        n_accept = 0
        for i, j, fi, fj in ranked:
            if len(new_us) >= budget:
                break
            ui = np.asarray(self._points[i].u, dtype=np.float64)
            uj = np.asarray(self._points[j].u, dtype=np.float64)
            edge_len = float(np.linalg.norm(uj - ui))
            if edge_len < 1e-15:
                continue
            r_edge = max(1e-4, min(r_floor, 0.2 * edge_len))
            edge_hits = 0
            for alpha in self._secant_alphas(fi, fj, t):
                if len(new_us) >= budget or edge_hits >= 2:
                    break
                n_try += 1
                alpha = float(min(0.85, max(0.08, alpha)))
                # Converge to T: predicted f must be closer to T than both ends'
                # worse endpoint, and stay inside soft outer shell.
                f_pred = float(fi + alpha * (fj - fi))
                if abs(f_pred - t) > f_cap + 1e-15:
                    continue
                if abs(f_pred - t) >= min(abs(fi - t), abs(fj - t)) - 1e-15:
                    # Not strictly between the better end and T — skip midpoints
                    # that sit mid-chord far from the level set.
                    if abs(f_pred - t) > abs(fi - t) and abs(f_pred - t) > abs(fj - t):
                        continue
                cand = clip_unit_cube(ui + alpha * (uj - ui))
                if float(np.linalg.norm(cand - ui)) < r_edge * 0.5:
                    continue
                if float(np.linalg.norm(cand - uj)) < r_edge * 0.5:
                    continue
                if not sep_index.far_enough(cand):
                    continue
                physical = map_u_to_physical(cand, self.vars)
                if self._selectionexp and not evaluate_selection(
                    self._selectionexp,
                    physical,
                    context=self._expression_context,
                ):
                    continue
                if self._unique_point_count() + len(new_us) >= self._max_points:
                    break
                sep_index.add(cand)
                new_us.append(cand)
                n_accept += 1
                edge_hits += 1
        self._logger.debug(
            "AdaptiveBridson root-correction detail\n"
            "    brackets raw / kept   -> %d / %d\n"
            "    trials / accepted     -> %d / %d\n"
            "    r_floor (blue noise)  -> %g  (toward-T only)",
            n_raw,
            len(ranked),
            n_try,
            n_accept,
            r_floor,
        )
        return new_us

    def _max_core_nn_gap(self, core_indices: Sequence[int]) -> float | None:
        """Max nearest-neighbor distance among cores in u-space (None if <2)."""
        cores = [int(i) for i in core_indices]
        if len(cores) < 2:
            return None
        pts = np.asarray(
            [self._points[i].u for i in cores], dtype=np.float64
        )
        try:
            from scipy.spatial import cKDTree

            tree = cKDTree(pts)
            dists, _ = tree.query(pts, k=2)
            # column 0 is self (0); column 1 is NN
            nn = np.asarray(dists[:, 1], dtype=np.float64)
            return float(np.max(nn))
        except Exception:
            max_gap = 0.0
            for a in range(len(pts)):
                best = float("inf")
                for b in range(len(pts)):
                    if a == b:
                        continue
                    d = float(np.linalg.norm(pts[a] - pts[b]))
                    if d < best:
                        best = d
                if best < float("inf") and best > max_gap:
                    max_gap = best
            return max_gap

    def _core_bounds(
        self, core_indices: Sequence[int]
    ) -> tuple[np.ndarray, np.ndarray] | None:
        """Axis-aligned u-space support bounds of the live core cloud."""
        cores = [int(i) for i in core_indices]
        if not cores:
            return None
        pts = np.asarray([self._points[i].u for i in cores], dtype=np.float64)
        return np.min(pts, axis=0), np.max(pts, axis=0)

    @staticmethod
    def _bounds_expanded(
        before: tuple[np.ndarray, np.ndarray] | None,
        after: tuple[np.ndarray, np.ndarray] | None,
        *,
        min_gain: float,
    ) -> bool:
        """True when new cores extend an explored contour end in u-space."""
        if before is None:
            return after is not None
        if after is None:
            return False
        lo_before, hi_before = before
        lo_after, hi_after = after
        gain = max(
            float(np.max(lo_before - lo_after)),
            float(np.max(hi_after - hi_before)),
        )
        return gain > max(float(min_gain), 1e-12)

    @staticmethod
    def _fill_made_structural_progress(
        *,
        core_gain: int,
        coverage_before: bool,
        best_improved: bool,
        extent_improved: bool,
        gap_improved: bool,
    ) -> bool:
        """Ignore interior core churn once coverage already exists."""
        return bool(
            best_improved
            or extent_improved
            or gap_improved
            or (not coverage_before and int(core_gain) > 0)
        )

    def _core_coverage_ok(self, cores: Sequence[int]) -> bool:
        """True when the core set covers the contour at the *current* r_g.

        Coverage = enough cores, one connected component, and max
        nearest-neighbor gap ≤ spacing factor·r_g.
        Until this holds, keep densifying and **appending** cores (do not shrink).
        """
        cores = [int(i) for i in cores]
        if len(cores) < int(self._min_cores_for_coverage):
            return False
        if len(self._core_components(cores)) != 1:
            return False
        gap = self._max_core_nn_gap(cores)
        if gap is None:
            return len(cores) >= int(self._min_cores_for_coverage)
        limit = float(self._core_spacing_factor) * float(self._radius)
        return gap <= limit + 1e-12

    def _effective_final_half_width(self) -> float:
        """Half-width for accepted contour export (final if set, else core)."""
        if self._final_half_width is not None:
            return float(self._final_half_width)
        return float(self._core_half_width)

    def _final_indices(self) -> list[int]:
        """Indices with |f−T| ≤ effective final half-width."""
        w = self._effective_final_half_width()
        t = self._target_value
        out: list[int] = []
        for i, p in enumerate(self._points):
            fi = _finite_f(p.f)
            if fi is None:
                continue
            if abs(fi - t) <= w + 1e-15:
                out.append(i)
        return out

    def _solidify_cores(self, cores: Sequence[int], *, reason: str) -> int:
        """Freeze current cores into the solidified set (only on r_g shrink)."""
        n_new = 0
        for i in cores:
            uuid = str(self._points[int(i)].uuid or "")
            if not uuid or uuid in self._solidified_core_uuids:
                continue
            self._solidified_core_uuids.add(uuid)
            n_new += 1
        if n_new or cores:
            self._logger.debug(
                "AdaptiveBridson solidify cores\n"
                "    reason                -> %s\n"
                "    cores this scale      -> %d  (+%d newly solidified)\n"
                "    solidified total      -> %d\n"
                "    r_g                   -> %g",
                reason,
                len(cores),
                n_new,
                len(self._solidified_core_uuids),
                self._radius,
            )
        return n_new

    def _contour_converged(
        self,
        *,
        t_min: float | None,
        t_max: float | None,
        fill_needed: bool,
    ) -> bool:
        """Primary exit: fill done and 2*r_g band (t_max−t_min) < threshold.

        Optional: if final_half_width is configured, require ≥1 final point.
        """
        if fill_needed:
            return False
        best = self._best_index()
        if best is None or not self._is_core_index(best):
            return False
        cores, _frontier, _outer, _classified_best = self._classify_bands()
        if not self._core_coverage_ok(cores):
            return False
        if not self._band_width_ok(t_min, t_max):
            return False
        if self._final_half_width_configured and not self._final_indices():
            return False
        return True

    def _should_advance_generation(
        self,
        *,
        t_min: float | None,
        t_max: float | None,
        fill_needed: bool,
        cores: Sequence[int],
    ) -> bool:
        """Next gen (r_g shrink) only after fill is done and band not yet thin.

        User flow:
          fill/bridge first → then if tmin/tmax OK exit, else propose next gen.
        """
        if fill_needed:
            return False
        if self._band_width_ok(t_min, t_max):
            return False  # should exit, not shrink
        if self._at_min_radius():
            return False
        if self._radius_shrink_mode == "every_generation":
            return True
        # Any finite best point is enough to refine the resolution.  In
        # particular, no-core bootstrap scales must shrink when strict
        # blue-noise packing leaves no legal same-scale candidates.
        if not cores and t_min is None:
            return self._best_index() is not None
        return True

    def _at_min_radius(self) -> bool:
        """True when r_g has reached the configured u-space resolution floor."""
        return self._radius <= self._min_radius + 1e-15

    def _next_radius(self) -> float:
        """Shrink r_g without crossing the configured minimum radius."""
        return max(
            float(self._min_radius),
            float(self._radius) * float(self._refinement_factor),
        )

    def _maybe_shrink_outer(
        self,
        *,
        n_core: int,
        open_brackets: int,
        raw_brackets: int = 0,
    ) -> bool:
        """Annealing: w_outer → w_core only when local nearest brackets = 0.

        ``open_brackets`` is the count of local nearest opposite-side pairs on
        the full control cloud. Outer-only pools hide real sign changes.
        """
        if n_core <= 0:
            self._logger.debug(
                "AdaptiveBridson outer anneal blocked (no cores)\n"
                "    raw full-cloud brackets      -> %d\n"
                "    local nearest brackets       -> %d",
                int(raw_brackets),
                int(open_brackets),
            )
            return False
        if int(open_brackets) > 0:
            self._logger.debug(
                "AdaptiveBridson outer anneal blocked (brackets remain)\n"
                "    raw full-cloud brackets      -> %d\n"
                "    local nearest brackets       -> %d\n"
                "    cores                        -> %d\n"
                "    outer half-w                 -> %g",
                int(raw_brackets),
                int(open_brackets),
                int(n_core),
                self._outer_half_width,
            )
            return False
        if self._outer_half_width <= self._core_half_width + 1e-15:
            return False
        new_w = max(
            self._core_half_width,
            self._outer_half_width * self._outer_shrink_factor,
        )
        if new_w >= self._outer_half_width - 1e-15:
            new_w = self._core_half_width
        old = self._outer_half_width
        self._outer_half_width = float(new_w)
        self._logger.debug(
            "AdaptiveBridson outer-band anneal  %g → %g\n"
            "    core half-w              -> %g\n"
            "    raw Delaunay / nearest    -> %d / %d  (must be 0 nearest)",
            old,
            self._outer_half_width,
            self._core_half_width,
            int(raw_brackets),
            int(open_brackets),
        )
        return True

    # ----------------------------------------------------------- local densify
    @staticmethod
    def _mst_bridge_pairs(
        positions: list[np.ndarray] | np.ndarray,
        *,
        min_dist: float,
        max_dist: float,
    ) -> list[tuple[int, int, float]]:
        """Kruskal MST edges with length in (min_dist, max_dist].

        Builds a **radius graph** via cKDTree (O(n log n + |E|)), not the
        dense O(n²) all-pairs loop — critical once live cores reach O(10³).
        """
        if isinstance(positions, np.ndarray):
            pts = np.asarray(positions, dtype=np.float64)
        else:
            n0 = len(positions)
            if n0 < 2 or max_dist <= 0:
                return []
            pts = np.asarray(positions, dtype=np.float64)
        n = int(pts.shape[0])
        if n < 2 or max_dist <= 0:
            return []

        # Sparse candidate edges: only pairs within max_dist.
        edges: list[tuple[float, int, int]] = []
        try:
            from scipy.spatial import cKDTree

            tree = cKDTree(pts)
            # query_pairs is far faster than nested Python loops for n≳500.
            pairs_idx = tree.query_pairs(r=float(max_dist), output_type="ndarray")
            if pairs_idx.size:
                for i, j in pairs_idx:
                    i = int(i)
                    j = int(j)
                    d = float(np.linalg.norm(pts[i] - pts[j]))
                    if d <= max_dist + 1e-15:
                        edges.append((d, i, j))
        except Exception:
            # Fallback (tiny n / no scipy): dense pairs.
            for i in range(n):
                pi = pts[i]
                for j in range(i + 1, n):
                    d = float(np.linalg.norm(pi - pts[j]))
                    if d <= max_dist + 1e-15:
                        edges.append((d, i, j))

        edges.sort(key=lambda t: (t[0], t[1], t[2]))
        parent = list(range(n))

        def find(a: int) -> int:
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a

        def unite(a: int, b: int) -> bool:
            ra, rb = find(a), find(b)
            if ra == rb:
                return False
            parent[rb] = ra
            return True

        pairs: list[tuple[int, int, float]] = []
        for d, i, j in edges:
            if unite(i, j) and d > min_dist + 1e-15:
                pairs.append((i, j, d))
        return pairs

    def _try_accept_candidate(
        self,
        cand: np.ndarray,
        *,
        sep_index: _SepIndex,
        new_us: list[np.ndarray],
        budget: int,
    ) -> bool:
        """Shared accept path for ball densify and gap bridges."""
        if len(new_us) >= budget:
            return False
        if self._unique_point_count() + len(new_us) >= self._max_points:
            return False
        cand = clip_unit_cube(np.asarray(cand, dtype=np.float64))
        if self._u_key(cand) in self._submitted_u_keys:
            return False
        if not sep_index.far_enough(cand):
            return False
        physical = map_u_to_physical(cand, self.vars)
        if self._selectionexp and not evaluate_selection(
            self._selectionexp,
            physical,
            context=self._expression_context,
        ):
            return False
        sep_index.add(cand)
        new_us.append(cand)
        return True

    def _bridge_core_gaps(
        self,
        positions: np.ndarray,
        *,
        r_window: float,
        r_sep: float,
        sep_index: _SepIndex,
        new_us: list[np.ndarray],
        budget: int,
        rng: np.random.Generator,
        core_local_ids: Sequence[int] | None = None,
    ) -> int:
        """Place points along MST segments between nearby live cores.

        ``positions`` is (n_core, d). If ``core_local_ids`` is set, only those
        local rows participate in the MST (thinned centers); otherwise all.
        """
        if not self._bridge_gaps or r_window <= 0 or r_sep <= 0:
            return 0
        if positions.shape[0] < 2:
            return 0
        if core_local_ids is None:
            pos = positions
        else:
            ids = [int(i) for i in core_local_ids]
            if len(ids) < 2:
                return 0
            pos = positions[np.asarray(ids, dtype=np.int64)]
        max_dist = float(self._bridge_span_factor) * float(r_window)
        pairs = self._mst_bridge_pairs(
            pos, min_dist=r_sep, max_dist=max_dist
        )
        n_before = len(new_us)
        for li, lj, dist in pairs:
            if len(new_us) >= budget:
                break
            if dist <= r_sep + 1e-15:
                continue
            n_slots = max(1, int(math.floor(dist / r_sep)) - 1)
            n_slots = min(n_slots, 8)
            u_a = pos[li]
            u_b = pos[lj]
            direction = u_b - u_a
            for s in range(1, n_slots + 1):
                if len(new_us) >= budget:
                    break
                t = s / (n_slots + 1)
                base = u_a + t * direction
                jitter_scale = 0.15 * r_sep
                noise = rng.normal(size=self._dim)
                dn = float(np.linalg.norm(direction))
                if dn > 1e-15:
                    unit = direction / dn
                    noise = noise - float(np.dot(noise, unit)) * unit
                nn = float(np.linalg.norm(noise))
                if nn > 1e-15:
                    noise = noise / nn * jitter_scale * float(rng.random())
                else:
                    noise = np.zeros(self._dim)
                cand = base + noise
                self._try_accept_candidate(
                    cand,
                    sep_index=sep_index,
                    new_us=new_us,
                    budget=budget,
                )
        return len(new_us) - n_before

    def _densify_local(
        self,
        core_indices: Sequence[int],
        *,
        r_window: float,
        r_sep: float,
        generation: int,
        extra_accepted: Sequence[np.ndarray] | None = None,
    ) -> list[np.ndarray]:
        """Propose new points: per-core r_window balls + inter-core gap bridges.

        Performance (late gens with O(10³) live cores):
          * spatial-hash exclusion (not O(N) scan per trial)
          * thinned proposal centers (spacing ~ r_window) — all cores still
            constrain Accepted, but only a covering subset opens windows
          * sparse KDTree MST for gap bridges

        ``extra_accepted``: already-proposed u-points this generation (e.g.
        root-correction) that must also exclude.
        """
        if not core_indices or r_window <= 0 or r_sep <= 0:
            return []
        t0 = time.perf_counter()
        rng = self._fill_rng(generation)
        cores = sorted(int(i) for i in core_indices)
        positions = np.asarray(
            [self._points[i].u for i in cores], dtype=np.float64
        )

        # Exclusion: every coordinate submitted over the whole run, including
        # sites pruned from the live control cloud, plus same-pass proposals.
        # Proposal centers remain the true core set.
        sep_index = self._submitted_sep_index(r_sep)
        if extra_accepted:
            for u in extra_accepted:
                sep_index.add(np.asarray(u, dtype=np.float64))

        # Prefer densifying around *every* core. Only thin if the core set is
        # huge (performance); still process cores first so none are dropped
        # in favor of later non-core centers.
        n_core = len(cores)
        if n_core <= 128:
            center_local = list(range(n_core))
        else:
            center_local = _thin_centers(
                positions,
                min_spacing=max(r_sep, 1e-15),
                order=list(range(n_core)),
            )
        new_us: list[np.ndarray] = []
        budget = self._max_new_per_generation

        # --- (1) densify near cores in a Bridson annulus ---
        # The caller uses [r_g, 2*r_g]; the history index enforces the inner
        # radius globally, including against pruned sites.
        r_local = min(float(r_window), max(float(r_sep) * 2.0, float(r_window) * 0.35))
        r_sep_local = min(float(r_sep), r_local * 0.5)
        for loc in center_local:
            if len(new_us) >= budget:
                break
            if self._unique_point_count() + len(new_us) >= self._max_points:
                break
            center = positions[loc]
            for _ in range(self._k_ref):
                if len(new_us) >= budget:
                    break
                if self._unique_point_count() + len(new_us) >= self._max_points:
                    break
                cand = sample_ball(center, r_local, rng, self._dim)
                dist_core = float(np.linalg.norm(cand - center))
                if dist_core > r_local + 1e-15:
                    continue
                # Together with the outer radius this is the standard Bridson
                # annulus [r_sep, r_window] around the active center.
                if dist_core < r_sep_local:
                    continue
                self._try_accept_candidate(
                    cand,
                    sep_index=sep_index,
                    new_us=new_us,
                    budget=budget,
                )

        n_ball = len(new_us)

        # --- (2) core–core MST bridge only along short gaps (anti-segmentation) ---
        # Cap bridge span so we do not draw long chords that fatten the band.
        n_bridge = self._bridge_core_gaps(
            positions,
            r_window=min(r_local * 1.5, float(r_window)),
            r_sep=r_sep_local,
            sep_index=sep_index,
            new_us=new_us,
            budget=budget,
            rng=rng,
            core_local_ids=center_local,
        )

        elapsed_ms = (time.perf_counter() - t0) * 1000.0
        self._logger.debug(
            "AdaptiveBridson densify (local-core Bridson + gap bridge)\n"
            "    cores (frozen)       -> %d  (proposal centers thinned → %d)\n"
            "    r_window / r_sep     -> %g / %g\n"
            "    bridge               -> %s  (span=%g → %d pts)\n"
            "    new proposals        -> %d  (ball=%d + bridge=%d / budget=%d, k_ref=%d)\n"
            "    densify wall         -> %.1f ms",
            len(cores),
            len(center_local),
            r_window,
            r_sep,
            "on" if self._bridge_gaps else "off",
            float(self._bridge_span_factor) * r_window,
            n_bridge,
            len(new_us),
            n_ball,
            n_bridge,
            budget,
            self._k_ref,
            elapsed_ms,
        )
        return new_us

    # --------------------------------------------------------------- submit
    def _next_uuid(self) -> str:
        idx = self._accepted_index
        self._accepted_index += 1
        return deterministic_sampler_uuid(
            prefix="abridson",
            seed=self._seed,
            sample_index=idx,
        )

    def _submit_generation(self, us: list[np.ndarray], *, generation: int) -> int:
        if not us:
            return 0
        self._require_redis("AdaptiveBridson")
        samples: list[Sample] = []
        # Submission boundary is the final invariant: no caller can bypass the
        # current-r_g blue-noise filter accidentally.
        for u in self._filter_blue_noise(us, min_distance=self._radius):
            if self._unique_point_count() >= self._max_points:
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
            self._submitted_u_keys.add(self._u_key(u))
            sample = self._build_sample(np.asarray(u, dtype=np.float64))
            sample.uuid = uuid
            samples.append(sample)
        return self._submit_sample_batch(samples)

    def _wait_generation(self, *, timeout: float) -> None:
        results = self.wait_for_generation(timeout=timeout)
        self.absorb_generation(results)

    def propose_generation(self) -> Sequence[Sample] | None:
        return None

    def absorb_generation(self, results: Sequence[Mapping[str, Any]]) -> None:
        from jarvishep2.feedback_return import (
            WIRE_LOGL_KEY,
            feedback_logl,
            is_unusable_logl,
            normalize_feedback_record,
        )

        for record in results:
            uuid = str(record.get("uuid", ""))
            index = self._uuid_to_index.get(uuid)
            if index is None:
                continue
            flat = normalize_feedback_record(record)
            payload = dict(self._points[index].x)
            for key, value in flat.items():
                if key == "uuid":
                    continue
                if key == WIRE_LOGL_KEY:
                    payload["LogL"] = value
                    payload["logL"] = value
                else:
                    payload[str(key)] = value
            target_expr = str(self._target_expression or "").strip()
            target_is_logl = target_expr in {"LogL", "logL"}
            if (
                target_is_logl
                and WIRE_LOGL_KEY in flat
                and is_unusable_logl(feedback_logl(flat))
            ):
                self._points[index].f = None
                self._failed_regions.append(
                    {"uuid": uuid, "u": list(self._points[index].u)}
                )
                continue
            try:
                self._points[index].f = self._eval_target(payload)
            except Exception:
                self._points[index].f = None
                self._failed_regions.append(
                    {"uuid": uuid, "u": list(self._points[index].u)}
                )

    def _prune_to_indices(self, keep: Sequence[int]) -> None:
        """Keep only listed indices as the next generation's seed cloud."""
        keep_set = sorted(set(int(i) for i in keep if 0 <= int(i) < len(self._points)))
        new_points = [self._points[i] for i in keep_set]
        self._points = new_points
        self._uuid_to_index = {p.uuid: i for i, p in enumerate(self._points) if p.uuid}

    # ---------------------------------------------------------------- finalize
    def finalize(self) -> dict[str, Any]:
        var_names = [str(v.name) for v in self.vars]
        finite = [
            (i, float(p.f))
            for i, p in enumerate(self._points)
            if _finite_f(p.f) is not None
        ]
        cores, frontier, outer, _best = self._classify_bands()
        final_idx = self._final_indices()
        w_final = self._effective_final_half_width()
        final_u = [list(self._points[i].u) for i in final_idx]
        final_f = [self._points[i].f for i in final_idx]
        final_x = [dict(self._points[i].x) for i in final_idx]

        # Primary point cloud: when final_half_width is configured, export only
        # accepted final contour points (not the whole densify core set).
        if self._final_half_width_configured:
            export_u, export_f, export_x = final_u, final_f, final_x
            n_export = len(final_idx)
        else:
            export_u = [list(p.u) for p in self._points]
            export_f = [p.f for p in self._points]
            export_x = [dict(p.x) for p in self._points]
            n_export = len(self._points)

        payload: dict[str, Any] = {
            "dim": self._dim,
            "algorithm": "outer_core_root_correction_bridson",
            "variable_names": var_names,
            "target_expression": self._target_expression,
            "target_value": self._target_value,
            "t_min": self._t_min,
            "t_max": self._t_max,
            "band_width": (
                None
                if self._t_min is None or self._t_max is None
                else float(self._t_max - self._t_min)
            ),
            "core_half_width": self._core_half_width,
            "outer_half_width": self._outer_half_width,
            "final_half_width": w_final,
            "final_half_width_configured": bool(self._final_half_width_configured),
            "threshold": self._threshold,
            "n_points_total": n_export,
            "n_unique_submitted": self._unique_point_count(),
            "n_control_points": len(self._points),
            "n_finite_f": len(finite),
            # generation == r_g shrink count; report scales completed + current.
            "n_generations": int(self._generation) + 1,
            "n_radius_refines": int(self._n_radius_refines),
            "fill_pass": int(self._fill_pass),
            "radius": self._radius,
            "min_radius": self._min_radius,
            "radius_shrink_mode": self._radius_shrink_mode,
            "converged": self._converged,
            "contour_status": "converged" if self._converged else "partial",
            "stop_reason": self._stop_reason,
            "failed_regions": list(self._failed_regions),
            "points_u": export_u,
            "points_f": export_f,
            "points_x": export_x,
            "live_core_count": len(cores),
            "n_frontier": len(frontier),
            "n_outer": len(outer),
            "open_brackets": int(self._open_brackets),
            "n_bracket_support": len(self._bracket_support_indices),
            "n_virtual_core_anchors": len(self._virtual_core_anchors),
            "n_final_points": len(final_idx),
            "final_points_u": final_u,
            "final_points_x": final_x,
            "final_points_f": final_f,
        }
        # Core cloud for diagnostics / plot hooks.
        if cores:
            payload["live_band_points_u"] = [
                list(self._points[i].u) for i in cores
            ]
            payload["n_live_band"] = len(cores)
        payload["frontier_points_u"] = [
            list(self._points[i].u) for i in frontier
        ]
        self._levelset = payload
        out_dir = str(self.config.get("task_result_dir") or os.getcwd())
        try:
            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, "levelset.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            self._logger.info(
                "AdaptiveBridson finished: reason=%s, converged=%s, "
                "points=%d, gen=%d, cores=%d, final=%d, "
                "half_width=%g%s, band=%s..%s (Δ=%s), wrote=%s",
                self._stop_reason,
                self._converged,
                payload["n_points_total"],
                payload["n_generations"],
                len(cores),
                len(final_idx),
                w_final,
                " configured" if self._final_half_width_configured else " =core",
                f"{self._t_min:g}" if self._t_min is not None else "n/a",
                f"{self._t_max:g}" if self._t_max is not None else "n/a",
                (
                    f"{self._t_max - self._t_min:g}"
                    if self._t_min is not None and self._t_max is not None
                    else "n/a"
                ),
                path,
            )
        except OSError as exc:
            self._logger.warning("failed to write levelset.json: %s", exc)
        return payload

    # ------------------------------------------------------------------- driver
    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Live-band loop; timeout applies independently to each generation."""
        if timeout is not None:
            generation_timeout = timeout
        self._require_redis("AdaptiveBridson.run_adaptive")
        self._ensure_seed_sequence()
        if self._graph is None:
            self._graph = self._make_graph()
        if self._target_fn is None:
            self._compile_target()

        if self._converged and self._levelset is not None and not self._pending_uuids:
            self._logger.info(
                "AdaptiveBridson already converged (checkpoint); skip run_adaptive"
            )
            return 0

        self._logger.info(
            "AdaptiveBridson start: target=%s=%g, r0=%g, min_r=%g, "
            "outer/core=%g/%g, shrink=%s, max_gen=%d, max_points=%d",
            self._target_expression,
            self._target_value,
            self._initial_radius,
            self._min_radius,
            self._outer_half_width,
            self._core_half_width,
            self._radius_shrink_mode,
            self._max_generations,
            self._max_points,
        )

        total_submitted = 0
        if not self._submitted_u_keys and self._points:
            self._submitted_u_keys = {self._u_key(p.u) for p in self._points}

        if not self._points:
            self._radius = self._initial_radius
            self._fill_pass = 0
            self._quiet_fill_passes = 0
            self._reset_endpoint_state()
            gen0 = self._generation_0_points()
            self._generation = 0  # scale index; no r_g shrink yet
            n = self._submit_generation(gen0, generation=0)
            total_submitted += n
            self._logger.debug(
                "AdaptiveBridson gen=0 (r_g=%g) initial Bridson submitted %d; "
                "waiting feedback…",
                self._radius,
                n,
            )
            self._wait_generation(timeout=generation_timeout)
            self.checkpoint_at_barrier(reason="als_generation_0")
        elif self._pending_uuids:
            self._logger.debug(
                "AdaptiveBridson resume: waiting %d pending feedback…",
                len(self._pending_uuids),
            )
            self._wait_generation(timeout=generation_timeout)

        while True:
            scale = self._prepare_scale_state()
            cores, frontier, outer = scale.cores, scale.frontier, scale.outer
            best, t_min, t_max = scale.best, scale.t_min, scale.t_max

            # Full-cloud straddles for root-correction / bridge fill.
            brackets = self._straddle_brackets(pool=None)
            self._open_brackets = len(brackets)
            nearest_brackets = self._nearest_straddle_brackets(
                max_distance=2.5 * float(self._radius)
            )
            local_pairs = nearest_brackets
            n_filt = len(nearest_brackets)
            local_brackets = nearest_brackets
            virtual_cores = self._bracket_root_anchors(local_brackets)
            self._bracket_support_indices = sorted(
                {
                    int(index)
                    for edge in local_brackets
                    for index in (edge[0], edge[1])
                }
            )
            self._virtual_core_anchors = [anchor.tolist() for anchor in virtual_cores]
            n_active_virtual = self._refresh_virtual_cores(virtual_cores, cores)
            coverage_ok = self._core_coverage_ok(cores)
            bridge_cls = self._classify_bridge_edges(
                local_pairs, cores, r_window=self._radius
            )
            n_A = len(bridge_cls["A"])
            n_B = len(bridge_cls["B"])
            n_C = len(bridge_cls["C"])
            n_D = len(bridge_cls["D"])
            endpoint_edges = self._actionable_endpoint_edges(
                local_pairs,
                cores,
                r_sep=self._radius,
            )
            n_active_endpoints = self._refresh_active_endpoints(cores)
            fill_needed = self._fill_needed(
                cores=cores,
                n_local_brackets=n_filt,
                coverage_ok=coverage_ok,
                densify_new=self._last_densify_new,
                densify_attempted=self._last_densify_attempted,
                n_bridge_A=n_A,
                n_bridge_B=n_B,
                n_bridge_C=n_C,
                n_endpoint_edges=len(endpoint_edges),
                n_active_endpoints=n_active_endpoints,
                n_active_virtual_cores=n_active_virtual,
                scale_quiet=(
                    self._quiet_fill_passes >= self._quiet_fill_limit
                ),
                force_scale_close=self._fill_pass >= self._max_fill_passes,
            )
            band_ok = self._band_width_ok(t_min, t_max)
            band_delta = (
                None
                if t_min is None or t_max is None
                else float(t_max) - float(t_min)
            )
            tau = (
                float(self._threshold)
                if self._threshold is not None
                else float(self._core_half_width)
            )

            best_f = (
                float(self._points[best].f)  # type: ignore[arg-type]
                if best is not None and _finite_f(self._points[best].f) is not None
                else float("nan")
            )
            gap = self._max_core_nn_gap(cores)
            gap_limit = float(self._core_spacing_factor) * float(self._radius)
            self._logger.debug(
                "AdaptiveBridson band check gen=%d fill_pass=%d  "
                "(gen = r_g shrink count; fill_pass = same-r_g densify)\n"
                "    best index / f         -> %s / %g\n"
                "    t_min / t_max / Δ      -> %s / %s / %s  (threshold=%g, band_ok=%s)\n"
                "    fill needed            -> %s\n"
                "    outer / core half-w    -> %g / %g\n"
                "    cores / frontier / all -> %d / %d / %d\n"
                "    brackets Delaunay/nearest-> %d / %d\n"
                "    bracket support / roots-> %d / %d  (active uncovered=%d; nearest-pairs=%d; edge≤2.5 r_g)\n"
                "    bridge A/B/C/D         -> %d / %d / %d / %d  "
                "(A=no-core, B=one-core, C=core-gap, D=done)\n"
                "    endpoints straddle/geom-> %d / %d  (persistent geom fronts)\n"
                "    coverage OK            -> %s  (components=%d, NN gap %s / limit %g, min_cores=%d)\n"
                "    points / r_g / refines -> %d / %g / %d",
                int(self._generation),
                int(self._fill_pass),
                str(best) if best is not None else "n/a",
                best_f,
                f"{t_min:g}" if t_min is not None else "n/a",
                f"{t_max:g}" if t_max is not None else "n/a",
                f"{band_delta:g}" if band_delta is not None else "n/a",
                tau,
                "yes" if band_ok else "no",
                "yes (root-corr/expand/densify/bridge)" if fill_needed else "no (fill done)",
                self._outer_half_width,
                self._core_half_width,
                len(cores),
                len(frontier),
                len(self._points),
                self._open_brackets,
                n_filt,
                len(self._bracket_support_indices),
                len(virtual_cores),
                n_active_virtual,
                len(nearest_brackets),
                n_A,
                n_B,
                n_C,
                n_D,
                len(endpoint_edges),
                n_active_endpoints,
                "yes" if coverage_ok else "no",
                len(self._core_components(cores)),
                f"{gap:g}" if gap is not None else "n/a",
                gap_limit,
                self._min_cores_for_coverage,
                len(self._points),
                self._radius,
                self._n_radius_refines,
            )

            decision = self._pre_fill_decision(
                best=best,
                t_min=t_min,
                t_max=t_max,
                fill_needed=fill_needed,
            )

            # Every terminal path is selected by the named decision above;
            # these branches retain their established diagnostics and cleanup.
            if decision.reason == "level-set not present in domain":
                self._stop_reason = decision.reason
                self._converged = False
                self._logger.warning(
                    "AdaptiveBridson stop: %s (no finite f)", self._stop_reason
                )
                self.finalize()
                return total_submitted

            # (3) Exit when fill done AND (t_max − t_min) < threshold.
            if decision.reason == "converged":
                self._converged = True
                self._stop_reason = decision.reason
                self._solidify_cores(cores, reason="exit_converged_tmin_tmax")
                self._logger.debug(
                    "AdaptiveBridson stop: converged (t_max-t_min < threshold)\n"
                    "    t_min / t_max / Δ      -> %g / %g / %g\n"
                    "    threshold              -> %g\n"
                    "    cores / r_g / gen      -> %d / %g / %d",
                    float(t_min),  # type: ignore[arg-type]
                    float(t_max),  # type: ignore[arg-type]
                    float(band_delta),  # type: ignore[arg-type]
                    tau,
                    len(cores),
                    self._radius,
                    int(self._generation),
                )
                self._prune_to_indices(cores if cores else outer)
                self.finalize()
                return total_submitted

            # The u-space resolution floor is a normal partial-completion
            # condition. Finish endpoint growth at this scale first; once fill
            # is quiet, do not create a finer generation below min_radius.
            if decision.reason == "min_radius":
                self._stop_reason = decision.reason
                self._converged = False
                self._solidify_cores(cores, reason="exit_min_radius")
                self._logger.debug(
                    "AdaptiveBridson stop: minimum radius reached\n"
                    "    r_g / min_radius      -> %g / %g\n"
                    "    t_min / t_max / Δ      -> %s / %s / %s\n"
                    "    threshold / cores      -> %g / %d",
                    self._radius,
                    self._min_radius,
                    f"{t_min:g}" if t_min is not None else "n/a",
                    f"{t_max:g}" if t_max is not None else "n/a",
                    f"{band_delta:g}" if band_delta is not None else "n/a",
                    tau,
                    len(cores),
                )
                self._prune_to_indices(cores if cores else outer)
                self.finalize()
                return total_submitted

            # generation == r_g shrink count; stop after enough scale refinements.
            if decision.reason == "max_generations":
                self._stop_reason = decision.reason
                self._converged = False
                self._logger.warning(
                    "AdaptiveBridson stop: max_generations=%d "
                    "(r_g shrinks done; cores=%d, band_ok=%s, fill_needed=%s, r_g=%g)",
                    self._max_generations,
                    len(cores),
                    band_ok,
                    fill_needed,
                    self._radius,
                )
                self._prune_to_indices(cores if cores else outer)
                self.finalize()
                return total_submitted

            if int(self._fill_pass) >= int(self._max_fill_passes):
                self._logger.warning(
                    "AdaptiveBridson scale quiet: max_fill_passes=%d at "
                    "gen=%d r_g=%g (cores=%d, band_ok=%s) — force t-band "
                    "decision / r_g shrink",
                    self._max_fill_passes,
                    int(self._generation),
                    self._radius,
                    len(cores),
                    band_ok,
                )

            if decision.reason == "max_points":
                self._stop_reason = decision.reason
                self._converged = False
                self._logger.warning(
                    "AdaptiveBridson stop: max_points=%d", self._max_points
                )
                self._prune_to_indices(cores if cores else outer)
                self.finalize()
                return total_submitted

            r_window = self._radius
            # Current r_g is the global blue-noise exclusion distance for all
            # proposal types at this scale.
            r_sep = self._radius
            budget = self._max_new_per_generation
            new_us: list[np.ndarray] = []

            # --- Fill phase (same gen): A/B/C, converge toward T, no gen++ ---
            # A: neither end core → secant toward T only
            # B: one end core → step along edge toward T (not outward off-band)
            # C: both cores, large u-gap → densify/bridge among *band* cores
            # densify balls: tight around in-band cores only
            n_densify = 0
            n_expand_b = 0
            n_endpoint_growth = 0
            n_direct_root = 0
            n_virtual_core = 0
            densify_attempted = False
            densify_centers: list[int] = []

            if fill_needed:
                # Edge/root proposals share the strict current-r_g blue-noise
                # exclusion over the complete submitted history.
                sep_fill = self._submitted_sep_index(r_sep)

                # Persistent geometric fronts get first claim on this pass's
                # budget. Sample a full hypersphere, wait for the entire batch
                # feedback barrier, then rebuild/rejudge the endpoint set.
                if self._active_endpoints:
                    n_endpoint_growth = self._active_endpoint_proposals(
                        sep_index=sep_fill,
                        new_us=new_us,
                        budget=budget,
                        rng=self._fill_rng(int(self._generation)),
                    )
                    self._logger.debug(
                        "AdaptiveBridson geometric endpoint growth\n"
                        "    active fronts          -> %d\n"
                        "    random annulus         -> [1, 2.5] r_g\n"
                        "    cap / front            -> %d\n"
                        "    proposals              -> %d\n"
                        "    blue-noise             -> r_g",
                        len(self._active_endpoints),
                        self._endpoint_omni_probes,
                        n_endpoint_growth,
                    )

                # A bracket pair already supplies a target-crossing secant.
                # Submit that deterministic root before spending evaluations
                # on random annulus fallback.
                if self._active_virtual_cores and len(new_us) < budget:
                    n_direct_root = self._virtual_core_direct_proposals(
                        sep_index=sep_fill,
                        new_us=new_us,
                        budget=budget,
                    )
                    self._logger.debug(
                        "AdaptiveBridson direct nearest-pair secants\n"
                        "    active uncovered roots -> %d\n"
                        "    direct proposals       -> %d  (before random fallback)",
                        len(self._active_virtual_cores),
                        n_direct_root,
                    )

                corr_budget = min(
                    max(0, budget - len(new_us)),
                    120 if not cores else 80,
                    max(16, self._max_new_per_generation // 4),
                )
                # Case A (+ any unclassified local pair as safety): root correction.
                edges_a = list(bridge_cls["A"])
                if not edges_a and n_filt > 0 and n_A + n_B + n_C == 0:
                    edges_a = list(local_pairs)
                if edges_a:
                    corr = self._root_correct_proposals(
                        edges_a,
                        sep_index=sep_fill,
                        r_sep=r_sep,
                        budget=corr_budget,
                        r_window=r_window,
                    )
                    new_us.extend(corr)
                    self._logger.debug(
                        "AdaptiveBridson fill A (no-core straddle → secant)\n"
                        "    edges A               -> %d\n"
                        "    proposals             -> %d",
                        len(edges_a),
                        len(corr),
                    )

                # Only roots whose direct secant was blocked may use the
                # random annulus. Cap this low-efficiency fallback per pass;
                # this is a proposal budget, not a cap on bracket discovery.
                if self._active_virtual_cores and len(new_us) < budget:
                    fallback_allowance = min(
                        256,
                        max(32, int(self._max_new_per_generation) // 16),
                    )
                    fallback_budget = min(
                        budget, len(new_us) + fallback_allowance
                    )
                    n_virtual_core = self._virtual_core_proposals(
                        sep_index=sep_fill,
                        new_us=new_us,
                        budget=fallback_budget,
                        rng=self._fill_rng(int(self._generation) + 104729),
                    )
                    self._logger.debug(
                        "AdaptiveBridson bracket discovery fallback\n"
                        "    support / detected roots -> %d / %d\n"
                        "    active uncovered roots   -> %d  (>2 r_g from core)\n"
                        "    guided / retrying roots  -> %d / %d\n"
                        "    annulus proposals        -> %d / cap %d  "
                        "(≤2/root, [1, 2.5] r_g)",
                        len(self._bracket_support_indices),
                        len(virtual_cores),
                        len(self._active_virtual_cores),
                        sum(
                            state.get("best_probe_dev") is not None
                            for state in self._active_virtual_cores.values()
                        ),
                        sum(
                            int(state.get("misses", 0)) > 0
                            for state in self._active_virtual_cores.values()
                        ),
                        n_virtual_core,
                        fallback_allowance,
                    )

                # Case B: expand from core along edge toward non-core endpoint.
                if bridge_cls["B"] and len(new_us) < budget:
                    core_set = {int(i) for i in cores}
                    n_expand_b = self._expand_from_core_along_edge(
                        bridge_cls["B"],
                        core_set,
                        sep_index=sep_fill,
                        r_sep=r_sep,
                        budget=budget,
                        new_us=new_us,
                    )
                    self._logger.debug(
                        "AdaptiveBridson fill B (core→toward T along edge)\n"
                        "    edges B               -> %d\n"
                        "    toward-T proposals    -> %d",
                        len(bridge_cls["B"]),
                        n_expand_b,
                    )

                # Isotropic densify ONLY around cores (expand the core set).
                # A/B endpoints must NOT become ball centers — that packs the
                # exclusion grid far off-contour and starves densify near cores.
                densify_centers = list(cores)
                if not densify_centers:
                    seed: list[int] = []
                    if best is not None:
                        seed.append(int(best))
                    seed.extend(int(i) for i in frontier)
                    if len(seed) < 3:
                        scored: list[tuple[float, int]] = []
                        for i, p in enumerate(self._points):
                            d = self._abs_dev(i)
                            if d is not None:
                                scored.append((d, i))
                        scored.sort()
                        for _d, i in scored[: max(3, len(scored))]:
                            if i not in seed:
                                seed.append(i)
                            if len(seed) >= 5:
                                break
                    densify_centers = seed
                    if densify_centers:
                        self._logger.debug(
                            "AdaptiveBridson bootstrap densify (no cores yet)\n"
                            "    seed centers          -> %d  (best+frontier/near-T)",
                            len(densify_centers),
                        )
                else:
                    ab_ends = self._bridge_endpoint_indices(bridge_cls)
                    self._logger.debug(
                        "AdaptiveBridson densify centers = in-band cores only\n"
                        "    n_core_centers        -> %d  (t_min/t_max gated)\n"
                        "    A/B ends (edge→T only)-> %d\n"
                        "    bridge A/B/C          -> %d / %d / %d\n"
                        "    densify window        -> Bridson annulus [r_g, 2*r_g]",
                        len(densify_centers),
                        len(ab_ends),
                        n_A,
                        n_B,
                        n_C,
                    )

                # Case C + densify: Bridson-annulus sampling around in-band
                # cores only, with global current-r_g exclusion.
                if densify_centers and len(new_us) < budget:
                    densify_attempted = True
                    # Standard Bridson proposal annulus [r_g, 2*r_g].
                    # A smaller window would be entirely rejected by the
                    # strict current-r_g history exclusion.
                    r_sep_use = r_sep
                    r_win_use = 2.0 * r_sep
                    old_budget = self._max_new_per_generation
                    self._max_new_per_generation = max(1, budget - len(new_us))
                    try:
                        dens = self._densify_local(
                            densify_centers,
                            r_window=r_win_use,
                            r_sep=r_sep_use,
                            generation=int(self._generation),
                            extra_accepted=list(new_us),
                        )
                    finally:
                        self._max_new_per_generation = old_budget
                    n_densify = len(dens)
                    new_us.extend(dens)

            # Final shared gate for root, endpoint, ball, and bridge proposals:
            # enforce current-r_g blue noise against complete history and peers.
            new_us = self._filter_blue_noise(
                new_us, min_distance=self._radius
            )

            if not new_us:
                # Count the attempted endpoint fan as one miss.  A front is
                # deliberately retained for several distinct fill-pass RNG
                # streams rather than being removed after a single blockage.
                n_active_now = self._refresh_active_endpoints(cores)
                n_active_virtual_now = self._refresh_virtual_cores(
                    virtual_cores, cores
                )
                if (
                    (n_active_now > 0 or n_active_virtual_now > 0)
                    and int(self._fill_pass) + 1 < int(self._max_fill_passes)
                ):
                    self._fill_pass = int(self._fill_pass) + 1
                    self._quiet_fill_passes += 1
                    self._logger.debug(
                        "AdaptiveBridson retry active endpoints at same r_g\n"
                        "    gen / fill_pass / r_g -> %d / %d / %g\n"
                        "    active fronts / roots  -> %d / %d\n"
                        "    quiet attempts         -> %d",
                        int(self._generation),
                        int(self._fill_pass),
                        self._radius,
                        n_active_now,
                        n_active_virtual_now,
                        self._quiet_fill_passes,
                    )
                    self.checkpoint_at_barrier(
                        reason=f"als_gen{int(self._generation)}_endpoint_retry{int(self._fill_pass)}"
                    )
                    continue
                # Fill quiet or no proposals: recompute fill/band, then exit or next gen.
                raw_full = self._straddle_brackets(pool=None)
                filt_full = self._nearest_straddle_brackets(
                    max_distance=2.5 * float(self._radius)
                )
                n_filt_now = len(filt_full)
                coverage_now = self._core_coverage_ok(cores)
                cls_now = self._classify_bridge_edges(
                    filt_full, cores, r_window=self._radius
                )
                endpoint_now = self._actionable_endpoint_edges(
                    filt_full,
                    cores,
                    r_sep=self._radius,
                )
                n_active_now = self._refresh_active_endpoints(cores)
                local_now = filt_full
                virtual_now = self._bracket_root_anchors(local_now)
                n_active_virtual_now = self._refresh_virtual_cores(
                    virtual_now, cores
                )
                fill_now = self._fill_needed(
                    cores=cores,
                    n_local_brackets=n_filt_now,
                    coverage_ok=coverage_now,
                    densify_new=0,
                    densify_attempted=True,
                    n_bridge_A=len(cls_now["A"]),
                    n_bridge_B=len(cls_now["B"]),
                    n_bridge_C=len(cls_now["C"]),
                    n_endpoint_edges=len(endpoint_now),
                    n_active_endpoints=n_active_now,
                    n_active_virtual_cores=n_active_virtual_now,
                    force_scale_close=True,
                )
                # Optional outer anneal only when no local nearest brackets.
                self._maybe_shrink_outer(
                    n_core=len(cores),
                    open_brackets=n_filt_now,
                    raw_brackets=len(raw_full),
                )
                if self._contour_converged(
                    t_min=t_min, t_max=t_max, fill_needed=fill_now
                ):
                    self._converged = True
                    self._stop_reason = "converged"
                    self._solidify_cores(cores, reason="exit_converged_tmin_tmax")
                    self._logger.debug(
                        "AdaptiveBridson stop: converged (t_max-t_min < threshold)\n"
                        "    t_min / t_max / Δ      -> %s / %s / %s\n"
                        "    threshold              -> %g",
                        f"{t_min:g}" if t_min is not None else "n/a",
                        f"{t_max:g}" if t_max is not None else "n/a",
                        (
                            f"{float(t_max) - float(t_min):g}"
                            if t_min is not None and t_max is not None
                            else "n/a"
                        ),
                        tau,
                    )
                    self._prune_to_indices(cores if cores else outer)
                    self.finalize()
                    return total_submitted
                if self._should_advance_generation(
                    t_min=t_min,
                    t_max=t_max,
                    fill_needed=fill_now,
                    cores=cores,
                ):
                    self._advance_generation(cores, t_min=t_min, t_max=t_max)
                    continue
                self._stop_reason = "no_refinement_candidates"
                self._converged = False
                self._logger.warning(
                    "AdaptiveBridson stop: no_refinement_candidates "
                    "at gen=%d fill_pass=%d (cores=%d, fill_needed=%s, band_ok=%s)",
                    int(self._generation),
                    int(self._fill_pass),
                    len(cores),
                    fill_now,
                    self._band_width_ok(t_min, t_max),
                )
                self._prune_to_indices(cores if cores else outer)
                self.finalize()
                return total_submitted

            registered_before = len(self._accumulated_core_uuids)
            best_dev_before = self._abs_dev(best) if best is not None else None
            core_bounds_before = self._core_bounds(cores)
            core_gap_before = self._max_core_nn_gap(cores)

            # Control cloud: keep *in-band* cores (t_min/t_max) + near-T support.
            # Do not let off-band densify junk displace cores. Bracket ends kept
            # only for edge fill continuity.
            keep: list[int] = []
            keep.extend(
                self._accumulated_core_indices(
                    t_min=t_min, t_max=t_max, require_in_band=True
                )
            )
            # Always retain currently absolute cores even if band not yet set.
            keep.extend(cores)
            for i, p in enumerate(self._points):
                if p.uuid and p.uuid in self._solidified_core_uuids:
                    keep.append(i)
                if p.uuid and any(
                    str(endpoint.get("anchor_uuid") or "") == p.uuid
                    for endpoint in self._active_endpoints.values()
                ):
                    keep.append(i)
            if frontier:
                keep.extend(frontier)
            if best is not None:
                keep.append(int(best))
            scored: list[tuple[float, int]] = []
            for i, _p in enumerate(self._points):
                d = self._abs_dev(i)
                if d is not None:
                    scored.append((d, i))
            scored.sort()
            for _d, i in scored[:40]:
                keep.append(i)
            filt_keep = self._nearest_straddle_brackets(
                max_distance=2.5 * float(self._radius)
            )
            for i, j, _fi, _fj in filt_keep:
                keep.append(int(i))
                keep.append(int(j))
            # Preserve local nearest opposite-side evidence independently of
            # the capped Delaunay bracket list. Otherwise pruning can erase
            # exactly the pair needed to repair a visible contour gap.
            self._prune_to_indices(keep)
            # Same-r_g densify: tag samples with current generation (scale), not fill_pass+1.
            n = self._submit_generation(new_us, generation=int(self._generation))
            total_submitted += n
            self._last_densify_new = int(n_densify)
            self._last_densify_attempted = bool(densify_attempted)
            self._fill_pass = int(self._fill_pass) + 1
            self._logger.debug(
                "AdaptiveBridson gen=%d fill_pass=%d submitted %d "
                "(direct-root=%d + virtual-core=%d + endpoint=%d + densify=%d); waiting feedback… "
                "[same-r_g densify @ r_g=%g, coverage=%s]",
                int(self._generation),
                int(self._fill_pass),
                n,
                n_direct_root,
                n_virtual_core,
                n_endpoint_growth,
                n_densify,
                r_window,
                "ok" if coverage_ok else "building",
            )
            self._wait_generation(timeout=generation_timeout)
            # Do NOT increment generation here — only r_g shrink advances gen.

            # After feedback: reclassify; continue fill if needed, else tmin/tmax
            # decides exit (next loop) or next gen (r_g shrink).
            cores_after, _fr, _ou, _ = self._classify_bands()
            best_a, t_min_a, t_max_a = self._radius_t_band()
            self._t_min, self._t_max = t_min_a, t_max_a
            self._register_cores(cores_after)
            core_gain = max(
                0, len(self._accumulated_core_uuids) - registered_before
            )
            cores_after = self._accumulated_core_indices(
                t_min=t_min_a, t_max=t_max_a, require_in_band=True
            )
            if not cores_after:
                cores_after = [
                    i
                    for i in self._accumulated_core_indices(
                        require_in_band=False
                    )
                    if (d := self._abs_dev(i)) is not None
                    and d <= float(self._core_half_width) + 1e-15
                ]
            self._live_core_indices = list(cores_after)
            brackets_after = self._straddle_brackets(pool=None)
            n_raw_after = len(brackets_after)
            filt_after = self._nearest_straddle_brackets(
                max_distance=2.5 * float(self._radius)
            )
            n_br_after = len(filt_after)
            nearest_after = filt_after
            local_after = filt_after
            virtual_after = self._bracket_root_anchors(local_after)
            self._bracket_support_indices = sorted(
                {
                    int(index)
                    for edge in local_after
                    for index in (edge[0], edge[1])
                }
            )
            self._virtual_core_anchors = [anchor.tolist() for anchor in virtual_after]
            n_active_virtual_after = self._refresh_virtual_cores(
                virtual_after, cores_after
            )
            coverage_after = self._core_coverage_ok(cores_after)
            cls_after = self._classify_bridge_edges(
                filt_after, cores_after, r_window=self._radius
            )
            endpoint_after = self._actionable_endpoint_edges(
                filt_after,
                cores_after,
                r_sep=r_sep,
            )
            n_active_after = self._refresh_active_endpoints(cores_after)
            best_dev_after = self._abs_dev(best_a) if best_a is not None else None
            best_improved = (
                best_dev_before is not None
                and best_dev_after is not None
                and best_dev_after
                < best_dev_before - max(1e-12, 0.02 * self._core_half_width)
            )
            core_bounds_after = self._core_bounds(cores_after)
            extent_improved = self._bounds_expanded(
                core_bounds_before,
                core_bounds_after,
                min_gain=0.25 * r_sep,
            )
            core_gap_after = self._max_core_nn_gap(cores_after)
            gap_improved = (
                core_gap_before is not None
                and core_gap_after is not None
                and core_gap_after < core_gap_before - 0.25 * r_sep
            )
            structural_progress = self._fill_made_structural_progress(
                core_gain=core_gain,
                coverage_before=coverage_ok,
                best_improved=best_improved,
                extent_improved=extent_improved,
                gap_improved=gap_improved,
            )
            if structural_progress:
                self._quiet_fill_passes = 0
            else:
                self._quiet_fill_passes += 1
            scale_quiet = self._quiet_fill_passes >= self._quiet_fill_limit
            fill_after = self._fill_needed(
                cores=cores_after,
                n_local_brackets=n_br_after,
                coverage_ok=coverage_after,
                densify_new=n_densify,
                densify_attempted=densify_attempted,
                n_bridge_A=len(cls_after["A"]),
                n_bridge_B=len(cls_after["B"]),
                n_bridge_C=len(cls_after["C"]),
                n_endpoint_edges=len(endpoint_after),
                n_active_endpoints=n_active_after,
                n_active_virtual_cores=n_active_virtual_after,
                scale_quiet=scale_quiet,
            )
            band_ok_after = self._band_width_ok(t_min_a, t_max_a)
            self._logger.debug(
                "AdaptiveBridson post-feedback (full control cloud)\n"
                "    Delaunay / nearest brackets  -> %d / %d\n"
                "    bracket support / roots      -> %d / %d (active=%d, guided=%d, retrying=%d, nearest-pairs=%d)\n"
                "    t_min / t_max / band_ok      -> %s / %s / %s\n"
                "    fill needed                  -> %s\n"
                "    cores / frontier / all       -> %d / %d / %d\n"
                "    core gain / structural       -> %d / %s\n"
                "    extent / gap / best progress -> %s / %s / %s\n"
                "    endpoints straddle / geom    -> %d / %d\n"
                "    quiet progress               -> %d of %d\n"
                "    gen / fill_pass / r_g        -> %d / %d / %g",
                n_raw_after,
                n_br_after,
                len(self._bracket_support_indices),
                len(virtual_after),
                n_active_virtual_after,
                sum(
                    state.get("best_probe_dev") is not None
                    for state in self._active_virtual_cores.values()
                ),
                sum(
                    int(state.get("misses", 0)) > 0
                    for state in self._active_virtual_cores.values()
                ),
                len(nearest_after),
                f"{t_min_a:g}" if t_min_a is not None else "n/a",
                f"{t_max_a:g}" if t_max_a is not None else "n/a",
                "yes" if band_ok_after else "no",
                "yes" if fill_after else "no",
                len(cores_after),
                len(_fr),
                len(self._points),
                core_gain,
                "yes" if structural_progress else "no",
                "yes" if extent_improved else "no",
                "yes" if gap_improved else "no",
                "yes" if best_improved else "no",
                len(endpoint_after),
                n_active_after,
                self._quiet_fill_passes,
                self._quiet_fill_limit,
                int(self._generation),
                int(self._fill_pass),
                self._radius,
            )
            self._maybe_shrink_outer(
                n_core=len(cores_after),
                open_brackets=n_br_after,
                raw_brackets=n_raw_after,
            )
            decision = self._post_fill_decision(
                t_min=t_min_a,
                t_max=t_max_a,
                fill_needed=fill_after,
                cores=cores_after,
            )
            # If fill still needed → next loop continues 补点 at same gen.
            if decision.action is LoopAction.CONTINUE and fill_after:
                self._logger.debug(
                    "AdaptiveBridson continue fill at gen=%d r_g=%g "
                    "(fill_pass=%d, band not judged for next gen yet)",
                    int(self._generation),
                    self._radius,
                    int(self._fill_pass),
                )
                self.checkpoint_at_barrier(
                    reason=f"als_gen{int(self._generation)}_fill{int(self._fill_pass)}"
                )
                continue
            # Fill done: if band thin → next loop will exit converged; else next gen.
            if decision.action is LoopAction.ADVANCE:
                self._advance_generation(cores_after, t_min=t_min_a, t_max=t_max_a)
            else:
                self._logger.debug(
                    "AdaptiveBridson fill done at gen=%d r_g=%g; "
                    "band_ok=%s (next loop exits if thin)",
                    int(self._generation),
                    self._radius,
                    band_ok_after,
                )
                self.checkpoint_at_barrier(
                    reason=f"als_gen{int(self._generation)}_fill{int(self._fill_pass)}"
                )

    def _pre_fill_decision(
        self,
        *,
        best: int | None,
        t_min: float | None,
        t_max: float | None,
        fill_needed: bool,
    ) -> LoopDecision:
        """Decide whether a scale can enter the same-radius fill phase."""
        if best is None:
            return LoopDecision.stop("level-set not present in domain")
        if self._contour_converged(
            t_min=t_min, t_max=t_max, fill_needed=fill_needed
        ):
            return LoopDecision.stop("converged")
        if not fill_needed and self._at_min_radius():
            return LoopDecision.stop("min_radius")
        if int(self._generation) >= int(self._max_generations):
            return LoopDecision.stop("max_generations")
        if self._unique_point_count() >= self._max_points:
            return LoopDecision.stop("max_points")
        return LoopDecision.continue_fill()

    def _prepare_scale_state(self) -> ScaleState:
        """Classify the live cloud and retain only scale-local core sites."""
        cores, frontier, outer, best = self._classify_bands()
        best_in_radius, t_min, t_max = self._radius_t_band()
        if best is None:
            best = best_in_radius
        self._t_min, self._t_max = t_min, t_max
        self._register_cores(cores)
        cores = self._accumulated_core_indices(
            t_min=t_min, t_max=t_max, require_in_band=True
        )
        if not cores:
            # Do not let a temporarily empty radius band discard valid cores.
            cores, frontier, outer, _ = self._classify_bands()
            self._register_cores(cores)
            cores = [
                index
                for index in self._accumulated_core_indices(
                    t_min=None, t_max=None, require_in_band=False
                )
                if (deviation := self._abs_dev(index)) is not None
                and deviation <= float(self._core_half_width) + 1e-15
            ]
        self._live_core_indices = list(cores)
        self._outer_indices = list(outer)
        self._frontier_indices = list(frontier)
        return ScaleState(
            cores=list(cores),
            frontier=list(frontier),
            outer=list(outer),
            best=best,
            t_min=t_min,
            t_max=t_max,
        )

    def _post_fill_decision(
        self,
        *,
        t_min: float | None,
        t_max: float | None,
        fill_needed: bool,
        cores: Sequence[int],
    ) -> LoopDecision:
        """Choose same-scale fill or the next radius after a feedback barrier."""
        if fill_needed:
            return LoopDecision.continue_fill()
        if self._should_advance_generation(
            t_min=t_min,
            t_max=t_max,
            fill_needed=False,
            cores=cores,
        ):
            return LoopDecision.advance()
        return LoopDecision.continue_fill()

    def _advance_generation(
        self,
        cores: Sequence[int],
        *,
        t_min: float | None,
        t_max: float | None,
    ) -> None:
        """Commit one radius shrink after an explicit ``ADVANCE`` decision."""
        self._solidify_cores(cores, reason=f"pre_shrink_r_g={self._radius:g}")
        old_r = self._radius
        self._radius = self._next_radius()
        self._n_radius_refines += 1
        self._generation = int(self._n_radius_refines)
        self._fill_pass = 0
        self._quiet_fill_passes = 0
        self._reset_endpoint_state()
        self._last_densify_new = 0
        self._last_densify_attempted = False
        threshold = (
            float(self._threshold)
            if self._threshold is not None
            else float(self._core_half_width)
        )
        self._logger.info(
            "AdaptiveBridson gen=%d: r_g %g→%g, band=%s..%s "
            "(Δ=%s, target<%g), cores=%d",
            int(self._generation),
            old_r,
            self._radius,
            f"{t_min:g}" if t_min is not None else "n/a",
            f"{t_max:g}" if t_max is not None else "n/a",
            (
                f"{float(t_max) - float(t_min):g}"
                if t_min is not None and t_max is not None
                else "n/a"
            ),
            threshold,
            len(cores),
        )
        self.checkpoint_at_barrier(reason=f"als_generation_{int(self._generation)}")

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
        state = self._feedback_export_state()
        state.update(
            {
                "points": [asdict(p) for p in self._points],
                "dim": self._dim,
                "accepted_index": self._accepted_index,
                "converged": self._converged,
                "levelset": self._levelset,
                "failed_regions": list(self._failed_regions),
                "stop_reason": self._stop_reason,
                "uuid_to_index": dict(self._uuid_to_index),
                "radius": self._radius,
                "min_radius": self._min_radius,
                "t_min": self._t_min,
                "t_max": self._t_max,
                "live_core_indices": list(self._live_core_indices),
                "outer_half_width": self._outer_half_width,
                "core_half_width": self._core_half_width,
                "final_half_width": self._final_half_width,
                "final_half_width_configured": self._final_half_width_configured,
                "open_brackets": self._open_brackets,
                "bracket_support_indices": list(self._bracket_support_indices),
                "virtual_core_anchors": list(self._virtual_core_anchors),
                "solidified_core_uuids": sorted(self._solidified_core_uuids),
                "accumulated_core_uuids": sorted(self._accumulated_core_uuids),
                "n_radius_refines": self._n_radius_refines,
                "fill_pass": self._fill_pass,
                "max_fill_passes": self._max_fill_passes,
                "submitted_u_keys": [
                    list(key) for key in sorted(self._submitted_u_keys)
                ],
                "quiet_fill_passes": self._quiet_fill_passes,
                "quiet_fill_limit": self._quiet_fill_limit,
                "active_endpoints": {
                    key: dict(value)
                    for key, value in sorted(self._active_endpoints.items())
                },
                "retired_endpoint_ids": sorted(self._retired_endpoint_ids),
                "active_virtual_cores": {
                    key: dict(value)
                    for key, value in sorted(self._active_virtual_cores.items())
                },
                "retired_virtual_core_ids": sorted(
                    self._retired_virtual_core_ids
                ),
                "endpoint_stall_limit": self._endpoint_stall_limit,
                "endpoint_omni_probes": self._endpoint_omni_probes,
                "algorithm": "outer_core_root_correction_bridson",
            }
        )
        return state

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
        self._accepted_index = int(state.get("accepted_index", 0) or 0)
        self._converged = bool(state.get("converged", False))
        self._levelset = state.get("levelset")
        self._failed_regions = list(state.get("failed_regions") or [])
        self._stop_reason = str(state.get("stop_reason") or "")
        self._uuid_to_index = {
            str(k): int(v) for k, v in dict(state.get("uuid_to_index") or {}).items()
        }
        if not self._uuid_to_index:
            self._uuid_to_index = {p.uuid: i for i, p in enumerate(self._points) if p.uuid}
        if state.get("min_radius") is not None:
            self._min_radius = min(
                float(state["min_radius"]), self._initial_radius
            )
        restored_radius = float(
            state.get("radius", self._initial_radius) or self._initial_radius
        )
        self._radius = max(self._min_radius, restored_radius)
        self._t_min = state.get("t_min")
        self._t_max = state.get("t_max")
        if self._t_min is not None:
            self._t_min = float(self._t_min)
        if self._t_max is not None:
            self._t_max = float(self._t_max)
        self._live_core_indices = [
            int(i) for i in (state.get("live_core_indices") or [])
        ]
        if state.get("outer_half_width") is not None:
            self._outer_half_width = float(state["outer_half_width"])
        if state.get("core_half_width") is not None:
            self._core_half_width = float(state["core_half_width"])
        if "final_half_width_configured" in state:
            self._final_half_width_configured = bool(
                state.get("final_half_width_configured")
            )
        if state.get("final_half_width") is not None:
            self._final_half_width = float(state["final_half_width"])
        elif not self._final_half_width_configured:
            self._final_half_width = None
        self._open_brackets = int(state.get("open_brackets", 0) or 0)
        self._bracket_support_indices = [
            int(index) for index in (state.get("bracket_support_indices") or [])
        ]
        self._virtual_core_anchors = [
            list(map(float, anchor))
            for anchor in (state.get("virtual_core_anchors") or [])
            if isinstance(anchor, (list, tuple))
        ]
        self._solidified_core_uuids = {
            str(u) for u in (state.get("solidified_core_uuids") or [])
        }
        self._accumulated_core_uuids = {
            str(u) for u in (state.get("accumulated_core_uuids") or [])
        }
        # Legacy checkpoints: seed accumulated cores from live indices if empty.
        if not self._accumulated_core_uuids and self._live_core_indices:
            for i in self._live_core_indices:
                if 0 <= int(i) < len(self._points) and self._points[int(i)].uuid:
                    self._accumulated_core_uuids.add(
                        str(self._points[int(i)].uuid)
                    )
        self._n_radius_refines = int(state.get("n_radius_refines", 0) or 0)
        self._fill_pass = int(state.get("fill_pass", 0) or 0)
        if state.get("max_fill_passes") is not None:
            self._max_fill_passes = max(1, int(state["max_fill_passes"]))
        raw_u_keys = state.get("submitted_u_keys") or []
        self._submitted_u_keys = {
            self._u_key(item)
            for item in raw_u_keys
            if isinstance(item, (list, tuple))
        }
        if not self._submitted_u_keys:
            self._submitted_u_keys = {self._u_key(p.u) for p in self._points}
        self._quiet_fill_passes = max(
            0, int(state.get("quiet_fill_passes", 0) or 0)
        )
        if state.get("quiet_fill_limit") is not None:
            self._quiet_fill_limit = max(1, int(state["quiet_fill_limit"]))
        raw_endpoints = state.get("active_endpoints") or {}
        self._active_endpoints = {
            str(key): dict(value)
            for key, value in dict(raw_endpoints).items()
            if isinstance(value, Mapping)
        }
        self._retired_endpoint_ids = {
            str(value) for value in (state.get("retired_endpoint_ids") or [])
        }
        raw_virtual = state.get("active_virtual_cores") or {}
        self._active_virtual_cores = {
            str(key): dict(value)
            for key, value in dict(raw_virtual).items()
            if isinstance(value, Mapping)
        }
        self._retired_virtual_core_ids = {
            str(value)
            for value in (state.get("retired_virtual_core_ids") or [])
        }
        if state.get("endpoint_stall_limit") is not None:
            self._endpoint_stall_limit = max(
                3, int(state["endpoint_stall_limit"])
            )
        if state.get("endpoint_omni_probes") is not None:
            self._endpoint_omni_probes = max(
                2 * self._dim, int(state["endpoint_omni_probes"])
            )
        # Keep generation aligned with shrink count after resume.
        self._generation = int(self._n_radius_refines)
        self._feedback_import_state(state)
        if self._target_expression:
            self._compile_target()
        self._graph = self._make_graph()
        self._logger.debug(
            "AdaptiveBridson restored runtime state\n"
            "    points / generation / pending -> %d / %s / %d\n"
            "    converged / stop_reason       -> %s / %s\n"
            "    radius / band                 -> %g / %s",
            len(self._points),
            getattr(self, "_generation", "?"),
            len(getattr(self, "_pending_uuids", ()) or ()),
            self._converged,
            self._stop_reason or "(none)",
            self._radius,
            (
                f"{self._t_min:g}..{self._t_max:g}"
                if self._t_min is not None and self._t_max is not None
                else "n/a"
            ),
        )


__all__ = [
    "AdaptiveBridsonSampler",
    "DelaunayGraph",
    "KNNGraph",
    "LevelSetPoint",
    "NeighborGraph",
    "clip_unit_cube",
    "sample_ball",
]
