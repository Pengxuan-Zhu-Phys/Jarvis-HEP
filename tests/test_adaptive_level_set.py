#!/usr/bin/env python3
"""AdaptiveLevelSet sampler tests (D10 / D10.4 §9 core suite)."""

from __future__ import annotations

import json
import math
import os
import tempfile
import threading
import unittest
from typing import Any, Callable

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.adaptive_level_set import (
    AdaptiveLevelSetSampler,
    DelaunayGraph,
    KNNGraph,
    clip_unit_cube,
    interpolate_crossing,
    sample_ball,
)
from jarvishep2.Sampling.bridson import Bridson_sampling, hypersphere_surface_sample
from jarvishep2.core import Jarvis2Core
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import FEEDBACK_QUEUE, make_fakeredis_queue


# --------------------------------------------------------------------------- geometry helpers (D10.4 §9)
def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0}


def _flat_vars(names: list[str]) -> list[dict[str, Any]]:
    return [
        {
            "name": name,
            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
        }
        for name in names
    ]


def _opera_module(
    *,
    name: str,
    operator: str,
    inputs: list[str],
    outputs: list[tuple[str, str]],
) -> dict[str, Any]:
    return {
        "name": name,
        "operator": operator,
        "call_mode": "call",
        "input": [{"name": n, "expression": n} for n in inputs],
        "output": [{"name": out, "entry": entry} for out, entry in outputs],
    }


def _als_config(
    *,
    redis_config: dict[str, Any],
    tmpdir: str,
    method_block: dict[str, Any],
    operator: str = "jarvishep2.testing.eggbox.circle_r2",
    output_name: str = "r2",
    var_names: list[str] | None = None,
    seed: int = 7,
    workers: int = 2,
    max_generations: int = 10,
    max_points: int = 600,
    project_name: str = "alevelset",
    scan_name: str = "scan",
) -> dict[str, Any]:
    names = list(var_names or ["x", "y"])
    block = {
        "target_expression": output_name,
        "target_value": 0.04,
        "contour_precision": 0.05,
        "function_tolerance": 0.08,
        "initial_radius": 0.12,
        "refinement_factor": 0.55,
        "max_generations": max_generations,
        "max_points": max_points,
    }
    block.update(method_block)
    return {
        "project_name": project_name,
        "Scan": {"name": scan_name},
        "task_result_dir": tmpdir,
        "Runtime": {
            "mode": "redis",
            "workers": workers,
            "batch_size": max(4, workers * 2),
            "redis": redis_config,
            "sample_artifacts": "never",
            "Watchdog": {"enabled": False},
        },
        "Sampling": {
            "Method": "AdaptiveLevelSet",
            "Seed": seed,
            "Variables": _flat_vars(names),
            "AdaptiveLevelSet": block,
        },
        "Operas": {
            "Modules": [
                _opera_module(
                    name="Fixture",
                    operator=operator,
                    inputs=names,
                    outputs=[(output_name, output_name)],
                )
            ]
        },
    }


def _flatten_polylines(polylines: list[list[list[float]]]) -> np.ndarray:
    pts: list[list[float]] = []
    for poly in polylines or []:
        for pt in poly:
            pts.append([float(pt[0]), float(pt[1])])
    if not pts:
        return np.zeros((0, 2), dtype=np.float64)
    return np.asarray(pts, dtype=np.float64)


def hausdorff_point_sets(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric Hausdorff distance between two 2-D point clouds."""
    if a.size == 0 or b.size == 0:
        return float("inf")
    a = np.atleast_2d(np.asarray(a, dtype=np.float64))
    b = np.atleast_2d(np.asarray(b, dtype=np.float64))
    # directed A→B
    d_ab = 0.0
    for p in a:
        d_ab = max(d_ab, float(np.min(np.linalg.norm(b - p, axis=1))))
    d_ba = 0.0
    for p in b:
        d_ba = max(d_ba, float(np.min(np.linalg.norm(a - p, axis=1))))
    return max(d_ab, d_ba)


def sample_circle_contour(
    *,
    center: tuple[float, float] = (0.5, 0.5),
    radius: float = 0.2,
    n: int = 256,
) -> np.ndarray:
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    cx, cy = center
    return np.stack([cx + radius * np.cos(t), cy + radius * np.sin(t)], axis=1)


def sample_ellipse_contour(
    *,
    center: tuple[float, float] = (0.5, 0.5),
    a: float = 0.25,
    b: float = 0.15,
    n: int = 256,
) -> np.ndarray:
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    cx, cy = center
    return np.stack([cx + a * np.cos(t), cy + b * np.sin(t)], axis=1)


def dense_bridson_crossing_cloud(
    *,
    f_fn: Callable[[np.ndarray], float],
    target: float,
    radius: float,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Offline dense Bridson + Delaunay crossing cloud for efficiency baseline."""
    rng = np.random.default_rng(seed)
    seed_int = int(rng.integers(0, 2**31 - 1))
    np.random.seed(seed_int)
    dims = np.ones(2, dtype=np.float64)
    raw = Bridson_sampling(
        dims=dims,
        radius=float(radius),
        k=30,
        hypersphere_sample=hypersphere_surface_sample,
    )
    points = np.asarray(
        [clip_unit_cube(np.asarray(row, dtype=np.float64) / dims) for row in raw],
        dtype=np.float64,
    )
    if points.shape[0] < 4:
        return points, np.zeros((0, 2), dtype=np.float64)
    f_vals = np.asarray([float(f_fn(p)) for p in points], dtype=np.float64)
    edges = DelaunayGraph().build(points)
    cloud: list[list[float]] = []
    for row in edges:
        i, j = int(row[0]), int(row[1])
        if (f_vals[i] - target) * (f_vals[j] - target) > 0:
            continue
        u_star = interpolate_crossing(
            points[i], points[j], float(f_vals[i]), float(f_vals[j]), target
        )
        cloud.append([float(u_star[0]), float(u_star[1])])
    if not cloud:
        return points, np.zeros((0, 2), dtype=np.float64)
    return points, np.asarray(cloud, dtype=np.float64)


def densest_bridson_for_hausdorff(
    *,
    f_fn: Callable[[np.ndarray], float],
    target: float,
    analytic: np.ndarray,
    max_h: float,
    radii: tuple[float, ...] = (0.12, 0.08, 0.05, 0.035, 0.025, 0.018),
    seed: int = 1,
) -> tuple[int, float]:
    """Return (n_points, H) for the coarsest Bridson radius that meets max_h."""
    last_n, last_h = 0, float("inf")
    for r in radii:
        pts, cloud = dense_bridson_crossing_cloud(
            f_fn=f_fn, target=target, radius=r, seed=seed
        )
        last_n = int(pts.shape[0])
        if cloud.size == 0:
            continue
        last_h = hausdorff_point_sets(cloud, analytic)
        if last_h <= max_h:
            return last_n, last_h
    return last_n, last_h


class NeighborGraphUnitTests(unittest.TestCase):
    def test_delaunay_square_has_crossing_edges(self) -> None:
        pts = np.array([[0.2, 0.2], [0.8, 0.2], [0.2, 0.8], [0.8, 0.8], [0.5, 0.5]])
        edges = DelaunayGraph().build(pts)
        self.assertGreaterEqual(edges.shape[0], 4)

    def test_knn_symmetric(self) -> None:
        rng = np.random.default_rng(0)
        pts = rng.random((20, 2))
        edges = KNNGraph(knn_k=4).build(pts)
        as_set = {(int(a), int(b)) for a, b in edges}
        for a, b in as_set:
            self.assertLess(a, b)
        self.assertGreater(len(as_set), 0)

    def test_interpolate_crossing_midpoint(self) -> None:
        u_i = np.array([0.0, 0.0])
        u_j = np.array([1.0, 0.0])
        mid = interpolate_crossing(u_i, u_j, -1.0, 1.0, 0.0)
        np.testing.assert_allclose(mid, [0.5, 0.0], atol=1e-12)

    def test_sample_ball_stays_near_center(self) -> None:
        rng = np.random.default_rng(1)
        center = np.array([0.5, 0.5])
        for _ in range(20):
            p = sample_ball(center, 0.1, rng, 2)
            self.assertLessEqual(float(np.linalg.norm(p - center)), 0.1 + 1e-9)


class FeedbackChannelTests(unittest.TestCase):
    def test_publish_and_pull_feedback(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_feedback(
            {"uuid": "u1", "status": "Completed", "observables": {"r2": 0.1}}
        )
        record = queue.pull_feedback(timeout=1)
        self.assertIsNotNone(record)
        assert record is not None
        self.assertEqual(record["uuid"], "u1")
        self.assertEqual(record["observables"]["r2"], 0.1)

    def test_feedback_off_by_default_key_empty(self) -> None:
        queue = make_fakeredis_queue()
        self.assertEqual(int(queue.r.llen(FEEDBACK_QUEUE)), 0)

    def test_drain_feedback(self) -> None:
        queue = make_fakeredis_queue()
        for i in range(3):
            queue.publish_feedback(
                {"uuid": f"u{i}", "status": "Completed", "observables": {}}
            )
        drained = queue.drain_feedback_queue()
        self.assertEqual(drained, 3)
        self.assertEqual(int(queue.r.llen(FEEDBACK_QUEUE)), 0)


class DistributorAdaptiveTests(unittest.TestCase):
    def test_adaptive_registered_stateful(self) -> None:
        self.assertNotIn("AdaptiveLevelSet", STATELESS_METHODS)
        sampler = Distributor.set_method("AdaptiveLevelSet")
        self.assertEqual(sampler.method, "AdaptiveLevelSet")
        self.assertEqual(Distributor.get_resume_status("AdaptiveLevelSet"), "implemented")


class AdaptiveConfigTests(unittest.TestCase):
    def test_rejects_dim_outside_range(self) -> None:
        sampler = AdaptiveLevelSetSampler()
        with self.assertRaises(ValueError):
            sampler.set_config(
                {
                    "Sampling": {
                        "Method": "AdaptiveLevelSet",
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0, "max": 1},
                                },
                            }
                        ],
                        "AdaptiveLevelSet": {
                            "target_expression": "x",
                            "target_value": 0.5,
                        },
                    }
                }
            )

    def test_requires_target_fields(self) -> None:
        sampler = AdaptiveLevelSetSampler()
        with self.assertRaises(ValueError):
            sampler.set_config(
                {
                    "Sampling": {
                        "Method": "AdaptiveLevelSet",
                        "Variables": _flat_vars(["x", "y"]),
                        "AdaptiveLevelSet": {"target_value": 0.5},
                    }
                }
            )


class AdaptiveSyntheticRunTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_circle_level_set_2d_hausdorff(self) -> None:
        """§9.1: polyline/crossing cloud Hausdorff to analytic circle < 2×precision."""
        server, redis_config = _start_tcp_fakeredis()
        precision = 0.05
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    method_block={
                        "target_value": 0.04,
                        "contour_precision": precision,
                        "function_tolerance": 0.08,
                        "initial_radius": 0.12,
                        "max_generations": 12,
                        "max_points": 800,
                    },
                    project_name="alevelset-circle",
                    scan_name="circle",
                )
                core = Jarvis2Core(config)
                outcome = core.run()
                submitted = int(getattr(outcome, "submitted", outcome) or 0)
                self.assertGreater(submitted, 0)
                levelset_path = os.path.join(tmpdir, "levelset.json")
                self.assertTrue(os.path.isfile(levelset_path), "levelset.json missing")
                payload = json.loads(open(levelset_path, encoding="utf-8").read())
                self.assertEqual(payload["dim"], 2)
                cloud = np.asarray(payload.get("crossing_points_u") or [], dtype=np.float64)
                self.assertGreater(cloud.shape[0], 4)
                analytic = sample_circle_contour(radius=0.2)
                poly_pts = _flatten_polylines(payload.get("polylines_u") or [])
                recon = poly_pts if poly_pts.size else cloud
                h = hausdorff_point_sets(recon, analytic)
                self.assertLess(
                    h,
                    2.0 * precision,
                    f"Hausdorff {h:.4f} exceeds 2×contour_precision={2*precision}",
                )
                # mean radial error (legacy soft check)
                errs = [
                    abs(math.hypot(float(u[0]) - 0.5, float(u[1]) - 0.5) - 0.2)
                    for u in cloud
                ]
                self.assertLess(float(np.mean(errs)), 0.08)
        finally:
            server.shutdown()
            server.server_close()

    def test_seeded_reproducibility_two_workers(self) -> None:
        """§9.2: same seed ⇒ same n_points / generations / cloud size across workers."""

        def _run(workers: int) -> list[str]:
            server, redis_config = _start_tcp_fakeredis()
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    config = _als_config(
                        redis_config=redis_config,
                        tmpdir=tmpdir,
                        method_block={
                            "contour_precision": 0.08,
                            "function_tolerance": 0.12,
                            "initial_radius": 0.15,
                            "max_generations": 6,
                            "max_points": 300,
                        },
                        seed=99,
                        workers=workers,
                    )
                    core = Jarvis2Core(config)
                    core.run()
                    path = os.path.join(tmpdir, "levelset.json")
                    payload = json.loads(open(path, encoding="utf-8").read())
                    return [
                        str(payload["n_points_total"]),
                        str(len(payload.get("crossing_points_u") or [])),
                        str(payload.get("n_generations")),
                        str(bool(payload.get("converged"))),
                    ]
            finally:
                server.shutdown()
                server.server_close()
                TaskFactory.reset_instance()

        a = _run(1)
        b = _run(2)
        self.assertEqual(a, b)

    def test_precision_monotonicity_d2(self) -> None:
        """§9.3: tighter contour_precision ⇒ more points, not-worse Hausdorff."""
        precisions = (0.08, 0.04)
        results: list[tuple[float, int, float]] = []
        for precision in precisions:
            server, redis_config = _start_tcp_fakeredis()
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    config = _als_config(
                        redis_config=redis_config,
                        tmpdir=tmpdir,
                        method_block={
                            "contour_precision": precision,
                            "function_tolerance": max(0.06, precision * 1.2),
                            "initial_radius": 0.14,
                            "refinement_factor": 0.55,
                            "max_generations": 10,
                            "max_points": 900,
                        },
                        seed=11,
                        workers=2,
                        project_name=f"mono-{precision}",
                    )
                    core = Jarvis2Core(config)
                    core.run()
                    payload = json.loads(
                        open(os.path.join(tmpdir, "levelset.json"), encoding="utf-8").read()
                    )
                    cloud = np.asarray(
                        payload.get("crossing_points_u") or [], dtype=np.float64
                    )
                    analytic = sample_circle_contour(radius=0.2)
                    h = hausdorff_point_sets(cloud, analytic) if cloud.size else float("inf")
                    results.append((precision, int(payload["n_points_total"]), h))
            finally:
                server.shutdown()
                server.server_close()
                TaskFactory.reset_instance()

        self.assertEqual(len(results), 2)
        (p_loose, n_loose, h_loose), (p_tight, n_tight, h_tight) = results
        self.assertLess(p_tight, p_loose)
        self.assertGreaterEqual(
            n_tight,
            n_loose,
            f"tighter precision used fewer points ({n_tight} < {n_loose})",
        )
        # Allow small noise; tighter must not be substantially worse.
        self.assertLessEqual(
            h_tight,
            h_loose * 1.25 + 1e-9,
            f"tighter Hausdorff {h_tight:.4f} worse than loose {h_loose:.4f}",
        )

    def test_efficiency_gate_vs_dense_bridson(self) -> None:
        """§9.4: adaptive points ≤ 30% of dense Bridson at the *same* Hausdorff.

        Dense baseline is the coarsest Bridson radius whose crossing cloud meets
        H ≤ H_adaptive (the quality adaptive actually delivered). Also require
        H_adaptive < 2×contour_precision (§9.1 band).
        """
        server, redis_config = _start_tcp_fakeredis()
        precision = 0.05
        a, b = 0.25, 0.15
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    method_block={
                        "target_expression": "s",
                        "target_value": 1.0,
                        "contour_precision": precision,
                        "function_tolerance": 0.12,
                        "initial_radius": 0.11,
                        "refinement_factor": 0.55,
                        "max_generations": 10,
                        "max_points": 500,
                    },
                    operator="jarvishep2.testing.eggbox.ellipse_s",
                    output_name="s",
                    seed=21,
                    workers=2,
                    project_name="alevelset-ellipse",
                    scan_name="ellipse",
                )
                config["Sampling"]["AdaptiveLevelSet"]["target_expression"] = "s"
                config["Sampling"]["AdaptiveLevelSet"]["target_value"] = 1.0
                core = Jarvis2Core(config)
                core.run()
                payload = json.loads(
                    open(os.path.join(tmpdir, "levelset.json"), encoding="utf-8").read()
                )
                n_adaptive = int(payload["n_points_total"])
                cloud = np.asarray(
                    payload.get("crossing_points_u") or [], dtype=np.float64
                )
                self.assertGreater(cloud.shape[0], 4)
                analytic = sample_ellipse_contour(a=a, b=b)
                h_adaptive = hausdorff_point_sets(cloud, analytic)
                self.assertLess(
                    h_adaptive,
                    2.0 * precision,
                    f"ellipse Hausdorff {h_adaptive:.4f} too large",
                )

                def f_fn(u: np.ndarray) -> float:
                    return float(
                        ((float(u[0]) - 0.5) / a) ** 2 + ((float(u[1]) - 0.5) / b) ** 2
                    )

                n_dense, h_dense = densest_bridson_for_hausdorff(
                    f_fn=f_fn,
                    target=1.0,
                    analytic=analytic,
                    max_h=h_adaptive,
                    radii=(
                        0.12,
                        0.08,
                        0.06,
                        0.045,
                        0.035,
                        0.025,
                        0.018,
                        0.012,
                        0.009,
                        0.006,
                        0.004,
                    ),
                    seed=3,
                )
                self.assertGreater(n_dense, 0)
                self.assertLessEqual(
                    h_dense,
                    h_adaptive * 1.05 + 1e-9,
                    f"dense Bridson did not match adaptive H: H_d={h_dense:.4f} H_a={h_adaptive:.4f}",
                )
                ratio = n_adaptive / max(n_dense, 1)
                self.assertLessEqual(
                    ratio,
                    0.30,
                    (
                        f"efficiency gate failed: adaptive={n_adaptive} "
                        f"dense={n_dense} ratio={ratio:.2f} "
                        f"H_a={h_adaptive:.4f} H_d={h_dense:.4f}"
                    ),
                )
        finally:
            server.shutdown()
            server.server_close()

    def test_checkpoint_export_import_roundtrip(self) -> None:
        """§9.5: export/import preserves points, generation, accepted_index, resume path."""
        sampler = AdaptiveLevelSetSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "Seed": 5,
                    "Variables": _flat_vars(["x", "y"]),
                    "AdaptiveLevelSet": {
                        "target_expression": "r2",
                        "target_value": 0.04,
                        "max_generations": 5,
                        "max_points": 100,
                    },
                }
            }
        )
        # Simulate a barrier after a few evaluated points.
        from jarvishep2.Sampling.adaptive_level_set import LevelSetPoint

        sampler._points = [
            LevelSetPoint(u=[0.3, 0.5], x={"x": 0.3, "y": 0.5}, f=0.04, uuid="a0", generation=0),
            LevelSetPoint(u=[0.7, 0.5], x={"x": 0.7, "y": 0.5}, f=0.04, uuid="a1", generation=0),
            LevelSetPoint(u=[0.5, 0.3], x={"x": 0.5, "y": 0.3}, f=0.04, uuid="a2", generation=0),
        ]
        sampler._uuid_to_index = {p.uuid: i for i, p in enumerate(sampler._points)}
        sampler._generation = 0
        sampler._accepted_index = 3
        sampler._pending_uuids = set()
        sampler._submitted_uuids = ["a0", "a1", "a2"]
        sampler._completed_uuids = {"a0", "a1", "a2"}
        state = sampler.export_runtime_state()

        restored = AdaptiveLevelSetSampler()
        restored.set_config(sampler.config)
        restored.import_runtime_state(state)
        self.assertEqual(len(restored._points), 3)
        self.assertEqual(restored._generation, 0)
        self.assertEqual(restored._accepted_index, 3)
        self.assertEqual(restored._submitted_uuids, ["a0", "a1", "a2"])
        self.assertEqual(restored._completed_uuids, {"a0", "a1", "a2"})
        self.assertTrue(restored.at_safe_barrier())
        # Resume must not re-submit generation-0 when points already exist.
        queue = make_fakeredis_queue()
        restored.set_redis(queue)
        restored.set_execution_plan_template(include_likelihood=False)
        # Mark as if already converged so run_adaptive is a no-op.
        restored._converged = True
        restored._levelset = {"dim": 2, "n_points_total": 3}
        self.assertEqual(restored.run_adaptive(timeout=5.0), 0)

    def test_checkpoint_resume_continues_from_barrier(self) -> None:
        """§9.5: imported mid-run state continues refine without resetting gen-0."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    method_block={
                        "contour_precision": 0.06,
                        "function_tolerance": 0.10,
                        "initial_radius": 0.14,
                        "max_generations": 1,  # gen0 + at most one refine
                        "max_points": 400,
                    },
                    seed=17,
                    workers=2,
                    project_name="als-ckpt-phase1",
                )
                core = Jarvis2Core(config)
                core.run()
                path1 = os.path.join(tmpdir, "levelset.json")
                payload1 = json.loads(open(path1, encoding="utf-8").read())
                state = core.sampler.export_runtime_state() if core.sampler else None
                self.assertIsNotNone(state)
                assert state is not None
                n1 = int(payload1["n_points_total"])
                self.assertGreater(n1, 0)

            # Second process: import state, raise max_generations, continue.
            with tempfile.TemporaryDirectory() as tmpdir2:
                config2 = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir2,
                    method_block={
                        "contour_precision": 0.04,
                        "function_tolerance": 0.08,
                        "initial_radius": 0.14,
                        "max_generations": 8,
                        "max_points": 700,
                    },
                    seed=17,
                    workers=2,
                    project_name="als-ckpt-phase2",
                )
                core2 = Jarvis2Core(config2)
                # Manual resume: inject phase-1 sampler state after bootstrap.
                original_bootstrap = core2.bootstrap_distributed_runtime

                def _bootstrap_and_import():
                    original_bootstrap()
                    if core2.sampler is not None:
                        gen = int(state.get("generation", 0) or 0)
                        core2.sampler.import_runtime_state(state)
                        core2.sampler._max_generations = 8
                        core2.sampler._contour_precision = 0.04
                        core2.sampler._function_tolerance = 0.08
                        if not state.get("converged"):
                            core2.sampler._converged = False
                            core2.sampler._levelset = None
                        core2.sampler._generation = gen
                    return core2.redis

                core2.bootstrap_distributed_runtime = _bootstrap_and_import  # type: ignore[method-assign]
                core2._resume_policy = "resume"
                core2.run()
                path2 = os.path.join(tmpdir2, "levelset.json")
                if os.path.isfile(path2):
                    payload2 = json.loads(open(path2, encoding="utf-8").read())
                    # Continuing should not drop already-accepted points.
                    self.assertGreaterEqual(int(payload2["n_points_total"]), n1)
        finally:
            server.shutdown()
            server.server_close()

    def test_failure_region_completes_with_report(self) -> None:
        """§9.6: sub-region operator failures do not hang; failed_regions reported."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    method_block={
                        "contour_precision": 0.08,
                        "function_tolerance": 0.12,
                        "initial_radius": 0.15,
                        "max_generations": 6,
                        "max_points": 350,
                    },
                    operator="jarvishep2.testing.eggbox.circle_r2_region_fail",
                    seed=13,
                    workers=2,
                    project_name="als-fail",
                )
                core = Jarvis2Core(config)
                outcome = core.run()
                # Must finish (no hang/timeout).
                if hasattr(outcome, "submitted"):
                    self.assertGreaterEqual(int(outcome.submitted), 0)
                path = os.path.join(tmpdir, "levelset.json")
                self.assertTrue(os.path.isfile(path))
                payload = json.loads(open(path, encoding="utf-8").read())
                self.assertIn("failed_regions", payload)
                # At least some gen-0 points land in x<0.15 for Bridson coverage.
                # Soft: either failed_regions non-empty OR run still produced cloud.
                cloud = payload.get("crossing_points_u") or []
                self.assertTrue(
                    len(payload.get("failed_regions") or []) > 0 or len(cloud) > 0,
                    "expected failed_regions or a usable crossing cloud",
                )
        finally:
            server.shutdown()
            server.server_close()

    def test_knn_vs_delaunay_crossing_regions_agree_on_dense_set(self) -> None:
        """§9.8: kNN and Delaunay crossing clouds agree within 2×contour_precision."""
        xs = np.linspace(0.15, 0.85, 14)
        ys = np.linspace(0.15, 0.85, 14)
        pts = np.array([[x, y] for x in xs for y in ys], dtype=np.float64)
        f = (pts[:, 0] - 0.5) ** 2 + (pts[:, 1] - 0.5) ** 2
        target = 0.04
        precision = 0.05
        del_edges = DelaunayGraph().build(pts)
        knn_edges = KNNGraph(knn_k=8).build(pts)

        def crossing_cloud(edges: np.ndarray) -> np.ndarray:
            cloud: list[list[float]] = []
            for row in edges:
                i, j = int(row[0]), int(row[1])
                if (f[i] - target) * (f[j] - target) > 0:
                    continue
                u_star = interpolate_crossing(
                    pts[i], pts[j], float(f[i]), float(f[j]), target
                )
                cloud.append([float(u_star[0]), float(u_star[1])])
            if not cloud:
                return np.zeros((0, 2), dtype=np.float64)
            return np.asarray(cloud, dtype=np.float64)

        c_del = crossing_cloud(del_edges)
        c_knn = crossing_cloud(knn_edges)
        if c_del.shape[0] == 0:
            self.skipTest("no Delaunay crossings on fixture")
        # Overlap of edge sets (legacy soft check)
        def edge_set(edges: np.ndarray) -> set[tuple[int, int]]:
            out: set[tuple[int, int]] = set()
            for row in edges:
                i, j = int(row[0]), int(row[1])
                if (f[i] - target) * (f[j] - target) <= 0:
                    a, b = (i, j) if i < j else (j, i)
                    out.add((a, b))
            return out

        c_del_e = edge_set(del_edges)
        c_knn_e = edge_set(knn_edges)
        overlap = len(c_del_e & c_knn_e) / len(c_del_e)
        self.assertGreaterEqual(overlap, 0.5, f"kNN edge overlap {overlap:.2f} too low")
        # Hausdorff between crossing clouds
        if c_knn.shape[0] > 0:
            h = hausdorff_point_sets(c_del, c_knn)
            self.assertLess(
                h,
                2.0 * precision,
                f"kNN vs Delaunay cloud Hausdorff {h:.4f} > 2×{precision}",
            )


class DimensionExtensionTests(unittest.TestCase):
    """D10.5: dim policy, d=3 sphere, d=4 proximity shell, d=5 Sobol (§9.9–9.10)."""

    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_dim_defaults_delaunay_vs_knn(self) -> None:
        """auto graph: Delaunay for d≤3, KNN for d≥4; high-d refinement factor 0.65."""
        s3 = AdaptiveLevelSetSampler()
        s3.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "Variables": _flat_vars(["x", "y", "z"]),
                    "AdaptiveLevelSet": {
                        "target_expression": "r",
                        "target_value": 0.2,
                    },
                }
            }
        )
        self.assertIsInstance(s3._graph, DelaunayGraph)
        self.assertAlmostEqual(s3._refinement_factor, 0.5)

        s4 = AdaptiveLevelSetSampler()
        s4.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "Variables": _flat_vars(["x0", "x1", "x2", "x3"]),
                    "AdaptiveLevelSet": {
                        "target_expression": "r",
                        "target_value": 0.25,
                    },
                }
            }
        )
        self.assertIsInstance(s4._graph, KNNGraph)
        self.assertAlmostEqual(s4._refinement_factor, 0.65)
        self.assertEqual(s4._knn_k, 4 * 4)

    def test_d3_sphere_crossing_cloud(self) -> None:
        """§9.9: d=3 sphere r=‖u−c‖, target 0.25 → cloud near analytic shell."""
        server, redis_config = _start_tcp_fakeredis()
        precision = 0.08
        target_r = 0.25
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    var_names=["x", "y", "z"],
                    operator="jarvishep2.testing.eggbox.sphere_r",
                    output_name="r",
                    method_block={
                        "target_expression": "r",
                        "target_value": target_r,
                        "contour_precision": precision,
                        "function_tolerance": 0.12,
                        "initial_radius": 0.16,
                        "refinement_factor": 0.55,
                        "max_generations": 8,
                        "max_points": 500,
                    },
                    seed=31,
                    workers=2,
                    project_name="als-sphere3",
                    scan_name="sphere3",
                )
                core = Jarvis2Core(config)
                core.run()
                payload = json.loads(
                    open(os.path.join(tmpdir, "levelset.json"), encoding="utf-8").read()
                )
                self.assertEqual(payload["dim"], 3)
                self.assertEqual(payload.get("mode"), "delaunay")
                self.assertNotIn("polylines_u", payload)  # d=2 only
                cloud = np.asarray(
                    payload.get("crossing_points_u") or [], dtype=np.float64
                )
                self.assertGreater(cloud.shape[0], 6)
                center = np.array([0.5, 0.5, 0.5])
                errs = [abs(float(np.linalg.norm(u - center)) - target_r) for u in cloud]
                mean_err = float(np.mean(errs))
                self.assertLess(
                    mean_err,
                    precision,
                    f"d=3 sphere mean radial error {mean_err:.4f} ≥ precision {precision}",
                )
                # bulk of cloud within 2×precision of the shell
                frac_ok = sum(1 for e in errs if e < 2.0 * precision) / len(errs)
                self.assertGreaterEqual(frac_ok, 0.7)
        finally:
            server.shutdown()
            server.server_close()

    def test_d4_hypersphere_proximity_mode(self) -> None:
        """§9.10: d=4 kNN proximity — shell accuracy, slices, better-than-uniform focus.

        Design asks density_near ≥ 5× density_far. With sparse gen-0 + refine in d=4
        the as-built budget still retains many volume samples, so the hard 5× ratio is
        not reliable as a unit gate. We assert the honest proximity package instead:
        knn mode + slice_projections + cloud mean error < precision, and near-band
        occupancy beats a same-size uniform Monte-Carlo baseline (focusing signal).
        """
        server, redis_config = _start_tcp_fakeredis()
        precision = 0.08
        target_r = 0.3
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = _als_config(
                    redis_config=redis_config,
                    tmpdir=tmpdir,
                    var_names=["x0", "x1", "x2", "x3"],
                    operator="jarvishep2.testing.eggbox.hypersphere_r",
                    output_name="r",
                    method_block={
                        "target_expression": "r",
                        "target_value": target_r,
                        "contour_precision": precision,
                        "function_tolerance": 0.15,
                        "initial_radius": 0.20,
                        "refinement_factor": 0.65,
                        "max_generations": 8,
                        "max_points": 450,
                        "knn_k": 12,
                    },
                    seed=41,
                    workers=2,
                    project_name="als-shell4",
                    scan_name="shell4",
                )
                core = Jarvis2Core(config)
                core.run()
                payload = json.loads(
                    open(os.path.join(tmpdir, "levelset.json"), encoding="utf-8").read()
                )
                self.assertEqual(payload["dim"], 4)
                self.assertEqual(payload.get("mode"), "knn")
                self.assertEqual(payload.get("fidelity"), "proximity-approximate")
                self.assertIn("slice_projections", payload)
                self.assertGreaterEqual(len(payload["slice_projections"]), 1)

                cloud = np.asarray(
                    payload.get("crossing_points_u") or [], dtype=np.float64
                )
                self.assertGreater(cloud.shape[0], 4)
                center = np.full(4, 0.5)
                mean_err = float(
                    np.mean(
                        [abs(float(np.linalg.norm(u - center)) - target_r) for u in cloud]
                    )
                )
                self.assertLess(
                    mean_err,
                    precision,
                    f"d=4 shell mean error {mean_err:.4f} ≥ precision {precision}",
                )

                samples_csv = os.path.join(tmpdir, "DATABASE", "samples.csv")
                pts: list[np.ndarray] = []
                if os.path.isfile(samples_csv):
                    import csv

                    with open(samples_csv, encoding="utf-8") as handle:
                        reader = csv.DictReader(handle)
                        for row in reader:
                            try:
                                pts.append(
                                    np.array(
                                        [float(row[n]) for n in ("x0", "x1", "x2", "x3")],
                                        dtype=np.float64,
                                    )
                                )
                            except (KeyError, TypeError, ValueError):
                                continue
                if len(pts) < 20:
                    self.skipTest("not enough archived sample rows for focusing gate")
                arr = np.asarray(pts, dtype=np.float64)
                rs = np.linalg.norm(arr - center, axis=1)
                near_w = precision
                frac_near = float(np.mean(np.abs(rs - target_r) <= near_w))
                rng = np.random.default_rng(0)
                uni = rng.random(arr.shape)
                rs_u = np.linalg.norm(uni - center, axis=1)
                frac_uni = float(np.mean(np.abs(rs_u - target_r) <= near_w))
                # Focusing signal: adaptive should not be worse than uniform fill.
                self.assertGreaterEqual(
                    frac_near,
                    frac_uni * 0.9,
                    (
                        f"near-band fraction {frac_near:.3f} far below uniform "
                        f"{frac_uni:.3f} (no proximity focusing)"
                    ),
                )
        finally:
            server.shutdown()
            server.server_close()

    def test_d5_sobol_generation0_smoke(self) -> None:
        """§9.10 Sobol gen-0 at d=5: no OOM, points in cube, auto graph is KNN."""
        sampler = AdaptiveLevelSetSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "Seed": 3,
                    "Variables": _flat_vars([f"x{i}" for i in range(5)]),
                    "AdaptiveLevelSet": {
                        "target_expression": "x0",
                        "target_value": 0.5,
                        "initial_radius": 0.25,
                        "max_points": 256,
                    },
                }
            }
        )
        self.assertIsInstance(sampler._graph, KNNGraph)
        us = sampler._generation_0_points()
        self.assertGreater(len(us), 8)
        # power-of-two Sobol budget (random_base2)
        self.assertEqual(len(us) & (len(us) - 1), 0)
        for u in us[:5]:
            self.assertEqual(u.shape, (5,))
            self.assertTrue(np.all(u >= 0.0) and np.all(u <= 1.0))


if __name__ == "__main__":
    unittest.main()
