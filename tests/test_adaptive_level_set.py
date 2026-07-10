#!/usr/bin/env python3
"""AdaptiveLevelSet sampler tests (D10)."""

from __future__ import annotations

import json
import math
import os
import tempfile
import threading
import time
import unittest

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
from jarvishep2.core import Jarvis2Core
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import FEEDBACK_QUEUE, make_fakeredis_queue
def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0}


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
            queue.publish_feedback({"uuid": f"u{i}", "status": "Completed", "observables": {}})
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
                                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
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
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
                            },
                            {
                                "name": "y",
                                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
                            },
                        ],
                        "AdaptiveLevelSet": {"target_value": 0.5},
                    }
                }
            )


class AdaptiveSyntheticRunTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_circle_level_set_2d_converges(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                config = {
                    "project_name": "alevelset-circle",
                    "Scan": {"name": "circle"},
                    "task_result_dir": tmpdir,
                    "Runtime": {
                        "mode": "redis",
                        "workers": 2,
                        "batch_size": 8,
                        "redis": redis_config,
                        "sample_artifacts": "never",
                        "Watchdog": {"enabled": False},
                    },
                    "Sampling": {
                        "Method": "AdaptiveLevelSet",
                        "Seed": 7,
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                            {
                                "name": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {"min": 0.0, "max": 1.0},
                                },
                            },
                        ],
                        "AdaptiveLevelSet": {
                            "target_expression": "r2",
                            "target_value": 0.04,
                            "contour_precision": 0.05,
                            "function_tolerance": 0.08,
                            "initial_radius": 0.12,
                            "refinement_factor": 0.55,
                            "max_generations": 12,
                            "max_points": 800,
                        },
                    },
                    "Operas": {
                        "Modules": [
                            {
                                "name": "Circle",
                                "operator": "jarvishep2.testing.eggbox.circle_r2",
                                "call_mode": "call",
                                "input": [
                                    {"name": "x", "expression": "x"},
                                    {"name": "y", "expression": "y"},
                                ],
                                "output": [{"name": "r2", "entry": "r2"}],
                            }
                        ]
                    },
                }
                core = Jarvis2Core(config)
                count = core.run()
                self.assertGreater(count, 0)
                levelset_path = os.path.join(tmpdir, "levelset.json")
                self.assertTrue(os.path.isfile(levelset_path), "levelset.json missing")
                payload = json.loads(open(levelset_path, encoding="utf-8").read())
                self.assertEqual(payload["dim"], 2)
                self.assertIn("polylines_u", payload)
                cloud = payload.get("crossing_points_u") or []
                self.assertGreater(len(cloud), 4)
                # mean distance of crossing points to circle radius 0.2
                errs = []
                for u in cloud:
                    r = math.hypot(float(u[0]) - 0.5, float(u[1]) - 0.5)
                    errs.append(abs(r - 0.2))
                mean_err = float(np.mean(errs))
                self.assertLess(
                    mean_err,
                    0.08,
                    f"crossing cloud mean radial error {mean_err} too large",
                )
        finally:
            server.shutdown()
            server.server_close()

    def test_seeded_reproducibility_two_workers(self) -> None:
        """Same seed ⇒ same accepted uuids regardless of worker count (barrier + priority)."""

        def _run(workers: int) -> list[str]:
            server, redis_config = _start_tcp_fakeredis()
            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    config = {
                        "task_result_dir": tmpdir,
                        "Runtime": {
                            "mode": "redis",
                            "workers": workers,
                            "batch_size": 4,
                            "redis": redis_config,
                            "sample_artifacts": "never",
                            "Watchdog": {"enabled": False},
                        },
                        "Sampling": {
                            "Method": "AdaptiveLevelSet",
                            "Seed": 99,
                            "Variables": [
                                {
                                    "name": "x",
                                    "distribution": {
                                        "type": "Flat",
                                        "parameters": {"min": 0.0, "max": 1.0},
                                    },
                                },
                                {
                                    "name": "y",
                                    "distribution": {
                                        "type": "Flat",
                                        "parameters": {"min": 0.0, "max": 1.0},
                                    },
                                },
                            ],
                            "AdaptiveLevelSet": {
                                "target_expression": "r2",
                                "target_value": 0.04,
                                "contour_precision": 0.08,
                                "function_tolerance": 0.12,
                                "initial_radius": 0.15,
                                "max_generations": 6,
                                "max_points": 300,
                            },
                        },
                        "Operas": {
                            "Modules": [
                                {
                                    "name": "Circle",
                                    "operator": "jarvishep2.testing.eggbox.circle_r2",
                                    "call_mode": "call",
                                    "input": [
                                        {"name": "x", "expression": "x"},
                                        {"name": "y", "expression": "y"},
                                    ],
                                    "output": [{"name": "r2", "entry": "r2"}],
                                }
                            ]
                        },
                    }
                    core = Jarvis2Core(config)
                    core.run()
                    path = os.path.join(tmpdir, "levelset.json")
                    payload = json.loads(open(path, encoding="utf-8").read())
                    # deterministic point uuids from sampler state not exported; use
                    # n_points_total + crossing cloud size as weak check + re-export via sampler
                    return [
                        str(payload["n_points_total"]),
                        str(len(payload.get("crossing_points_u") or [])),
                        str(payload.get("n_generations")),
                    ]
            finally:
                server.shutdown()
                server.server_close()
                TaskFactory.reset_instance()

        a = _run(1)
        b = _run(2)
        self.assertEqual(a, b)

    def test_knn_vs_delaunay_crossing_regions_agree_on_dense_set(self) -> None:
        """§9.8: on a dense evaluated set, kNN finds similar crossings as Delaunay."""
        # Analytic f on a grid-like point set around a circle
        xs = np.linspace(0.15, 0.85, 12)
        ys = np.linspace(0.15, 0.85, 12)
        pts = np.array([[x, y] for x in xs for y in ys], dtype=np.float64)
        f = (pts[:, 0] - 0.5) ** 2 + (pts[:, 1] - 0.5) ** 2
        target = 0.04
        del_edges = DelaunayGraph().build(pts)
        knn_edges = KNNGraph(knn_k=8).build(pts)

        def crossings(edges: np.ndarray) -> set[tuple[int, int]]:
            out: set[tuple[int, int]] = set()
            for row in edges:
                i, j = int(row[0]), int(row[1])
                if (f[i] - target) * (f[j] - target) <= 0:
                    a, b = (i, j) if i < j else (j, i)
                    out.add((a, b))
            return out

        c_del = crossings(del_edges)
        c_knn = crossings(knn_edges)
        # kNN should recover a substantial fraction of Delaunay crossings
        if not c_del:
            self.skipTest("no Delaunay crossings on fixture")
        overlap = len(c_del & c_knn) / len(c_del)
        self.assertGreaterEqual(overlap, 0.5, f"kNN overlap {overlap:.2f} too low")


class GraphModeSmokeTests(unittest.TestCase):
    def test_d5_sobol_generation0_smoke(self) -> None:
        """Sobol gen-0 path at d=5 does not OOM and produces points."""
        sampler = AdaptiveLevelSetSampler()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "AdaptiveLevelSet",
                    "Seed": 3,
                    "Variables": [
                        {
                            "name": f"x{i}",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0.0, "max": 1.0},
                            },
                        }
                        for i in range(5)
                    ],
                    "AdaptiveLevelSet": {
                        "target_expression": "x0",
                        "target_value": 0.5,
                        "initial_radius": 0.25,
                        "max_points": 256,
                    },
                }
            }
        )
        us = sampler._generation_0_points()
        self.assertGreater(len(us), 8)
        for u in us[:5]:
            self.assertEqual(u.shape, (5,))
            self.assertTrue(np.all(u >= 0.0) and np.all(u <= 1.0))


if __name__ == "__main__":
    unittest.main()
