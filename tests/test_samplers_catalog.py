#!/usr/bin/env python3
"""Sampler catalog gates — Bridson Poisson-disk / blue-noise path."""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from jarvishep2.Sampling.bridson import Bridson, Bridson_sampling, hypersphere_surface_sample
from jarvishep2.Sampling.grid import Grid, grid_sampling
from jarvishep2.Sampling.randoms import RandomS
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue

from test_worker_mvp import _start_tcp_fakeredis


def _normalize_bridson_records(records: list[dict]) -> list[dict[str, float]]:
    normalized = []
    for row in records:
        normalized.append(
            {
                "x": round(float(row["x"]), 8),
                "y": round(float(row["y"]), 8),
                "z": round(float(row["z"]), 8),
                "LogL_Z": round(float(row["LogL_Z"]), 8),
            }
        )
    return sorted(normalized, key=lambda item: (item["x"], item["y"]))


TESTS_ROOT = os.path.dirname(__file__)
BRIDSON_YAML = os.path.join(TESTS_ROOT, "parity_project", "bridson_opera.yaml")
RANDOM_YAML = os.path.join(TESTS_ROOT, "parity_project", "random_opera.yaml")
GRID_YAML = os.path.join(TESTS_ROOT, "parity_project", "grid_opera.yaml")
CSV_YAML = os.path.join(TESTS_ROOT, "parity_project", "csv_opera.yaml")


def _stop_factory_workers() -> None:
    factory = TaskFactory._instance
    if factory is not None:
        try:
            factory.shutdown(wait=False)
        except Exception:
            pass
    TaskFactory.reset_instance()


class BridsonAlgorithmTests(unittest.TestCase):
    def test_bridson_sampling_seeded_point_count(self) -> None:
        np.random.seed(7)
        points = Bridson_sampling(
            dims=np.array([1.0, 1.0]),
            radius=0.35,
            k=30,
            hypersphere_sample=hypersphere_surface_sample,
        )
        self.assertGreaterEqual(points.shape[0], 4)
        self.assertEqual(points.shape[1], 2)
        np.random.seed(7)
        repeat = Bridson_sampling(
            dims=np.array([1.0, 1.0]),
            radius=0.35,
            k=30,
            hypersphere_sample=hypersphere_surface_sample,
        )
        np.testing.assert_array_equal(points, repeat)

    def test_neighborhood_axis_fix_does_not_raise_in_2d(self) -> None:
        np.random.seed(0)
        points = Bridson_sampling(
            dims=np.array([1.0, 1.0]),
            radius=0.2,
            k=10,
            hypersphere_sample=hypersphere_surface_sample,
        )
        self.assertGreater(points.shape[0], 0)


class DistributorDispatchTests(unittest.TestCase):
    def test_distributor_resolves_stateless_samplers(self) -> None:
        self.assertEqual(Distributor.set_method("Bridson").method, "Bridson")
        self.assertEqual(Distributor.set_method("Random").method, "Random")
        self.assertEqual(Distributor.set_method("Grid").method, "Grid")
        self.assertEqual(Distributor.set_method("CSV").method, "CSV")


class GridAlgorithmTests(unittest.TestCase):
    def test_grid_sampling_cartesian_product_size(self) -> None:
        points = grid_sampling(np.array([2, 3], dtype=np.int64))
        self.assertEqual(points.shape, (6, 2))

    def test_grid_sampling_clips_endpoints_for_transforms(self) -> None:
        points = grid_sampling(np.array([3], dtype=np.int64))
        eps = np.finfo(np.float64).eps
        self.assertGreater(float(points[0, 0]), 0.0)
        self.assertLessEqual(float(points[-1, 0]), 1.0 - eps)


def _grid_test_config(**overrides: object) -> dict:
    cfg = {
        "Runtime": {"mode": "redis", "workers": 1},
        "Sampling": {
            "Method": "Grid",
            "Seed": 0,
            "Variables": [
                {
                    "name": "x",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0, "max": 1, "num": 2},
                    },
                },
                {
                    "name": "y",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0, "max": 1, "num": 2},
                    },
                },
            ],
        },
    }
    cfg.update(overrides)
    return cfg


class GridSamplerUnitTests(unittest.TestCase):
    def test_grid_requires_num_per_variable(self) -> None:
        sampler = Grid()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis"},
                "Sampling": {
                    "Method": "Grid",
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        },
                    ],
                },
            }
        )
        with self.assertRaisesRegex(ValueError, "parameters.num"):
            sampler.initialize()

    def test_checkpoint_roundtrip_restores_grid_cursor(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        sampler.initialize()
        state = sampler.export_runtime_state()
        expected = sampler.propose_next()
        self.assertIsNotNone(expected)
        restored = Grid()
        restored.set_config(sampler.config)
        restored.import_runtime_state(state)
        actual = restored.propose_next()
        self.assertIsNotNone(actual)
        self.assertEqual(actual.uuid, expected.uuid)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)

    def test_repropose_unfinished_requeues_pending_uuids(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        queue = make_fakeredis_queue()
        queue.connect()
        sampler.set_redis(queue)
        sampler.set_execution_plan_template(include_likelihood=False)
        sampler.initialize()
        first = sampler.propose_next()
        second = sampler.propose_next()
        assert first is not None and second is not None
        sampler._submit(first)
        sampler._submit(second)
        sampler._submitted_uuids.extend([first.uuid, second.uuid])
        queue.drain_task_queue()
        sampler.mark_completed(first.uuid)
        sampler.set_resume_repropose_hint(True)

        requeued = sampler.repropose_unfinished()
        self.assertEqual(requeued, [second.uuid])
        self.assertEqual(int(queue.r.llen("hep:task_queue")), 1)


class BridsonSamplerUnitTests(unittest.TestCase):

    def test_checkpoint_roundtrip_restores_grid_cursor(self) -> None:
        sampler = Bridson()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 1},
                "Sampling": {
                    "Method": "Bridson",
                    "Radius": 0.35,
                    "MaxAttempt": 30,
                    "Seed": 42,
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1, "length": 1},
                            },
                        },
                        {
                            "name": "y",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1, "length": 1},
                            },
                        },
                    ],
                },
            }
        )
        sampler.initialize()
        state = sampler.export_runtime_state()
        expected = sampler.propose_next()
        self.assertIsNotNone(expected)
        restored = Bridson()
        restored.set_config(sampler.config)
        restored.import_runtime_state(state)
        actual = restored.propose_next()
        self.assertIsNotNone(actual)
        self.assertEqual(actual.uuid, expected.uuid)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)


class StatelessDistributedRunTests(unittest.TestCase):
    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def _run_task_yaml(self, yaml_path: str, *, expected_count: int | None = None) -> int:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core()
                core.load_task_yaml(yaml_path)
                core.config["task_result_dir"] = tmpdir
                core.config["Runtime"]["redis"] = redis_config
                core.config["Runtime"]["Watchdog"] = {"enabled": False}
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                count = core.run(write_run_summary=False)
                self.assertGreater(count, 0)
                if expected_count is not None:
                    self.assertEqual(count, expected_count)

                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                records = _normalize_bridson_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(len(records), count)
                for row in records:
                    self.assertIn("x", row)
                    self.assertIn("y", row)
                    self.assertIn("LogL_Z", row)
                return count
        finally:
            server.shutdown()
            server.server_close()

    def test_bridson_yaml_end_to_end_via_core_run(self) -> None:
        self._run_task_yaml(BRIDSON_YAML)

    def test_random_yaml_end_to_end_via_core_run(self) -> None:
        self._run_task_yaml(RANDOM_YAML, expected_count=6)

    def test_grid_yaml_end_to_end_via_core_run(self) -> None:
        self._run_task_yaml(GRID_YAML, expected_count=9)

    def test_csv_yaml_end_to_end_via_core_run(self) -> None:
        self._run_task_yaml(CSV_YAML, expected_count=10)

    def test_random_seeded_proposals_are_reproducible(self) -> None:
        cfg = {
            "Runtime": {"mode": "redis", "workers": 1},
            "Sampling": {
                "Method": "Random",
                "Point number": 4,
                "Seed": 99,
                "Variables": [
                    {
                        "name": "x",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0, "max": 1, "length": 1},
                        },
                    },
                    {
                        "name": "y",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0, "max": 1, "length": 1},
                        },
                    },
                ],
            },
        }
        first = RandomS()
        first.set_config(cfg)
        coords_a = [first.propose_next().u_coords.tolist() for _ in range(4)]
        second = RandomS()
        second.set_config(cfg)
        coords_b = [second.propose_next().u_coords.tolist() for _ in range(4)]
        self.assertEqual(coords_a, coords_b)


if __name__ == "__main__":
    unittest.main()