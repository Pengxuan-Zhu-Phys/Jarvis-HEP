#!/usr/bin/env python3
"""Sampler catalog gates — Bridson Poisson-disk / blue-noise path."""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from jarvishep2.Sampling.bridson import Bridson, Bridson_sampling, hypersphere_surface_sample
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor
from jarvishep2.factory import TaskFactory

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


class BridsonSamplerUnitTests(unittest.TestCase):
    def test_distributor_resolves_bridson(self) -> None:
        sampler = Distributor.set_method("Bridson")
        self.assertEqual(sampler.method, "Bridson")

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


class BridsonDistributedRunTests(unittest.TestCase):
    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def test_bridson_yaml_end_to_end_via_core_run(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core()
                core.load_task_yaml(BRIDSON_YAML)
                core.config["task_result_dir"] = tmpdir
                core.config["Runtime"]["redis"] = redis_config
                core.config["Runtime"]["Watchdog"] = {"enabled": False}
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                count = core.run(write_run_summary=False)
                self.assertGreater(count, 0)

                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                records = _normalize_bridson_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(len(records), count)
                for row in records:
                    self.assertIn("x", row)
                    self.assertIn("y", row)
                    self.assertIn("LogL_Z", row)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()