#!/usr/bin/env python3
"""Sampler catalog gates — Bridson Poisson-disk / blue-noise path."""

from __future__ import annotations

import os
import pickle
import tempfile
import unittest
from unittest import mock

import numpy as np

from jarvishep2.Sampling import sampling_utils
from jarvishep2.Sampling.bridson import (
    Bridson,
    Bridson_sampling,
    hypersphere_surface_sample,
)
from jarvishep2.Sampling.csv_sampler import CSVSampler
from jarvishep2.log_kv import format_duration
from jarvishep2.Sampling.grid import Grid, grid_sampling
from jarvishep2.Sampling.randoms import RandomS
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue

from test_worker_mvp import _start_tcp_fakeredis


def _restore_snapshot(restored, source, state):
    # These unit tests supply their point-set fixture directly; all.pickle
    # loading and cross-file identity checks live in test_sample_checkpoint.
    from jarvishep2.sampling.completion_sampler import PointSetSampler
    if isinstance(source, PointSetSampler):
        restored._P = source._P.copy()
        restored._point_uuids = tuple(source._point_uuid(i) for i in range(len(source._P)))
        restored._dataset_id = state["dataset_id"]
    restored.restore_checkpoint(state)


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
V1_BRIDSON_OPERAS_YAML = os.path.join(
    TESTS_ROOT,
    "parity_project",
    "bridson_operas_v1.yaml",
)


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


class SelectionExpressionCacheTests(unittest.TestCase):
    def test_selection_expression_compiles_once_for_repeated_points(self) -> None:
        sampling_utils._SELECTION_EXPRESSIONS.clear_cache()
        for index in range(100):
            accepted = sampling_utils.evaluate_selection(
                "x + y < 2",
                {"x": index / 100.0, "y": 0.25},
            )
            self.assertIsInstance(accepted, bool)

        cache = sampling_utils._SELECTION_EXPRESSIONS.cache_info()
        self.assertEqual(cache.misses, 1)
        self.assertEqual(cache.hits, 99)


class DistributorDispatchTests(unittest.TestCase):
    def test_distributor_resolves_stateless_samplers(self) -> None:
        self.assertEqual(Distributor.set_method("Bridson").method, "Bridson")
        self.assertEqual(Distributor.set_method("Random").method, "Random")
        self.assertEqual(Distributor.set_method("Grid").method, "Grid")
        self.assertEqual(Distributor.set_method("CSV").method, "CSV")

    def test_unknown_method_lists_available(self) -> None:
        with self.assertRaises(NotImplementedError) as raised:
            Distributor.set_method("Dynestyy")
        message = str(raised.exception)
        self.assertIn("Dynestyy", message)
        self.assertIn("Available:", message)
        self.assertIn("Bridson", message)
        self.assertIn("Dynesty", message)

    def test_register_new_sampler_without_editing_set_method(self) -> None:
        from jarvishep2.Sampling.sampler import SamplingVirtial
        from jarvishep2.distributor import STATELESS_METHODS

        class _UnitDummy(SamplingVirtial):
            def __init__(self) -> None:
                super().__init__()
                self.method = "UnitDummy"

        Distributor.register(
            "UnitDummy",
            lambda: _UnitDummy(),
            stateless=True,
            resume="implemented",
            override=True,
        )
        try:
            sampler = Distributor.set_method("UnitDummy")
            self.assertEqual(sampler.method, "UnitDummy")
            self.assertIn("UnitDummy", STATELESS_METHODS)
            self.assertEqual(Distributor.get_resume_status("UnitDummy"), "implemented")
        finally:
            Distributor.unregister("UnitDummy")
        self.assertNotIn("UnitDummy", STATELESS_METHODS)
        with self.assertRaises(NotImplementedError):
            Distributor.set_method("UnitDummy")


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
            "Bounds": {"seed": 0},
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
        state = sampler.checkpoint_state()
        expected = sampler.propose_next()
        self.assertIsNotNone(expected)
        restored = Grid()
        restored.set_config(sampler.config)
        _restore_snapshot(restored, sampler, state)
        actual = restored.propose_next()
        self.assertIsNotNone(actual)
        self.assertEqual(actual.uuid, expected.uuid)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)

    def test_stateless_samples_have_contiguous_logical_indices(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        sampler.initialize()
        first = sampler.propose_next()
        second = sampler.propose_next()
        assert first is not None and second is not None
        self.assertEqual((first.sample_index, second.sample_index), (0, 1))
        self.assertEqual(first.uuid, sampler._uuid_for_accepted_index(0))
        self.assertEqual(second.uuid, sampler._uuid_for_accepted_index(1))

    def test_grid_submit_progress_logs_permille_heartbeats(self) -> None:
        sampler = Grid()
        config = _grid_test_config()
        for variable in config["Sampling"]["Variables"]:
            variable["distribution"]["parameters"]["num"] = 10
        sampler.set_config(config)
        info_messages: list[str] = []
        warning_messages: list[str] = []

        class _Capture:
            def info(self, msg, *args, **kwargs):  # noqa: ANN001
                info_messages.append(str(msg) % args if args else str(msg))

            def warning(self, msg, *args, **kwargs):  # noqa: ANN001
                warning_messages.append(str(msg) % args if args else str(msg))

        sampler._logger = _Capture()  # type: ignore[assignment]
        sampler.initialize()
        total = int(sampler.info["NSamples"])
        for _ in range(total // 2):
            self.assertIsNotNone(sampler.propose_next())

        progress_lines = [
            line
            for line in info_messages + warning_messages
            if "‰ of" in line and "grid points sampled" in line
        ]
        self.assertTrue(progress_lines, "expected Grid submit progress heartbeats")
        self.assertTrue(any(line.startswith("0‰ of 0/100") for line in progress_lines))
        self.assertGreaterEqual(len(progress_lines), 2)
        for line in warning_messages:
            if "‰ of" not in line:
                continue
            permille = int(line.split("‰", 1)[0])
            self.assertGreater(permille, 0)
            self.assertEqual(permille % 10, 0)

class BridsonSamplerUnitTests(unittest.TestCase):
    def test_format_duration_matches_v1_shape(self) -> None:
        self.assertRegex(format_duration(0.192), r"^\d{2}:\d{2}:\d{2}\.\d{3}$")
        self.assertEqual(format_duration(0.0)[:8], "00:00:00")

    def test_bridson_seed_zero_resume_restores_coordinates(self) -> None:
        """seed 0 must reseed (not skip); otherwise uuid matches but physics drifts."""
        cfg = {
            "EnvReqs": {"V2": {"workers": 1, "checkpoint": {"heartbeat": 30}}},
            "Sampling": {
                "Method": "Bridson",
                "Bounds": {"seed": 0, "radius": 0.35, "max_attempt": 12},
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
        first = Bridson()
        first.set_config(cfg)
        first.initialize()
        self.assertIsNotNone(first.propose_next())
        state = first.checkpoint_state()
        expected = first.propose_next()
        self.assertIsNotNone(expected)

        restored = Bridson()
        restored.set_config(cfg)
        _restore_snapshot(restored, first, state)
        actual = restored.propose_next()
        self.assertIsNotNone(actual)
        assert expected is not None and actual is not None
        self.assertEqual(actual.uuid, expected.uuid)
        self.assertEqual(actual.sample_index, expected.sample_index)
        np.testing.assert_allclose(actual.u_coords, expected.u_coords, atol=1e-12)

    def test_random_seed_zero_is_deterministic_across_fresh_inits(self) -> None:
        cfg = {
            "EnvReqs": {"V2": {"workers": 1}},
            "Sampling": {
                "Method": "Random",
                "Bounds": {"seed": 0, "point_number": 5},
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
        a = RandomS()
        a.set_config(cfg)
        ua = [a.propose_next().u_coords.copy() for _ in range(3)]
        b = RandomS()
        b.set_config(cfg)
        ub = [b.propose_next().u_coords.copy() for _ in range(3)]
        for pa, pb in zip(ua, ub):
            np.testing.assert_allclose(pa, pb, atol=1e-12)

    def test_batch_size_uses_normalized_runtime_default(self) -> None:
        """Missing Runtime.batch_size must keep FixedSetSampler/runtime default (256), not max_worker."""
        sampler = Bridson()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 2},
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"radius": 0.35, "max_attempt": 30, "max_worker": 2, "seed": 1},
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
        self.assertEqual(sampler._batch_size, 256)
        self.assertEqual(sampler._max_inflight, 2)

    def test_max_inflight_backpressure_limits_outstanding_tasks(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        sampler._max_inflight = 2
        pending = []
        peaks = []

        class Queue:
            def read_proc_board(self, role):
                return {}

            def push_task(self, task):
                pending.append(task)
                peaks.append(len(pending))
                assert task["sample_index"] in sampler._unfinished

            def push_many_tasks(self, tasks):
                for task in tasks:
                    self.push_task(task)

            def poll_feedback_batch(self, limit):
                records = []
                while pending and len(records) < limit:
                    task = pending.pop()
                    records.append({"uuid": task["uuid"],
                                    "sample_index": task["sample_index"],
                                    "terminal_status": "Completed"})
                return records

        sampler.set_redis(Queue())
        total = sampler.run_distributed(timeout=2)
        self.assertEqual(total, 4)
        self.assertLessEqual(max(peaks), 2)
        self.assertFalse(sampler._unfinished)
        self.assertEqual(sampler._completed_count, total)

    def test_submit_progress_heartbeat_logs_permille(self) -> None:
        """V1-style ``N‰ of i/N samples submitted`` heartbeat while proposing."""
        sampler = Bridson()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 1},
                "Sampling": {
                    "Method": "Bridson",
                    # Smaller radius → denser grid → multiple ‰ heartbeats.
                    "Bounds": {"radius": 0.08, "max_attempt": 30, "seed": 7},
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
        info_messages: list[str] = []
        warning_messages: list[str] = []

        class _Capture:
            def info(self, msg, *args, **kwargs):  # noqa: ANN001
                text = str(msg) % args if args else str(msg)
                info_messages.append(text)

            def warning(self, msg, *args, **kwargs):  # noqa: ANN001
                text = str(msg) % args if args else str(msg)
                warning_messages.append(text)

        sampler._logger = _Capture()  # type: ignore[assignment]
        sampler.initialize()
        total = int(sampler.info["NSamples"])
        self.assertGreaterEqual(total, 20)
        # Drain enough of the grid so several ‰ milestones fire.
        for _ in range(max(1, total // 2)):
            if sampler.propose_next() is None:
                break

        progress_lines = [
            line
            for line in info_messages + warning_messages
            if "‰ of" in line and "samples submitted" in line
        ]
        self.assertTrue(progress_lines, "expected Bridson submit progress heartbeats")
        self.assertTrue(any(line.startswith("0‰ of") for line in progress_lines))
        self.assertGreaterEqual(len(progress_lines), 2)
        # 1% milestones use warning level when they appear.
        warning_progress = [line for line in warning_messages if "‰ of" in line]
        for line in warning_progress:
            permille = int(line.split("‰", 1)[0])
            self.assertEqual(permille % 10, 0)
            self.assertGreater(permille, 0)

    def test_checkpoint_roundtrip_restores_grid_cursor(self) -> None:
        sampler = Bridson()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 1},
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"radius": 0.35, "max_attempt": 30, "seed": 42},
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
        state = sampler.checkpoint_state()
        expected = sampler.propose_next()
        self.assertIsNotNone(expected)
        restored = Bridson()
        restored.set_config(sampler.config)
        _restore_snapshot(restored, sampler, state)
        actual = restored.propose_next()
        self.assertIsNotNone(actual)
        self.assertEqual(actual.uuid, expected.uuid)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)


class RandomSamplerProgressTests(unittest.TestCase):
    @staticmethod
    def _config(*, points: int = 1000) -> dict:
        return {
            "Runtime": {"mode": "redis", "workers": 1},
            "Sampling": {
                "Method": "Random",
                "Bounds": {"point_number": points, "seed": 7},
                "Variables": [
                    {
                        "name": "x",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0, "max": 1, "length": 1},
                        },
                    },
                ],
            },
        }

    def test_logs_each_permille_and_warns_at_each_percent(self) -> None:
        sampler = RandomS()
        sampler.set_config(self._config())
        info_messages: list[str] = []
        warning_messages: list[str] = []

        class _Capture:
            def info(self, msg, *args, **kwargs):  # noqa: ANN001
                info_messages.append(str(msg) % args if args else str(msg))

            def warning(self, msg, *args, **kwargs):  # noqa: ANN001
                warning_messages.append(str(msg) % args if args else str(msg))

        sampler._logger = _Capture()  # type: ignore[assignment]
        for _ in range(10):
            self.assertIsNotNone(sampler.propose_next())

        self.assertTrue(any(line.startswith("0‰ of 0/1000") for line in info_messages))
        for permille in range(1, 10):
            self.assertTrue(
                any(line.startswith(f"{permille}‰ of {permille}/1000") for line in info_messages),
                f"missing info heartbeat for {permille}‰",
            )
        self.assertTrue(any(line.startswith("10‰ of 10/1000") for line in warning_messages))
        self.assertFalse(any(line.startswith("10‰") for line in info_messages))

class StatelessDistributedRunTests(unittest.TestCase):
    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def _run_task_yaml(
        self,
        yaml_path: str,
        *,
        expected_count: int | None = None,
        render_plots: bool = False,
    ) -> int:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core()
                core.load_task_yaml(yaml_path)
                core.config["task_result_dir"] = tmpdir
                core.config["Runtime"]["redis"] = redis_config
                core.config["Runtime"]["Watchdog"] = {"enabled": False}
                # The JarvisPLOT render takes ~26 s of each run; only the
                # Bridson run checks it.
                core.config["skip_render_plots"] = not render_plots
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                count = core.run(write_run_summary=False)
                self.assertGreater(count, 0)
                if expected_count is not None:
                    self.assertEqual(count, expected_count)
                if render_plots:
                    self.assertEqual(core.info["plot_scenes"].get("rendered"), "true")

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
        self._run_task_yaml(BRIDSON_YAML, render_plots=True)

    def test_v1_shaped_bridson_operas_yaml_end_to_end(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core()
                core.load_task_yaml(V1_BRIDSON_OPERAS_YAML)
                self.assertEqual(core.config["Runtime"]["mode"], "redis")
                core.config["task_root"] = tmpdir
                core.config["project_root"] = tmpdir
                core.config["task_result_dir"] = tmpdir
                core.config["Runtime"]["redis"] = redis_config
                core.config["Runtime"]["Watchdog"] = {"enabled": False}
                core.config["skip_render_plots"] = True
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                count = core.run(write_run_summary=False)
                self.assertGreater(count, 0)

                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                records = SimpleHDF5Writer(db_path).read_records()
                self.assertEqual(len(records), count)
                for row in records:
                    self.assertIn("xx", row)
                    self.assertIn("yy", row)
                    self.assertIn("z", row)
                    self.assertIn("LogL_Z", row)
                    self.assertIn("LogL", row)
        finally:
            server.shutdown()
            server.server_close()

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
                "Bounds": {"point_number": 4, "seed": 99},
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


class CheckpointScalingTests(unittest.TestCase):
    @staticmethod
    def _random_config(points: int = 20_000) -> dict:
        return {
            "Runtime": {"mode": "redis", "workers": 1, "batch_size": 8},
            "Sampling": {
                "Method": "Random",
                "Bounds": {"point_number": points, "seed": 123},
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

    def test_random_export_import_restores_the_seeded_generator(self) -> None:
        sampler = RandomS()
        sampler.set_config(self._random_config(points=10))
        sampler.propose_next()
        sampler.propose_next()
        state = sampler.checkpoint_state()
        expected = sampler.propose_next()
        assert expected is not None
        restored = RandomS()
        restored.set_config(self._random_config(points=10))
        _restore_snapshot(restored, sampler, state)
        actual = restored.propose_next()
        assert actual is not None
        self.assertEqual(actual.sample_index, expected.sample_index)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)

    def test_random_checkpoint_payload_is_constant_after_twenty_thousand_samples(self) -> None:
        """D21.12 P3/P4: no UUID or coordinate history enters the payload."""
        sampler = RandomS()
        sampler.set_config(self._random_config())
        sampler.initialize()
        before = pickle.dumps(sampler.checkpoint_state(), protocol=pickle.HIGHEST_PROTOCOL)
        while (sample := sampler.propose_next()) is not None:
            sampler.accept_completion({"uuid": sample.uuid,
                                       "sample_index": sample.sample_index,
                                       "terminal_status": "Completed"})
        after_state = sampler.checkpoint_state()
        after = pickle.dumps(after_state, protocol=pickle.HIGHEST_PROTOCOL)
        self.assertLess(abs(len(after) - len(before)), 1024)
        self.assertNotIn("u_by_uuid", after_state)
        self.assertNotIn("uuid_by_accepted_index", after_state)
        self.assertNotIn("submitted_uuids", after_state)
        self.assertFalse(hasattr(sampler, "_u_by_uuid"))
        self.assertFalse(hasattr(sampler, "_uuid_by_accepted_index"))


if __name__ == "__main__":
    unittest.main()
