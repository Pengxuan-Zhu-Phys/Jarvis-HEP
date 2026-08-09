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

    def test_resume_fast_forwards_durable_prefix_without_submitting_physics_work(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        sampler.initialize()
        self.assertEqual(sampler.advance_to_persisted_prefix(2), 2)
        next_sample = sampler.propose_next()
        assert next_sample is not None
        self.assertEqual(next_sample.sample_index, 2)

    def test_bridson_prefix_replay_suppresses_submit_progress(self) -> None:
        """D21.14: local resume replay must not log 'samples submited' ‰ lines."""
        from jarvishep2.Sampling.bridson import Bridson

        sampler = Bridson()
        sampler.set_config(
            {
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"radius": 0.4, "max_attempt": 5, "seed": 7},
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0.0, "max": 1.0, "length": 1.0},
                            },
                        },
                        {
                            "name": "y",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0.0, "max": 1.0, "length": 1.0},
                            },
                        },
                    ],
                },
                "Runtime": {},
            }
        )
        sampler.initialize()
        target = min(5, int(sampler.info.get("NSamples") or 0))
        self.assertGreater(target, 0)
        messages: list[str] = []

        def _capture(msg: object, *args: object, **_kwargs: object) -> None:
            text = str(msg)
            if args:
                try:
                    text = text % args
                except Exception:
                    text = f"{text} {args}"
            messages.append(text)

        with (
            mock.patch.object(sampler._logger, "info", side_effect=_capture),
            mock.patch.object(sampler._logger, "warning", side_effect=_capture),
        ):
            advanced = sampler.advance_to_persisted_prefix(target)
        self.assertEqual(advanced, target)
        self.assertFalse(
            any("submited" in msg or "submitted" in msg for msg in messages),
            messages,
        )
        self.assertEqual(sampler.barinfo, {})
        # After replay, normal proposes may emit progress again.
        messages.clear()
        with (
            mock.patch.object(sampler._logger, "info", side_effect=_capture),
            mock.patch.object(sampler._logger, "warning", side_effect=_capture),
        ):
            sample = sampler.propose_next()
        self.assertIsNotNone(sample)
        self.assertTrue(any("submited" in msg for msg in messages), messages)

    def test_restore_point_advances_only_after_batch_is_persisted(self) -> None:
        sampler = Grid()
        sampler.set_config(_grid_test_config())
        queue = make_fakeredis_queue()
        queue.connect()
        sampler.set_redis(queue)
        sampler.set_execution_plan_template(include_likelihood=False)
        sampler.initialize()
        sampler.configure_checkpoint(checkpoint_file="", save_callback=lambda **_: None)
        try:
            first_batch = [sampler.propose_next(), sampler.propose_next()]
            first_batch = [sample for sample in first_batch if sample is not None]
            sampler.request_checkpoint_capture()
            sampler.flush_batch(first_batch)
            self.assertEqual(sampler.checkpoint_runtime_state(safe=False)["index"], 0)

            # A delayed sample index prevents promotion across its gap (D21.12 A8).
            sampler.set_persisted_prefix(1)
            self.assertEqual(sampler.checkpoint_runtime_state(safe=False)["index"], 0)
            sampler.set_persisted_prefix(2)
            self.assertEqual(sampler.checkpoint_runtime_state(safe=False)["index"], 2)

            third = sampler.propose_next()
            assert third is not None
            sampler.request_checkpoint_capture()
            sampler.flush_batch([third])
            self.assertEqual(sampler.checkpoint_runtime_state(safe=False)["index"], 2)
        finally:
            sampler.shutdown_checkpointing()


class BridsonSamplerUnitTests(unittest.TestCase):
    def test_format_duration_matches_v1_shape(self) -> None:
        self.assertRegex(format_duration(0.192), r"^\d{2}:\d{2}:\d{2}\.\d{3}$")
        self.assertEqual(format_duration(0.0)[:8], "00:00:00")

    def test_bridson_seed_zero_resume_restores_coordinates(self) -> None:
        """seed 0 must reseed (not skip); otherwise uuid matches but physics drifts."""
        cfg = {
            "EnvReqs": {"V2": {"workers": 1, "checkpoint_heartbeat_sec": 10}},
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
        state = first.export_runtime_state()
        expected = first.propose_next()
        self.assertIsNotNone(expected)

        restored = Bridson()
        restored.set_config(cfg)
        restored.import_runtime_state(state)
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
        """run_distributed must not exceed max_inflight without worker completions."""
        sampler = Bridson()
        sampler.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 2, "batch_size": 8},
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"radius": 0.12, "max_attempt": 30, "max_worker": 2, "seed": 3},
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
        self.assertEqual(sampler._max_inflight, 2)

        pushed_events: list[int] = []
        finish_after = {"n": 0}

        class _FakeRedis:
            def fetch_sample_stats(self):
                # Finish one sample whenever more than max_inflight would block forever.
                finished = max(0, finish_after["n"])
                return {"completed": finished, "failed": 0, "running": 0}

            def get_queue_lengths(self):
                return {"task_queue_length": 0}

        def _submit(sample):  # noqa: ANN001
            pushed_events.append(1)
            # Simulate a slow worker: allow progress only after 2 are outstanding.
            if len(pushed_events) - finish_after["n"] >= sampler._max_inflight:
                finish_after["n"] += 1

        def _submit_group(samples):  # noqa: ANN001
            for sample in samples:
                _submit(sample)

        sampler.redis = _FakeRedis()  # type: ignore[assignment]
        sampler.runtime_mode = "redis"
        sampler._submit = _submit  # type: ignore[method-assign]
        sampler._submit_group = _submit_group  # type: ignore[method-assign]

        # Track peak outstanding = submitted - finished_at_that_time
        peaks: list[int] = []
        original_wait = sampler._wait_for_inflight_room

        def _wait(pushed, *, base_finished, max_inflight):  # noqa: ANN001
            peaks.append(sampler._inflight(pushed, base_finished=base_finished))
            original_wait(pushed, base_finished=base_finished, max_inflight=max_inflight)

        sampler._wait_for_inflight_room = _wait  # type: ignore[method-assign]
        total = sampler.run_distributed()
        self.assertGreater(total, 0)
        self.assertEqual(total, len(pushed_events))
        # Never wait while already above the cap (inflight recorded before wait).
        self.assertTrue(all(p <= sampler._max_inflight for p in peaks))

    def test_submit_progress_heartbeat_logs_permille(self) -> None:
        """V1-style ``N‰ of i/N samples submited`` heartbeat while proposing."""
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
            if "‰ of" in line and "samples submited" in line
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
        state = sampler.export_runtime_state()
        expected = sampler.propose_next()
        assert expected is not None
        restored = RandomS()
        restored.set_config(self._random_config(points=10))
        restored.import_runtime_state(state)
        actual = restored.propose_next()
        assert actual is not None
        self.assertEqual(actual.sample_index, expected.sample_index)
        np.testing.assert_array_equal(actual.u_coords, expected.u_coords)

    def test_stateless_checkpoint_fields_are_all_explicitly_classified(self) -> None:
        random = RandomS()
        random.set_config(self._random_config(points=10))
        random.assert_checkpoint_attribute_contract()

        grid = Grid()
        grid.set_config(_grid_test_config())
        grid.initialize()
        grid.assert_checkpoint_attribute_contract()

        bridson = Bridson()
        bridson.set_config(
            {
                "Runtime": {"mode": "redis", "workers": 1},
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"radius": 0.3, "max_attempt": 8},
                    "Variables": [
                        {
                            "name": name,
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1, "length": 1},
                            },
                        }
                        for name in ("x", "y")
                    ],
                },
            }
        )
        bridson.assert_checkpoint_attribute_contract()

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "points.csv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("uuid,x\nu0,0.2\n")
            csv_sampler = CSVSampler()
            csv_sampler.set_config(
                {
                    "Runtime": {"mode": "redis", "workers": 1},
                    "Sampling": {"Method": "CSV", "Bounds": {"path": path}},
                }
            )
            csv_sampler.assert_checkpoint_attribute_contract()

    def test_heartbeat_only_requests_a_capture_until_the_next_batch_boundary(self) -> None:
        sampler = RandomS()
        sampler.set_config(self._random_config(points=10))
        calls = 0
        original = sampler.export_runtime_state

        def counted_export() -> dict:
            nonlocal calls
            calls += 1
            return original()

        sampler.export_runtime_state = counted_export  # type: ignore[method-assign]
        sampler.configure_checkpoint(checkpoint_file="", save_callback=lambda **_: None)
        try:
            calls = 0
            sampler.request_checkpoint_capture()
            self.assertEqual(calls, 0)
            self.assertTrue(sampler.capture_checkpoint_at_batch_boundary())
            self.assertEqual(calls, 1)
            self.assertFalse(sampler.capture_checkpoint_at_batch_boundary())
            self.assertEqual(calls, 1)
        finally:
            sampler.shutdown_checkpointing()

    def test_random_checkpoint_payload_is_constant_after_twenty_thousand_samples(self) -> None:
        """D21.12 P3/P4: no UUID or coordinate history enters the payload."""
        sampler = RandomS()
        sampler.set_config(self._random_config())
        sampler.initialize()
        before = pickle.dumps(sampler.export_runtime_state(), protocol=pickle.HIGHEST_PROTOCOL)
        while sampler.propose_next() is not None:
            pass
        after_state = sampler.export_runtime_state()
        after = pickle.dumps(after_state, protocol=pickle.HIGHEST_PROTOCOL)
        self.assertLess(abs(len(after) - len(before)), 1024)
        self.assertNotIn("u_by_uuid", after_state)
        self.assertNotIn("uuid_by_accepted_index", after_state)
        self.assertNotIn("submitted_uuids", after_state)
        self.assertFalse(hasattr(sampler, "_u_by_uuid"))
        self.assertFalse(hasattr(sampler, "_uuid_by_accepted_index"))


if __name__ == "__main__":
    unittest.main()
