#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
import sys
import threading
import tempfile
import unittest
from copy import deepcopy

import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.bridson import Bridson  # noqa: E402
from jarvishep.Sampling.csv_sampler import CSVSampler  # noqa: E402
from jarvishep.Sampling.grid import Grid  # noqa: E402
from jarvishep.Sampling.randoms import RandomS  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


def _build_run_spec(config: dict, scan_name: str, task_root: str) -> dict:
    return {
        "raw_yaml_text": None,
        "normalized_config": deepcopy(config),
        "scan_name": scan_name,
        "task_root": task_root,
        "task_result_dir": os.path.join(task_root, "RESULTS"),
        "logs_dir": os.path.join(task_root, "LOGS"),
        "images_dir": os.path.join(task_root, "IMAGES"),
        "worker_parallel": 1,
        "sampler_method": config.get("Sampling", {}).get("Method"),
        "workflow": {},
        "workflow_layers": {},
    }


def _pending_info(base_dir: str, uuid: str, params: dict, *, extra: dict | None = None) -> dict:
    info = {
        "uuid": uuid,
        "params": deepcopy(params),
        "observables": {**deepcopy(params), "uuid": uuid},
        "save_dir": os.path.join(base_dir, uuid),
        "run_log": os.path.join(base_dir, uuid, "Sample_running.log"),
        "logger": None,
        "handlers": {},
        "status": "Running",
    }
    if extra:
        info.update(deepcopy(extra))
    return info


class TestDirectSamplerCheckpoint(unittest.TestCase):
    def _setup_sample_context(self, sampler, tempdir: str, method_cfg: dict) -> None:
        sampler.logger = _NoopLogger()
        sampler.info = {
            "sample": {
                "task_result_dir": tempdir,
                "sample_dirs": os.path.join(tempdir, "SAMPLE"),
                "archive_samples": False,
            }
        }
        sampler.set_config(method_cfg)

    def _checkpoint_roundtrip(self, sampler, checkpoint_root: str, run_spec: dict, pending_info: dict):
        sampler.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        sampler.set_runtime_checkpoint_context(run_spec=run_spec, factory_blueprint={})
        sampler.future_to_sample = {object(): type("Pending", (), {"info": pending_info})()}
        self.assertTrue(sampler.persist_runtime_checkpoint(force=True, reason="unit-test"))
        self.assertTrue(os.path.exists(os.path.join(checkpoint_root, "state.pkl")))
        self.assertFalse(os.path.exists(os.path.join(checkpoint_root, "runtime_state.pkl")))

        restored = sampler.__class__()
        restored.logger = _NoopLogger()
        restored.info = deepcopy(sampler.info)
        restored.set_config(run_spec["normalized_config"])
        restored.configure_runtime_checkpointing(checkpoint_root, logger=_NoopLogger())
        self.assertTrue(restored.restore_runtime_checkpoint_if_available())
        return restored

    def test_grid_checkpoint_roundtrip_restores_grid_points_and_pending_samples(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-grid-") as tmp:
            sampler = Grid()
            config = {
                "Scan": {"sample_directory": {"limit": 4, "width": 3}},
                "Sampling": {
                    "Method": "Grid",
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "num": 3}},
                        }
                    ],
                },
            }
            self._setup_sample_context(sampler, tmp, config)
            sampler._P = np.array([[0.1], [0.2], [0.3]], dtype=np.float64)
            sampler._index = 2
            sampler.bucket_alloc.set_state({"base_path": sampler.bucket_alloc.base_path, "bucket": 3, "count": 2, "limit": 4, "width": 3})
            pending = _pending_info(tmp, "grid-pending", {"x": 0.2})
            run_spec = _build_run_spec(config, "grid-resume", tmp)

            restored = self._checkpoint_roundtrip(sampler, os.path.join(tmp, "checkpoints"), run_spec, pending)

            np.testing.assert_allclose(restored._P, sampler._P)
            self.assertEqual(restored._index, 2)
            self.assertEqual(restored._runtime_pending_samples[0]["uuid"], "grid-pending")
            self.assertEqual(restored.bucket_alloc.get_state()["bucket"], 3)
            self.assertEqual(restored.bucket_alloc.get_state()["count"], 2)

    def test_random_checkpoint_roundtrip_restores_rng_and_pending_samples(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-random-") as tmp:
            sampler = RandomS()
            config = {
                "Scan": {"sample_directory": {"limit": 4, "width": 3}},
                "Sampling": {
                    "Method": "Random",
                    "Point number": 5,
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                        }
                    ],
                },
            }
            self._setup_sample_context(sampler, tmp, config)
            sampler._index = 0
            sampler.bucket_alloc.set_state({"base_path": sampler.bucket_alloc.base_path, "bucket": 2, "count": 1, "limit": 4, "width": 3})
            np.random.seed(1234)
            first = sampler.next_sample()
            pending = _pending_info(tmp, "random-pending", first)
            run_spec = _build_run_spec(config, "random-resume", tmp)

            restored = self._checkpoint_roundtrip(sampler, os.path.join(tmp, "checkpoints"), run_spec, pending)
            self.assertEqual(restored._runtime_pending_samples[0]["uuid"], "random-pending")

            actual = restored.next_sample()
            control = RandomS()
            self._setup_sample_context(control, tmp, config)
            control._index = 1
            control.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            control.set_runtime_checkpoint_context(run_spec=run_spec, factory_blueprint={})
            self.assertTrue(control.restore_runtime_checkpoint_if_available())
            expected = control.next_sample()
            self.assertEqual(actual, expected)
            self.assertEqual(restored._index, sampler._index + 1)

    def test_bridson_checkpoint_roundtrip_restores_barinfo_and_pending_samples(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-bridson-") as tmp:
            sampler = Bridson()
            config = {
                "Scan": {"sample_directory": {"limit": 4, "width": 3}},
                "Sampling": {
                    "Method": "Bridson",
                    "Radius": 0.1,
                    "MaxAttempt": 5,
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "length": 1.0}},
                        },
                        {
                            "name": "y",
                            "description": "y",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "length": 1.0}},
                        },
                    ],
                },
            }
            self._setup_sample_context(sampler, tmp, config)
            sampler._P = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64)
            sampler._index = 1
            sampler.barinfo = {"total": 2, "t0": 1.0, "permille": 500}
            sampler.bucket_alloc.set_state({"base_path": sampler.bucket_alloc.base_path, "bucket": 4, "count": 3, "limit": 4, "width": 3})
            pending = _pending_info(tmp, "bridson-pending", {"x": 0.1, "y": 0.2})
            run_spec = _build_run_spec(config, "bridson-resume", tmp)

            restored = self._checkpoint_roundtrip(sampler, os.path.join(tmp, "checkpoints"), run_spec, pending)

            np.testing.assert_allclose(restored._P, sampler._P)
            self.assertEqual(restored._index, 1)
            self.assertEqual(restored.barinfo["permille"], 500)
            self.assertEqual(restored._runtime_pending_samples[0]["uuid"], "bridson-pending")
            self.assertEqual(restored.bucket_alloc.get_state()["bucket"], 4)

    def test_bridson_checkpoint_roundtrip_sanitizes_unserializable_pending_handles(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-bridson-lock-") as tmp:
            sampler = Bridson()
            config = {
                "Scan": {"sample_directory": {"limit": 4, "width": 3}},
                "Sampling": {
                    "Method": "Bridson",
                    "Radius": 0.1,
                    "MaxAttempt": 5,
                    "Variables": [
                        {
                            "name": "x",
                            "description": "x",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "length": 1.0}},
                        },
                        {
                            "name": "y",
                            "description": "y",
                            "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0, "length": 1.0}},
                        },
                    ],
                },
            }
            self._setup_sample_context(sampler, tmp, config)
            sampler._P = np.array([[0.1, 0.2], [0.3, 0.4]], dtype=np.float64)
            sampler._index = 1
            sampler.barinfo = {"total": 2, "t0": 1.0, "permille": 500}
            sampler.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            sampler.set_runtime_checkpoint_context(
                run_spec=_build_run_spec(config, "bridson-lock-resume", tmp),
                factory_blueprint={},
            )
            pending_info = {
                "uuid": "bridson-lock-pending",
                "params": {"x": 0.1, "y": 0.2},
                "observables": {"x": 0.1, "y": 0.2, "uuid": "bridson-lock-pending"},
                "save_dir": os.path.join(tmp, "pending", "bridson-lock-pending"),
                "run_log": os.path.join(tmp, "pending", "bridson-lock-pending", "Sample_running.log"),
                "logger": object(),
                "handlers": {"lock": threading.Lock()},
                "status": "Running",
            }
            sampler.future_to_sample = {object(): type("PendingSample", (), {"info": pending_info})()}

            self.assertTrue(sampler.persist_runtime_checkpoint(force=True, reason="unit-test"))

            restored = Bridson()
            self._setup_sample_context(restored, tmp, config)
            restored.configure_runtime_checkpointing(os.path.join(tmp, "checkpoints"), logger=_NoopLogger())
            restored.set_runtime_checkpoint_context(
                run_spec=_build_run_spec(config, "bridson-lock-resume", tmp),
                factory_blueprint={},
            )
            self.assertTrue(restored.restore_runtime_checkpoint_if_available())
            self.assertEqual(restored._runtime_pending_samples[0]["uuid"], "bridson-lock-pending")
            self.assertEqual(restored._runtime_pending_samples[0]["handlers"], {})
            self.assertIsNone(restored._runtime_pending_samples[0]["logger"])

    def test_csv_checkpoint_roundtrip_restores_cursor_and_uuid_tracking(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-") as tmp:
            csv_path = os.path.join(tmp, "points.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["uuid", "x"])
                writer.writeheader()
                writer.writerow({"uuid": "u1", "x": "1"})
                writer.writerow({"uuid": "u2", "x": "2"})
                writer.writerow({"uuid": "u3", "x": "3"})

            sampler = CSVSampler()
            config = {
                "Scan": {"sample_directory": {"limit": 4, "width": 3}},
                "Sampling": {
                    "Method": "CSV",
                    "CSV": {"path": csv_path, "uuid_column": "uuid"},
                },
            }
            self._setup_sample_context(sampler, tmp, config)
            sampler.set_bucket_alloc()
            sampler._runtime_csv_cursor = 2
            sampler._runtime_seen_source_uuid = {"u1", "u2"}
            sampler.bucket_alloc.set_state({"base_path": sampler.bucket_alloc.base_path, "bucket": 2, "count": 1, "limit": 4, "width": 3})
            pending = _pending_info(tmp, "csv-pending", {"x": 2}, extra={"csv_row_index": 2, "csv_source_uuid": "u2"})
            run_spec = _build_run_spec(config, "csv-resume", tmp)

            restored = self._checkpoint_roundtrip(sampler, os.path.join(tmp, "checkpoints"), run_spec, pending)

            self.assertEqual(restored._runtime_csv_cursor, 2)
            self.assertEqual(restored._runtime_seen_source_uuid, {"u1", "u2"})
            self.assertEqual(restored._runtime_pending_samples[0]["csv_row_index"], 2)
            self.assertEqual(restored._runtime_pending_samples[0]["csv_source_uuid"], "u2")
            self.assertEqual(restored.next_sample(), {"x": 3})


if __name__ == "__main__":
    unittest.main()
