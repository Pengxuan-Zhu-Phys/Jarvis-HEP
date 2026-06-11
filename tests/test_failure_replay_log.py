#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from datetime import datetime


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.likelihood import LogLikelihood  # noqa: E402
from jarvishep.Module.module import Module  # noqa: E402
from jarvishep.factory import WorkerFactory  # noqa: E402
from jarvishep.moduleManager import ModuleManager  # noqa: E402
from jarvishep.sample import Sample, materialize_failure_artifacts  # noqa: E402
from jarvishep.sample_logger import BufferedSampleLogger, SampleLogger  # noqa: E402


class _LazySampleConfigMixin:
    @staticmethod
    def _lazy_sample_cfg(tmpdir: str) -> dict:
        return {
            "sample_dirs": tmpdir,
            "task_result_dir": tmpdir,
            "sample_artifacts": "auto",
            "workflow_has_calculator": False,
            "workflow_references_sdir": False,
        }


class BufferedSampleLoggerTests(_LazySampleConfigMixin, unittest.TestCase):
    def test_replay_preserves_structured_and_raw_format(self):
        fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            logger = BufferedSampleLogger(
                extra={"module": "Sample@test-uuid"},
                time_provider=lambda: fixed_dt,
            )
            child = logger.bind(module="Sample@test-uuid (Likelihood)")
            child.info("hello world")
            child.bind(raw=True).info("stdout line")
            child.warning("warn")
            logger.replay_to(log_path)

            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(),
                    (
                        "\n·•· Sample@test-uuid (Likelihood) \n"
                        "\t-> 01-02 03:04:05.678 - [INFO] >>> \n"
                        "hello world\n"
                        "stdout line\n"
                        "\n·•· Sample@test-uuid (Likelihood) \n"
                        "\t-> 01-02 03:04:05.678 - [WARNING] >>> \n"
                        "warn"
                    ),
                )

    def test_buffer_is_bounded(self):
        logger = BufferedSampleLogger(max_events=2, extra={"module": "Sample@test"})
        logger.info("one")
        logger.info("two")
        logger.info("three")
        self.assertEqual(len(logger.events), 2)
        self.assertEqual(logger.events[0].message, "two")
        self.assertEqual(logger.events[1].message, "three")


class FailureReplayLogTests(_LazySampleConfigMixin, unittest.TestCase):
    def test_failed_lazy_sample_replays_buffered_events(self):
        fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0, "z": 100.0})
            cfg = self._lazy_sample_cfg(tmpdir)
            sample.set_config(cfg)
            buffered = sample.info["logger"]
            self.assertIsInstance(buffered, BufferedSampleLogger)
            buffered._time_provider = lambda: fixed_dt  # type: ignore[attr-defined]

            sample.start()
            likelihood = LogLikelihood([{"name": "LogL_Z", "expression": "LogGauss(z, 100, 10)"}])
            likelihood.calculate(sample.observables, sample.info)
            sample.info["status"] = "Failed"

            save_dir = materialize_failure_artifacts(sample.info, error="workflow module failure")
            self.assertIsNotNone(save_dir)

            log_path = os.path.join(save_dir, "Sample_running.log")
            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()

            self.assertIn("Sample ->", text)
            self.assertIn("Evaluating   LogL_Z:", text)
            self.assertIn("Sample failed -> workflow module failure", text)
            self.assertIn("·•·", text)
            self.assertNotIn("Sample created into the Disk", text)

    def test_successful_lazy_sample_discards_buffer_without_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config(self._lazy_sample_cfg(tmpdir))
            sample.start()
            sample.info["logger"].info("transient event")
            sample.close()

            self.assertFalse(sample.info.get("_materialized"))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, sample.uuid)))

    def test_materialized_samples_keep_direct_file_logging(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            cfg = self._lazy_sample_cfg(tmpdir)
            cfg["sample_artifacts"] = "always"
            sample.set_config(cfg)

            self.assertIsInstance(sample.info["logger"], SampleLogger)
            self.assertTrue(os.path.exists(sample.info["run_log"]))
            with open(sample.info["run_log"], "r", encoding="utf-8") as handle:
                self.assertIn("Sample created into the Disk", handle.read())


class _FailOperasModule(Module):
    def __init__(self):
        super().__init__("FailOperas", selection=None)
        self.name = "FailOperas"
        self.type = "Operas"
        self.calls = 0

    def execute(self, observables, sample_info):
        self.calls += 1
        sample_info["logger"].info("operas reached failure point")
        raise RuntimeError("boom")


class _InMemoryDatabase:
    def add_data(self, row):
        return None


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class FailureReplayWorkflowTests(unittest.TestCase):
    def setUp(self):
        WorkerFactory._instance = None
        ModuleManager._instance = None

    def tearDown(self):
        WorkerFactory._instance = None
        ModuleManager._instance = None

    def test_factory_failure_replays_workflow_logs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            manager = ModuleManager()
            manager.set_logger(_NoopLogger())
            manager.set_config({"Sampling": {"ModuleFailurePolicy": "fail-fast"}})
            manager.workflow = {2: ["FailOperas"]}
            manager.module_pools = {"FailOperas": _FailOperasModule()}
            manager._database = _InMemoryDatabase()

            factory = WorkerFactory(manager)
            factory.configure(max_workers=1)

            sample = Sample({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            sample.start()

            with self.assertRaises(RuntimeError):
                factory._execute_workflow_tracked(sample.info)

            self.assertTrue(sample.info.get("_materialized"))
            log_path = os.path.join(sample.info["save_dir"], "Sample_running.log")
            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("operas reached failure point", text)
            self.assertIn("Sample failed -> boom", text)


if __name__ == "__main__":
    unittest.main()