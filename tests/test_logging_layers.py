#!/usr/bin/env python3
"""Unit tests for jarvishep2 two-layer logging (WP-D0.3)."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from datetime import datetime

from jarvishep2.logging import (
    BufferedSampleLogger,
    SampleLogger,
    get_jarvis_logger,
    setup_jarvis_logging,
    shutdown_jarvis_logging,
)
from jarvishep2.sample import Sample, materialize_failure_artifacts


class NoLoguruGuardTests(unittest.TestCase):
    def test_hot_path_has_no_loguru_imports(self):
        root = Path(__file__).resolve().parents[1] / "jarvishep2"
        hits = []
        for path in root.rglob("*.py"):
            text = path.read_text(encoding="utf-8")
            if "import loguru" in text or "from loguru" in text:
                hits.append(str(path.relative_to(root)))
        self.assertEqual(hits, [])


class TopLevelLoggingTests(unittest.TestCase):
    def tearDown(self):
        shutdown_jarvis_logging()

    def test_contract_format_includes_timestamp_level_name_and_context(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = setup_jarvis_logging(
                log_dir=tmpdir,
                role="worker",
                console=False,
                use_queue=True,
            )
            logger = get_jarvis_logger("worker").bind(worker_id="w-1")
            logger.info("start sample", extra={"sample_uuid": "abc-123"})

            shutdown_jarvis_logging()
            time.sleep(0.05)

            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()

            self.assertIn("INFO", text)
            self.assertIn("jarvis_hep.worker", text)
            self.assertIn("start sample", text)
            self.assertIn("worker_id=w-1", text)
            self.assertIn("sample_uuid=abc-123", text)
            self.assertRegex(text, r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}")

    def test_two_layers_keep_summary_separate_from_sample_detail(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = setup_jarvis_logging(
                log_dir=tmpdir,
                role="worker",
                console=False,
                use_queue=True,
            )
            top = get_jarvis_logger("worker").bind(worker_id="w-1")
            top.info("start sample", extra={"sample_uuid": "abc-123"})

            sample_log_path = os.path.join(tmpdir, "samples", "abc-123.log")
            sample_logger = SampleLogger.open(sample_log_path, module="Sample@abc-123")
            sample_logger.info("detailed calculator trace")
            sample_logger.close()

            shutdown_jarvis_logging()
            time.sleep(0.05)

            with open(log_path, "r", encoding="utf-8") as handle:
                top_text = handle.read()
            with open(sample_log_path, "r", encoding="utf-8") as handle:
                sample_text = handle.read()

            self.assertIn("start sample", top_text)
            self.assertNotIn("detailed calculator trace", top_text)
            self.assertIn("detailed calculator trace", sample_text)
            self.assertNotIn("start sample", sample_text)

    def test_bind_returns_new_adapter_without_mutating_parent(self):
        shutdown_jarvis_logging()
        parent = get_jarvis_logger("factory")
        child = parent.bind(worker_id="w-2")
        self.assertIsNot(child, parent)
        self.assertEqual(parent.extra, {})
        self.assertEqual(child.extra, {"worker_id": "w-2"})


class SampleLoggingReplayTests(unittest.TestCase):
    def test_child_logger_binds_via_logger_name_without_materialization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            child = sample.child_logger(module=f"{sample.info['logger_name']} (Likelihood)")
            assert child is not None
            child.info("bound without materialize")
            self.assertFalse(sample.info["_materialized"])
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

    def test_failure_materialization_replays_buffered_events(self):
        fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0, "z": 100.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            buffered = sample.info["logger"]
            self.assertIsInstance(buffered, BufferedSampleLogger)
            buffered._time_provider = lambda: fixed_dt  # type: ignore[attr-defined]

            sample.start()
            sample.info["logger"].info("Evaluating   LogL_Z: value computed")
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
            sample = Sample.from_params({"x": 1.0})
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
            sample.info["logger"].info("transient event")
            sample.close()

            self.assertFalse(sample.info.get("_materialized"))
            self.assertFalse(os.path.exists(os.path.join(tmpdir, sample.uuid)))

    def test_no_logger_on_task_dict_wire(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            wire = sample.to_task_dict()
            self.assertNotIn("logger", wire)

            rebuilt = Sample.from_task_dict(wire)
            self.assertIsNone(rebuilt._logger)
            self.assertEqual(rebuilt.info, {})


class QueueLoggingPerformanceTests(unittest.TestCase):
    def tearDown(self):
        shutdown_jarvis_logging()

    def test_queue_handler_logging_completes_quickly(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            setup_jarvis_logging(
                log_dir=tmpdir,
                role="perf",
                console=False,
                use_queue=True,
            )
            logger = get_jarvis_logger("perf").bind(phase="bulk")
            started = time.perf_counter()
            for idx in range(10_000):
                logger.info("bulk event", extra={"idx": idx})
            elapsed = time.perf_counter() - started
            shutdown_jarvis_logging()
            self.assertLess(elapsed, 5.0)


class SubprocessLoggingSmokeTests(unittest.TestCase):
    def test_setup_is_spawn_safe_without_loguru(self):
        project_root = str(Path(__file__).resolve().parents[1])
        script = f"""
import os, sys, tempfile
sys.path.insert(0, {project_root!r})
from jarvishep2.logging import setup_jarvis_logging, get_jarvis_logger, shutdown_jarvis_logging
tmpdir = tempfile.mkdtemp()
setup_jarvis_logging(log_dir=tmpdir, role="spawn", console=False, use_queue=True)
get_jarvis_logger("spawn").info("spawn ok", extra={{"pid": os.getpid()}})
shutdown_jarvis_logging()
print("OK")
"""
        proc = subprocess.run(
            [sys.executable, "-c", script],
            cwd=project_root,
            text=True,
            capture_output=True,
            timeout=30,
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        self.assertIn("OK", proc.stdout)


if __name__ == "__main__":
    unittest.main()