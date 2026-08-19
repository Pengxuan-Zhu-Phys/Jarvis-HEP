#!/usr/bin/env python3
"""Unit tests for jarvishep2 two-layer logging (WP-D0.3)."""

from __future__ import annotations

from contextlib import redirect_stderr
import io
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
    JARVIS_HEP_LOG_DOMAIN,
    SampleLogger,
    component_log_path,
    get_jarvis_logger,
    resolve_record_component,
    scan_logs_dir,
    setup_jarvis_logging,
    shutdown_jarvis_logging,
)
import logging
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
            logger = get_jarvis_logger("worker", worker_id=1).bind(worker_id="w-1")
            logger.info("start sample", extra={"sample_uuid": "abc-123"})
            logger.debug("retained debug detail")

            shutdown_jarvis_logging()
            time.sleep(0.05)

            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()

            self.assertIn("INFO", text)
            self.assertIn("·•·", text)
            self.assertIn("Jarvis-HEP.Worker.01", text)
            self.assertIn("start sample", text)
            self.assertIn("worker_id=w-1", text)
            self.assertIn("sample_uuid=abc-123", text)
            self.assertIn("retained debug detail", text)
            # V1-style timestamp: MM-DD HH:mm:ss.SSS
            self.assertRegex(text, r"\d{2}-\d{2} \d{2}:\d{2}:\d{2}\.\d{3}")

    def test_default_console_level_is_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            setup_jarvis_logging(log_dir=tmpdir, role="worker", console=True, use_queue=False)
            logger = logging.getLogger(JARVIS_HEP_LOG_DOMAIN)
            console_handlers = [
                handler for handler in logger.handlers if isinstance(handler, logging.StreamHandler)
            ]
            self.assertTrue(console_handlers)
            self.assertEqual(console_handlers[0].level, logging.WARNING)

    def test_v1_style_formatter_raw_and_module_bind(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = setup_jarvis_logging(
                log_dir=tmpdir,
                role="core",
                console=False,
                use_queue=False,
            )
            logger = get_jarvis_logger("core").bind(module="Jarvis-HEP")
            logger.warning("structured line")
            logger.bind(raw=True).warning("RAW BANNER\nline2")
            shutdown_jarvis_logging()

            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("·•· Jarvis-HEP", text)
            self.assertIn("[WARNING]", text)
            self.assertIn("structured line", text)
            self.assertIn("RAW BANNER\nline2", text)
            self.assertNotIn("·•·", text.split("RAW BANNER", 1)[1][:20])

    def test_scan_log_path_layout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            logs = scan_logs_dir(tmpdir, "myscan")
            self.assertEqual(logs, os.path.join(os.path.abspath(tmpdir), "logs", "myscan"))
            core_path = component_log_path(logs, "core")
            self.assertTrue(core_path.endswith(os.path.join("myscan", "core.log")))
            worker_path = component_log_path(logs, "worker", worker_id=3)
            self.assertTrue(worker_path.endswith(os.path.join("myscan", "worker-03.log")))
            archiver_path = component_log_path(logs, "archiver")
            self.assertTrue(archiver_path.endswith(os.path.join("myscan", "archiver.log")))

            path = setup_jarvis_logging(
                scan_logs_dir=logs,
                component="core",
                role="core",
                console=False,
                use_queue=False,
            )
            self.assertEqual(path, os.path.abspath(core_path))
            self.assertTrue(os.path.isfile(path))
            logger = get_jarvis_logger("core", module="Jarvis-HEP")
            logger.info("control line")
            shutdown_jarvis_logging()
            with open(path, encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("·•· Jarvis-HEP", text)
            self.assertIn("control line", text)

    def test_worker_component_log_under_scan_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            logs = scan_logs_dir(tmpdir, "egg")
            path = setup_jarvis_logging(
                scan_logs_dir=logs,
                component="worker",
                worker_id=1,
                console=False,
                use_queue=False,
            )
            self.assertEqual(
                path,
                os.path.abspath(os.path.join(logs, "worker-01.log")),
            )
            log = get_jarvis_logger("worker", worker_id=1)
            log.info("worker ready")
            shutdown_jarvis_logging()
            with open(path, encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("·•· Jarvis-HEP.Worker.01", text)
            self.assertIn("worker ready", text)

    def test_build_worker_config_always_sets_scan_logs_dir(self):
        """Worker blueprint must stamp logs/<scan>/ so processes never use cwd jarvis_worker_PID."""
        from jarvishep2.worker_config import build_worker_config

        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "outputs", "ScanA")
            os.makedirs(out)
            wc = build_worker_config(
                {
                    "scan_name": "ScanA",
                    "task_root": tmpdir,
                    "Sampling": {"Method": "Bridson"},
                    "Operas": {"Modules": []},
                },
                task_result_dir=out,
            )
            expected = scan_logs_dir(tmpdir, "ScanA")
            self.assertEqual(wc.get("logs_dir"), expected)
            self.assertTrue(wc["logs_dir"].endswith(os.path.join("logs", "ScanA")))

    def test_resolve_record_component_by_logger_name_and_module(self) -> None:
        rec = logging.LogRecord(
            name="jarvis_hep.factory.watchdog",
            level=logging.INFO,
            pathname=__file__,
            lineno=1,
            msg="x",
            args=(),
            exc_info=None,
        )
        self.assertEqual(resolve_record_component(rec), "factory")
        rec2 = logging.LogRecord(
            name="jarvis_hep.sampler.dynesty",
            level=logging.INFO,
            pathname=__file__,
            lineno=1,
            msg="y",
            args=(),
            exc_info=None,
        )
        rec2.jarvis_module = "MultiNest.Inner"  # type: ignore[attr-defined]
        self.assertEqual(resolve_record_component(rec2), "sampler")
        rec3 = logging.LogRecord(
            name="jarvis_hep.core",
            level=logging.INFO,
            pathname=__file__,
            lineno=1,
            msg="z",
            args=(),
            exc_info=None,
        )
        self.assertEqual(resolve_record_component(rec3), "core")

    def test_control_multi_sink_splits_factory_and_sampler(self) -> None:
        """Control process multi-sink: core/factory/sampler/archiver/datarecorder."""
        with tempfile.TemporaryDirectory() as tmpdir:
            logs = scan_logs_dir(tmpdir, "split")
            setup_jarvis_logging(
                role="core",
                component="core",
                scan_logs_dir=logs,
                console=False,
                use_queue=False,
                multi_sink=True,
            )
            get_jarvis_logger("core").info("control-only line")
            get_jarvis_logger("factory").info("factory-only line")
            get_jarvis_logger("sampler.dynesty").info("sampler-only line")
            get_jarvis_logger(
                "sampler.multinest.inner",
                module="Jarvis-HEP.Sampler.MultiNest.Inner",
            ).info("progress-only line")
            get_jarvis_logger("archiver").info("archiver-only line")
            get_jarvis_logger("datarecorder").info(
                "DATABASE progress ->\n\trows written\t-> 10"
            )
            shutdown_jarvis_logging()

            core_text = Path(os.path.join(logs, "core.log")).read_text(encoding="utf-8")
            factory_text = Path(os.path.join(logs, "factory.log")).read_text(
                encoding="utf-8"
            )
            sampler_text = Path(os.path.join(logs, "sampler.log")).read_text(
                encoding="utf-8"
            )
            archiver_text = Path(os.path.join(logs, "archiver.log")).read_text(
                encoding="utf-8"
            )
            data_text = Path(os.path.join(logs, "datarecorder.log")).read_text(
                encoding="utf-8"
            )

            self.assertIn("control-only line", core_text)
            self.assertNotIn("factory-only line", core_text)
            self.assertNotIn("sampler-only line", core_text)
            self.assertNotIn("progress-only line", core_text)
            self.assertNotIn("DATABASE progress", core_text)

            self.assertIn("factory-only line", factory_text)
            self.assertNotIn("control-only line", factory_text)

            self.assertIn("sampler-only line", sampler_text)
            self.assertIn("progress-only line", sampler_text)
            self.assertNotIn("control-only line", sampler_text)

            self.assertIn("archiver-only line", archiver_text)
            self.assertNotIn("control-only line", archiver_text)
            self.assertNotIn("DATABASE progress", archiver_text)

            self.assertIn("DATABASE progress", data_text)
            self.assertIn("Jarvis-HEP.DataRecorder", data_text)
            self.assertNotIn("control-only line", data_text)
            self.assertIn("Jarvis-HEP.Factory", factory_text)
            self.assertIn("Jarvis-HEP.Sampler", sampler_text)
            self.assertIn("Jarvis-HEP.Archiver", archiver_text)

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

    def test_sample_terminal_default_threshold_is_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            terminal = io.StringIO()
            with redirect_stderr(terminal):
                setup_jarvis_logging(
                    log_dir=tmpdir,
                    role="worker",
                    console=True,
                    console_level="WARNING",
                    use_queue=False,
                )
                sample_logger = SampleLogger.open(
                    os.path.join(tmpdir, "sample.log"),
                    module="Sample@abc-123",
                    extra={
                        "to_console": True,
                        "sample_console_level": "ERROR",
                    },
                )
                sample_logger.info("sample info hidden")
                sample_logger.warning("sample warning hidden")
                sample_logger.error("sample error visible")
                sample_logger.close()
                shutdown_jarvis_logging()

            text = terminal.getvalue()
            self.assertNotIn("sample info hidden", text)
            self.assertNotIn("sample warning hidden", text)
            self.assertIn("sample error visible", text)

    def test_sample_terminal_check_threshold_includes_debug(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            terminal = io.StringIO()
            with redirect_stderr(terminal):
                setup_jarvis_logging(
                    log_dir=tmpdir,
                    role="worker",
                    console=True,
                    console_level="DEBUG",
                    use_queue=False,
                )
                sample_logger = SampleLogger.open(
                    os.path.join(tmpdir, "sample.log"),
                    module="Sample@check-123",
                    extra={
                        "to_console": True,
                        "sample_console_level": "DEBUG",
                    },
                )
                sample_logger.debug("sample debug visible")
                sample_logger.close()
                shutdown_jarvis_logging()

            self.assertIn("sample debug visible", terminal.getvalue())

    def test_sample_terminal_raw_event_does_not_repeat_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            terminal = io.StringIO()
            with redirect_stderr(terminal):
                setup_jarvis_logging(
                    log_dir=tmpdir,
                    role="worker",
                    console=True,
                    console_level="DEBUG",
                    use_queue=False,
                )
                sample_logger = SampleLogger.open(
                    os.path.join(tmpdir, "sample.log"),
                    module="Sample@raw-123",
                    extra={
                        "to_console": True,
                        "sample_console_level": "DEBUG",
                    },
                )
                sample_logger.info("structured sample line")
                sample_logger.bind(raw=True).info("stdout line")
                sample_logger.close()
                shutdown_jarvis_logging()

            text = terminal.getvalue()
            self.assertIn("structured sample line", text)
            self.assertIn("stdout line", text)
            self.assertIn("·•· Sample@raw-123", text)
            raw_tail = text.split("structured sample line", 1)[1]
            self.assertIn("stdout line", raw_tail)
            self.assertNotIn("·•· Sample@raw-123", raw_tail)
            self.assertNotIn("stdout line\n\n", text)

    def test_bind_returns_new_adapter_without_mutating_parent(self):
        shutdown_jarvis_logging()
        parent = get_jarvis_logger("factory")  # default ·•· label "Jarvis-HEP.Factory"

        child = parent.bind(worker_id="w-2")
        self.assertIsNot(child, parent)
        self.assertEqual(parent.extra, {"jarvis_module": "Jarvis-HEP.Factory"})
        self.assertEqual(
            child.extra,
            {"jarvis_module": "Jarvis-HEP.Factory", "worker_id": "w-2"},
        )
        # V1-style module= is remapped onto jarvis_module (LogRecord-safe).
        labeled = parent.bind(module="Jarvis-HEP")
        self.assertEqual(labeled.extra.get("jarvis_module"), "Jarvis-HEP")
        self.assertNotIn("module", labeled.extra)
        self.assertEqual(parent.extra.get("jarvis_module"), "Jarvis-HEP.Factory")


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
