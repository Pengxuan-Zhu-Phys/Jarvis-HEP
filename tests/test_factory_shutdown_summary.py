#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import os
import sys
import time
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.factory import WorkerFactory  # noqa: E402


class _CaptureLogger:
    def __init__(self):
        self.records = []

    def info(self, msg, *args, **kwargs):
        self.records.append(("INFO", str(msg)))

    def warning(self, msg, *args, **kwargs):
        self.records.append(("WARNING", str(msg)))

    def error(self, msg, *args, **kwargs):
        self.records.append(("ERROR", str(msg)))


class _FakeModuleManager:
    def execute_workflow(self, sample_info):
        return float(sample_info.get("value", 1.0))

    def _module_failure_policy(self):
        return "fail-fast"


class _RaisingModuleManager:
    def __init__(self, policy="fail-fast"):
        self._policy = str(policy)

    def execute_workflow(self, sample_info):
        raise RuntimeError("boom")

    def _module_failure_policy(self):
        return self._policy


class FactoryShutdownSummaryTests(unittest.TestCase):
    def setUp(self):
        WorkerFactory._instance = None

    def tearDown(self):
        try:
            factory = WorkerFactory._instance
            if factory is not None and getattr(factory, "executor", None) is not None:
                factory.shutdown(wait=True, cancel_futures=True)
        finally:
            WorkerFactory._instance = None

    def test_shutdown_emits_summary_log(self):
        logger = _CaptureLogger()
        manager = _FakeModuleManager()

        factory = WorkerFactory()
        factory.configure(module_manager=manager, max_workers=2)
        factory.set_logger(logger)

        futures = []
        for idx in range(3):
            futures.append(factory.submit_task({"uuid": f"u-{idx}", "value": idx + 1.0}))

        for fut in futures:
            self.assertIsInstance(fut, concurrent.futures.Future)
            self.assertIsNotNone(fut.result(timeout=2.0))

        factory.shutdown(wait=True, cancel_futures=False)

        summaries = [
            msg for level, msg in logger.records
            if level == "WARNING" and "Factory shutdown summary" in msg
        ]
        self.assertTrue(summaries, "Expected shutdown summary log in factory warnings")
        summary = summaries[-1]
        self.assertRegex(summary, r"submitted\s+\|\s+3\b")
        self.assertRegex(summary, r"ok\s+\|\s+3\b")
        self.assertRegex(summary, r"failed\s+\|\s+0\b")

    def test_warning_status_interval_stays_at_1000(self):
        logger = _CaptureLogger()
        manager = _FakeModuleManager()

        factory = WorkerFactory()
        factory.configure(module_manager=manager, max_workers=2)
        factory.set_logger(logger)

        factory._update_status_interval(1)
        self.assertEqual(factory._status_interval, 1000)
        factory._update_status_interval(1000)
        self.assertEqual(factory._status_interval, 1000)
        factory._update_status_interval(100_000)
        self.assertEqual(factory._status_interval, 1000)
        factory._update_status_interval(5_000_000)
        factory._update_status_interval(1)
        self.assertEqual(factory._status_interval, 1000)

        fut = factory.submit_task({"uuid": "pressure-001", "value": 1.0})
        self.assertIsNotNone(fut.result(timeout=2.0))

        factory.shutdown(wait=True, cancel_futures=False)

    def test_progress_logs_emit_info_every_100_and_warning_every_1000(self):
        logger = _CaptureLogger()
        manager = _FakeModuleManager()

        factory = WorkerFactory()
        factory.configure(module_manager=manager, max_workers=2)
        factory.set_logger(logger)

        for count in range(100, 1001, 100):
            with factory._status_lock:
                factory.task_count = count
            factory.print_status(task_count=count)

        submitted_logs = [
            (level, msg)
            for level, msg in logger.records
            if msg.startswith("Submitted ")
        ]
        info_logs = [msg for level, msg in submitted_logs if level == "INFO"]
        warning_logs = [msg for level, msg in submitted_logs if level == "WARNING"]

        self.assertEqual(len(info_logs), 10)
        self.assertEqual(len(warning_logs), 1)
        self.assertIn("Submitted 100 tasks, time for last 100 tasks:", info_logs[0])
        self.assertIn("Submitted 1000 tasks, time for last 1000 tasks:", warning_logs[0])

    def test_continue_policy_downgrades_sample_exception_without_killing_factory(self):
        logger = _CaptureLogger()
        manager = _RaisingModuleManager(policy="continue")

        factory = WorkerFactory()
        factory.configure(module_manager=manager, max_workers=2)
        factory.set_logger(logger)

        sample_info = {"uuid": "u-fail", "observables": {"seed": 1}}
        fut = factory.submit_task(sample_info)

        self.assertEqual(fut.result(timeout=2.0), 1.0)
        self.assertEqual(sample_info["status"], "Failed")

        factory.shutdown(wait=True, cancel_futures=False)

        self.assertEqual(factory.completed_ok, 0)
        self.assertEqual(factory.completed_failed, 1)
        errors = [
            msg for level, msg in logger.records
            if level == "ERROR" and "Workflow exception downgraded by ModuleFailurePolicy=continue" in msg
        ]
        self.assertTrue(errors)

    def test_low_throughput_mode_emits_minutely_summary(self):
        logger = _CaptureLogger()
        manager = _FakeModuleManager()

        factory = WorkerFactory()
        factory.configure(module_manager=manager, max_workers=2)
        factory.set_logger(logger)

        now_mono = time.monotonic()
        factory._factory_start_ts = now_mono - 301.0
        factory._last_low_throughput_log_ts = now_mono

        with factory._status_lock:
            factory.task_count = 120
            factory._last_status_count = 100
            factory.completed_ok = 100
            factory.completed_failed = 0
            factory.last_time_mark = time.time() - 61.0

        factory.print_status(task_count=120)
        summary_count_1 = len(
            [
                msg
                for level, msg in logger.records
                if level == "WARNING" and "Factory low-throughput summary ->" in msg
            ]
        )
        self.assertEqual(summary_count_1, 1)

        with factory._status_lock:
            factory.last_time_mark = time.time() - 1.0
        factory.print_status(task_count=121)
        summary_count_2 = len(
            [
                msg
                for level, msg in logger.records
                if level == "WARNING" and "Factory low-throughput summary ->" in msg
            ]
        )
        self.assertEqual(summary_count_2, 1)


if __name__ == "__main__":
    unittest.main()
