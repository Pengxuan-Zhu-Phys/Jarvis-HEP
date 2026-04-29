#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import sys
import tempfile
import time
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.factory import WorkerFactory  # noqa: E402
from jarvishep.monitoring.run_summary import (  # noqa: E402
    RunSummaryCollector,
    RunSummaryRenderer,
    _utc_iso_from_epoch,
)


class _CaptureLogger:
    def __init__(self):
        self.records = []

    def info(self, msg, *args, **kwargs):
        self.records.append(("INFO", str(msg)))

    def warning(self, msg, *args, **kwargs):
        self.records.append(("WARNING", str(msg)))

    def error(self, msg, *args, **kwargs):
        self.records.append(("ERROR", str(msg)))


class _FakeFactory:
    def get_run_metrics(self):
        return {
            "submitted": 10,
            "ok": 8,
            "failed": 2,
            "configured_workers": 4,
            "peak_active_workers": 3,
            "mean_active_workers": 1.75,
            "total_point_eval_sec": 30.0,
            "completed_durations_sec": [1.0, 2.0, 3.0, 4.0],
            "retry_count": 3,
        }


class _FakeResourceSampler:
    def __init__(self):
        self.started = False
        self.stopped = False

    def start(self):
        self.started = True

    def stop(self):
        self.stopped = True

    def summary(self):
        return {
            "cpu_percent_mean": 12.5,
            "memory_rss_mb_peak": 256.0,
            "open_file_peak": 48,
        }


class _SleepyModuleManager:
    def execute_workflow(self, sample_info):
        time.sleep(float(sample_info.get("sleep", 0.01)))
        if sample_info.get("fail"):
            raise RuntimeError("boom")
        return float(sample_info.get("value", 1.0))


class RunSummaryCollectorTests(unittest.TestCase):
    def test_collector_builds_summary_and_writes_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            resource_sampler = _FakeResourceSampler()
            collector = RunSummaryCollector(
                output_dir=tmpdir,
                project_name="demo-project",
                sampler_name="Random",
                configured_workers=4,
                run_label="demo-scan",
                run_id="demo-run-001",
                resource_sampler=resource_sampler,
            )
            collector.attach_factory(_FakeFactory())
            collector.start()
            collector.record_external_command(duration_sec=4.5, ok=True)
            collector.record_external_command(duration_sec=0.5, ok=False, timed_out=True)
            collector.record_skipped_subtasks(2)

            collector._start_epoch = 100.0
            collector._start_time_iso = _utc_iso_from_epoch(100.0)
            collector._end_epoch = 160.0
            collector._end_time_iso = _utc_iso_from_epoch(160.0)
            collector._resource_summary = resource_sampler.summary()
            collector._finished = True

            summary = collector.build_summary()

            self.assertEqual(summary["run_id"], "demo-run-001")
            self.assertEqual(summary["project_name"], "demo-project")
            self.assertEqual(summary["sampler_name"], "Random")
            self.assertAlmostEqual(summary["wall_time_sec"], 60.0)
            self.assertEqual(summary["total_points_submitted"], 10)
            self.assertEqual(summary["total_points_finished"], 8)
            self.assertEqual(summary["total_points_failed"], 2)
            self.assertAlmostEqual(summary["success_rate"], 0.8)
            self.assertAlmostEqual(summary["throughput_points_per_min"], 8.0)
            self.assertEqual(summary["configured_workers"], 4)
            self.assertEqual(summary["peak_active_workers"], 3)
            self.assertAlmostEqual(summary["mean_active_workers"], 1.75)
            self.assertAlmostEqual(summary["avg_point_eval_sec"], 2.5)
            self.assertAlmostEqual(summary["median_point_eval_sec"], 2.5)
            self.assertAlmostEqual(summary["time_in_external_tools_sec"], 5.0)
            self.assertAlmostEqual(summary["time_in_framework_sec"], 25.0)
            self.assertAlmostEqual(summary["framework_overhead_fraction"], 25.0 / 30.0)
            self.assertEqual(summary["retry_count"], 3)
            self.assertEqual(summary["crashed_subtasks"], 1)
            self.assertEqual(summary["skipped_subtasks"], 2)
            self.assertAlmostEqual(summary["cpu_percent_mean"], 12.5)
            self.assertAlmostEqual(summary["memory_rss_mb_peak"], 256.0)
            self.assertEqual(summary["open_file_peak"], 48)

            renderer = RunSummaryRenderer()
            rendered = renderer.render(summary)
            self.assertIn("[Run Overview]", rendered)
            self.assertIn("[Execution Breakdown]", rendered)
            self.assertIn(
                "Tracked cumulative timing across completed points; not equal to wall-clock time.",
                rendered,
            )
            self.assertIn("Internal Workflow Time, Total (s)", rendered)
            self.assertIn("Internal Time Fraction", rendered)
            self.assertIn("Points Completed", rendered)
            self.assertIn("Mean CPU Usage (%)", rendered)
            self.assertIn("| Points Submitted", rendered)
            self.assertRegex(rendered, r"\| Run ID\s+\|\s+demo-run-001\s+\|")
            self.assertRegex(rendered, r"\|\s+80\.00%\s+\|")
            self.assertRegex(rendered, r"\| Skipped Tasks\s+\|\s+2\s+\|")

            outputs = renderer.write_outputs(summary, tmpdir, rendered_text=rendered)
            self.assertTrue(os.path.exists(outputs["json"]))
            self.assertTrue(os.path.exists(outputs["csv"]))
            self.assertTrue(os.path.exists(outputs["txt"]))

            with open(outputs["json"], "r", encoding="utf-8") as f1:
                payload = json.load(f1)
            self.assertEqual(payload["run_id"], "demo-run-001")
            self.assertEqual(payload["total_points_finished"], 8)

            with open(outputs["csv"], "r", encoding="utf-8", newline="") as f1:
                rows = list(csv.DictReader(f1))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["run_id"], "demo-run-001")
            self.assertEqual(rows[0]["configured_workers"], "4")

            with open(outputs["txt"], "r", encoding="utf-8") as f1:
                txt_payload = f1.read()
            self.assertEqual(txt_payload, rendered)


class WorkerFactoryRunMetricsTests(unittest.TestCase):
    def setUp(self):
        WorkerFactory._instance = None

    def tearDown(self):
        try:
            factory = WorkerFactory._instance
            if factory is not None and getattr(factory, "executor", None) is not None:
                factory.shutdown(wait=True, cancel_futures=True)
        finally:
            WorkerFactory._instance = None

    def test_factory_tracks_eval_durations_and_parallelism(self):
        logger = _CaptureLogger()
        factory = WorkerFactory()
        factory.configure(module_manager=_SleepyModuleManager(), max_workers=2)
        factory.set_logger(logger)
        factory.log_executor = None

        futures = [
            factory.submit_task({"uuid": "ok-1", "value": 1.0, "sleep": 0.05}),
            factory.submit_task({"uuid": "ok-2", "value": 2.0, "sleep": 0.05}),
            factory.submit_task({"uuid": "bad-1", "fail": True, "sleep": 0.02, "NAttempt": 2}),
        ]

        self.assertEqual(futures[0].result(timeout=2.0), 1.0)
        self.assertEqual(futures[1].result(timeout=2.0), 2.0)
        with self.assertRaises(RuntimeError):
            futures[2].result(timeout=2.0)

        metrics = factory.get_run_metrics()
        self.assertEqual(metrics["submitted"], 3)
        self.assertEqual(metrics["ok"], 2)
        self.assertEqual(metrics["failed"], 1)
        self.assertEqual(metrics["configured_workers"], 2)
        self.assertGreaterEqual(metrics["peak_active_workers"], 1)
        self.assertLessEqual(metrics["peak_active_workers"], 2)
        self.assertGreater(metrics["mean_active_workers"], 0.0)
        self.assertGreater(metrics["total_point_eval_sec"], 0.0)
        self.assertEqual(len(metrics["completed_durations_sec"]), 2)
        self.assertEqual(metrics["retry_count"], 1)

        factory.shutdown(wait=True, cancel_futures=False)


if __name__ == "__main__":
    unittest.main()
