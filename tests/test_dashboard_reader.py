#!/usr/bin/env python3
"""WP-D5.2 dashboard reader + run_summary tests."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from unittest import mock

from jarvishep2.client import run_monitor
from jarvishep2.dashboard import SnapshotReader, attach_reader, format_monitor_view
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import (
    CALC_STATUS,
    chain_feedback_queue,
    calc_status_busy_field,
    calc_status_free_field,
)
from jarvishep2.monitoring.run_summary import (
    RUN_SUMMARY_FIELD_ORDER,
    RunSummaryRenderer,
    build_run_summary,
    validate_run_summary,
)
from jarvishep2.redis_queue import make_fakeredis_queue


class DashboardReaderTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_snapshot_reader_maps_redis_keys(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        queue.push_task({"uuid": "u1", "u_coords": [0.1], "execution_plan": []})
        reader = attach_reader(redis=queue)
        view = reader.read()
        self.assertTrue(view.has_active_scan())
        self.assertEqual(view.queues["task_queue_length"], 1)
        self.assertGreaterEqual(view.op_counts.get("task", 0), 1)

    def test_snapshot_reader_reads_only_explicit_feedback_shards(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        assert queue.r is not None
        queue.r.rpush(chain_feedback_queue(0), "known-chain-result")
        queue.r.rpush(chain_feedback_queue(9), "unlisted-chain-result")
        view = attach_reader(redis=queue, feedback_chain_ids=[0, 1]).read()
        self.assertEqual(view.queues["feedback_shards"], {"0": 1, "1": 0})
        self.assertEqual(view.queues["feedback_queue_length"], 0)

    def test_worker_projection_keeps_file_operator_evidence(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        queue.heartbeat(
            "3",
            status="busy",
            pid=110,
            ts=1.0,
            board_ttl_sec=30,
            file_operation_pid=220,
            file_operation_pgid=220,
            file_operation_mode="process",
        )
        worker = attach_reader(redis=queue, owner_ids=["3"]).read().workers[0]
        self.assertEqual(worker["file_operation_pid"], 220)
        self.assertEqual(worker["file_operation_pgid"], 220)
        self.assertEqual(worker["file_operation_mode"], "process")
        self.assertEqual(worker["board_ttl_sec"], 30)

    def test_snapshot_reader_is_read_only(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        writes = {"count": 0}
        assert queue.r is not None
        for method_name in ("set", "incr", "hset", "rpush", "delete", "hincrby", "lpush"):
            original = getattr(queue.r, method_name)

            def _guard(name: str, original_fn):
                def wrapped(*args, **kwargs):
                    writes["count"] += 1
                    return original_fn(*args, **kwargs)

                return wrapped

            setattr(queue.r, method_name, _guard(method_name, original))
        attach_reader(redis=queue).read()
        self.assertEqual(writes["count"], 0)

    def test_snapshot_reader_uses_occupancy_when_want_is_off(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        queue.register_calc_pool("DemoCalc", 2)
        pack = queue.acquire_calc("DemoCalc", timeout=1)
        self.assertEqual(pack, "001")
        view = attach_reader(redis=queue).read()
        self.assertEqual(
            view.occupancy["DemoCalc"],
            {"free": 1, "busy": 1, "slots": 2},
        )
        self.assertEqual(view.calculators["DemoCalc"]["busy"], 1)
        status = queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 2)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 0)
        text = format_monitor_view(view)
        self.assertIn("calculator_occupancy:", text)
        self.assertIn("DemoCalc", text)

    def test_snapshot_reader_skips_calc_status_when_names_are_given(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        queue.register_calc_pool("DemoCalc", 2)
        queue.acquire_calc("DemoCalc", timeout=1)
        original = queue.r.hgetall

        def _guard(key, *args, **kwargs):
            if str(key) == CALC_STATUS:
                raise AssertionError("must not HGETALL hep:calculator:status")
            return original(key, *args, **kwargs)

        queue.r.hgetall = _guard  # type: ignore[method-assign]
        view = attach_reader(
            redis=queue,
            calculator_names=["DemoCalc"],
            calculator_slots={"DemoCalc": 2},
        ).read()
        self.assertEqual(view.occupancy["DemoCalc"]["busy"], 1)
        self.assertEqual(view.occupancy["DemoCalc"]["free"], 1)

    def test_run_monitor_exits_when_no_active_scan(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        code = run_monitor(redis=queue)
        self.assertEqual(code, 1)

    def test_validate_run_summary_requires_frozen_fields(self) -> None:
        with self.assertRaises(ValueError):
            validate_run_summary({"run_id": "x"})
        summary = build_run_summary(factory_metrics={"submitted": 1, "ok": 1, "failed": 0})
        validate_run_summary(summary)
        self.assertEqual(list(summary.keys()), list(RUN_SUMMARY_FIELD_ORDER))

    def test_build_run_summary_writes_schema_ordered_files(self) -> None:
        summary = build_run_summary(
            factory_metrics={"submitted": 4, "ok": 3, "failed": 1},
            project_name="eggbox",
            sampler_name="SamplingVirtial",
            run_id="run-001",
            start_epoch=100.0,
            end_epoch=160.0,
            configured_workers=2,
        )
        # 3 finished / 60s wall → 0.05 samples/sec, 20 sec amortized / sample
        self.assertAlmostEqual(float(summary["samples_per_sec"]), 0.05)
        self.assertAlmostEqual(float(summary["samples_per_min"]), 3.0)
        self.assertAlmostEqual(float(summary["avg_sample_sec"]), 20.0)
        self.assertAlmostEqual(float(summary["avg_point_eval_sec"]), 20.0)
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = RunSummaryRenderer().write_outputs(summary, tmpdir)
            self.assertTrue(os.path.exists(paths["json"]))
            self.assertTrue(os.path.exists(paths["csv"]))
            self.assertTrue(os.path.exists(paths["txt"]))
            with open(paths["txt"], encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("[Scan Performance]", text)
            self.assertIn("samples / sec", text)


class ClientMonitorTests(unittest.TestCase):
    def test_client_main_without_reference_lists_scan_choices(self) -> None:
        """Legacy ``Jarvis --monitor`` lists scan choices before attaching."""
        from jarvishep2.client import main

        queue = make_fakeredis_queue()
        queue.connect()
        with mock.patch("jarvishep2.redis_queue.RedisQueue", return_value=queue):
            with mock.patch.object(queue, "connect", return_value=None):
                with mock.patch.object(queue, "close", return_value=None):
                    # Documented legacy form is top-level --monitor, not bare YAML.
                    code = main(["--monitor"])
        self.assertEqual(code, 0)


if __name__ == "__main__":
    unittest.main()
