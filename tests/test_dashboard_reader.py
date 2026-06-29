#!/usr/bin/env python3
"""WP-D5.2 dashboard reader + run_summary tests."""

from __future__ import annotations

import os
import sys
import tempfile
import unittest
from unittest import mock

from jarvishep2.client import run_monitor
from jarvishep2.dashboard import SnapshotReader, attach_reader
from jarvishep2.factory import TaskFactory
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

    def test_snapshot_reader_is_read_only(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        writes = {"count": 0}
        assert queue.r is not None
        for method_name in ("set", "incr", "hset", "rpush", "delete"):
            original = getattr(queue.r, method_name)

            def _guard(name: str, original_fn):
                def wrapped(*args, **kwargs):
                    writes["count"] += 1
                    return original_fn(*args, **kwargs)

                return wrapped

            setattr(queue.r, method_name, _guard(method_name, original))
        attach_reader(redis=queue).read()
        self.assertEqual(writes["count"], 0)

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
        with tempfile.TemporaryDirectory() as tmpdir:
            paths = RunSummaryRenderer().write_outputs(summary, tmpdir)
            self.assertTrue(os.path.exists(paths["json"]))
            self.assertTrue(os.path.exists(paths["csv"]))
            self.assertTrue(os.path.exists(paths["txt"]))


class ClientMonitorTests(unittest.TestCase):
    def test_client_main_no_scan_exit_code(self) -> None:
        from jarvishep2.client import main

        queue = make_fakeredis_queue()
        queue.connect()
        with mock.patch("jarvishep2.client.RedisQueue", return_value=queue):
            with mock.patch.object(queue, "connect", return_value=None):
                with mock.patch.object(queue, "close", return_value=None):
                    code = main(["task.yaml", "--monitor"])
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()