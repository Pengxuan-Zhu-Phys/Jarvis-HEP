#!/usr/bin/env python3
"""WP-D5.1 op_count-gated monitor snapshot tests."""

from __future__ import annotations

import time
import unittest
from typing import Any

from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue


class MonitorSnapshotTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_get_monitor_snapshot_does_not_touch_redis(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory._snapshot = {
            "timestamp": time.time(),
            "workers": [],
            "workers_alive": 0,
            "workers_total": 0,
            "task_queue_length": 3,
            "archive_queue_length": 1,
            "sample_stats": {"completed": 2},
            "calculator_status": {},
            "op_counts": {"task": 1, "worker": 0, "calculator": 0, "sample": 2},
        }

        def _forbidden_redis_call(*_args: Any, **_kwargs: Any) -> None:
            raise AssertionError("get_monitor_snapshot must not call Redis")

        assert queue.r is not None
        queue.r.get = _forbidden_redis_call  # type: ignore[method-assign]
        queue.r.hgetall = _forbidden_redis_call  # type: ignore[method-assign]
        queue.r.llen = _forbidden_redis_call  # type: ignore[method-assign]

        snapshot = factory.get_monitor_snapshot()
        self.assertEqual(snapshot["task_queue_length"], 3)
        self.assertEqual(snapshot["sample_stats"]["completed"], 2)

    def test_collect_latest_status_gates_hgetall_on_idle_ticks(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        assert queue.r is not None

        factory._snapshot = factory._collect_latest_status()
        hgetall = queue.r.hgetall
        calls = {"count": 0}

        def _counting_hgetall(*args: Any, **kwargs: Any) -> Any:
            calls["count"] += 1
            return hgetall(*args, **kwargs)

        queue.r.hgetall = _counting_hgetall  # type: ignore[method-assign]
        factory._collect_latest_status()
        self.assertEqual(calls["count"], 0)

    def test_collect_latest_status_refreshes_sample_stats_on_op_count_bump(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue

        factory._snapshot = factory._collect_latest_status()
        queue.submit_result({"uuid": "sample-1", "status": "Completed"})
        updated = factory._collect_latest_status()

        self.assertEqual(updated["sample_stats"].get("completed"), 1)
        self.assertGreaterEqual(updated["op_counts"]["sample"], 1)

    def test_factory_collect_status_issues_no_redis_writes(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory._snapshot = factory._collect_latest_status()
        assert queue.r is not None

        write_methods = (
            "set",
            "mset",
            "incr",
            "incrby",
            "hset",
            "hincrby",
            "rpush",
            "lpush",
            "delete",
            "hdel",
        )
        originals = {name: getattr(queue.r, name) for name in write_methods}
        writes = {"count": 0}

        def _guard(name: str, original: Any) -> Any:
            def guarded(*args: Any, **kwargs: Any) -> Any:
                writes["count"] += 1
                return original(*args, **kwargs)

            return guarded

        for method_name in write_methods:
            setattr(queue.r, method_name, _guard(method_name, originals[method_name]))

        for _ in range(5):
            factory._snapshot = factory._collect_latest_status()
        self.assertEqual(writes["count"], 0)

    def test_get_run_metrics_projects_sample_stats(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        queue.push_task(
            {
                "uuid": "task-1",
                "u_coords": [0.0],
                "execution_plan": [],
            }
        )
        queue.submit_result({"uuid": "task-1", "status": "Completed", "observables": {"x": 1.0}})

        metrics = factory.get_run_metrics()
        self.assertGreaterEqual(metrics["submitted"], 1)
        self.assertEqual(metrics["ok"], 1)


if __name__ == "__main__":
    unittest.main()