#!/usr/bin/env python3
"""WP-D5.1 / D9.4 op_count-gated monitor snapshot tests (explicit TaskFactory)."""

from __future__ import annotations

import time
import unittest
from types import SimpleNamespace
from typing import Any

from jarvishep2.dashboard import attach_reader
from jarvishep2.factory import TaskFactory, _MonitorLoop, _Watchdog
from jarvishep2.redis_queue import (
    INFLIGHT,
    PROC_ARCHIVER,
    PROC_CHILDREN,
    PROC_CORE,
    PROC_REDIS,
    PROC_WORKER,
    TASK_QUEUE,
    WORKER_STATUS,
    encode_payload,
    make_fakeredis_queue,
)
from jarvishep2.runtime.worker import Worker


class MonitorSnapshotTests(unittest.TestCase):
    def setUp(self) -> None:
        self.factory = TaskFactory()

    def tearDown(self) -> None:
        self.factory.shutdown(wait=False)

    def test_collaborators_are_composed(self) -> None:
        self.assertIsInstance(self.factory._monitor, _MonitorLoop)
        self.assertIsInstance(self.factory._watchdog, _Watchdog)

    def test_get_monitor_snapshot_does_not_touch_redis(self) -> None:
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory._monitor._snapshot = {
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
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue
        assert queue.r is not None

        factory._monitor._snapshot = factory._collect_latest_status()
        hgetall = queue.r.hgetall
        calls = {"count": 0}

        def _counting_hgetall(*args: Any, **kwargs: Any) -> Any:
            calls["count"] += 1
            return hgetall(*args, **kwargs)

        queue.r.hgetall = _counting_hgetall  # type: ignore[method-assign]
        factory._collect_latest_status()
        self.assertEqual(calls["count"], 0)

    def test_collect_latest_status_refreshes_sample_stats_on_op_count_bump(self) -> None:
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue

        factory._monitor._snapshot = factory._collect_latest_status()
        queue.submit_result({"uuid": "sample-1", "status": "Completed"})
        updated = factory._collect_latest_status()

        self.assertEqual(updated["sample_stats"].get("completed"), 1)
        self.assertGreaterEqual(updated["op_counts"]["sample"], 1)

    def test_factory_collect_status_issues_no_redis_writes(self) -> None:
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory._monitor._snapshot = factory._collect_latest_status()
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
            factory._monitor._snapshot = factory._collect_latest_status()
        self.assertEqual(writes["count"], 0)

    def test_get_run_metrics_projects_sample_stats(self) -> None:
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue
        queue.push_task(
            {
                "uuid": "task-1",
                "u_coords": [0.0],
                "execution_plan": [],
            }
        )
        queue.submit_result(
            {"uuid": "task-1", "status": "Completed", "observables": {"x": 1.0}}
        )

        metrics = factory.get_run_metrics()
        self.assertGreaterEqual(metrics["submitted"], 1)
        self.assertEqual(metrics["ok"], 1)
        # D9.4 honest None for untracked gauges
        self.assertIsNone(metrics["mean_active_workers"])
        self.assertIsNone(metrics["total_point_eval_sec"])

    def test_resume_scan_deletes_out_of_range_proc_and_inflight(self) -> None:
        queue = make_fakeredis_queue()
        queue.r.hset(PROC_CORE, mapping={"role": "core"})
        queue.r.hset(PROC_ARCHIVER, mapping={"role": "archiver"})
        queue.r.hset(PROC_REDIS, mapping={"role": "redis"})
        queue.r.hset(PROC_WORKER.format(id="0"), mapping={"role": "worker"})
        queue.r.hset(PROC_WORKER.format(id="99"), mapping={"role": "worker"})
        queue.r.hset(PROC_CHILDREN.format(id="77"), mapping={"role": "children"})
        queue.r.rpush(INFLIGHT.format(worker="0"), "in-range")
        queue.r.rpush(INFLIGHT.format(worker="99"), "leftover")
        queue.r.rpush(TASK_QUEUE, "transport")
        queue.r.hset(WORKER_STATUS.format(id="3"), mapping={"status": "busy"})

        result = queue.reconcile_resume_ephemeral(completed=12, failed=3)

        self.assertGreaterEqual(int(result["deleted_keys"]), 1)
        self.assertEqual(queue.r.hgetall(PROC_CORE), {})
        self.assertEqual(queue.read_proc_board("archiver"), {})
        self.assertEqual(queue.read_proc_board("redis"), {})
        self.assertEqual(queue.read_proc_board("worker", owner_id="0"), {})
        self.assertEqual(queue.read_proc_board("worker", owner_id="99"), {})
        self.assertEqual(queue.read_proc_board("children", owner_id="77"), {})
        self.assertEqual(int(queue.r.llen(INFLIGHT.format(worker="0"))), 0)
        self.assertEqual(int(queue.r.llen(INFLIGHT.format(worker="99"))), 0)
        self.assertEqual(int(queue.r.llen(TASK_QUEUE)), 0)
        self.assertEqual(queue.r.hgetall(WORKER_STATUS.format(id="3")), {})
        stats = queue.fetch_sample_stats()
        self.assertEqual(int(stats["completed"]), 12)
        self.assertEqual(int(stats["failed"]), 3)
        self.assertEqual(int(stats["running"]), 0)

    def test_snapshot_raw_contains_proc_core(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_proc_board("core", role="core", pid=42, scan_mode="running")
        queue.publish_proc_board("archiver", role="archiver", status="running")
        queue.publish_proc_board("redis", role="redis", status="running")
        queue.publish_proc_board("worker", owner_id="0", status="idle", pid=7)

        def _forbidden_scan(*_args: Any, **_kwargs: Any):
            raise AssertionError("snapshot_raw must not SCAN worker/proc keys")

        queue.r.scan_iter = _forbidden_scan  # type: ignore[method-assign]

        snapshot = queue.snapshot_raw()
        self.assertEqual(int(snapshot["proc_core"]["pid"]), 42)
        self.assertEqual(snapshot["proc_core"]["scan_mode"], "running")
        self.assertEqual(snapshot["proc_archiver"]["role"], "archiver")
        self.assertEqual(snapshot["proc_redis"]["role"], "redis")
        self.assertNotIn("proc_workers", snapshot)

        with_workers = queue.snapshot_raw(owner_ids=["0", "1"])
        self.assertEqual(with_workers["proc_workers"]["0"]["status"], "idle")
        self.assertEqual(with_workers["proc_workers"]["1"], {})
        self.assertNotIn("current_task", with_workers["proc_workers"]["0"])
        self.assertNotIn("u_coords", with_workers["proc_workers"]["0"])

    def test_snapshot_reader_issues_no_redis_writes(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_proc_board("core", role="core", scan_mode="running")
        queue.heartbeat("0", status="idle", pid=11, ts=1.0)
        assert queue.r is not None
        write_methods = ("set", "mset", "incr", "incrby", "hset", "hincrby", "rpush", "lpush", "delete", "hdel")
        originals = {name: getattr(queue.r, name) for name in write_methods}
        writes = {"count": 0}

        def _guard(original: Any):
            def guarded(*args: Any, **kwargs: Any) -> Any:
                writes["count"] += 1
                return original(*args, **kwargs)

            return guarded

        for method_name in write_methods:
            setattr(queue.r, method_name, _guard(originals[method_name]))

        view = attach_reader(redis=queue, owner_ids=["0"]).read()
        self.assertEqual(writes["count"], 0)
        self.assertEqual(view.proc_core.get("scan_mode"), "running")
        self.assertEqual(len(view.workers), 1)
        self.assertEqual(view.workers[0]["status"], "idle")
        self.assertNotIn("current_task", view.workers[0])
        self.assertNotIn("u_coords", view.workers[0])
        self.assertNotIn("execution_plan", view.workers[0])

    def test_heartbeat_status_has_no_current_task_blob(self) -> None:
        queue = make_fakeredis_queue()
        queue.r.hset(
            WORKER_STATUS.format(id="0"),
            mapping={"current_task": "stale-blob"},
        )
        queue.heartbeat(
            "0",
            status="busy",
            pid=1,
            ts=1.0,
            current_sample="uuid-1",
            current_task="blob",
            u_coords=[0.1],
            execution_plan=[],
        )
        status = queue.r.hgetall(WORKER_STATUS.format(id="0"))
        self.assertNotIn("current_task", status)
        self.assertNotIn("u_coords", status)
        self.assertNotIn("execution_plan", status)
        self.assertEqual(status["status"], "busy")
        board = queue.read_proc_board("worker", owner_id="0")
        self.assertNotIn("current_task", board)
        self.assertNotIn("u_coords", board)
        self.assertEqual(board["current_uuid"], "uuid-1")
        self.assertIsNone(queue.decode_heartbeat_task(status))

        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        worker._redis = queue
        worker._current_task = {
            "uuid": "uuid-1",
            "u_coords": [0.1],
            "execution_plan": [],
        }
        worker._current_sample_uuid = "uuid-1"
        worker._heartbeat("busy")
        stored = queue.r.hgetall(WORKER_STATUS.format(id="0"))
        self.assertNotIn("current_task", stored)
        self.assertIsNone(queue.decode_heartbeat_task(stored))

    def test_requeue_from_inflight_without_current_task(self) -> None:
        queue = make_fakeredis_queue()
        payload = {"uuid": "alive", "u_coords": [0.1], "execution_plan": []}
        queue.r.lpush(INFLIGHT.format(worker="0"), encode_payload(payload, codec="json"))
        queue.heartbeat("0", status="busy", pid=1, current_sample="alive", current_task="ignored")
        heartbeat = queue.r.hgetall(WORKER_STATUS.format(id="0"))
        self.assertNotIn("current_task", heartbeat)
        self.assertIsNone(queue.decode_heartbeat_task(heartbeat))

        watchdog = _Watchdog(SimpleNamespace(redis=queue))
        self.assertTrue(watchdog.requeue_in_flight_task(heartbeat, worker_id="0"))
        queued = [
            queue._decode_task_payload(raw)
            for raw in queue.r.lrange(TASK_QUEUE, 0, -1)
        ]
        self.assertEqual([item["uuid"] for item in queued], ["alive"])
        self.assertEqual(int(queued[0]["_retry_count"]), 1)
        self.assertEqual(int(queue.r.llen(INFLIGHT.format(worker="0"))), 0)

    def test_fetch_workers_redis_reads_proc_boards(self) -> None:
        factory = self.factory
        queue = make_fakeredis_queue()
        factory.redis = queue
        worker = SimpleNamespace(worker_id=0, pid=None, is_alive=lambda: False)
        factory.workers = [worker]  # type: ignore[assignment]
        queue.heartbeat("0", status="idle", pid=9, ts=1.0, current_task="blob")
        queue.publish_proc_board(
            "worker", owner_id="0", status="busy", pid=9, current_uuid="u1"
        )
        rows = factory._fetch_workers_redis()
        self.assertEqual(rows["0"]["status"], "busy")
        self.assertEqual(rows["0"]["current_uuid"], "u1")
        self.assertNotIn("current_task", rows["0"])


if __name__ == "__main__":
    unittest.main()
