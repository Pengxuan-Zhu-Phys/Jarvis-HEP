#!/usr/bin/env python3
"""Unit tests for jarvishep2.redis_queue (WP-D0.2/D0.4, fakeredis)."""

from __future__ import annotations

import threading
import time
import unittest
from uuid import UUID, uuid4

import numpy as np

from jarvishep2.redis_queue import (
    ARCHIVE_QUEUE,
    CALC_BUSY_PACKS,
    CALC_FREE,
    CALC_STATUS,
    OP_COUNT,
    RESULTS,
    SAMPLE_STATS,
    TASK_QUEUE,
    WORKER_STATUS,
    RedisQueue,
    TaskValidationError,
    calc_busy_packs_key,
    calc_free_list_key,
    calc_status_busy_field,
    calc_status_free_field,
    decode_payload,
    encode_payload,
    make_fakeredis_queue,
)


def _minimal_task(**overrides) -> dict:
    payload = {"uuid": str(uuid4()), "u_coords": [0.1]}
    payload.update(overrides)
    return payload


class RedisQueueKeyNamespaceTests(unittest.TestCase):
    def test_exact_key_strings_match_design(self):
        self.assertEqual(TASK_QUEUE, "hep:task_queue")
        self.assertEqual(CALC_FREE, "calc:free:{name}")
        self.assertEqual(CALC_BUSY_PACKS, "calc:busy:{name}")
        self.assertEqual(ARCHIVE_QUEUE, "hep:archive_queue")
        self.assertEqual(WORKER_STATUS, "hep:worker:status:{id}")
        self.assertEqual(CALC_STATUS, "hep:calculator:status")
        self.assertEqual(SAMPLE_STATS, "hep:sample:stats")
        self.assertEqual(RESULTS, "hep:results:{uuid}")
        self.assertEqual(OP_COUNT, "hep:{kind}:op_count")

    def test_calc_free_list_key_returns_pool_key_not_status_field(self):
        self.assertEqual(calc_free_list_key("DemoCalc"), "calc:free:DemoCalc")
        self.assertEqual(calc_status_free_field("DemoCalc"), "DemoCalc:free")
        self.assertEqual(calc_status_busy_field("DemoCalc"), "DemoCalc:busy")
        self.assertNotEqual(calc_free_list_key("DemoCalc"), calc_status_free_field("DemoCalc"))


class RedisQueueTests(unittest.TestCase):
    def setUp(self):
        self.queue = make_fakeredis_queue(codec="json")

    def test_task_round_trip_with_numpy_and_uuid(self):
        task = {
            "uuid": str(uuid4()),
            "u_coords": np.array([0.1, 0.2], dtype=np.float64),
            "meta_id": UUID("12345678-1234-5678-1234-567812345678"),
        }
        self.queue.push_task(task)
        pulled = self.queue.pull_task(timeout=1)
        self.assertIsNotNone(pulled)
        assert pulled is not None
        self.assertEqual(pulled["uuid"], task["uuid"])
        self.assertEqual(pulled["u_coords"], [0.1, 0.2])
        self.assertEqual(pulled["meta_id"], str(task["meta_id"]))

    def test_push_task_rejects_malformed_payload(self):
        with self.assertRaises(TaskValidationError):
            self.queue.push_task({"u_coords": [1.0]})
        with self.assertRaises(TaskValidationError):
            self.queue.push_task({"uuid": "abc"})
        with self.assertRaises(TaskValidationError):
            self.queue.push_many_tasks([{"uuid": "only-uuid"}])

    def test_push_many_tasks_increments_op_count_once(self):
        tasks = [_minimal_task() for _ in range(5)]
        self.queue.push_many_tasks(tasks)
        self.assertEqual(self.queue.get_op_count("task"), 5)
        self.assertEqual(self.queue.r.llen(TASK_QUEUE), 5)

    def test_push_task_increments_op_count(self):
        self.queue.push_task(_minimal_task(uuid="a"))
        self.assertEqual(self.queue.get_op_count("task"), 1)
        self.queue.push_task(_minimal_task(uuid="b"))
        self.assertEqual(self.queue.get_op_count("task"), 2)

    def test_calc_pool_cap_and_distinct_pack_ids(self):
        self.queue.register_calc_pool("DemoCalc", 2)
        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 2)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 0)
        self.assertEqual(self.queue.r.llen(calc_free_list_key("DemoCalc")), 2)

        pack_a = self.queue.acquire_calc("DemoCalc", timeout=1)
        pack_b = self.queue.acquire_calc("DemoCalc", timeout=1)
        pack_c = self.queue.acquire_calc("DemoCalc", timeout=1)

        self.assertIsNotNone(pack_a)
        self.assertIsNotNone(pack_b)
        self.assertIsNone(pack_c)
        self.assertNotEqual(pack_a, pack_b)

        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 0)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 2)
        self.assertEqual(self.queue.get_op_count("calculator"), 2)

        self.queue.release_calc("DemoCalc", pack_a)
        pack_d = self.queue.acquire_calc("DemoCalc", timeout=1)
        self.assertIsNotNone(pack_d)

        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 0)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 2)

    def test_release_calc_requires_known_pack_id(self):
        self.queue.register_calc_pool("DemoCalc", 1)
        pack = self.queue.acquire_calc("DemoCalc", timeout=1)
        assert pack is not None
        with self.assertRaises(ValueError):
            self.queue.release_calc("DemoCalc", "not-a-real-pack")
        self.queue.release_calc("DemoCalc", pack)

    def test_concurrent_acquire_release_keeps_counts_consistent(self):
        self.queue.register_calc_pool("DemoCalc", 2)
        errors: list[BaseException] = []

        def worker() -> None:
            try:
                for _ in range(20):
                    pack = self.queue.acquire_calc("DemoCalc", timeout=2)
                    if pack is None:
                        continue
                    time.sleep(0.001)
                    self.queue.release_calc("DemoCalc", pack)
            except BaseException as exc:  # pragma: no cover - surfaced via errors list
                errors.append(exc)

        threads = [threading.Thread(target=worker) for _ in range(4)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()

        self.assertEqual(errors, [])
        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 2)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 0)
        self.assertEqual(self.queue.r.llen(calc_free_list_key("DemoCalc")), 2)
        self.assertEqual(self.queue.r.hlen(calc_busy_packs_key("DemoCalc")), 0)

    def test_submit_result_updates_archive_queue_and_stats(self):
        info = {"uuid": "done-1", "status": "Completed", "LogL": -1.5}
        self.queue.submit_result(info)
        self.assertEqual(self.queue.get_op_count("sample"), 1)
        self.assertEqual(int(self.queue.r.hget(SAMPLE_STATS, "completed")), 1)

        pulled = self.queue.pull_result(timeout=1)
        self.assertEqual(pulled, info)
        self.assertEqual(self.queue.r.llen(ARCHIVE_QUEUE), 0)

    def test_submit_result_rejects_missing_uuid(self):
        with self.assertRaises(TaskValidationError):
            self.queue.submit_result({"status": "Completed"})

    def test_heartbeat_preserves_numeric_types(self):
        self.queue.heartbeat("worker-1", status="idle", pid=1234, load=0.75)
        stored = self.queue.r.hgetall(WORKER_STATUS.format(id="worker-1"))
        self.assertEqual(stored["status"], "idle")
        self.assertEqual(int(stored["pid"]), 1234)
        self.assertAlmostEqual(float(stored["load"]), 0.75)

        snapshot = self.queue.snapshot_raw()
        self.assertEqual(snapshot["op_counts"]["worker"], 1)

    def test_snapshot_raw(self):
        self.queue.push_task(_minimal_task(uuid="t1"))
        snapshot = self.queue.snapshot_raw()
        self.assertEqual(snapshot["task_queue_length"], 1)
        self.assertEqual(snapshot["op_counts"]["task"], 1)

    def test_codec_round_trip_json(self):
        payload = {"arr": np.array([1, 2, 3])}
        encoded = encode_payload(payload, codec="json")
        decoded = decode_payload(encoded, codec="json")
        self.assertEqual(decoded["arr"], [1, 2, 3])

    def test_incr_op_rejects_invalid_kind(self):
        with self.assertRaises(ValueError):
            self.queue.incr_op("invalid")


class RedisQueueMsgpackTests(unittest.TestCase):
    def test_msgpack_codec_when_available(self):
        try:
            import msgpack  # noqa: F401
        except ImportError:
            self.skipTest("msgpack not installed")

        queue = make_fakeredis_queue(codec="msgpack")
        task = {"uuid": "x", "u_coords": np.array([1.0])}
        queue.push_task(task)
        pulled = queue.pull_task(timeout=1)
        self.assertEqual(pulled["uuid"], "x")
        self.assertEqual(pulled["u_coords"], [1.0])


if __name__ == "__main__":
    unittest.main()