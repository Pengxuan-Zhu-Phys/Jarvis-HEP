#!/usr/bin/env python3
"""D26.1 proc-board mixin, dual-write heartbeat, and fresh-reset deletes."""

from __future__ import annotations

import unittest

from jarvishep2.redis_queue import (
    ARCHIVER_BOARD_TTL_SEC,
    INFLIGHT,
    PROC_ARCHIVER,
    PROC_BOARD_TTL_SEC,
    PROC_CHILDREN,
    PROC_CORE,
    PROC_WORKER,
    WORKER_STATUS,
    make_fakeredis_queue,
)


class ProcBoardMixinTests(unittest.TestCase):
    def setUp(self) -> None:
        self.queue = make_fakeredis_queue(codec="json")

    def test_publish_overlay_last_write_wins(self) -> None:
        self.queue.publish_proc_board("core", role="core", pid=1, host="a")
        self.queue.publish_proc_board("core", pid=2, status="running")
        board = self.queue.read_proc_board("core")
        self.assertEqual(board["role"], "core")
        self.assertEqual(board["host"], "a")
        self.assertEqual(int(board["pid"]), 2)
        self.assertEqual(board["status"], "running")

    def test_unknown_role_raises(self) -> None:
        with self.assertRaises(ValueError):
            self.queue.publish_proc_board("monitor", pid=1)
        with self.assertRaises(ValueError):
            self.queue.read_proc_board("factory")

    def test_worker_and_children_require_owner_id(self) -> None:
        with self.assertRaises(ValueError):
            self.queue.publish_proc_board("worker", status="idle")
        with self.assertRaises(ValueError):
            self.queue.publish_proc_board("children", owner_id="")
        with self.assertRaises(ValueError):
            self.queue.read_proc_board("worker")

    def test_missing_key_reads_as_empty_dict(self) -> None:
        self.assertEqual(self.queue.read_proc_board("core"), {})
        self.assertEqual(
            self.queue.read_proc_board("worker", owner_id="0"), {}
        )

    def test_publish_sets_expire_ttl(self) -> None:
        self.queue.publish_proc_board("core", role="core")
        core_ttl = int(self.queue.r.ttl(PROC_CORE))
        self.assertGreater(core_ttl, 0)
        self.assertLessEqual(core_ttl, PROC_BOARD_TTL_SEC)

        self.queue.publish_proc_board("archiver", role="archiver")
        archiver_ttl = int(self.queue.r.ttl(PROC_ARCHIVER))
        self.assertGreater(archiver_ttl, 0)
        self.assertLessEqual(archiver_ttl, ARCHIVER_BOARD_TTL_SEC)
        self.assertGreater(archiver_ttl, PROC_BOARD_TTL_SEC)

        self.queue.publish_proc_board(
            "worker", owner_id="0", role="worker", ttl_sec=15
        )
        worker_ttl = int(self.queue.r.ttl(PROC_WORKER.format(id="0")))
        self.assertGreater(worker_ttl, 0)
        self.assertLessEqual(worker_ttl, 15)

    def test_reset_run_ephemeral_keys_deletes_proc_and_inflight(self) -> None:
        self.queue.r.hset(PROC_CORE, mapping={"role": "core"})
        self.queue.r.hset(PROC_WORKER.format(id="0"), mapping={"role": "worker"})
        self.queue.r.hset(PROC_CHILDREN.format(id="0"), mapping={"role": "children"})
        self.queue.r.rpush(INFLIGHT.format(worker="0"), "task-0")
        self.queue.r.hset(PROC_WORKER.format(id="99"), mapping={"role": "worker"})
        self.queue.reset_run_ephemeral_keys(worker_ids=4)
        self.assertEqual(self.queue.r.hgetall(PROC_CORE), {})
        self.assertEqual(self.queue.read_proc_board("worker", owner_id="0"), {})
        self.assertEqual(self.queue.read_proc_board("children", owner_id="0"), {})
        self.assertEqual(int(self.queue.r.llen(INFLIGHT.format(worker="0"))), 0)
        leftover = self.queue.read_proc_board("worker", owner_id="99")
        self.assertEqual(leftover.get("role"), "worker")

    def test_heartbeat_dual_writes_status_and_worker_board(self) -> None:
        self.queue.heartbeat(
            "0",
            status="idle",
            pid=1,
            ts=1.0,
            current_task="blob",
            held_calc_packs={"A": "001"},
        )
        board = self.queue.read_proc_board("worker", owner_id="0")
        self.assertTrue(board)
        self.assertEqual(board["role"], "worker")
        self.assertEqual(board["status"], "idle")
        self.assertEqual(int(board["pid"]), 1)
        self.assertEqual(board["current_uuid"], "")
        self.assertNotIn("current_task", board)
        self.assertNotIn("held_calc_packs", board)
        self.assertEqual(int(board["held_calc_n"]), 1)
        status = self.queue.r.hgetall(WORKER_STATUS.format(id="0"))
        self.assertEqual(status["status"], "idle")
        self.assertEqual(status["current_task"], "blob")
        self.assertTrue(status.get("held_calc_packs"))

    def test_heartbeat_renews_children_board_ttl(self) -> None:
        self.queue.publish_children_board(
            "0",
            file_operation_pid=1,
            file_operation_pgid=1,
            calc_pgids=[],
            ttl_sec=10,
        )
        key = PROC_CHILDREN.format(id="0")
        self.assertLessEqual(int(self.queue.r.ttl(key)), 10)
        self.queue.heartbeat("0", status="idle", pid=1, ts=1.0, board_ttl_sec=25)
        ttl = int(self.queue.r.ttl(key))
        self.assertGreater(ttl, 10)
        self.assertLessEqual(ttl, 25)

    def test_two_heartbeats_incr_op_count_once_each(self) -> None:
        self.queue.heartbeat("0", status="idle", pid=1, ts=1.0)
        self.queue.heartbeat("0", status="busy", pid=1, ts=2.0)
        self.assertEqual(self.queue.get_op_count("worker"), 2)

    def test_publish_proc_board_does_not_incr_op_count(self) -> None:
        before = self.queue.get_op_count("worker")
        self.queue.publish_proc_board("core", role="core", pid=1)
        self.queue.publish_proc_board("worker", owner_id="0", status="idle")
        self.assertEqual(self.queue.get_op_count("worker"), before)

    def test_drop_proc_board_deletes_key(self) -> None:
        self.queue.publish_proc_board("worker", owner_id="0", role="worker")
        self.assertTrue(self.queue.read_proc_board("worker", owner_id="0"))
        self.queue.drop_proc_board("worker", owner_id="0")
        self.assertEqual(self.queue.read_proc_board("worker", owner_id="0"), {})

    def test_touch_proc_board_missing_key_returns_false(self) -> None:
        self.assertFalse(self.queue.touch_proc_board("core"))
        self.queue.publish_proc_board("core", role="core")
        self.assertTrue(self.queue.touch_proc_board("core"))


if __name__ == "__main__":
    unittest.main()
