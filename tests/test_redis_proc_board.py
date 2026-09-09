#!/usr/bin/env python3
"""D26.1 proc-board mixin, dual-write heartbeat, and fresh-reset deletes."""

from __future__ import annotations

import json
import threading
import unittest
from unittest import mock

from jarvishep2.redis_queue import (
    ARCHIVER_BOARD_TTL_SEC,
    INFLIGHT,
    PROC_ARCHIVER,
    PROC_BOARD_TTL_SEC,
    PROC_CHILDREN,
    PROC_CORE,
    PROC_WORKER,
    SAMPLE_STATS,
    TASK_QUEUE,
    WORKER_STATUS,
    RedisQueue,
    encode_payload,
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

    def test_reconcile_resume_scan_deletes_out_of_range_proc_and_inflight(self) -> None:
        self.queue.r.hset(PROC_CORE, mapping={"role": "core"})
        self.queue.r.hset(PROC_WORKER.format(id="99"), mapping={"role": "worker"})
        self.queue.r.hset(PROC_CHILDREN.format(id="77"), mapping={"role": "children"})
        self.queue.r.rpush(INFLIGHT.format(worker="99"), "leftover")
        self.queue.reconcile_resume_ephemeral(completed=4, failed=1)
        self.assertEqual(self.queue.r.hgetall(PROC_CORE), {})
        self.assertEqual(self.queue.read_proc_board("worker", owner_id="99"), {})
        self.assertEqual(self.queue.read_proc_board("children", owner_id="77"), {})
        self.assertEqual(int(self.queue.r.llen(INFLIGHT.format(worker="99"))), 0)
        self.assertEqual(int(self.queue.fetch_sample_stats()["completed"]), 4)

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
        self.assertNotIn("u_coords", board)
        self.assertEqual(int(board["held_calc_n"]), 1)
        status = self.queue.r.hgetall(WORKER_STATUS.format(id="0"))
        self.assertEqual(status["status"], "idle")
        self.assertNotIn("current_task", status)
        self.assertTrue(status.get("held_calc_packs"))

    def test_heartbeat_overlays_children_pgids_without_extra_incr(self) -> None:
        before = self.queue.get_op_count("worker")
        self.queue.heartbeat(
            "0",
            status="busy",
            pid=1,
            ts=1.0,
            file_operation_pid=10,
            file_operation_pgid=10,
            calc_pgids=[20, 21],
            board_ttl_sec=25,
        )
        self.assertEqual(self.queue.get_op_count("worker"), before + 1)
        board = self.queue.read_children_board("0")
        self.assertEqual(int(board["file_operation_pid"]), 10)
        self.assertEqual(int(board["file_operation_pgid"]), 10)
        self.assertEqual(json.loads(board["calc_pgids"]), [20, 21])
        self.assertEqual(board["updated_reason"], "heartbeat")
        self.assertGreater(int(self.queue.r.ttl(PROC_CHILDREN.format(id="0"))), 10)

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


def _task(uuid: str, **extra: object) -> dict:
    payload = {"uuid": uuid, "u_coords": [0.1]}
    payload.update(extra)
    return payload


class InflightOwnershipTests(unittest.TestCase):
    def setUp(self) -> None:
        self.queue = make_fakeredis_queue(codec="json")
        self.inflight = INFLIGHT.format(worker="0")

    def _plant_inflight(self, *tasks: dict) -> None:
        # LPUSH so the first argument is the remaining occupancy head.
        for task in tasks:
            self.queue.r.lpush(self.inflight, encode_payload(task, codec="json"))

    def test_occupancy_llen_3_bounces_to_original_head(self) -> None:
        t1 = _task("t1")
        t2 = _task("t2")
        t3 = _task("t3")
        self._plant_inflight(t1, t2, t3)
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 3)
        self.queue.r.hset(SAMPLE_STATS, mapping={"running": 0})

        got = self.queue._run_occupancy(self.inflight)
        self.assertEqual(got["uuid"], "t1")
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 1)
        remaining = self.queue.get_inflight_task("0")
        self.assertEqual(remaining["uuid"], "t1")
        bounced = [
            self.queue._decode_task_payload(raw)["uuid"]
            for raw in self.queue.r.lrange(TASK_QUEUE, 0, -1)
        ]
        self.assertEqual(bounced, ["t2", "t3"])
        self.assertEqual(int(self.queue.r.hget(SAMPLE_STATS, "running") or 0), 0)

    def test_occupancy_n_start_1_increments_running_once(self) -> None:
        self.queue.push_task(_task("only"))
        got = self.queue.pull_task_to_inflight("0", timeout=0)
        self.assertEqual(got["uuid"], "only")
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 1)
        self.assertEqual(int(self.queue.r.hget(SAMPLE_STATS, "running") or 0), 1)

    def test_ack_mismatch_leaves_list_and_never_deletes(self) -> None:
        self._plant_inflight(_task("keep-me"), _task("extra"))
        delete = mock.Mock(side_effect=AssertionError("ACK must not DEL"))
        self.queue.r.delete = delete  # type: ignore[method-assign]
        self.assertFalse(self.queue.ack_inflight_task("0", "wrong-uuid"))
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 2)
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "extra")
        delete.assert_not_called()
        self.assertTrue(self.queue.ack_inflight_task("0", "extra"))
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 1)
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "keep-me")
        delete.assert_not_called()

    def test_factory_reclaim_two_inflight_retries_each(self) -> None:
        from types import SimpleNamespace

        from jarvishep2.runtime.factory import _Watchdog

        t1 = _task("r1", _retry_count=0)
        t2 = _task("r2", _retry_count=1)
        self.queue.push_task(_task("t3"))
        self.queue.push_task(_task("t4"))
        self._plant_inflight(t1, t2)
        watchdog = _Watchdog(SimpleNamespace(redis=self.queue))
        watchdog.max_sample_retries = 3
        self.assertTrue(watchdog.requeue_in_flight_task({}, worker_id="0"))
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 0)
        queued = [
            self.queue._decode_task_payload(raw)
            for raw in self.queue.r.lrange(TASK_QUEUE, 0, -1)
        ]
        self.assertEqual([item["uuid"] for item in queued], ["r1", "r2", "t3", "t4"])
        by_uuid = {item["uuid"]: int(item.get("_retry_count") or 0) for item in queued}
        self.assertEqual(by_uuid["r1"], 1)
        self.assertEqual(by_uuid["r2"], 2)

    def test_reclaim_lpush_batch_is_one_transaction(self) -> None:
        from types import SimpleNamespace

        from jarvishep2.runtime.factory import _Watchdog

        self._plant_inflight(_task("r1"), _task("r2"))
        executes = {"n": 0}
        real_pipeline = self.queue.r.pipeline

        def spy_pipeline(*args, **kwargs):
            pipe = real_pipeline(*args, **kwargs)
            real_execute = pipe.execute

            def counting_execute(*a, **k):
                executes["n"] += 1
                return real_execute(*a, **k)

            pipe.execute = counting_execute  # type: ignore[method-assign]
            return pipe

        self.queue.r.pipeline = spy_pipeline  # type: ignore[method-assign]
        watchdog = _Watchdog(SimpleNamespace(redis=self.queue))
        self.assertTrue(watchdog.requeue_in_flight_task({}, worker_id="0"))
        self.assertEqual(executes["n"], 1)
        queued = [
            self.queue._decode_task_payload(raw)["uuid"]
            for raw in self.queue.r.lrange(TASK_QUEUE, 0, -1)
        ]
        self.assertEqual(queued, ["r1", "r2"])

    def test_empty_reclaim_does_not_fallback_to_heartbeat(self) -> None:
        from types import SimpleNamespace

        from jarvishep2.runtime.factory import _Watchdog

        heartbeat = {
            "current_task": encode_payload(_task("already-acked"), codec="json"),
        }
        watchdog = _Watchdog(SimpleNamespace(redis=self.queue))
        self.assertFalse(watchdog.requeue_in_flight_task(heartbeat, worker_id="0"))
        self.assertEqual(int(self.queue.r.llen(TASK_QUEUE)), 0)

    def test_heartbeat_fallback_when_reclaim_missing(self) -> None:
        from types import SimpleNamespace

        from jarvishep2.runtime.factory import _Watchdog

        class _NoReclaim:
            def __init__(self, queue):
                self.queue = queue

            def decode_heartbeat_task(self, heartbeat):
                return self.queue.decode_heartbeat_task(heartbeat)

            def lpush_task(self, task):
                self.queue.lpush_task(task)

            def push_task(self, task):
                self.queue.push_task(task)

            def submit_result(self, info):
                self.queue.submit_result(info)

        heartbeat = {
            "current_task": encode_payload(_task("hb-only"), codec="json"),
        }
        watchdog = _Watchdog(SimpleNamespace(redis=_NoReclaim(self.queue)))
        self.assertTrue(watchdog.requeue_in_flight_task(heartbeat, worker_id="0"))
        queued = [
            self.queue._decode_task_payload(raw)
            for raw in self.queue.r.lrange(TASK_QUEUE, 0, -1)
        ]
        self.assertEqual([item["uuid"] for item in queued], ["hb-only"])
        self.assertEqual(int(queued[0]["_retry_count"]), 1)

    def test_reclaim_drains_all_inflight_items(self) -> None:
        self._plant_inflight(_task("a"), _task("b"))
        payloads = self.queue.reclaim_inflight_task("0")
        uuids = {item["uuid"] for item in payloads}
        self.assertEqual(uuids, {"a", "b"})
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 0)
        self.assertEqual(self.queue.reclaim_inflight_task("0"), [])

    def test_get_inflight_task_does_not_consume(self) -> None:
        self._plant_inflight(_task("head"))
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "head")
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "head")
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 1)

    def test_list_inflight_worker_ids_scans_keys(self) -> None:
        self.queue.r.lpush(INFLIGHT.format(worker="3"), encode_payload(_task("x"), codec="json"))
        self.queue.r.lpush(INFLIGHT.format(worker="9"), encode_payload(_task("y"), codec="json"))
        self.assertEqual(sorted(self.queue.list_inflight_worker_ids()), ["3", "9"])

    def test_fakeredis_without_blmove_uses_steal_lua(self) -> None:
        self.queue.push_task(_task("stolen"))
        self.queue.r.blmove = mock.Mock(  # type: ignore[method-assign]
            side_effect=Exception("unknown command 'BLMOVE'")
        )
        eval_scripts: list[str] = []
        real_eval = self.queue.r.eval

        def spy_eval(script, nkeys, *args):
            eval_scripts.append(str(script))
            return real_eval(script, nkeys, *args)

        self.queue.r.eval = spy_eval  # type: ignore[method-assign]
        got = self.queue.pull_task_to_inflight("0", timeout=1)
        self.assertEqual(got["uuid"], "stolen")
        self.assertTrue(any("occupied" in script for script in eval_scripts))
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "stolen")
        occupied = self.queue.pull_task_to_inflight("0", timeout=1)
        self.assertIsNone(occupied)

    def test_production_pull_raises_without_blmove(self) -> None:
        queue = RedisQueue({"host": "127.0.0.1", "port": 1, "db": 0})
        blocking = mock.MagicMock()
        blocking.blmove.side_effect = Exception("unknown command 'BLMOVE'")
        queue.r = blocking
        queue.r_ctrl = mock.MagicMock()
        queue._allow_steal_inflight = False
        with self.assertRaises(RuntimeError) as ctx:
            queue.pull_task_to_inflight("0", timeout=1)
        self.assertIn("BLMOVE", str(ctx.exception))
        self.assertIn("6.2", str(ctx.exception))

    def test_blmove_timeout_with_llen_gt_1_runs_occupancy(self) -> None:
        self._plant_inflight(_task("t1"), _task("t2"), _task("t3"))
        self.queue.r.blmove = mock.Mock(return_value=None)  # type: ignore[method-assign]
        got = self.queue.pull_task_to_inflight("0", timeout=1)
        self.assertEqual(got["uuid"], "t1")
        self.assertEqual(int(self.queue.r.llen(self.inflight)), 1)
        self.assertEqual(self.queue.get_inflight_task("0")["uuid"], "t1")

    def test_ack_and_reclaim_do_not_both_yield(self) -> None:
        self._plant_inflight(_task("one"))
        acked: list[bool] = []
        reclaimed: list[list] = []

        def do_ack() -> None:
            acked.append(self.queue.ack_inflight_task("0", "one"))

        def do_reclaim() -> None:
            reclaimed.append(self.queue.reclaim_inflight_task("0"))

        ack_thread = threading.Thread(target=do_ack)
        reclaim_thread = threading.Thread(target=do_reclaim)
        ack_thread.start()
        reclaim_thread.start()
        ack_thread.join()
        reclaim_thread.join()
        ack_got = bool(acked and acked[0])
        rec_items = reclaimed[0] if reclaimed else []
        rec_got = bool(rec_items)
        self.assertTrue(ack_got or rec_got)
        self.assertFalse(ack_got and rec_got)
        if rec_got:
            self.assertEqual(rec_items[0]["uuid"], "one")

    def test_ack_then_reclaim_and_reclaim_then_ack_are_exclusive(self) -> None:
        self._plant_inflight(_task("seq-a"))
        self.assertTrue(self.queue.ack_inflight_task("0", "seq-a"))
        self.assertEqual(self.queue.reclaim_inflight_task("0"), [])

        self._plant_inflight(_task("seq-b"))
        payloads = self.queue.reclaim_inflight_task("0")
        self.assertEqual([item["uuid"] for item in payloads], ["seq-b"])
        self.assertFalse(self.queue.ack_inflight_task("0", "seq-b"))

    def test_require_blmove_raises_on_missing_command(self) -> None:
        client = mock.MagicMock()
        client.execute_command.return_value = []
        client.info.return_value = {"redis_version": "6.0.16"}
        queue = RedisQueue(client=client)
        with self.assertRaises(RuntimeError) as ctx:
            queue.require_blmove()
        message = str(ctx.exception)
        self.assertIn("BLMOVE", message)
        self.assertIn("6.2", message)
        self.assertIn("Ubuntu 22.04", message)

    def test_connect_skips_require_blmove_for_injected_client(self) -> None:
        queue = make_fakeredis_queue()
        with mock.patch.object(RedisQueue, "require_blmove") as require:
            queue.connect()
            require.assert_not_called()

    def test_blmove_uses_blocking_client_occupancy_uses_ctrl(self) -> None:
        import fakeredis

        server = fakeredis.FakeServer()
        blocking = fakeredis.FakeStrictRedis(server=server, decode_responses=True)
        control = fakeredis.FakeStrictRedis(server=server, decode_responses=True)
        queue = RedisQueue({"codec": "json"}, client=blocking)
        queue.r_ctrl = control

        class _CtrlBoom:
            def blmove(self, *_a, **_k):
                raise AssertionError("control client must not BLMOVE")

            def eval(self, *args, **kwargs):
                return control.eval(*args, **kwargs)

            def llen(self, *args, **kwargs):
                return control.llen(*args, **kwargs)

            def lrange(self, *args, **kwargs):
                return control.lrange(*args, **kwargs)

            def delete(self, *args, **kwargs):
                return control.delete(*args, **kwargs)

            def lpop(self, *args, **kwargs):
                return control.lpop(*args, **kwargs)

            def lpush(self, *args, **kwargs):
                return control.lpush(*args, **kwargs)

            def lindex(self, *args, **kwargs):
                return control.lindex(*args, **kwargs)

            def hincrby(self, *args, **kwargs):
                return control.hincrby(*args, **kwargs)

        queue.r_ctrl = _CtrlBoom()  # type: ignore[assignment]
        queue.push_task(_task("via-blmove"))
        got = queue.pull_task_to_inflight("0", timeout=0)
        self.assertEqual(got["uuid"], "via-blmove")
        self.assertEqual(queue.get_inflight_task("0")["uuid"], "via-blmove")


if __name__ == "__main__":
    unittest.main()
