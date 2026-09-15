"""INF-15a / want API: default-off calculator inflight telemetry."""

from __future__ import annotations

import json
import threading
import time
import unittest
from types import SimpleNamespace
from unittest import mock

import numpy as np

from jarvishep2.redis_queue import (
    CALC_STATUS,
    CONTROL_LOCK,
    INFLIGHT,
    MONITOR_CALC_BUSY,
    MONITOR_SAMPLE_RUNNING,
    MONITOR_WANT,
    MONITOR_WANT_TTL_SEC,
    SAMPLE_STATS,
    WORKER_STATUS,
    calc_busy_packs_key,
    calc_free_list_key,
    calc_shared_free_list_key,
    calc_status_busy_field,
    calc_status_free_field,
    make_fakeredis_queue,
)
from jarvishep2.monitor.collector import Collector, open_collector, want_channels_for_page
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.runtime.sample_executor import SampleExecutor
from jarvishep2.sample import Sample
from jarvishep2.worker import Worker


class MonitorWantApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.queue = make_fakeredis_queue()

    def test_set_monitor_want_is_last_writer_wins_hash_plus_expire(self) -> None:
        self.queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=11, ttl_sec=5)
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "ch:calc"), "1")
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "ch:sample"), "0")
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "pid"), "11")
        ttl = int(self.queue.r.ttl(MONITOR_WANT) or 0)
        self.assertGreater(ttl, 0)
        self.assertLessEqual(ttl, MONITOR_WANT_TTL_SEC)
        self.queue.set_monitor_want(ch_calc=False, ch_sample=True, pid=22, ttl_sec=5)
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "ch:calc"), "0")
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "ch:sample"), "1")
        self.assertEqual(self.queue.r.hget(MONITOR_WANT, "pid"), "22")
        self.assertIsNone(self.queue.r.get(CONTROL_LOCK))

    def test_clear_monitor_want_is_idempotent(self) -> None:
        self.queue.clear_monitor_want()
        self.queue.set_monitor_want(ch_calc=True, ch_sample=False)
        self.queue.clear_monitor_want()
        self.assertEqual(int(self.queue.r.exists(MONITOR_WANT) or 0), 0)
        self.queue.clear_monitor_want()
        self.assertEqual(int(self.queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_set_monitor_want_never_set_nx_control_lock(self) -> None:
        self.queue.claim_control_lock("core-owner", ttl_sec=30)
        self.queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=9)
        self.assertEqual(self.queue.get_control_lock_owner(), "core-owner")


class MonitorInflightHotPathTests(unittest.TestCase):
    def setUp(self) -> None:
        self.queue = make_fakeredis_queue()
        self.queue.register_calc_pool("DemoCalc", 2)

    def test_want_off_owns_pack_without_status_op_count_or_sidecar(self) -> None:
        pack, want_calc = self.queue._acquire_calc("DemoCalc", timeout=1, worker_id="3")
        self.assertEqual(pack, "001")
        self.assertFalse(want_calc)
        busy = self.queue.r.hgetall(calc_busy_packs_key("DemoCalc"))
        self.assertEqual(busy.get("001") or busy.get(b"001"), "running")
        self.assertEqual(self.queue.r.hlen(calc_busy_packs_key("DemoCalc")), 1)
        self.assertEqual(int(self.queue.r.exists(MONITOR_CALC_BUSY.format(name="DemoCalc")) or 0), 0)
        self.assertEqual(self.queue.get_op_count("calculator"), 0)
        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 2)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 0)
        occupancy = self.queue.fetch_calc_occupancy(["DemoCalc"], slots={"DemoCalc": 2})
        self.assertEqual(occupancy["DemoCalc"], {"free": 1, "busy": 1, "slots": 2})

    def test_want_calc_writes_sidecar_and_status(self) -> None:
        self.queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=1)
        pack, want_calc = self.queue._acquire_calc("DemoCalc", timeout=1, worker_id="3")
        self.assertEqual(pack, "001")
        self.assertTrue(want_calc)
        self.assertEqual(
            self.queue.fetch_monitor_calc_busy("DemoCalc"),
            {"001": "3"},
        )
        self.assertEqual(self.queue.get_op_count("calculator"), 1)
        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("DemoCalc")]), 1)
        self.assertEqual(int(status[calc_status_busy_field("DemoCalc")]), 1)
        self.queue.release_calc("DemoCalc", pack)
        self.assertEqual(self.queue.fetch_monitor_calc_busy("DemoCalc"), {})
        self.assertEqual(self.queue.r.llen(calc_free_list_key("DemoCalc")), 2)

    def test_want_sample_without_calc_skips_calc_sidecar(self) -> None:
        overlay = '{"uuid":"abc","step":"DemoCalc","t0":1}'
        self.queue.set_monitor_want(ch_calc=False, ch_sample=True, pid=1)
        pack, want_calc = self.queue._acquire_calc(
            "DemoCalc", timeout=1, worker_id="7", overlay_json=overlay
        )
        self.assertEqual(pack, "001")
        self.assertFalse(want_calc)
        self.assertEqual(self.queue.fetch_monitor_calc_busy("DemoCalc"), {})
        running = self.queue.fetch_monitor_sample_running()
        self.assertEqual(running["7"]["uuid"], "abc")
        self.assertEqual(self.queue.get_op_count("calculator"), 0)

    def test_want_ttl_expiry_makes_the_next_claim_fail_closed(self) -> None:
        self.queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=1, ttl_sec=1)
        time.sleep(1.1)

        pack, want_calc = self.queue._acquire_calc("DemoCalc", timeout=1, worker_id="3")

        self.assertEqual(pack, "001")
        self.assertFalse(want_calc)
        self.assertEqual(self.queue.fetch_monitor_calc_busy("DemoCalc"), {})
        self.assertEqual(self.queue.get_op_count("calculator"), 0)
        self.assertEqual(self.queue.r.hlen(calc_busy_packs_key("DemoCalc")), 1)

    def test_double_release_still_raises_when_lua_returns_pair(self) -> None:
        pack = self.queue.acquire_calc("DemoCalc", timeout=1)
        self.queue.release_calc("DemoCalc", pack)
        with self.assertRaises(ValueError):
            self.queue.release_calc("DemoCalc", pack)
        self.assertEqual(
            set(self.queue.r.lrange(calc_free_list_key("DemoCalc"), 0, -1)),
            {"001", "002"},
        )

    def test_junk_token_does_not_hincrby_status_when_want_off(self) -> None:
        key = calc_free_list_key("EggBox")
        self.queue.r.delete(key)
        self.queue.r.rpush(key, "ready", "001")
        self.queue.r.hset(
            CALC_STATUS,
            mapping={
                calc_status_free_field("EggBox"): 2,
                calc_status_busy_field("EggBox"): 0,
            },
        )
        pack = self.queue.acquire_calc("EggBox", timeout=2)
        self.assertEqual(pack, "001")
        status = self.queue.r.hgetall(CALC_STATUS)
        self.assertEqual(int(status[calc_status_free_field("EggBox")]), 2)

    def test_reset_drops_monitor_keys(self) -> None:
        self.queue.set_monitor_want(ch_calc=True, ch_sample=True, pid=1)
        self.queue.r.hset(MONITOR_CALC_BUSY.format(name="DemoCalc"), "001", "3")
        self.queue.r.hset(MONITOR_SAMPLE_RUNNING, "3", '{"uuid":"x"}')
        self.queue.reset_run_ephemeral_keys(calculator_names=["DemoCalc"], worker_ids=[0])
        self.assertEqual(int(self.queue.r.exists(MONITOR_WANT) or 0), 0)
        self.assertEqual(int(self.queue.r.exists(MONITOR_CALC_BUSY.format(name="DemoCalc")) or 0), 0)
        self.assertEqual(int(self.queue.r.exists(MONITOR_SAMPLE_RUNNING) or 0), 0)


class MonitorSharedPoolWantOffTests(unittest.TestCase):
    def test_want_off_still_updates_packmode_and_affinity(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_shared_calc_pool("Prep", 1, modes=["fast", "full"])
        first = queue.acquire_shared_calc("Prep", "fast", modes=["fast", "full"], timeout=1)
        self.assertEqual(first, ("001", None))
        self.assertEqual(queue.get_op_count("calculator"), 0)
        self.assertTrue(queue.release_shared_calc("Prep", "001", "fast"))
        self.assertEqual(
            list(queue.r.lrange(calc_shared_free_list_key("Prep", "fast"), 0, -1)),
            ["001"],
        )
        self.assertEqual(queue.r.hget("calc:packmode:Prep", "001"), "fast")
        warm = queue.acquire_shared_calc("Prep", "fast", modes=["fast", "full"], timeout=1)
        self.assertEqual(warm, ("001", "fast"))


class MonitorAlwaysOnWantOffTests(unittest.TestCase):
    def test_inflight_and_sample_stats_still_move_with_no_want_key(self) -> None:
        queue = make_fakeredis_queue()
        queue.push_task({"uuid": "always-on-sample", "u_coords": [0.0]})

        task = queue.pull_task_to_inflight("3", timeout=0)

        self.assertEqual(task["uuid"], "always-on-sample")
        self.assertEqual(queue.r.llen(INFLIGHT.format(worker="3")), 1)
        self.assertEqual(int(queue.r.hget(SAMPLE_STATS, "running") or 0), 1)
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)


class _PassCalc:
    config = {"name": "DemoCalc"}

    def should_run(self, _info):
        return True

    def acquire_pack_id(self, _pack):
        return None

    def prepare_runtime(self, _info):
        return None

    def execute(self, _info, runtime_prepared=True):
        return {}


def _executor_worker(queue, *, worker_id: int = 3) -> SimpleNamespace:
    lock = threading.Lock()
    heartbeats: list[str | None] = []
    worker = SimpleNamespace(
        worker_id=worker_id,
        worker_config={"calc_acquire_timeout": 1},
        _redis=queue,
        _calculators={"DemoCalc": _PassCalc()},
        _shared_mode_sets={},
        _held_calc_packs={},
        _current_sample_uuid="abc-uuid",
        _observables_lock=None,
        _heartbeat_calls=heartbeats,
        _heartbeat=lambda status=None: heartbeats.append(status),
        _hb_lock=lambda: lock,
        _merge_calculator_observables=lambda *_a, **_k: None,
    )
    return worker


class SampleExecutorInflightTests(unittest.TestCase):
    def test_want_off_skips_full_heartbeat_but_publishes_held_packs(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 1)
        log = _install_write_log(queue)
        worker = _executor_worker(queue)
        SampleExecutor(worker)._run_calculator_step(
            "DemoCalc",
            Sample(uuid="abc-uuid", u_coords=np.array([0.0])),
        )
        self.assertEqual(worker._heartbeat_calls, [])
        raw = queue.r.hget(WORKER_STATUS.format(id="3"), "held_calc_packs")
        self.assertIsNotNone(raw)
        self.assertEqual(queue.fetch_monitor_calc_busy("DemoCalc"), {})
        self.assertEqual(queue.get_op_count("calculator"), 0)
        self.assertEqual(queue.r.llen(calc_free_list_key("DemoCalc")), 1)
        held_updates = [
            op
            for op, key in log
            if op == "hset" and key == WORKER_STATUS.format(id="3")
        ]
        self.assertEqual(len(held_updates), 2)

    def test_want_calc_beats_and_writes_sidecar(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 1)
        queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=1)
        worker = _executor_worker(queue)
        SampleExecutor(worker)._run_calculator_step(
            "DemoCalc",
            Sample(uuid="abc-uuid", u_coords=np.array([0.0])),
        )
        self.assertEqual(worker._heartbeat_calls, ["busy", "busy"])
        self.assertEqual(queue.fetch_monitor_calc_busy("DemoCalc"), {})
        self.assertEqual(queue.get_op_count("calculator"), 2)


class SampleOverlayTests(unittest.TestCase):
    def test_overlay_writes_only_when_ch_sample(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_sample_overlay("3", {"uuid": "u", "step": "", "t0": 1})
        self.assertEqual(queue.fetch_monitor_sample_running(), {})
        queue.set_monitor_want(ch_calc=False, ch_sample=True, pid=1)
        queue.publish_sample_overlay("3", {"uuid": "u", "step": "", "t0": 1})
        self.assertEqual(queue.fetch_monitor_sample_running()["3"]["uuid"], "u")
        queue.clear_sample_overlay("3")
        self.assertEqual(queue.fetch_monitor_sample_running(), {})

    def test_overlay_skipped_on_calc_page(self) -> None:
        queue = make_fakeredis_queue()
        queue.set_monitor_want(ch_calc=True, ch_sample=False, pid=1)
        queue.publish_sample_overlay("3", {"uuid": "u", "step": "", "t0": 1})
        self.assertEqual(queue.fetch_monitor_sample_running(), {})

    def test_overlay_failure_is_fail_open_and_does_not_strand_a_pack(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 1)
        worker = SimpleNamespace(
            worker_id=3,
            _redis=queue,
            _current_sample_uuid="sample-with-sidecar-error",
        )
        with mock.patch.object(
            queue,
            "publish_sample_overlay",
            side_effect=ConnectionError("monitor sidecar unavailable"),
        ):
            Worker._publish_sample_overlay(worker, "start")

        pack = queue.acquire_calc("DemoCalc", timeout=1)
        self.assertEqual(pack, "001")
        queue.release_calc("DemoCalc", pack)
        self.assertEqual(queue.r.llen(calc_free_list_key("DemoCalc")), 1)


def _install_write_log(queue) -> list[tuple[str, str | None]]:
    log: list[tuple[str, str | None]] = []

    def wrap(name: str, original):
        def wrapped(*args, **kwargs):
            key = None
            if args:
                if name == "eval" and len(args) >= 3:
                    key = args[2]
                else:
                    key = args[0]
            log.append((name, None if key is None else str(key)))
            return original(*args, **kwargs)

        return wrapped

    seen: set[int] = set()
    for client in (queue.r, queue.r_ctrl):
        if client is None or id(client) in seen:
            continue
        seen.add(id(client))
        for method_name in (
            "set",
            "hset",
            "incr",
            "rpush",
            "lpush",
            "delete",
            "hincrby",
            "expire",
            "eval",
        ):
            if hasattr(client, method_name):
                setattr(client, method_name, wrap(method_name, getattr(client, method_name)))
    return log


def _install_tui_isolation_guard(queue) -> None:
    """Fail if Collector reaches outside the monitor-key Redis allowance."""

    def text(value) -> str:
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return str(value)

    def allowed_monitor_key(value) -> bool:
        return text(value) == MONITOR_WANT

    def guard(name: str, original):
        def wrapped(*args, **kwargs):
            if name == "set":
                raise AssertionError("TUI must never SET, including SET NX")
            if name == "eval":
                if len(args) < 3 or not allowed_monitor_key(args[2]):
                    raise AssertionError("TUI EVAL may only target hep:monitor:want")
            elif name in {"delete", "expire", "hset", "hdel", "hincrby"}:
                keys = args if name == "delete" else args[:1]
                if not keys or not all(allowed_monitor_key(key) for key in keys):
                    raise AssertionError(f"TUI {name} escaped hep:monitor:want")
            elif name in {"incr", "incrby", "rpush", "lpush"}:
                raise AssertionError(f"TUI must never issue {name}")
            elif name in {"scan", "scan_iter", "keys", "lrange"}:
                raise AssertionError(f"TUI must never issue {name}")
            elif name == "get" and args and text(args[0]).startswith("hep:results:"):
                raise AssertionError("TUI must never GET result payloads")
            return original(*args, **kwargs)

        return wrapped

    def guarded_pipeline(original):
        def wrapped(*args, **kwargs):
            pipeline = original(*args, **kwargs)
            for method_name in ("hset", "hdel", "expire", "delete", "incr", "incrby"):
                if hasattr(pipeline, method_name):
                    setattr(
                        pipeline,
                        method_name,
                        guard(method_name, getattr(pipeline, method_name)),
                    )
            return pipeline

        return wrapped

    seen: set[int] = set()
    for client in (queue.r, queue.r_ctrl):
        if client is None or id(client) in seen:
            continue
        seen.add(id(client))
        for method_name in (
            "set",
            "hset",
            "hdel",
            "delete",
            "expire",
            "hincrby",
            "incr",
            "incrby",
            "rpush",
            "lpush",
            "eval",
            "scan",
            "scan_iter",
            "keys",
            "lrange",
            "get",
        ):
            if hasattr(client, method_name):
                setattr(client, method_name, guard(method_name, getattr(client, method_name)))
        if hasattr(client, "pipeline"):
            setattr(client, "pipeline", guarded_pipeline(client.pipeline))


class CollectorWantTests(unittest.TestCase):
    def test_overview_rate_uses_completed_delta_and_monotonic_time(self) -> None:
        queue = make_fakeredis_queue()
        clock = iter((100.0, 110.0))
        collector = Collector(
            queue,
            pid=9,
            monotonic=lambda: next(clock),
            sampler_metadata={
                "method": "Random",
                "family": "simple",
                "config": {"point_number": 100},
            },
        )
        queue.publish_proc_board(
            "core",
            sampler_status=json.dumps(
                {
                    "schema": 1,
                    "method": "Random",
                    "state": "running",
                    "progress": {"kind": "finite", "current": 10, "target": 100},
                    "metrics": {},
                }
            ),
        )
        queue.r.hset(
            SAMPLE_STATS,
            mapping={"completed": 4, "running": 0, "failed": 0},
        )
        first = collector.tick("overview")
        queue.r.hset(SAMPLE_STATS, "completed", 9)
        queue.publish_proc_board(
            "core",
            sampler_status=json.dumps(
                {
                    "schema": 1,
                    "method": "Random",
                    "state": "running",
                    "progress": {"kind": "finite", "current": 15, "target": 100},
                    "metrics": {},
                }
            ),
        )
        second = collector.tick("overview")

        self.assertEqual(first.sample_rate, "—")
        self.assertEqual(second.sample_rate, "30.0 / min")
        self.assertEqual(first.sampler_display.progress.eta, "—")
        self.assertEqual(second.sampler_display.progress.eta, "0:02:50")

    def test_page_channel_map(self) -> None:
        self.assertEqual(want_channels_for_page("calculators"), (True, False))
        self.assertEqual(want_channels_for_page("samples"), (False, True))
        for slug in ("overview", "workers", "factory", "sampler", "host", "splash"):
            self.assertIsNone(want_channels_for_page(slug))

    def test_want_only_on_pages_5_and_6_and_del_once(self) -> None:
        queue = make_fakeredis_queue()
        log = _install_write_log(queue)
        collector = Collector(queue, pid=9)
        collector.tick("overview")
        collector.tick("workers")
        collector.tick("factory")
        collector.tick("sampler")
        collector.tick("host")
        self.assertEqual(log, [])
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

        collector.tick("calculators")
        self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:calc"), "1")
        self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:sample"), "0")
        self.assertTrue(log)
        self.assertTrue(all(key == MONITOR_WANT for _op, key in log if key))
        self.assertNotIn(CONTROL_LOCK, [key for _op, key in log])

        log.clear()
        collector.tick("overview")
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)
        deletes = [op for op, key in log if op == "delete" and key == MONITOR_WANT]
        self.assertEqual(len(deletes), 1)

        log.clear()
        collector.tick("overview")
        collector.tick("host")
        self.assertEqual(log, [])

        collector.tick("samples")
        self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:sample"), "1")
        self.assertEqual(queue.r.hget(MONITOR_WANT, "ch:calc"), "0")
        collector.close()
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_occupancy_and_sidecar_reads_do_not_write(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 2)
        log = _install_write_log(queue)
        collector = Collector(
            queue,
            pid=9,
            pool_names=["DemoCalc"],
            slots={"DemoCalc": 2},
            selected_calc="DemoCalc",
        )
        frame = collector.tick("calculators")
        self.assertEqual(frame.occupancy["DemoCalc"]["slots"], 2)
        self.assertEqual(frame.calc_busy, {})
        writes = {key for _op, key in log if key}
        self.assertEqual(writes, {MONITOR_WANT})

    def test_tick_failure_keeps_the_last_good_frame_and_marks_it_stale(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 2)
        collector = Collector(
            queue,
            pid=9,
            pool_names=["DemoCalc"],
            slots={"DemoCalc": 2},
            selected_calc="DemoCalc",
        )
        good = collector.tick("calculators")
        with mock.patch.object(
            queue,
            "fetch_monitor_calc_busy",
            side_effect=ConnectionError("read timeout"),
        ):
            stale = collector.tick("calculators")

        self.assertTrue(stale.stale)
        self.assertEqual(stale.error, "read timeout")
        self.assertEqual(stale.page, "calculators")
        self.assertEqual(stale.occupancy, good.occupancy)
        self.assertFalse(collector.last_frame.stale)
        self.assertEqual(collector.last_frame, good)

    def test_failed_want_write_is_nonfatal_and_fail_closed(self) -> None:
        queue = make_fakeredis_queue()
        collector = Collector(queue, pid=9)
        with mock.patch.object(
            queue,
            "set_monitor_want",
            side_effect=ConnectionError("want eval unavailable"),
        ):
            frame = collector.tick("calculators")

        self.assertTrue(frame.stale)
        self.assertFalse(frame.want_on)
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_standard_pages_project_the_existing_read_only_snapshot(self) -> None:
        queue = make_fakeredis_queue()
        queue.publish_proc_board("core", role="core", pid=44001, status="running")
        queue.publish_proc_board("archiver", role="archiver", pid=44002, status="idle")
        queue.publish_proc_board("redis", role="redis", pid=44003, status="ready")
        queue.publish_proc_board(
            "worker",
            owner_id="3",
            pid=44004,
            status="busy",
            current_sample="active-sample",
            held_calc_n=1,
        )
        queue.push_task({"uuid": "queued-sample", "u_coords": [0.0]})
        collector = Collector(
            queue,
            owner_ids=["3"],
            pid=9,
            host_snapshotter=lambda: {
                "available": True,
                "cpu_percent": 25.0,
                "memory_used": 8 * 1024**3,
                "memory_total": 32 * 1024**3,
                "swap_used": 0,
                "swap_total": 8 * 1024**3,
                "load": (1.0, 0.8, 0.6),
                "processes": [{"role": "worker", "pid": 44004, "alive": True}],
            },
        )

        overview = collector.tick("overview")
        workers = collector.tick("workers")
        factory = collector.tick("factory")
        sampler = collector.tick("sampler")
        host = collector.tick("host")

        self.assertEqual(overview.queues["task_queue_length"], 1)
        self.assertEqual(workers.workers[0]["worker_id"], "3")
        self.assertEqual(workers.workers[0]["current_uuid"], "active-sample")
        self.assertEqual(factory.proc_core["status"], "running")
        self.assertEqual(factory.proc_archiver["status"], "idle")
        self.assertEqual(factory.proc_redis["status"], "ready")
        self.assertEqual(sampler.queues["task_queue_length"], 1)
        self.assertEqual(host.host["cpu_percent"], 25.0)
        self.assertEqual(host.host["processes"][0]["pid"], 44004)
        self.assertEqual(int(queue.r.exists(MONITOR_WANT) or 0), 0)

    def test_live_collector_obeys_the_tui_redis_isolation_contract(self) -> None:
        queue = make_fakeredis_queue()
        queue.register_calc_pool("DemoCalc", 1)
        queue.r.hset(SAMPLE_STATS, mapping={"completed": 2, "running": 1})
        queue.r.rpush(INFLIGHT.format(worker="3"), '{"uuid":"work-only"}')
        queue.claim_control_lock("core-owner", ttl_sec=30)
        _install_tui_isolation_guard(queue)
        collector = Collector(
            queue,
            owner_ids=["3"],
            pid=9,
            pool_names=["DemoCalc"],
            slots={"DemoCalc": 1},
        )

        for page in (
            "overview",
            "workers",
            "factory",
            "sampler",
            "host",
            "calculators",
            "samples",
        ):
            collector.tick(page)
        collector.close()

        self.assertEqual(queue.get_control_lock_owner(), "core-owner")

    def test_open_collector_uses_runtime_pool_catalogue(self) -> None:
        queue = make_fakeredis_queue()
        choice = ScanChoice(
            reference="R1",
            name="live-scan",
            control_pid=44001,
            process_count=1,
            pids=(44001,),
            simulated=False,
        )
        metadata = {
            "redis": {"host": "127.0.0.1", "port": 6381, "db": 2},
            "calculator_pools": {"Slow": 3},
            "calculator_shared": {"Prep": {"modes": ["fast", "full"], "n": 2}},
        }
        with mock.patch("jarvishep2.process_cleanup.list_active_scans", return_value=[]), mock.patch(
            "jarvishep2.process_cleanup.resolve_scan_reference", return_value=object()
        ), mock.patch(
            "jarvishep2.process_cleanup.runtime_metadata_for_scan", return_value=metadata
        ), mock.patch("jarvishep2.monitor.collector.RedisQueue", return_value=queue), mock.patch.object(
            queue, "connect", return_value=None
        ):
            collector = open_collector(choice)

        assert collector is not None
        self.assertEqual(collector._pool_names, ["Slow", "Prep"])
        self.assertEqual(collector._slots, {"Slow": 3, "Prep": 2})
        collector.close()


if __name__ == "__main__":
    unittest.main()
