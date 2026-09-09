#!/usr/bin/env python3
"""WP-D6.1 Worker failure recovery tests."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.slow

import inspect
import logging
import os
import signal
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from types import SimpleNamespace
from typing import Any
from unittest import mock

from fakeredis import TcpFakeServer

import jarvishep2.factory as factory_module
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import (
    INFLIGHT,
    RedisQueue,
    calc_status_busy_field,
    calc_status_free_field,
    make_fakeredis_queue,
)

from test_layer_concurrency import _slow_calc_module, SLOW_A_SCRIPT
from test_worker_calculator import _start_tcp_fakeredis, _worker_config
from test_worker_mvp import (
    BENCHMARK_OPERA_MODULE,
    LIKELIHOOD_EXPRESSIONS,
    SAMPLING_VARIABLES,
)


def _opera_worker_config(tmpdir: str, *, delay_sec: float = 0.0) -> dict[str, Any]:
    config = {
        "sample_config": {
            "task_result_dir": tmpdir,
            "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
            "sample_artifacts": "auto",
            "workflow_has_calculator": False,
            "workflow_references_sdir": False,
        },
        "mapper": {
            "type": "flat",
            "variables": SAMPLING_VARIABLES,
        },
        "opera_modules": [BENCHMARK_OPERA_MODULE],
        "likelihood_expressions": LIKELIHOOD_EXPRESSIONS,
        "pull_timeout": 1,
        "handoff_to_staging": False,
    }
    if delay_sec > 0:
        config["test_process_delay_sec"] = delay_sec
    return config


def _task(uuid: str, *, x: float, y: float, shift: float) -> dict[str, Any]:
    return {
        "uuid": uuid,
        "u_coords": [x, y, shift],
        "execution_plan": [
            {"name": "TrivialEggbox", "type": "opera", "layer": 0},
            {"name": "LogL_Z", "type": "likelihood", "layer": 1},
        ],
    }


def _wait_until(predicate, *, timeout: float = 30.0, poll: float = 0.05) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(poll)
    raise TimeoutError("condition not met before timeout")


def _worker_is_busy(factory: TaskFactory, worker_id: int = 0) -> bool:
    if factory.redis is None:
        return False
    heartbeat = factory.redis.fetch_worker_status([str(worker_id)]).get(str(worker_id), {})
    return str(heartbeat.get("status") or "").lower() == "busy"


def _worker_holds_calc_slot(
    factory: TaskFactory,
    calc_name: str,
    *,
    worker_id: int = 0,
) -> bool:
    if factory.redis is None or not _worker_is_busy(factory, worker_id=worker_id):
        return False
    heartbeat = factory.redis.fetch_worker_status([str(worker_id)]).get(str(worker_id), {})
    held = factory.redis.decode_heartbeat_held_packs(heartbeat)
    return calc_name in held


def _drain_archive_uuids(redis: RedisQueue, expected: int, *, timeout: float = 45.0) -> list[str]:
    seen: list[str] = []
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and len(seen) < expected:
        result = redis.pull_result(timeout=1)
        if result is None:
            continue
        seen.append(str(result.get("uuid")))
    return seen


class WorkerFailureTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_sigkill_worker_requeues_and_completes_scan(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        expected_uuids = [f"sample-{index}" for index in range(3)]
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_opera_worker_config(tmpdir, delay_sec=2.0))
                factory.start_watchdog(
                    stale_sec=30.0,
                    poll_interval_sec=0.1,
                    max_sample_retries=3,
                )
                assert factory.redis is not None
                for index, sample_uuid in enumerate(expected_uuids):
                    factory.redis.push_task(
                        _task(
                            sample_uuid,
                            x=0.1 * (index + 1),
                            y=0.2 * (index + 1),
                            shift=0.3 * (index + 1),
                        )
                    )

                _wait_until(lambda: _worker_is_busy(factory), timeout=15.0)
                worker = factory.workers[0]
                self.assertIsNotNone(worker.pid)
                os.kill(worker.pid, signal.SIGKILL)
                worker.join(timeout=5.0)

                _wait_until(lambda: any(item.is_alive() for item in factory.workers), timeout=10.0)
                self.assertGreaterEqual(factory._respawn_count, 1)

                completed = _drain_archive_uuids(factory.redis, len(expected_uuids))
                self.assertEqual(sorted(completed), sorted(expected_uuids))
                self.assertEqual(len(completed), len(set(completed)))
                factory.shutdown()
        finally:
            server.shutdown()
            server.server_close()

    def test_sigkill_worker_releases_held_calc_slot(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        slow_module = _slow_calc_module("SlowA", SLOW_A_SCRIPT, "a")
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                worker_cfg = _worker_config(tmpdir)
                worker_cfg["calculator_modules"] = [slow_module]
                worker_cfg["calculator_pools"] = {"SlowA": 1}
                worker_cfg["likelihood_expressions"] = []
                worker_cfg["handoff_to_staging"] = False
                worker_cfg["pull_timeout"] = 1

                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **worker_cfg)
                factory.start_watchdog(
                    stale_sec=30.0,
                    poll_interval_sec=0.1,
                    max_sample_retries=2,
                )
                assert factory.redis is not None
                factory.redis.push_task(
                    {
                        "uuid": "slow-sample-1",
                        "u_coords": [0.1, 0.2, 0.3],
                        "execution_plan": [
                            {"name": "SlowA", "type": "calculator", "layer": 0},
                        ],
                    }
                )

                _wait_until(
                    lambda: _worker_holds_calc_slot(factory, "SlowA"),
                    timeout=15.0,
                )
                busy_worker = next(
                    worker
                    for worker in factory.workers
                    if worker.is_alive() and worker.pid is not None
                )
                os.kill(busy_worker.pid, signal.SIGKILL)
                busy_worker.join(timeout=5.0)

                _wait_until(
                    lambda: factory.redis is not None
                    and int(
                        factory.redis.fetch_calculator_status().get(
                            calc_status_busy_field("SlowA"),
                            0,
                        )
                        or 0
                    )
                    == 0,
                    timeout=10.0,
                )
                status = factory.redis.fetch_calculator_status()
                self.assertEqual(int(status.get(calc_status_busy_field("SlowA"), 0) or 0), 0)
                self.assertGreaterEqual(
                    int(status.get(calc_status_free_field("SlowA"), 0) or 0),
                    1,
                )
                factory.shutdown()
        finally:
            server.shutdown()
            server.server_close()


class IdleWatchdogTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def _alive_stub(self, worker_id: int = 0, *, spawned_at: float | None = None):
        return SimpleNamespace(
            worker_id=worker_id,
            pid=4242,
            is_alive=lambda: True,
            _spawned_at=time.time() if spawned_at is None else spawned_at,
        )

    def _inspect(self, factory: TaskFactory) -> list[str]:
        reasons: list[str] = []

        def _record(worker, *, reason: str) -> None:
            reasons.append(reason)

        factory._watchdog.handle_worker_failure = _record  # type: ignore[method-assign]
        factory._watchdog.inspect_workers()
        return reasons

    def test_idle_stale_heartbeat_is_recovered(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub()]
        factory._watchdog.stale_sec = 1.0
        old = time.time() - 10.0
        queue.heartbeat("0", status="idle", last_heartbeat=old, ts=old)
        self.assertEqual(self._inspect(factory), ["stale_heartbeat"])

    def test_idle_fresh_empty_inflight_is_not_killed(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub()]
        factory._watchdog.stale_sec = 1.0
        now = time.time()
        queue.heartbeat("0", status="idle", last_heartbeat=now, ts=now)
        self.assertEqual(self._inspect(factory), [])

    def test_idle_fresh_with_inflight_is_inflight_without_busy(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub()]
        factory._watchdog.stale_sec = 30.0
        now = time.time()
        queue.heartbeat("0", status="idle", last_heartbeat=now, ts=now)
        assert queue.r is not None
        queue.r.rpush(INFLIGHT.format(worker="0"), '{"uuid":"orphan-1"}')
        self.assertEqual(self._inspect(factory), ["inflight_without_busy"])

    def test_starting_fresh_with_inflight_is_inflight_without_busy(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub()]
        now = time.time()
        queue.heartbeat("0", status="starting", last_heartbeat=now, ts=now)
        assert queue.r is not None
        queue.r.rpush(INFLIGHT.format(worker="0"), '{"uuid":"orphan-2"}')
        self.assertEqual(self._inspect(factory), ["inflight_without_busy"])

    def test_missing_heartbeat_past_spawned_at_is_stale(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub(spawned_at=time.time() - 10.0)]
        factory._watchdog.stale_sec = 1.0
        self.assertEqual(self._inspect(factory), ["stale_heartbeat"])

    def test_missing_heartbeat_within_spawn_grace_is_not_killed(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub(spawned_at=time.time())]
        factory._watchdog.stale_sec = 30.0
        self.assertEqual(self._inspect(factory), [])

    def test_busy_fresh_heartbeat_is_not_killed(self) -> None:
        factory = TaskFactory({})
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory.workers = [self._alive_stub()]
        factory._watchdog.stale_sec = 1.0
        now = time.time()
        queue.heartbeat("0", status="busy", last_heartbeat=now, ts=now)
        self.assertEqual(self._inspect(factory), [])


class LongCalculatorHeartbeatTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_long_calculator_survives_short_stale_threshold(self) -> None:
        """The periodic heartbeat thread must keep a Worker alive while
        module.execute() runs longer than Watchdog.stale_sec; a false
        recovery would hand its PackID directory to a second owner."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                script = os.path.join(tmpdir, "very_slow.py")
                with open(script, "w", encoding="utf-8") as handle:
                    handle.write(
                        "#!/usr/bin/env python3\n"
                        "import json, sys, time\n"
                        "time.sleep(3.0)\n"
                        "open(sys.argv[1], 'w').write(json.dumps({'a': 1}))\n"
                    )
                slow_module = _slow_calc_module("VerySlow", script, "a")
                worker_cfg = _worker_config(tmpdir)
                worker_cfg["calculator_modules"] = [slow_module]
                worker_cfg["calculator_pools"] = {"VerySlow": 1}
                worker_cfg["likelihood_expressions"] = []
                worker_cfg["handoff_to_staging"] = False
                worker_cfg["pull_timeout"] = 1
                worker_cfg["heartbeat_interval_sec"] = 0.25

                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **worker_cfg)
                factory.start_watchdog(
                    stale_sec=1.0,
                    poll_interval_sec=0.2,
                    max_sample_retries=2,
                )
                assert factory.redis is not None
                factory.redis.push_task(
                    {
                        "uuid": "long-calc-1",
                        "u_coords": [0.1, 0.2, 0.3],
                        "execution_plan": [
                            {"name": "VerySlow", "type": "calculator", "layer": 0},
                        ],
                    }
                )

                completed = _drain_archive_uuids(factory.redis, 1, timeout=30.0)
                self.assertEqual(completed, ["long-calc-1"])
                self.assertEqual(
                    factory._respawn_count,
                    0,
                    "watchdog must not falsely recover a heartbeating Worker",
                )
                factory.shutdown()
        finally:
            server.shutdown()
            server.server_close()


class WorkerFailureOrderingTests(unittest.TestCase):
    def test_kill_precedes_slot_sweep(self) -> None:
        """Stable PackIDs alias shadow directories, so a stale-heartbeat Worker
        must be force-stopped BEFORE its held slots return to the free list."""
        order: list[str] = []

        class _StubWorker:
            def __init__(self, worker_id: int, *_args: Any, **_kwargs: Any) -> None:
                self.worker_id = worker_id
                self.pid = 4242

            def start(self) -> None:
                order.append("respawn")

            def is_alive(self) -> bool:
                return False

        dead_worker = SimpleNamespace(worker_id=7, pid=1234)
        redis_stub = SimpleNamespace(
            decode_heartbeat_subprocess_pids=lambda heartbeat: [999999],
            decode_heartbeat_held_packs=lambda heartbeat: {"SlowA": "001"},
            sweep_held_calc_slots=lambda held: (order.append("sweep"), 1)[1],
        )
        fake_factory = SimpleNamespace(
            _recovery_lock=threading.Lock(),
            _last_recovered_pid={},
            redis=redis_stub,
            _force_stop_worker=lambda worker: order.append("stop"),
            _kill_orphan_process_groups=lambda pids: (order.append("killpg"), 0)[1],
            _worker_heartbeat=lambda worker_id: {},
            _requeue_in_flight_task=lambda heartbeat: False,
            workers=[dead_worker],
            _redis_connection_config={},
            _worker_spawn_template={},
            _respawn_count=0,
            _peak_workers_alive=0,
            _alive_workers=lambda: [],
            _logger=logging.getLogger("test.worker_failure"),
        )

        with mock.patch.object(factory_module, "Worker", _StubWorker):
            TaskFactory._handle_worker_failure(
                fake_factory, dead_worker, reason="stale_heartbeat"
            )

        self.assertEqual(order, ["stop", "killpg", "sweep", "respawn"])
        replacement = fake_factory.workers[0]
        self.assertTrue(hasattr(replacement, "_spawned_at"))
        self.assertGreater(replacement._spawned_at, 0)

    def test_kill_orphan_process_groups_reaps_setsid_child(self) -> None:
        """A child in its own session survives its parent's SIGKILL; the
        watchdog must be able to reap it via the heartbeat-recorded PID."""
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(60)"],
            start_new_session=True,
        )
        try:
            killed = TaskFactory._kill_orphan_process_groups([child.pid])
            self.assertEqual(killed, 1)
            self.assertIsNotNone(child.wait(timeout=5.0))
        finally:
            if child.poll() is None:  # pragma: no cover - cleanup on failure
                child.kill()

    def test_kill_orphan_process_groups_skips_non_leaders_and_dead_pids(self) -> None:
        # A child WITHOUT start_new_session shares our process group, so its
        # pid is not a group id and must be skipped (pid-recycling guard).
        child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(30)"])
        try:
            self.assertEqual(TaskFactory._kill_orphan_process_groups([child.pid]), 0)
            self.assertIsNone(child.poll())
        finally:
            child.kill()
            child.wait(timeout=5.0)
        self.assertEqual(TaskFactory._kill_orphan_process_groups([999999999]), 0)
        self.assertEqual(TaskFactory._kill_orphan_process_groups([]), 0)


def _wait_ps_command_prefix(pid: int, prefix: str, *, timeout: float = 5.0) -> str:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        completed = subprocess.run(
            ["ps", "-ax", "-o", "pid=,command="],
            capture_output=True,
            text=True,
            check=False,
        )
        for line in (completed.stdout or "").splitlines():
            parts = line.strip().split(None, 1)
            if len(parts) != 2:
                continue
            try:
                found = int(parts[0])
            except ValueError:
                continue
            if found == int(pid) and parts[1].startswith(prefix):
                return parts[1]
        time.sleep(0.05)
    raise TimeoutError(f"ps command= for pid {pid} never started with {prefix!r}")


class ChildrenBoardOrphanTests(unittest.TestCase):
    def test_kill_orphan_comments_forbid_proc_comm(self) -> None:
        from jarvishep2.runtime.factory import _Watchdog

        source = inspect.getsource(_Watchdog.kill_orphan_from_children_board)
        source += inspect.getsource(_Watchdog.kill_orphan_process_groups)
        self.assertIn("/proc", source)
        self.assertIn("comm", source)
        self.assertIn("16", source)
        self.assertIn("ps", source)
        self.assertIn("command=", source)
        self.assertNotIn('"/proc/', source)
        self.assertNotIn("'/proc/", source)
        self.assertNotIn("os.kill(", source)
        self.assertIn("os.killpg(", source)
        self.assertNotIn("_signal_process_tree", source)

    def test_ps_timeout_skips_file_operation_killpg(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            start_new_session=True,
        )
        try:
            with (
                mock.patch(
                    "jarvishep2.runtime.factory.subprocess.run",
                    side_effect=subprocess.TimeoutExpired("ps", 2),
                ),
                mock.patch.object(factory_module.os, "killpg") as killpg,
            ):
                killed = TaskFactory._kill_orphan_from_children_board(
                    {
                        "file_operation_pid": child.pid,
                        "file_operation_pgid": child.pid,
                        "calc_pgids": [],
                    }
                )
            self.assertEqual(killed, 0)
            killpg.assert_not_called()
            self.assertIsNone(child.poll())
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_ps_timeout_still_kills_calculator_session_leader(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            start_new_session=True,
        )
        try:
            with mock.patch(
                "jarvishep2.runtime.factory.subprocess.run",
                side_effect=subprocess.TimeoutExpired("ps", 2),
            ):
                killed = TaskFactory._kill_orphan_from_children_board(
                    {
                        "file_operation_pgid": "",
                        "calc_pgids": [child.pid],
                    }
                )
            self.assertEqual(killed, 1)
            self.assertIsNotNone(child.wait(timeout=5.0))
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_empty_pgid_is_killpg_noop(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            start_new_session=True,
        )
        try:
            with (
                mock.patch.object(factory_module.os, "kill") as kill_one,
                mock.patch.object(factory_module.os, "killpg") as killpg,
            ):
                killed = TaskFactory._kill_orphan_from_children_board(
                    {
                        "file_operation_pid": child.pid,
                        "file_operation_pgid": "",
                        "calc_pgids": [],
                    }
                )
            self.assertEqual(killed, 0)
            kill_one.assert_not_called()
            killpg.assert_not_called()
            self.assertIsNone(child.poll())
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_empty_pgid_without_live_child_returns_zero(self) -> None:
        self.assertEqual(
            TaskFactory._kill_orphan_from_children_board(
                {
                    "file_operation_pid": 123,
                    "file_operation_pgid": "",
                    "calc_pgids": [],
                }
            ),
            0,
        )

    def test_file_operation_without_title_is_not_killed(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            start_new_session=True,
        )
        try:
            killed = TaskFactory._kill_orphan_from_children_board(
                {
                    "file_operation_pid": child.pid,
                    "file_operation_pgid": child.pid,
                    "calc_pgids": [],
                }
            )
            self.assertEqual(killed, 0)
            self.assertIsNone(child.poll())
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_setproctitle_file_operation_session_leader_is_killed(self) -> None:
        child = subprocess.Popen(
            [
                sys.executable,
                "-c",
                "import time, setproctitle; "
                "setproctitle.setproctitle('Jarvis-FileOperation:scan'); "
                "time.sleep(60)",
            ],
            start_new_session=True,
        )
        try:
            _wait_ps_command_prefix(child.pid, "Jarvis-FileOperation")
            killed = TaskFactory._kill_orphan_from_children_board(
                {
                    "file_operation_pid": child.pid,
                    "file_operation_pgid": child.pid,
                    "calc_pgids": [],
                }
            )
            self.assertEqual(killed, 1)
            self.assertIsNotNone(child.wait(timeout=5.0))
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_calculator_session_leader_killed_without_title_check(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"],
            start_new_session=True,
        )
        try:
            killed = TaskFactory._kill_orphan_from_children_board(
                {
                    "file_operation_pgid": "",
                    "calc_pgids": [child.pid],
                }
            )
            self.assertEqual(killed, 1)
            self.assertIsNotNone(child.wait(timeout=5.0))
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_calculator_non_leader_is_not_killed(self) -> None:
        child = subprocess.Popen(
            [sys.executable, "-c", "import time; time.sleep(30)"]
        )
        try:
            with mock.patch.object(factory_module.os, "kill") as kill_one:
                killed = TaskFactory._kill_orphan_from_children_board(
                    {
                        "file_operation_pgid": "",
                        "calc_pgids": [child.pid],
                    }
                )
            self.assertEqual(killed, 0)
            kill_one.assert_not_called()
            self.assertIsNone(child.poll())
        finally:
            if child.poll() is None:
                child.kill()
                child.wait(timeout=5.0)

    def test_handle_worker_failure_prefers_children_board_pgids(self) -> None:
        order: list[str] = []

        class _StubWorker:
            def __init__(self, worker_id: int, *_args: Any, **_kwargs: Any) -> None:
                self.worker_id = worker_id
                self.pid = 4242

            def start(self) -> None:
                order.append("respawn")

            def is_alive(self) -> bool:
                return False

        dead_worker = SimpleNamespace(worker_id=3, pid=2222)
        redis_stub = SimpleNamespace(
            read_children_board=lambda worker_id: {
                "file_operation_pgid": "",
                "calc_pgids": [111],
            },
            decode_heartbeat_subprocess_pids=lambda heartbeat: (order.append("hb-pids"), [999])[1],
            decode_heartbeat_held_packs=lambda heartbeat: {},
            sweep_held_calc_slots=lambda held: (order.append("sweep"), 0)[1],
        )
        fake_factory = SimpleNamespace(
            _recovery_lock=threading.Lock(),
            _last_recovered_pid={},
            redis=redis_stub,
            _force_stop_worker=lambda worker: order.append("stop"),
            _kill_orphan_from_children_board=lambda board: (
                order.append("killpg-board"),
                0,
            )[1],
            _kill_orphan_process_groups=lambda pids: (order.append("killpg-hb"), 0)[1],
            _worker_heartbeat=lambda worker_id: {},
            _requeue_in_flight_task=lambda heartbeat: False,
            workers=[dead_worker],
            _redis_connection_config={},
            _worker_spawn_template={},
            _respawn_count=0,
            _peak_workers_alive=0,
            _alive_workers=lambda: [],
            _logger=logging.getLogger("test.worker_failure"),
        )

        with mock.patch.object(factory_module, "Worker", _StubWorker):
            TaskFactory._handle_worker_failure(
                fake_factory, dead_worker, reason="process_exit"
            )

        self.assertEqual(order, ["stop", "killpg-board", "sweep", "respawn"])
        self.assertNotIn("killpg-hb", order)
        self.assertNotIn("hb-pids", order)


if __name__ == "__main__":
    unittest.main()