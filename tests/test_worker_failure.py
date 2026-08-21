#!/usr/bin/env python3
"""WP-D6.1 Worker failure recovery tests."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.slow

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
    RedisQueue,
    calc_status_busy_field,
    calc_status_free_field,
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


if __name__ == "__main__":
    unittest.main()