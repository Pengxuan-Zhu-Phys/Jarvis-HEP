"""Deterministic watchdog races; kept in the default fast gate."""

from __future__ import annotations

import threading
import time
from types import SimpleNamespace
from unittest import mock

import pytest

import jarvishep2.runtime.factory as factory_module
from jarvishep2.redis_queue import INFLIGHT, PROC_WORKER, WORKER_STATUS, make_fakeredis_queue
from jarvishep2.runtime.factory import TaskFactory
from jarvishep2.runtime.worker import Worker


@pytest.fixture
def rig(monkeypatch):
    factory = TaskFactory({})
    factory.redis = make_fakeredis_queue()
    worker = SimpleNamespace(worker_id=0, pid=4242, _spawned_at=time.time(), is_alive=lambda: True)
    factory.workers = [worker]
    now = [100.0]
    monkeypatch.setattr(factory_module, "time", SimpleNamespace(
        time=time.time, monotonic=lambda: now[0], sleep=time.sleep,
    ))
    factory.redis.heartbeat("0", status="idle", pid=worker.pid, ts=time.time())
    return factory, worker, now


def plant(factory, uuid="a", retry=0):
    key = INFLIGHT.format(worker="0")
    factory.redis.r.delete(key)
    factory.redis.r.rpush(key, factory.redis.encode_task_for_heartbeat({"uuid": uuid, "_retry_count": retry, "u_coords": [0.1]}))


def inspect(rig):
    factory, _, _ = rig
    with mock.patch.object(factory._watchdog, "handle_worker_failure") as recover:
        factory._inspect_workers()
    return [call.kwargs["reason"] for call in recover.call_args_list]


def test_different_tasks_never_accumulate_idle_time(rig):
    factory, _, now = rig
    for index in range(40):
        plant(factory, f"task-{index}")
        now[0] += 2.0
        assert inspect(rig) == []


def test_same_task_recovers_only_after_full_grace(rig):
    factory, _, now = rig
    plant(factory)
    assert inspect(rig) == []
    now[0] += 29.9
    assert inspect(rig) == []
    now[0] += 0.1
    assert inspect(rig) == ["inflight_without_busy"]


@pytest.mark.parametrize("change", ["retry", "pid", "spawn"])
def test_new_attempt_or_process_gets_its_own_grace(rig, change):
    factory, worker, now = rig
    plant(factory)
    assert inspect(rig) == []
    now[0] += 29.0
    if change == "retry":
        plant(factory, retry=1)
    elif change == "pid":
        worker.pid += 1
        factory.redis.heartbeat("0", status="idle", pid=worker.pid, ts=time.time())
    else:
        worker._spawned_at += 1.0
    assert inspect(rig) == []
    now[0] += 2.0
    assert inspect(rig) == []


@pytest.mark.parametrize("change", ["busy", "empty", "unknown", "read_error"])
def test_interrupted_observations_reset_grace(rig, change):
    factory, _, now = rig
    plant(factory)
    assert inspect(rig) == []
    now[0] += 29.0
    if change == "busy":
        factory.redis.heartbeat("0", status="busy", ts=time.time())
    elif change == "empty":
        factory.redis.r.delete(INFLIGHT.format(worker="0"))
    elif change == "unknown":
        plant(factory, uuid="")
    if change == "read_error":
        with mock.patch.object(factory.redis, "_watchdog_snapshot", side_effect=TimeoutError("offline")):
            assert inspect(rig) == []
    else:
        assert inspect(rig) == []
    assert 0 not in factory._watchdog._inflight_idle_since
    factory.redis.heartbeat("0", status="idle", ts=time.time())
    plant(factory)
    now[0] += 2.0
    assert inspect(rig) == []


def test_stale_heartbeat_is_not_masked_by_changing_inflight(rig):
    factory, _, now = rig
    plant(factory)
    assert inspect(rig) == []
    plant(factory, "b")
    now[0] += 2.0
    factory.redis.heartbeat("0", status="idle", ts=time.time() - 61.0)
    assert inspect(rig) == ["stale_heartbeat"]


def test_prior_process_heartbeat_cannot_mark_replacement_stale(rig):
    factory, _, _ = rig
    factory.redis.heartbeat("0", status="idle", pid=1, ts=time.time() - 61.0)
    plant(factory)
    assert inspect(rig) == []
    assert 0 not in factory._watchdog._inflight_idle_since


def test_elapsed_time_ignores_wall_clock_jump(rig, monkeypatch):
    factory, _, now = rig
    plant(factory)
    assert inspect(rig) == []
    wall = time.time()
    monkeypatch.setattr(factory_module.time, "time", lambda: wall - 3600.0)
    now[0] += 30.0
    assert inspect(rig) == ["inflight_without_busy"]


@pytest.mark.parametrize("change", ["busy", "empty", "task", "retry", "replacement", "read_error"])
def test_final_check_cancels_recovery_without_touching_queue(rig, change):
    factory, worker, now = rig
    plant(factory)
    assert inspect(rig) == []
    now[0] += 31.0
    # The first snapshot authorizes recovery; state changes while waiting for
    # the recovery lock, before the last snapshot and any SIGKILL/reclaim.
    class ChangeOnLock:
        def __enter__(self):
            if change == "busy":
                factory.redis.heartbeat("0", status="busy", ts=time.time())
            elif change == "empty":
                factory.redis.r.delete(INFLIGHT.format(worker="0"))
            elif change == "task":
                plant(factory, "next")
            elif change == "retry":
                plant(factory, retry=1)
            elif change == "replacement":
                factory.workers = [SimpleNamespace(**vars(worker))]
            elif change == "read_error":
                factory.redis._watchdog_snapshot = mock.Mock(side_effect=TimeoutError("offline"))

        def __exit__(self, *args):
            pass

    factory._recovery_lock = ChangeOnLock()
    with (
        mock.patch.object(factory, "_force_stop_worker") as kill,
        mock.patch.object(factory, "_requeue_in_flight_task") as reclaim,
    ):
        factory._inspect_workers()
    kill.assert_not_called()
    reclaim.assert_not_called()
    assert factory._respawn_count == 0


def test_confirmed_stall_requeues_once_and_logs_actual_identity(rig, caplog):
    factory, worker, now = rig
    factory._scan_mode = "stopping"  # Exercise recovery without spawning a real process.
    plant(factory)
    assert inspect(rig) == []
    now[0] += 30.0
    with mock.patch.object(factory, "_force_stop_worker") as kill:
        factory._inspect_workers()
        factory._handle_worker_failure(worker, reason="inflight_without_busy")
    kill.assert_called_once_with(worker)
    task = factory.redis.pull_task(timeout=0)
    assert task["uuid"] == "a"
    assert task["_retry_count"] == 1
    assert factory.redis.get_inflight_task("0") is None
    assert "inflight_uuid=a" in caplog.text
    assert "requeued_uuids=['a']" in caplog.text


def test_snapshot_is_read_only_and_uses_one_transaction(rig):
    factory, worker, _ = rig
    plant(factory)
    ctrl = factory.redis._ctrl()
    before = {key: ctrl.dump(key) for key in ctrl.keys("*")}
    with mock.patch.object(ctrl, "pipeline", wraps=ctrl.pipeline) as pipeline:
        snapshot = factory.redis._watchdog_snapshot("0")
    pipeline.assert_called_once_with(transaction=True)
    assert snapshot["task"]["uuid"] == "a"
    assert snapshot["heartbeat"]["pid"] == str(worker.pid)
    assert snapshot["inflight_count"] == 1
    assert {key: ctrl.dump(key) for key in ctrl.keys("*")} == before
    # Board expiry falls back to the legacy heartbeat from the SAME transaction.
    ctrl.delete(PROC_WORKER.format(id="0"))
    assert factory.redis._watchdog_snapshot("0")["heartbeat"] == ctrl.hgetall(WORKER_STATUS.format(id="0"))


def test_heartbeat_serializes_snapshot_and_write_without_holding_state_lock():
    worker = Worker(0, {"host": "127.0.0.1", "port": 1}, {})
    worker._hb_lock()  # Child-side initialization before starting threads.
    idle_entered = threading.Event()
    release_idle = threading.Event()
    busy_attempted = threading.Event()
    busy_published = threading.Event()
    writes = []
    failures = []

    def publish(wid, **fields):
        try:
            # Calculator ownership updates may take this lock during Redis I/O.
            assert worker._hb_lock().acquire(blocking=False)
            worker._hb_lock().release()
            if fields["status"] == "idle":
                idle_entered.set()
                assert release_idle.wait(3.0)
            writes.append(fields)
            if fields["status"] == "busy":
                busy_published.set()
        except BaseException as exc:
            failures.append(exc)

    def busy():
        busy_attempted.set()
        with worker._hb_lock():
            worker._current_sample_uuid = "next"
        worker._heartbeat("busy")

    worker._redis = SimpleNamespace(heartbeat=publish)
    idle_thread = threading.Thread(target=worker._heartbeat, args=("idle",))
    busy_thread = threading.Thread(target=busy)
    idle_thread.start()
    try:
        assert idle_entered.wait(3.0)
        busy_thread.start()
        assert busy_attempted.wait(3.0)
        # An already prepared idle write cannot be overtaken by busy.
        assert not busy_published.wait(0.05)
    finally:
        release_idle.set()
        idle_thread.join(3.0)
        if busy_thread.ident is not None:
            busy_thread.join(3.0)
    assert not idle_thread.is_alive() and not busy_thread.is_alive()
    assert not failures
    assert [row["status"] for row in writes] == ["idle", "busy"]
    assert writes[0]["current_sample"] is None
    assert writes[1]["current_sample"] == "next"


def test_default_factory_policy_matches_normalized_config():
    factory = TaskFactory({})
    assert factory._watchdog.stale_sec == 60.0
    assert factory._watchdog.poll_interval_sec == 2.0
    assert factory._watchdog.inflight_idle_grace_sec == 30.0
    with mock.patch.object(factory._watchdog, "start") as start:
        factory.start_watchdog()
    assert start.call_args.kwargs["stale_sec"] == 60.0
    assert start.call_args.kwargs["poll_interval_sec"] == 2.0
    assert start.call_args.kwargs["inflight_idle_grace_sec"] == 30.0
