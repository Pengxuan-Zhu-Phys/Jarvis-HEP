#!/usr/bin/env python3
"""TaskFactory — Worker lifecycle manager for Jarvis-HEP V2.

D9.4: composition of private ``_MonitorLoop`` / ``_Watchdog`` collaborators;
explicit construction preferred over the deprecated singleton shell.
"""

from __future__ import annotations

import copy
import os
import signal
import threading
import time
from collections.abc import Callable
from typing import Any

from jarvishep2.calculator_pools import register_calculator_pools
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.worker import Worker


class _MonitorLoop:
    """Op_count-gated snapshot collector + background updater (D5.1 / D9.4)."""

    def __init__(self, factory: TaskFactory) -> None:
        self._factory = factory
        self._snapshot: dict[str, Any] = {}
        self._snapshot_lock = threading.RLock()
        self._last_op_counts: dict[str, int] = {}
        self._thread: threading.Thread | None = None
        self._running = False
        self._logger = get_jarvis_logger("factory.monitor")

    def get_snapshot(self) -> dict[str, Any]:
        with self._snapshot_lock:
            return copy.deepcopy(self._snapshot)

    def clear(self) -> None:
        self._last_op_counts.clear()
        with self._snapshot_lock:
            self._snapshot = {}

    def _worker_status(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for worker in self._factory.workers:
            rows.append(
                {
                    "worker_id": worker.worker_id,
                    "pid": worker.pid,
                    "alive": worker.is_alive(),
                    "name": worker.name,
                }
            )
        return rows

    def _fetch_workers_redis(self) -> dict[str, dict[str, Any]]:
        redis = self._factory.redis
        if redis is None:
            return {}
        worker_ids = [str(worker.worker_id) for worker in self._factory.workers]
        return redis.fetch_worker_status(worker_ids)

    @staticmethod
    def _carry_forward_section(
        prev: dict[str, Any],
        key: str,
        *,
        refresh: bool,
        supplier: Callable[[], Any],
    ) -> Any:
        if refresh:
            return supplier()
        value = prev.get(key)
        if value is None:
            return supplier()
        return copy.deepcopy(value)

    def _subsystem_refresh_needed(
        self,
        kind: str,
        op_counts: dict[str, int],
        prev: dict[str, Any],
        section_key: str,
    ) -> bool:
        current = int(op_counts.get(kind, 0))
        last_seen = self._last_op_counts.get(kind, -1)
        return current > last_seen or section_key not in prev

    def collect_latest_status(self) -> dict[str, Any]:
        """Collect queue/stats snapshot from Redis with op_count gating."""
        factory = self._factory
        workers_alive = sum(1 for worker in factory.workers if worker.is_alive())
        factory._peak_workers_alive = max(factory._peak_workers_alive, workers_alive)
        snap: dict[str, Any] = {
            "timestamp": time.time(),
            "workers": self._worker_status(),
            "workers_alive": workers_alive,
            "workers_total": len(factory.workers),
        }
        redis = factory.redis
        if redis is None:
            return snap

        with self._snapshot_lock:
            prev = self._snapshot

        op_counts = redis.get_all_op_counts()
        snap["op_counts"] = op_counts

        lengths = redis.get_queue_lengths()
        snap["task_queue_length"] = lengths["task_queue_length"]
        snap["archive_queue_length"] = lengths["archive_queue_length"]

        snap["worker_heartbeats"] = self._carry_forward_section(
            prev,
            "worker_heartbeats",
            refresh=self._subsystem_refresh_needed(
                "worker", op_counts, prev, "worker_heartbeats"
            ),
            supplier=self._fetch_workers_redis,
        )
        snap["calculator_status"] = self._carry_forward_section(
            prev,
            "calculator_status",
            refresh=self._subsystem_refresh_needed(
                "calculator", op_counts, prev, "calculator_status"
            ),
            supplier=redis.fetch_calculator_status,
        )
        snap["sample_stats"] = self._carry_forward_section(
            prev,
            "sample_stats",
            refresh=self._subsystem_refresh_needed(
                "sample", op_counts, prev, "sample_stats"
            ),
            supplier=redis.fetch_sample_stats,
        )

        self._last_op_counts = dict(op_counts)
        return snap

    def start(self, *, update_hz: float = 120.0) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._running = True
        interval = 1.0 / max(1.0, float(update_hz))

        def _loop() -> None:
            while self._running:
                try:
                    latest = self.collect_latest_status()
                    with self._snapshot_lock:
                        self._snapshot = latest
                except Exception as exc:
                    self._logger.warning("monitor snapshot update failed -> %s", exc)
                time.sleep(interval)

        self._thread = threading.Thread(
            target=_loop,
            name="Jarvis2-FactoryMonitor",
            daemon=True,
        )
        self._thread.start()

    def stop(self, *, join_timeout: float = 2.0) -> None:
        self._running = False
        if self._thread is not None:
            self._thread.join(timeout=join_timeout)
            if self._thread.is_alive():
                self._logger.warning("monitor thread did not stop within timeout")
            self._thread = None


class _Watchdog:
    """Worker liveness watchdog: stale heartbeat / process exit → recover (D6.1 / D9.4)."""

    def __init__(self, factory: TaskFactory) -> None:
        self._factory = factory
        self._thread: threading.Thread | None = None
        self._running = False
        self.stale_sec = 30.0
        self.poll_interval_sec = 1.0
        self.max_sample_retries = 3
        self._logger = get_jarvis_logger("factory.watchdog")

    def start(
        self,
        *,
        enabled: bool = True,
        stale_sec: float = 30.0,
        poll_interval_sec: float = 1.0,
        max_sample_retries: int = 3,
    ) -> None:
        if not enabled:
            return
        if self._thread is not None and self._thread.is_alive():
            return
        self.stale_sec = max(1.0, float(stale_sec))
        self.poll_interval_sec = max(0.1, float(poll_interval_sec))
        self.max_sample_retries = max(0, int(max_sample_retries))
        self._running = True

        def _loop() -> None:
            while self._running:
                try:
                    self.inspect_workers()
                except Exception as exc:
                    self._logger.warning("watchdog inspection failed -> %s", exc)
                time.sleep(self.poll_interval_sec)

        self._thread = threading.Thread(
            target=_loop,
            name="Jarvis2-FactoryWatchdog",
            daemon=True,
        )
        self._thread.start()

    def stop(self, *, join_timeout: float = 2.0) -> None:
        self._running = False
        if self._thread is not None:
            self._thread.join(timeout=join_timeout)
            if self._thread.is_alive():
                self._logger.warning("watchdog thread did not stop within timeout")
            self._thread = None

    @staticmethod
    def heartbeat_timestamp(heartbeat: dict[str, Any]) -> float:
        for key in ("last_heartbeat", "ts"):
            raw = heartbeat.get(key)
            if raw is None or raw == "":
                continue
            try:
                return float(raw)
            except (TypeError, ValueError):
                continue
        return 0.0

    def worker_heartbeat(self, worker_id: int) -> dict[str, Any]:
        redis = self._factory.redis
        if redis is None:
            return {}
        rows = redis.fetch_worker_status([str(worker_id)])
        return dict(rows.get(str(worker_id)) or {})

    def inspect_workers(self) -> None:
        for worker in list(self._factory.workers):
            if not worker.is_alive():
                self.handle_worker_failure(worker, reason="process_exit")
                continue
            heartbeat = self.worker_heartbeat(worker.worker_id)
            status = str(heartbeat.get("status") or "").strip().lower()
            if status not in {"busy", "starting"}:
                continue
            last_seen = self.heartbeat_timestamp(heartbeat)
            if last_seen <= 0:
                continue
            if (time.time() - last_seen) <= self.stale_sec:
                continue
            self.handle_worker_failure(worker, reason="stale_heartbeat")

    @staticmethod
    def force_stop_worker(worker: Worker) -> None:
        if not worker.is_alive():
            return
        if worker.pid is not None:
            try:
                os.kill(worker.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        worker.join(timeout=2.0)
        if worker.is_alive():
            worker.terminate()
            worker.join(timeout=2.0)

    def requeue_in_flight_task(self, heartbeat: dict[str, Any]) -> bool:
        redis = self._factory.redis
        if redis is None:
            return False
        task = redis.decode_heartbeat_task(heartbeat)
        if task is None:
            return False
        retry_count = int(task.get("_retry_count", 0) or 0)
        if retry_count >= self.max_sample_retries:
            sample_uuid = str(task.get("uuid") or heartbeat.get("current_sample") or "")
            if sample_uuid:
                redis.submit_result(
                    {
                        "uuid": sample_uuid,
                        "status": "Failed",
                        "observables": {},
                        "error": "worker_failure_retries_exhausted",
                    }
                )
            return False
        task["_retry_count"] = retry_count + 1
        redis.push_task(task)
        return True

    @staticmethod
    def kill_orphan_process_groups(pids: list[int]) -> int:
        """SIGKILL calculator process groups orphaned by a dead Worker.

        Children run with ``start_new_session=True``, so killing the Worker
        does not kill them; left alone they would keep writing into a PackID
        shadow directory that is about to be handed to a new owner. Only
        session leaders are targeted (``getpgid(pid) == pid``) so a recycled
        OS pid can never match an unrelated process.
        """
        killed = 0
        for pid in pids:
            try:
                if os.getpgid(pid) != pid:
                    continue
                os.killpg(pid, signal.SIGKILL)
                killed += 1
            except (ProcessLookupError, PermissionError, OSError):
                continue
        return killed

    def handle_worker_failure(self, worker: Worker, *, reason: str) -> None:
        # Lifecycle recovery stays on TaskFactory so duck-typed unit stubs
        # (``TaskFactory._handle_worker_failure(fake, worker, ...)``) keep working.
        self._factory._handle_worker_failure(worker, reason=reason)


class TaskFactory:
    """Worker lifecycle + monitor + watchdog for Jarvis-HEP V2.

    Prefer constructing an instance and holding it on ``Jarvis2Core.factory``.
    ``get_instance`` remains as a deprecated compatibility shell for older tests
    and monitor attach.

    The Factory does **not** execute tasks, hold calculators, or own Sample
    objects. During a normal run it only manages Worker processes and serves
    read-only monitor snapshots (invariant #6 — one Sample per Worker).

    Internals (D9.4): ``_MonitorLoop`` owns snapshot collection; ``_Watchdog``
    owns dead-Worker recovery. Public methods stay on this facade.
    """

    _instance: TaskFactory | None = None
    _lock = threading.Lock()

    def __init__(self, redis_config: dict[str, Any] | None = None) -> None:
        self.redis_config = dict(redis_config or {})
        self.redis: RedisQueue | None = None
        self.workers: list[Worker] = []
        self._run_started_at: float | None = None
        self._peak_workers_alive = 0
        self._worker_spawn_template: dict[str, Any] = {}
        self._redis_connection_config: dict[str, Any] = {}
        self._recovery_lock = threading.Lock()
        self._last_recovered_pid: dict[int, int | None] = {}
        self._respawn_count = 0
        self._logger = get_jarvis_logger("factory")
        self._monitor = _MonitorLoop(self)
        self._watchdog = _Watchdog(self)

    @classmethod
    def get_instance(cls, redis_config: dict[str, Any] | None = None) -> TaskFactory:
        """Deprecated process-local singleton shell.

        Prefer ``TaskFactory(redis_config)`` held by the caller (e.g. Jarvis2Core).
        """
        with cls._lock:
            if cls._instance is None:
                cls._instance = cls(redis_config)
            elif redis_config and not cls._instance.redis_config:
                # Only fill empty config; never silently overwrite an active factory.
                cls._instance.redis_config.update(redis_config)
            return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset the deprecated singleton shell (tests / monitor attach)."""
        with cls._lock:
            cls._instance = None

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        """Build the control-process Redis client."""
        self.redis = RedisQueue(self.redis_config, client=client)
        self.redis.connect()
        return self.redis

    def _alive_workers(self) -> list[Worker]:
        return [worker for worker in self.workers if worker.is_alive()]

    def start_workers(self, n: int, **worker_kwargs: Any) -> list[Worker]:
        """Spawn ``n`` Worker processes via the spawn context.

        Each :class:`Worker` receives only picklable connection settings derived
        from the control-process queue; child Workers open their own Redis clients
        in ``run()``.

        Raises:
            RuntimeError: If Redis is not initialized or workers are already running.
        """
        if n <= 0:
            return []
        if self.redis is None:
            raise RuntimeError(
                "TaskFactory.init_redis() must be called before start_workers()"
            )

        alive = self._alive_workers()
        if alive:
            raise RuntimeError(
                f"start_workers refused: {len(alive)} worker process(es) already running"
            )
        self.workers = [worker for worker in self.workers if worker.is_alive()]

        redis_config = self.redis.connection_config()
        shared_config = dict(worker_kwargs)
        self._redis_connection_config = dict(redis_config)
        self._worker_spawn_template = copy.deepcopy(shared_config)
        register_calculator_pools(self.redis, shared_config)
        started: list[Worker] = []
        for worker_id in range(n):
            config = copy.deepcopy(shared_config)
            worker = Worker(worker_id, redis_config, config)
            worker.start()
            started.append(worker)
        self.workers.extend(started)
        if self._run_started_at is None:
            self._run_started_at = time.time()
        self._peak_workers_alive = max(
            self._peak_workers_alive, len(self._alive_workers())
        )
        self._logger.info("started %d worker process(es)", len(started))
        return started

    def get_run_metrics(self) -> dict[str, Any]:
        """Project run counters for run_summary (WP-D5.2). Honest ``None`` for untracked."""
        submitted = 0
        ok = 0
        failed = 0
        if self.redis is not None:
            stats = self.redis.fetch_sample_stats()
            lengths = self.redis.get_queue_lengths()
            ok = int(stats.get("completed", 0) or 0)
            failed = int(stats.get("failed", 0) or 0)
            running = int(stats.get("running", 0) or 0)
            queued = int(lengths.get("task_queue_length", 0) or 0)
            submitted = max(
                int(self.redis.get_op_count("task")),
                ok + failed + running + queued,
            )
        alive = sum(1 for worker in self.workers if worker.is_alive())
        return {
            "submitted": submitted,
            "ok": ok,
            "failed": failed,
            "configured_workers": len(self.workers),
            "peak_active_workers": max(self._peak_workers_alive, alive),
            # Instantaneous alive count is not a true mean; leave unset until tracked.
            "mean_active_workers": None,
            "total_point_eval_sec": None,
            "completed_durations_sec": None,
            "retry_count": None,
            "run_started_at": self._run_started_at,
        }

    def stop_all_workers(self, *, graceful: bool = True, join_timeout: float = 30.0) -> None:
        """Stop all tracked Workers; SIGTERM first, then force-terminate."""
        if not self.workers:
            return
        if graceful:
            self.request_worker_shutdown()
            for worker in self.workers:
                if worker.is_alive():
                    worker.join(timeout=join_timeout)
        for worker in self.workers:
            if worker.is_alive():
                worker.terminate()
                worker.join(timeout=5.0)
        self.workers.clear()

    def request_worker_shutdown(self) -> None:
        """Signal live Workers to stop after the current Sample."""
        for worker in self.workers:
            if worker.is_alive() and worker.pid is not None:
                try:
                    os.kill(worker.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    def get_monitor_snapshot(self) -> dict[str, Any]:
        """Return an in-memory deepcopy of the latest monitor snapshot."""
        return self._monitor.get_snapshot()

    # Compatibility aliases for pre-D9.4 tests that poked private fields.
    @property
    def _snapshot(self) -> dict[str, Any]:
        return self._monitor._snapshot

    @_snapshot.setter
    def _snapshot(self, value: dict[str, Any]) -> None:
        self._monitor._snapshot = value

    @property
    def _updater_thread(self) -> threading.Thread | None:
        return self._monitor._thread

    @property
    def _running(self) -> bool:
        return self._monitor._running

    @property
    def _watchdog_thread(self) -> threading.Thread | None:
        return self._watchdog._thread

    @property
    def _watchdog_running(self) -> bool:
        return self._watchdog._running

    @property
    def _watchdog_stale_sec(self) -> float:
        return self._watchdog.stale_sec

    @_watchdog_stale_sec.setter
    def _watchdog_stale_sec(self, value: float) -> None:
        self._watchdog.stale_sec = float(value)

    @property
    def _watchdog_poll_interval_sec(self) -> float:
        return self._watchdog.poll_interval_sec

    @_watchdog_poll_interval_sec.setter
    def _watchdog_poll_interval_sec(self, value: float) -> None:
        self._watchdog.poll_interval_sec = float(value)

    @property
    def _max_sample_retries(self) -> int:
        return self._watchdog.max_sample_retries

    @_max_sample_retries.setter
    def _max_sample_retries(self, value: int) -> None:
        self._watchdog.max_sample_retries = int(value)

    # --- thin facades kept for tests / call sites that reach into private methods ---
    def _collect_latest_status(self) -> dict[str, Any]:
        return self._monitor.collect_latest_status()

    def _worker_status(self) -> list[dict[str, Any]]:
        return self._monitor._worker_status()

    def _fetch_workers_redis(self) -> dict[str, dict[str, Any]]:
        return self._monitor._fetch_workers_redis()

    def _start_snapshot_updater(self, *, update_hz: float = 120.0) -> None:
        self._monitor.start(update_hz=update_hz)

    def start_monitor(self, *, update_hz: float = 120.0) -> None:
        """Launch the background snapshot updater thread."""
        self._monitor.start(update_hz=update_hz)

    def start_watchdog(
        self,
        *,
        enabled: bool = True,
        stale_sec: float = 30.0,
        poll_interval_sec: float = 1.0,
        max_sample_retries: int = 3,
    ) -> None:
        """Launch the Worker watchdog (WP-D6.1)."""
        self._watchdog.start(
            enabled=enabled,
            stale_sec=stale_sec,
            poll_interval_sec=poll_interval_sec,
            max_sample_retries=max_sample_retries,
        )

    @staticmethod
    def _heartbeat_timestamp(heartbeat: dict[str, Any]) -> float:
        return _Watchdog.heartbeat_timestamp(heartbeat)

    def _worker_heartbeat(self, worker_id: int) -> dict[str, Any]:
        watchdog = getattr(self, "_watchdog", None)
        if watchdog is not None:
            return watchdog.worker_heartbeat(worker_id)
        return {}

    def _inspect_workers(self) -> None:
        self._watchdog.inspect_workers()

    def _force_stop_worker(self, worker: Worker) -> None:
        _Watchdog.force_stop_worker(worker)

    def _requeue_in_flight_task(self, heartbeat: dict[str, Any]) -> bool:
        # Prefer composed watchdog when present; duck-typed stubs may override.
        watchdog = getattr(self, "_watchdog", None)
        if watchdog is not None:
            return watchdog.requeue_in_flight_task(heartbeat)
        return False

    @staticmethod
    def _kill_orphan_process_groups(pids: list[int]) -> int:
        return _Watchdog.kill_orphan_process_groups(pids)

    def _handle_worker_failure(self, worker: Worker, *, reason: str) -> None:
        """Kill → sweep slots → requeue task → respawn (D6.1 recovery).

        Implemented on the factory so unit tests can call
        ``TaskFactory._handle_worker_failure(stub, worker, ...)`` with a
        SimpleNamespace that supplies the same private hooks.
        """
        worker_id = int(worker.worker_id)
        dead_pid = worker.pid
        with self._recovery_lock:
            if self._last_recovered_pid.get(worker_id) == dead_pid:
                return
            if self.redis is None:
                return

            # Kill first, sweep after: with stable PackIDs a slot IS its shadow
            # directory, so returning a stale-heartbeat Worker's slot while the
            # Worker is still alive would let a new owner write into the same
            # directory concurrently.
            self._force_stop_worker(worker)

            heartbeat = self._worker_heartbeat(worker_id)
            orphan_pids = self.redis.decode_heartbeat_subprocess_pids(heartbeat)
            orphans_killed = self._kill_orphan_process_groups(orphan_pids)
            held_packs = self.redis.decode_heartbeat_held_packs(heartbeat)
            released = self.redis.sweep_held_calc_slots(held_packs)
            requeued = self._requeue_in_flight_task(heartbeat)
            self.workers = [
                item for item in self.workers if item.worker_id != worker_id
            ]
            replacement = Worker(
                worker_id,
                self._redis_connection_config,
                copy.deepcopy(self._worker_spawn_template),
            )
            replacement.start()
            self.workers.append(replacement)
            self._last_recovered_pid[worker_id] = dead_pid
            self._respawn_count += 1
            self._peak_workers_alive = max(
                self._peak_workers_alive, len(self._alive_workers())
            )
            logger = getattr(self, "_logger", None)
            if logger is not None:
                logger.warning(
                    "recovered dead worker %d (reason=%s, orphans_killed=%d, "
                    "released_slots=%d, requeued=%s, new_pid=%s)",
                    worker_id,
                    reason,
                    orphans_killed,
                    released,
                    requeued,
                    replacement.pid,
                )

    def shutdown(self, *, wait: bool = True) -> None:
        """Stop monitor, signal Workers, join, and close the Redis connection."""
        self._watchdog.stop()
        self._monitor.stop()
        self.request_worker_shutdown()
        self.stop_all_workers(graceful=wait)
        if self.redis is not None:
            self.redis.close()
            self.redis = None
        self._monitor.clear()
        self._run_started_at = None
        self._peak_workers_alive = 0
        self._worker_spawn_template.clear()
        self._redis_connection_config.clear()
        self._last_recovered_pid.clear()
        self._respawn_count = 0
        self._logger.info("TaskFactory shutdown complete")


__all__ = ["TaskFactory"]
