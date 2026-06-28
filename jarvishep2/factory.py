#!/usr/bin/env python3
"""TaskFactory — Worker lifecycle manager for Jarvis-HEP V2."""

from __future__ import annotations

import copy
import threading
import time
from typing import Any

from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.worker import Worker


class TaskFactory:
    """Process-local singleton that spawns and monitors Workers.

    The Factory does **not** execute tasks, hold calculators, or own Sample
    objects. During a normal run it only manages Worker processes and serves
    read-only monitor snapshots (invariant #6 — one Sample per Worker).
    """

    _instance: TaskFactory | None = None
    _lock = threading.Lock()

    def __init__(self, redis_config: dict[str, Any] | None = None) -> None:
        self.redis_config = dict(redis_config or {})
        self.redis: RedisQueue | None = None
        self.workers: list[Worker] = []
        self._snapshot: dict[str, Any] = {}
        self._snapshot_lock = threading.RLock()
        self.last_op_counts = {
            "worker": 0,
            "calculator": 0,
            "sample": 0,
            "task": 0,
        }
        self._updater_thread: threading.Thread | None = None
        self._running = False
        self._logger = get_jarvis_logger("factory")

    @classmethod
    def get_instance(cls, redis_config: dict[str, Any] | None = None) -> TaskFactory:
        """Return the process-local TaskFactory singleton."""
        with cls._lock:
            if cls._instance is None:
                cls._instance = cls(redis_config)
            elif redis_config:
                cls._instance.redis_config.update(redis_config)
            return cls._instance

    @classmethod
    def reset_instance(cls) -> None:
        """Reset the singleton (test helper)."""
        with cls._lock:
            cls._instance = None

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        """Build the control-process Redis client."""
        self.redis = RedisQueue(self.redis_config, client=client)
        self.redis.connect()
        return self.redis

    def start_workers(self, n: int, **worker_kwargs: Any) -> list[Worker]:
        """Spawn ``n`` Worker processes via the spawn context.

        ``worker_kwargs`` are merged into the picklable ``worker_config`` dict
        passed to each :class:`Worker`. Workers connect to Redis independently
        inside their child process and ``blpop`` from ``hep:task_queue``.
        """
        if n <= 0:
            return []
        if self.redis is None:
            raise RuntimeError("TaskFactory.init_redis() must be called before start_workers()")

        shared_config = dict(worker_kwargs)
        started: list[Worker] = []
        for worker_id in range(n):
            config = copy.deepcopy(shared_config)
            worker = Worker(worker_id, self.redis_config, config)
            worker.start()
            started.append(worker)
        self.workers.extend(started)
        self._logger.info("started %d worker process(es)", len(started))
        return started

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
                import os
                import signal

                try:
                    os.kill(worker.pid, signal.SIGTERM)
                except ProcessLookupError:
                    pass

    def get_monitor_snapshot(self) -> dict[str, Any]:
        """Return an in-memory deepcopy of the latest monitor snapshot."""
        with self._snapshot_lock:
            return copy.deepcopy(self._snapshot)

    def _worker_status(self) -> list[dict[str, Any]]:
        rows: list[dict[str, Any]] = []
        for worker in self.workers:
            rows.append(
                {
                    "worker_id": worker.worker_id,
                    "pid": worker.pid,
                    "alive": worker.is_alive(),
                    "name": worker.name,
                }
            )
        return rows

    def _collect_latest_status(self) -> dict[str, Any]:
        """Collect queue/stats snapshot from Redis (D5.1 adds op_count gating)."""
        snap: dict[str, Any] = {
            "timestamp": time.time(),
            "workers": self._worker_status(),
            "workers_alive": sum(1 for worker in self.workers if worker.is_alive()),
            "workers_total": len(self.workers),
        }
        if self.redis is None:
            return snap
        raw = self.redis.snapshot_raw()
        snap["task_queue_length"] = raw.get("task_queue_length", 0)
        snap["archive_queue_length"] = raw.get("archive_queue_length", 0)
        snap["sample_stats"] = raw.get("sample_stats", {})
        snap["calculator_status"] = raw.get("calculator_status", {})
        snap["op_counts"] = raw.get("op_counts", {})
        return snap

    def start_monitor(self, *, update_hz: float = 2.0) -> None:
        """Launch the background snapshot updater thread."""
        if self._updater_thread is not None and self._updater_thread.is_alive():
            return
        self._running = True
        interval = 1.0 / max(0.1, float(update_hz))

        def _loop() -> None:
            while self._running:
                try:
                    latest = self._collect_latest_status()
                    with self._snapshot_lock:
                        self._snapshot = latest
                except Exception as exc:
                    self._logger.warning("monitor snapshot update failed -> %s", exc)
                time.sleep(interval)

        self._updater_thread = threading.Thread(
            target=_loop,
            name="Jarvis2-FactoryMonitor",
            daemon=True,
        )
        self._updater_thread.start()

    def shutdown(self, *, wait: bool = True) -> None:
        """Stop monitor, signal Workers, join, and release Redis handle."""
        self._running = False
        if self._updater_thread is not None:
            self._updater_thread.join(timeout=2.0)
            self._updater_thread = None
        self.request_worker_shutdown()
        self.stop_all_workers(graceful=wait)
        self.redis = None
        self._logger.info("TaskFactory shutdown complete")


__all__ = ["TaskFactory"]