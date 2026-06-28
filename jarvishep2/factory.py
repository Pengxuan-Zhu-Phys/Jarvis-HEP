#!/usr/bin/env python3
"""TaskFactory — Worker lifecycle manager for Jarvis-HEP V2."""

from __future__ import annotations

import copy
import threading
import time
from collections.abc import Callable
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
        self._last_op_counts: dict[str, int] = {}
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

    def _alive_workers(self) -> list[Worker]:
        return [worker for worker in self.workers if worker.is_alive()]

    def start_workers(self, n: int, **worker_kwargs: Any) -> list[Worker]:
        """Spawn ``n`` Worker processes via the spawn context.

        The live control-process :attr:`redis` queue is passed into each
        :class:`Worker`; only its picklable connection settings cross the spawn
        boundary. Child Workers open their own Redis clients in ``run()``.

        Raises:
            RuntimeError: If Redis is not initialized or workers are already running.
        """
        if n <= 0:
            return []
        if self.redis is None:
            raise RuntimeError("TaskFactory.init_redis() must be called before start_workers()")

        alive = self._alive_workers()
        if alive:
            raise RuntimeError(
                f"start_workers refused: {len(alive)} worker process(es) already running"
            )
        self.workers = [worker for worker in self.workers if worker.is_alive()]

        shared_config = dict(worker_kwargs)
        started: list[Worker] = []
        for worker_id in range(n):
            config = copy.deepcopy(shared_config)
            worker = Worker(worker_id, self.redis, config)
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

    def _fetch_workers_redis(self) -> dict[str, dict[str, Any]]:
        if self.redis is None:
            return {}
        worker_ids = [str(worker.worker_id) for worker in self.workers]
        return self.redis.fetch_worker_status(worker_ids)

    def _carry_forward_section(
        self,
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

    def _collect_latest_status(self) -> dict[str, Any]:
        """Collect queue/stats snapshot from Redis with op_count gating."""
        snap: dict[str, Any] = {
            "timestamp": time.time(),
            "workers": self._worker_status(),
            "workers_alive": sum(1 for worker in self.workers if worker.is_alive()),
            "workers_total": len(self.workers),
        }
        if self.redis is None:
            return snap

        with self._snapshot_lock:
            prev = self._snapshot

        op_counts = self.redis.get_all_op_counts()
        snap["op_counts"] = op_counts

        lengths = self.redis.get_queue_lengths()
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
            supplier=self.redis.fetch_calculator_status,
        )
        snap["sample_stats"] = self._carry_forward_section(
            prev,
            "sample_stats",
            refresh=self._subsystem_refresh_needed("sample", op_counts, prev, "sample_stats"),
            supplier=self.redis.fetch_sample_stats,
        )

        self._last_op_counts = dict(op_counts)
        return snap

    def _start_snapshot_updater(self, *, update_hz: float = 120.0) -> None:
        """Run the background snapshot updater at ~100–120 Hz."""
        if self._updater_thread is not None and self._updater_thread.is_alive():
            return
        self._running = True
        interval = 1.0 / max(1.0, float(update_hz))

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

    def start_monitor(self, *, update_hz: float = 120.0) -> None:
        """Launch the background snapshot updater thread."""
        self._start_snapshot_updater(update_hz=update_hz)

    def shutdown(self, *, wait: bool = True) -> None:
        """Stop monitor, signal Workers, join, and close the Redis connection."""
        self._running = False
        if self._updater_thread is not None:
            self._updater_thread.join(timeout=2.0)
            self._updater_thread = None
        self.request_worker_shutdown()
        self.stop_all_workers(graceful=wait)
        if self.redis is not None:
            self.redis.close()
            self.redis = None
        self._last_op_counts.clear()
        self._logger.info("TaskFactory shutdown complete")


__all__ = ["TaskFactory"]