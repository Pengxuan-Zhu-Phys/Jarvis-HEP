#!/usr/bin/env python3
"""Jarvis2 control-process orchestrator for the distributed runtime."""

from __future__ import annotations

import os
import time
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.archiver import SimpleArchiver
from jarvishep2.factory import TaskFactory
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.redis_queue import RedisQueue, make_fakeredis_queue
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample


class Jarvis2Core:
    """Minimal distributed run orchestrator for the D1.1 opera-only MVP."""

    def __init__(self, config: Mapping[str, Any] | None = None) -> None:
        self.config = dict(config or {})
        self.runtime = get_runtime_block(self.config)
        self.info: dict[str, Any] = {}
        self.redis: RedisQueue | None = None
        self.factory: TaskFactory | None = None
        self.archiver: SimpleArchiver | None = None
        self.sampler: Any = None
        self._logger = get_jarvis_logger("core")

    def init_logger(self) -> None:
        setup_jarvis_logging(role="core")

    def is_redis_runtime(self) -> bool:
        """Return True when the distributed Redis path should be used."""
        return str(self.runtime.get("mode", "auto")).strip().lower() == "redis"

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        redis_config = dict(self.runtime.get("redis") or {})
        if client is not None:
            self.redis = RedisQueue(redis_config, client=client)
        elif not redis_config:
            self.redis = make_fakeredis_queue()
        else:
            self.redis = RedisQueue(redis_config)
        self.redis.connect()
        return self.redis

    def init_factory(self, worker_config: Mapping[str, Any] | None = None) -> TaskFactory | None:
        if not self.is_redis_runtime():
            self._logger.info("Runtime.mode != redis; skipping TaskFactory bring-up")
            return None

        redis_config = dict(self.runtime.get("redis") or {})
        self.factory = TaskFactory.get_instance(redis_config)
        if self.redis is not None:
            self.factory.redis = self.redis
        else:
            self.factory.init_redis()

        workers = int(self.runtime.get("workers", 1) or 1)
        if workers <= 0:
            workers = 1

        merged_config = dict(worker_config or {})
        self.factory.start_workers(workers, **merged_config)
        self.factory.start_monitor(update_hz=2.0)
        self._logger.info("TaskFactory started with %d worker(s)", workers)
        return self.factory

    def init_archiver(self, db_path: str | None = None) -> SimpleArchiver:
        if self.redis is None:
            raise RuntimeError("init_redis() must run before init_archiver()")
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        os.makedirs(database_dir, exist_ok=True)
        resolved_db_path = db_path or os.path.join(database_dir, "samples.hdf5")
        self.archiver = SimpleArchiver(self.redis, resolved_db_path)
        self.archiver.start()
        return self.archiver

    def set_sampler(self, sampler: Any) -> None:
        self.sampler = sampler
        if self.redis is not None and hasattr(sampler, "set_redis"):
            sampler.set_redis(self.redis)

    def submit_samples(self, samples: Sequence[Sample]) -> None:
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        if hasattr(self.sampler, "_submit_group"):
            self.sampler._submit_group(list(samples))
        else:
            for sample in samples:
                self.sampler._submit(sample)

    def wait_for_results(
        self,
        expected: int,
        *,
        timeout: float = 30.0,
        poll_interval: float = 0.1,
    ) -> None:
        if self.archiver is None:
            raise RuntimeError("archiver is not configured")
        deadline = time.monotonic() + max(0.1, float(timeout))
        while time.monotonic() < deadline:
            self.archiver.drain(idle_timeout=poll_interval)
            if self.archiver.records_written >= expected:
                return
            time.sleep(poll_interval)
        raise TimeoutError(
            f"timed out waiting for {expected} archived results; got {self.archiver.records_written}"
        )

    def shutdown(self, *, wait: bool = True) -> None:
        if self.archiver is not None:
            self.archiver.stop(wait=wait, drain=True)
            self.archiver = None
        if self.factory is not None:
            self.factory.shutdown(wait=wait)
            self.factory = None
        self.redis = None


__all__ = ["Jarvis2Core"]