#!/usr/bin/env python3
"""Simplified Archiver for D1.1 — drain hep:archive_queue and write DATABASE."""

from __future__ import annotations

import threading
import time
from typing import Any

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, delete_paths, normalize_delete_method
from jarvishep2.redis_queue import RedisQueue


class SimpleArchiver:
    """Background thread that persists Worker results to HDF5."""

    def __init__(
        self,
        redis_queue: RedisQueue,
        db_path: str,
        *,
        poll_timeout: float = 1.0,
        delete_method: str = DEFAULT_DELETE_METHOD,
    ) -> None:
        self.redis = redis_queue
        self.writer = SimpleHDF5Writer(db_path)
        self.poll_timeout = max(0.05, float(poll_timeout))
        self.delete_method = normalize_delete_method(delete_method)
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        self.records_written = 0

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop_event.clear()
        self._thread = threading.Thread(
            target=self._run_loop,
            name="Jarvis2-SimpleArchiver",
            daemon=True,
        )
        self._thread.start()

    def _run_loop(self) -> None:
        timeout = max(1, int(round(self.poll_timeout)))
        while not self._stop_event.is_set():
            result = self.redis.pull_result(timeout=timeout)
            if result is None:
                continue
            observables = result.get("observables", {})
            if isinstance(observables, dict) and observables:
                self.writer.add_record(observables)
                self.records_written += 1

    def drain(self, *, idle_timeout: float = 2.0) -> int:
        """Drain remaining archive-queue items synchronously."""
        timeout = max(1, int(round(self.poll_timeout)))
        idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
        drained = 0
        while time.monotonic() < idle_deadline:
            result = self.redis.pull_result(timeout=timeout)
            if result is None:
                continue
            idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
            observables = result.get("observables", {})
            if isinstance(observables, dict) and observables:
                self.writer.add_record(observables)
                self.records_written += 1
                drained += 1
        return drained

    def cleanup_staging(self, paths: list[str] | tuple[str, ...]) -> None:
        """Delete archived staging directories using the configured backend."""
        delete_paths(list(paths), method=self.delete_method, missing_ok=True)

    def stop(self, *, wait: bool = True, drain: bool = True) -> None:
        if drain:
            self.drain()
        self._stop_event.set()
        thread = self._thread
        if thread is not None and wait:
            thread.join(timeout=5.0)
        self._thread = None


__all__ = ["SimpleArchiver"]