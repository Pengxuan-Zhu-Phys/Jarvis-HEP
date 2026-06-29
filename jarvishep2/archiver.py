#!/usr/bin/env python3
"""Archiver for Layer-2 persistence: staging handoff + DATABASE writes (WP-D4)."""

from __future__ import annotations

import os
import threading
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.archive_handoff import (
    archive_staging_to_sample,
    list_product_names,
    normalize_move_strategy,
)
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, delete_paths, normalize_delete_method
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_config import ARCHIVER_DEFAULTS

Process = get_spawn_context().Process


class ArchiveProcessor:
    """Consume archive-queue payloads: move staging → SAMPLE and write DATABASE rows."""

    def __init__(
        self,
        writer: SimpleHDF5Writer,
        *,
        sample_root: str,
        delete_method: str = DEFAULT_DELETE_METHOD,
        strategy: str = "move",
        delete_after_archive: bool = True,
        batch_size: int = 1,
        flush_interval_sec: float = 1.0,
    ) -> None:
        self.writer = writer
        self.sample_root = os.path.abspath(str(sample_root))
        self.delete_method = normalize_delete_method(delete_method)
        self.strategy = normalize_move_strategy(strategy)
        self.delete_after_archive = bool(delete_after_archive)
        self.batch_size = max(1, int(batch_size))
        self.flush_interval_sec = max(0.05, float(flush_interval_sec))
        self.records_written = 0
        self.acked_uuids: set[str] = set()
        self._batch: list[dict[str, Any]] = []
        self._last_flush = time.monotonic()
        self._lock = threading.Lock()
        os.makedirs(self.sample_root, exist_ok=True)

    def ingest(self, result: Mapping[str, Any]) -> int:
        """Queue one archive payload; flush when batch/interval thresholds are met."""
        with self._lock:
            self._batch.append(dict(result))
            if len(self._batch) >= self.batch_size:
                return self._flush_batch_locked()
            if time.monotonic() - self._last_flush >= self.flush_interval_sec:
                return self._flush_batch_locked()
            return 0

    def flush_batch(self) -> int:
        """Persist all queued payloads."""
        with self._lock:
            return self._flush_batch_locked()

    def _flush_batch_locked(self) -> int:
        if not self._batch:
            self._last_flush = time.monotonic()
            return 0
        written = 0
        for result in self._batch:
            if self._archive_one(result):
                written += 1
        self._batch.clear()
        self._last_flush = time.monotonic()
        return written

    def _archive_one(self, result: Mapping[str, Any]) -> bool:
        uuid = str(result.get("uuid", "")).strip()
        if not uuid:
            return False

        staging_path = str(result.get("staging_path") or "").strip()
        save_dir = str(result.get("save_dir") or "").strip()
        destination = os.path.join(self.sample_root, uuid)
        if os.path.isdir(destination):
            if uuid in self.acked_uuids:
                return True
            observables = result.get("observables", {})
            if isinstance(observables, Mapping) and observables:
                record = dict(observables)
                record.setdefault("product_list", list_product_names(destination))
                self.writer.add_record(record)
                self.records_written += 1
                self.acked_uuids.add(uuid)
                return True
            self.acked_uuids.add(uuid)
            return True
        if staging_path and os.path.isdir(staging_path):
            destination = archive_staging_to_sample(
                staging_path,
                self.sample_root,
                uuid,
                strategy=self.strategy,
            )
            if self.delete_after_archive and os.path.lexists(staging_path):
                delete_paths([staging_path], method=self.delete_method, missing_ok=True)
        elif save_dir and os.path.isdir(save_dir):
            destination = archive_staging_to_sample(
                save_dir,
                self.sample_root,
                uuid,
                strategy=self.strategy,
            )
        elif not os.path.isdir(destination):
            return False

        observables = result.get("observables", {})
        if isinstance(observables, Mapping) and observables:
            record = dict(observables)
            if os.path.isdir(destination):
                record.setdefault("product_list", list_product_names(destination))
            self.writer.add_record(record)
            self.records_written += 1
            self.acked_uuids.add(uuid)
            return True
        return False

    def persistence_state(self) -> dict[str, Any]:
        acked = sorted(self.acked_uuids)
        return {
            "acked_uuids": acked,
            "acked_uuids_highwater": len(acked),
            "next_sample_index": len(acked),
        }


class SimpleArchiver:
    """Background thread that drains ``hep:archive_queue``."""

    def __init__(
        self,
        redis_queue: RedisQueue,
        db_path: str,
        *,
        sample_root: str,
        poll_timeout: float = 1.0,
        delete_method: str = DEFAULT_DELETE_METHOD,
        archiver_config: Mapping[str, Any] | None = None,
    ) -> None:
        cfg = dict(ARCHIVER_DEFAULTS)
        if isinstance(archiver_config, Mapping):
            cfg.update(archiver_config)
        self.redis = redis_queue
        self.poll_timeout = max(0.05, float(poll_timeout))
        self.processor = ArchiveProcessor(
            SimpleHDF5Writer(db_path),
            sample_root=sample_root,
            delete_method=delete_method,
            strategy=str(cfg.get("strategy", "move")),
            delete_after_archive=bool(cfg.get("delete_after_archive", True)),
            batch_size=int(cfg.get("batch_size", 1)),
            flush_interval_sec=float(cfg.get("flush_interval_sec", 1.0)),
        )
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None

    @property
    def records_written(self) -> int:
        return self.processor.records_written

    def persistence_state(self) -> dict[str, Any]:
        return self.processor.persistence_state()

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
                if self.processor._batch:
                    self.processor.flush_batch()
                continue
            self.processor.ingest(result)

    def drain(self, *, idle_timeout: float = 2.0) -> int:
        """Drain remaining archive-queue items synchronously."""
        timeout = max(1, int(round(self.poll_timeout)))
        idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
        drained = 0
        while time.monotonic() < idle_deadline:
            result = self.redis.pull_result(timeout=timeout)
            if result is None:
                flushed = self.processor.flush_batch()
                drained += flushed
                continue
            idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
            drained += self.processor.ingest(result)
        drained += self.processor.flush_batch()
        return drained

    def cleanup_staging(self, paths: list[str] | tuple[str, ...]) -> None:
        """Delete archived staging directories using the configured backend."""
        delete_paths(list(paths), method=self.processor.delete_method, missing_ok=True)

    def stop(self, *, wait: bool = True, drain: bool = True) -> None:
        if drain:
            self.drain()
        self._stop_event.set()
        thread = self._thread
        if thread is not None and wait:
            thread.join(timeout=5.0)
        self._thread = None


class ArchiverProcess(Process):
    """Spawned Archiver process for Layer-2 persistence (WP-D4.1)."""

    def __init__(
        self,
        redis_config: Mapping[str, Any],
        *,
        db_path: str,
        sample_root: str,
        delete_method: str = DEFAULT_DELETE_METHOD,
        archiver_config: Mapping[str, Any] | None = None,
        poll_timeout: float = 1.0,
        name: str = "HEP2-Archiver",
    ) -> None:
        super().__init__(name=name, daemon=False)
        self.redis_config = dict(redis_config)
        self.db_path = str(db_path)
        self.sample_root = str(sample_root)
        self.delete_method = str(delete_method)
        self.archiver_config = dict(archiver_config or {})
        self.poll_timeout = max(0.05, float(poll_timeout))
        self._stop_event = get_spawn_context().Event()
        self.records_written = get_spawn_context().Value("i", 0)

    def run(self) -> None:
        redis = RedisQueue(self.redis_config)
        redis.connect()
        processor = ArchiveProcessor(
            SimpleHDF5Writer(self.db_path),
            sample_root=self.sample_root,
            delete_method=self.delete_method,
            strategy=str(self.archiver_config.get("strategy", "move")),
            delete_after_archive=bool(self.archiver_config.get("delete_after_archive", True)),
            batch_size=int(self.archiver_config.get("batch_size", 1)),
            flush_interval_sec=float(self.archiver_config.get("flush_interval_sec", 1.0)),
        )
        timeout = max(1, int(round(self.poll_timeout)))
        while not self._stop_event.is_set():
            result = redis.pull_result(timeout=timeout)
            if result is None:
                flushed = processor.flush_batch()
                with self.records_written.get_lock():
                    self.records_written.value += flushed
                continue
            written = processor.ingest(result)
            with self.records_written.get_lock():
                self.records_written.value += written
        processor.flush_batch()
        with self.records_written.get_lock():
            self.records_written.value = processor.records_written
        redis.close()

    def stop(self, *, wait: bool = True, timeout: float = 5.0) -> None:
        self._stop_event.set()
        if wait:
            self.join(timeout=max(0.1, float(timeout)))


__all__ = ["ArchiveProcessor", "ArchiverProcess", "SimpleArchiver"]