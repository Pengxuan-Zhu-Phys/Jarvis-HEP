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
from jarvishep2.database import SimpleHDF5Writer, StreamingHDF5Writer
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, delete_paths, normalize_delete_method
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_config import ARCHIVER_DEFAULTS
from jarvishep2.sample_bucket import pack_bucket_dir

Process = get_spawn_context().Process


# Keep the default console useful on long scans: DATABASE writes can reach
# tens of thousands of rows, so per-50-row INFO updates overwhelm the sampler
# milestones.  The forced final update is still emitted regardless of this
# interval.
DATARECORDER_PROGRESS_INTERVAL = 500


def _fmt_int(n: Any) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return str(n)


def _short_path(path: str, *, sample_root: str = "") -> str:
    """Prefer SAMPLE-relative or basename for compact logs."""
    text = str(path or "").strip()
    if not text:
        return ""
    root = str(sample_root or "").rstrip(os.sep)
    if root and text.startswith(root + os.sep):
        rel = text[len(root) + 1 :]
        return rel or os.path.basename(text)
    # Collapse .../SAMPLE/001369 → SAMPLE/001369 when possible.
    parts = text.replace("\\", "/").split("/")
    if "SAMPLE" in parts:
        idx = parts.index("SAMPLE")
        return "/".join(parts[idx:])
    return os.path.basename(text) or text


def _format_arrow_block(title: str, rows: list[tuple[str, Any]]) -> str:
    """Multi-line ``title ->`` block with tab-aligned key/value lines.

    No leading blank line; title has no ``[Archiver]`` prefix (module label
    already comes from ``·•· Jarvis-HEP.Archiver``).
    """
    text = str(title or "").strip()
    if text.startswith("[Archiver]"):
        text = text[len("[Archiver]") :].strip()
    label_w = max((len(str(k)) for k, _ in rows), default=0)
    lines = [f"{text} ->"]
    for key, value in rows:
        lines.append(f"\t{str(key):<{label_w}}\t-> {value}")
    return "\n".join(lines)


def _log_archiver_kv(
    logger: Any,
    title: str,
    rows: list[tuple[str, Any]],
    *,
    level: str = "info",
) -> None:
    """Emit Archiver status as ``title ->`` with tab-aligned fields."""
    if logger is None:
        return
    emit = getattr(logger, level, None) or getattr(logger, "info", None)
    if emit is None:
        return
    block = _format_arrow_block(title, rows)
    try:
        emit(block)
    except Exception:
        try:
            flat = " ".join(f"{k}={v}" for k, v in rows)
            emit("%s", f"{title} -> {flat}")
        except Exception:
            pass


class ArchiveProcessor:
    """Consume archive-queue payloads: move staging → SAMPLE and write DATABASE rows."""

    def __init__(
        self,
        writer: SimpleHDF5Writer | StreamingHDF5Writer,
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
        self._last_flushed: list[dict[str, Any]] = []
        self._last_flush = time.monotonic()
        self._lock = threading.Lock()
        os.makedirs(self.sample_root, exist_ok=True)

    @classmethod
    def from_config(
        cls,
        writer: SimpleHDF5Writer,
        *,
        sample_root: str,
        delete_method: str = DEFAULT_DELETE_METHOD,
        archiver_config: Mapping[str, Any] | None = None,
    ) -> ArchiveProcessor:
        """Build a processor from the Runtime.Archiver-style config map."""
        cfg = dict(ARCHIVER_DEFAULTS)
        if isinstance(archiver_config, Mapping):
            cfg.update(archiver_config)
        return cls(
            writer,
            sample_root=sample_root,
            delete_method=delete_method,
            strategy=str(cfg.get("strategy", "move")),
            delete_after_archive=bool(cfg.get("delete_after_archive", True)),
            batch_size=int(cfg.get("batch_size", 1)),
            flush_interval_sec=float(cfg.get("flush_interval_sec", 1.0)),
        )

    def has_pending(self) -> bool:
        """Return True when a partial batch is waiting to flush."""
        with self._lock:
            return bool(self._batch)

    def _flush_due_locked(self) -> bool:
        """Require both a full batch and the configured persistence interval."""
        return (
            len(self._batch) >= self.batch_size
            and time.monotonic() - self._last_flush >= self.flush_interval_sec
        )

    def ingest(self, result: Mapping[str, Any]) -> int:
        """Queue one archive payload; persist only when both thresholds are met."""
        with self._lock:
            self._batch.append(dict(result))
            if self._flush_due_locked():
                return self._flush_batch_locked()
            return 0

    def flush_batch(self, *, force: bool = True) -> int:
        """Persist queued payloads; shutdown/manual drains may force a tail batch."""
        with self._lock:
            return self._flush_batch_locked(force=force)

    def take_last_flushed(self) -> list[dict[str, Any]]:
        """Return and clear payloads successfully written by the last flush."""
        with self._lock:
            items = list(self._last_flushed)
            self._last_flushed.clear()
            return items

    def _flush_batch_locked(self, *, force: bool = False) -> int:
        self._last_flushed = []
        if not self._batch:
            self._last_flush = time.monotonic()
            return 0
        if not force and not self._flush_due_locked():
            return 0

        batched_writer = isinstance(self.writer, StreamingHDF5Writer)
        previous_records_written = self.records_written
        newly_acked: list[str] = []
        if batched_writer:
            self.writer.begin_batch()
        written = 0
        try:
            for result in self._batch:
                uuid = str(result.get("uuid", "")).strip()
                was_acked = bool(uuid and uuid in self.acked_uuids)
                if self._archive_one(result):
                    written += 1
                    self._last_flushed.append(dict(result))
                    if uuid and not was_acked:
                        newly_acked.append(uuid)
            if batched_writer:
                self.writer.commit_batch()
        except BaseException:
            if batched_writer:
                self.writer.abort_batch()
            for uuid in newly_acked:
                self.acked_uuids.discard(uuid)
            self.records_written = previous_records_written
            self._last_flushed = []
            raise
        self._batch.clear()
        self._last_flush = time.monotonic()
        return written

    def _archive_one(self, result: Mapping[str, Any]) -> bool:
        uuid = str(result.get("uuid", "")).strip()
        if not uuid:
            return False
        if uuid in self.acked_uuids:
            return True

        staging_path = str(result.get("staging_path") or "").strip()
        save_dir = str(result.get("save_dir") or "").strip()
        bucket_dir = str(result.get("bucket_dir") or "").strip()
        # Prefer the Worker-final path (SAMPLE/<bucket>/<uuid> under direct handoff).
        destination = ""
        if save_dir and os.path.isdir(save_dir):
            destination = save_dir
        elif bucket_dir:
            candidate = os.path.join(bucket_dir, uuid)
            if os.path.isdir(candidate):
                destination = candidate
        if not destination:
            destination = os.path.join(self.sample_root, uuid)

        if staging_path and os.path.isdir(staging_path) and not os.path.isdir(destination):
            # Legacy staging hop: move into SAMPLE/<uuid>/ (or keep under bucket if set).
            if bucket_dir:
                os.makedirs(bucket_dir, exist_ok=True)
                destination = archive_staging_to_sample(
                    staging_path,
                    bucket_dir,
                    uuid,
                    strategy=self.strategy,
                )
            else:
                destination = archive_staging_to_sample(
                    staging_path,
                    self.sample_root,
                    uuid,
                    strategy=self.strategy,
                )
            if self.delete_after_archive and os.path.lexists(staging_path):
                delete_paths([staging_path], method=self.delete_method, missing_ok=True)

        observables = result.get("observables", {})
        # Always persist observables when present — even if the on-disk sample
        # tree was already packed/pruned — so wait_for_results can complete.
        if isinstance(observables, Mapping) and observables:
            record = dict(observables)
            if destination and os.path.isdir(destination):
                record.setdefault("product_list", list_product_names(destination))
            if bucket_dir:
                record.setdefault("bucket_dir", bucket_dir)
            if result.get("bucket_id") is not None:
                record.setdefault("bucket_id", result.get("bucket_id"))
            self.writer.add_record(record)
            self.records_written += 1
            self.acked_uuids.add(uuid)
            return True

        # No observables: still ack if the sample dir exists (artifact-only sample).
        if destination and os.path.isdir(destination):
            self.acked_uuids.add(uuid)
            self.records_written += 1
            return True
        # Nothing to write and no dir — count as archived failure for the queue
        # item (already popped). Prefer writing a minimal stub so drains finish.
        self.writer.add_record({"uuid": uuid, "status": str(result.get("status") or "Completed")})
        self.records_written += 1
        self.acked_uuids.add(uuid)
        return True

    def persistence_state(self) -> dict[str, Any]:
        acked = sorted(self.acked_uuids)
        return {
            "acked_uuids": acked,
            "acked_uuids_highwater": len(acked),
            "next_sample_index": len(acked),
        }



class SimpleArchiver:
    """Background thread that drains ``hep:archive_queue`` + ready SAMPLE buckets."""

    def __init__(
        self,
        redis_queue: RedisQueue,
        db_path: str,
        *,
        sample_root: str,
        poll_timeout: float = 1.0,
        delete_method: str = DEFAULT_DELETE_METHOD,
        archiver_config: Mapping[str, Any] | None = None,
        logger: Any | None = None,
    ) -> None:
        self.redis = redis_queue
        self.sample_root = os.path.abspath(str(sample_root))
        self.poll_timeout = max(0.05, float(poll_timeout))
        cfg = dict(ARCHIVER_DEFAULTS)
        if isinstance(archiver_config, Mapping):
            cfg.update(archiver_config)
        self.pack_buckets = bool(cfg.get("pack_buckets", True))
        self.buckets_packed = 0
        self.processor = ArchiveProcessor.from_config(
            StreamingHDF5Writer(db_path),
            sample_root=sample_root,
            delete_method=delete_method,
            archiver_config=archiver_config,
        )
        self._stop_event = threading.Event()
        self._thread: threading.Thread | None = None
        parent = os.path.dirname(self.sample_root)
        self._manifest_jsonl = os.path.join(parent, "DATABASE", "archive_manifest.jsonl")
        # Own logger identity: ``·•· Jarvis-HEP.Archiver`` (pack / drain / lifecycle).
        self._logger = logger or get_jarvis_logger("archiver")
        # DATABASE / samples.hdf5 row counter → separate DataRecorder sink.
        self._datarecorder_logger = get_jarvis_logger("datarecorder")
        self._last_progress_written = -1
        self._progress_interval = DATARECORDER_PROGRESS_INTERVAL

    @property
    def records_written(self) -> int:
        return self.processor.records_written

    def persistence_state(self) -> dict[str, Any]:
        state = self.processor.persistence_state()
        state["buckets_packed"] = int(self.buckets_packed)
        return state

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
        _log_archiver_kv(
            self._logger,
            "loop started",
            [
                ("sample_root", _short_path(self.sample_root) or self.sample_root),
                ("pack_buckets", self.pack_buckets),
            ],
            level="debug",
        )

    def _note_last_flushed(self) -> None:
        """Tell Redis each freshly written sample is archived (may enqueue pack)."""
        for item in self.processor.take_last_flushed():
            bucket_id = item.get("bucket_id")
            if bucket_id is None:
                continue
            try:
                self.redis.note_sample_archived(bucket_id)
            except Exception as exc:
                try:
                    self._logger.warning(
                        "note_sample_archived failed for bucket %s -> %s",
                        bucket_id,
                        exc,
                    )
                except Exception:
                    pass

    def _maybe_log_progress(self, *, force: bool = False) -> None:
        written = int(self.records_written)
        if not force and written == self._last_progress_written:
            return
        if (
            not force
            and written > 0
            and written - self._last_progress_written < self._progress_interval
        ):
            return
        self._last_progress_written = written
        # HDF5 / DATABASE flush progress — not pack lifecycle (that stays on Archiver).
        _log_archiver_kv(
            self._datarecorder_logger,
            "DATABASE progress",
            [
                ("rows written", _fmt_int(written)),
                ("buckets packed", _fmt_int(self.buckets_packed)),
            ],
        )

    def _ingest_result(self, result: Mapping[str, Any]) -> int:
        written = self.processor.ingest(result)
        if written:
            self._note_last_flushed()
            self._maybe_log_progress()
        return written

    def _flush_and_note(self, *, force: bool = False) -> int:
        written = self.processor.flush_batch(force=force)
        if written:
            self._note_last_flushed()
            self._maybe_log_progress()
        return written

    def _pack_ready_buckets(self) -> int:
        """Tar sealed buckets that Redis marked ready (archived == assigned)."""
        if not self.pack_buckets:
            return 0
        packed = 0
        while True:
            ready = self.redis.pull_ready_bucket(timeout=0)
            if ready is None:
                break
            if not ready.get("pack", True):
                continue
            bucket_dir = str(ready.get("bucket_dir") or "")
            bucket_id = ready.get("bucket_id")
            sample_root = str(ready.get("sample_root") or self.sample_root)
            if not bucket_dir:
                continue
            try:
                pack_bucket_dir(
                    bucket_dir,
                    sample_root=sample_root,
                    manifest_jsonl_path=self._manifest_jsonl,
                    prune=True,
                )
                if bucket_id is not None:
                    self.redis.mark_bucket_packed(bucket_id)
                packed += 1
                self.buckets_packed += 1
                _log_archiver_kv(
                    self._logger,
                    "SAMPLE bucket packed",
                    [
                        ("bucket id", bucket_id),
                        (
                            "path",
                            _short_path(bucket_dir, sample_root=sample_root),
                        ),
                        ("total packed", _fmt_int(self.buckets_packed)),
                    ],
                    level="debug",
                )
            except Exception as exc:
                # Leave packing flag set; operator can inspect / re-seal later.
                _log_archiver_kv(
                    self._logger,
                    "SAMPLE bucket pack failed",
                    [
                        ("bucket id", bucket_id),
                        (
                            "path",
                            _short_path(bucket_dir, sample_root=sample_root),
                        ),
                        ("error", exc),
                    ],
                    level="warning",
                )
                continue
        return packed

    def _run_loop(self) -> None:
        # Thread shares the control process title; only Process-mode Archiver
        # renames the OS process. Keep this loop name for debugger visibility.
        timeout = max(1, int(round(self.poll_timeout)))
        while not self._stop_event.is_set():
            result = self.redis.pull_result(timeout=timeout)
            if result is None:
                if self.processor.has_pending():
                    self._flush_and_note()
                self._pack_ready_buckets()
                continue
            self._ingest_result(result)
            # Pack only after any flush from this ingest has noted archived counts.
            self._pack_ready_buckets()

    def drain(self, *, idle_timeout: float = 2.0) -> int:
        """Drain remaining archive-queue items synchronously."""
        timeout = max(1, int(round(self.poll_timeout)))
        idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
        drained = 0
        while time.monotonic() < idle_deadline:
            result = self.redis.pull_result(timeout=timeout)
            if result is None:
                flushed = self._flush_and_note(force=True)
                drained += flushed
                drained += self._pack_ready_buckets()
                continue
            idle_deadline = time.monotonic() + max(0.1, float(idle_timeout))
            drained += self._ingest_result(result)
            drained += self._pack_ready_buckets()
        drained += self._flush_and_note(force=True)
        drained += self._pack_ready_buckets()
        self._maybe_log_progress(force=True)
        _log_archiver_kv(
            self._logger,
            "drain complete",
            [
                ("flushed this drain", _fmt_int(drained)),
                ("rows written (total)", _fmt_int(self.records_written)),
                ("buckets packed (total)", _fmt_int(self.buckets_packed)),
            ],
            level="debug",
        )
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
        if thread is None or not thread.is_alive():
            self.processor.writer.close()
        _log_archiver_kv(
            self._logger,
            "stopped",
            [
                ("rows written", _fmt_int(self.records_written)),
                ("buckets packed", _fmt_int(self.buckets_packed)),
            ],
        )


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
        name: str = "Jarvis2-Archiver",
        scan_name: str | None = None,
        log_dir: str | None = None,
        log_path: str | None = None,
    ) -> None:
        super().__init__(name=name, daemon=False)
        self.redis_config = dict(redis_config)
        self.db_path = str(db_path)
        self.sample_root = str(sample_root)
        self.delete_method = str(delete_method)
        self.archiver_config = dict(archiver_config or {})
        self.poll_timeout = max(0.05, float(poll_timeout))
        self.scan_name = str(scan_name or "").strip()
        self.log_dir = str(log_dir or "").strip() or None
        self.log_path = str(log_path or "").strip() or None
        self._stop_event = get_spawn_context().Event()
        self.records_written = get_spawn_context().Value("i", 0)

    def run(self) -> None:
        from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
        from jarvishep2.proc_title import archiver_title, set_process_title

        set_process_title(archiver_title(scan_name=self.scan_name or None))
        # Child process: own logging sinks (do not share control QueueListener).
        # Optional logging policy stamped from control (same CLI as Workers).
        silence = bool(self.archiver_config.get("log_silence", False))
        console_level = str(self.archiver_config.get("console_level") or "INFO")
        log_level = str(self.archiver_config.get("log_level") or "INFO")
        setup_kwargs: dict[str, Any] = {
            "role": "archiver",
            "component": "archiver",
            "level": log_level,
            "console": not silence,
            "console_level": console_level,
            "silence": silence,
            "use_queue": True,
        }
        # Prefer scan multi-sink (archiver.log + datarecorder.log). Only force a
        # single log_path when no scan logs dir is available.
        if self.log_dir:
            setup_kwargs["scan_logs_dir"] = self.log_dir
            setup_kwargs["multi_sink"] = True
        elif self.log_path:
            setup_kwargs["log_path"] = self.log_path
        setup_jarvis_logging(**setup_kwargs)
        logger = get_jarvis_logger("archiver")
        _log_archiver_kv(
            logger,
            "process started",
            [
                ("pid", os.getpid()),
                ("scan", self.scan_name or "?"),
                ("sample_root", _short_path(self.sample_root) or self.sample_root),
                ("database", _short_path(self.db_path) or self.db_path),
            ],
        )

        redis = RedisQueue(self.redis_config)
        redis.connect()
        # Reuse SimpleArchiver loop so process mode also packs sealed buckets.
        archiver = SimpleArchiver(
            redis,
            self.db_path,
            sample_root=self.sample_root,
            poll_timeout=self.poll_timeout,
            delete_method=self.delete_method,
            archiver_config=self.archiver_config,
            logger=logger,
        )
        archiver.start()
        while not self._stop_event.is_set():
            time.sleep(0.1)
            with self.records_written.get_lock():
                self.records_written.value = int(archiver.records_written)
        archiver.stop(wait=True, drain=True)
        with self.records_written.get_lock():
            self.records_written.value = int(archiver.records_written)
        try:
            logger.debug(
                "Archiver process exiting pid=%s records_written=%s",
                os.getpid(),
                int(self.records_written.value),
            )
        except Exception:
            pass
        redis.close()

    def drain(self, *, idle_timeout: float = 2.0) -> int:
        """Best-effort: process mode drains on stop; expose no-op for core finalize."""
        # Active draining happens inside the child SimpleArchiver before exit.
        _ = idle_timeout
        try:
            return int(self.records_written.value)  # type: ignore[attr-defined]
        except Exception:
            return 0

    def persistence_state(self) -> dict[str, Any]:
        """Control-side checkpoint stub (child owns the real ack set).

        Process-mode Archiver does not share acked_uuids with the parent; resume
        re-proposes unfinished sampler work. Return an empty persistence bag so
        ``save_runtime_checkpoint`` / tests share the SimpleArchiver API.
        """
        return {"acked_uuids": [], "records_written": int(self.drain())}

    def stop(self, *, wait: bool = True, timeout: float = 5.0, drain: bool = True) -> None:
        """Stop the child Archiver.

        ``drain`` is accepted for API parity with :class:`SimpleArchiver`. The
        child always drains its in-process SimpleArchiver on exit
        (``archiver.stop(wait=True, drain=True)`` in :meth:`run`); the flag is
        therefore informational for callers and does not change parent behavior.
        """
        _ = drain
        self._stop_event.set()
        if wait:
            self.join(timeout=max(0.1, float(timeout)))


__all__ = ["ArchiveProcessor", "ArchiverProcess", "SimpleArchiver"]
