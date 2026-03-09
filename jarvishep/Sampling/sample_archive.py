#!/usr/bin/env python3
from __future__ import annotations

import json
import multiprocessing as mp
import os
import queue
import time
import shutil
import tarfile

from jarvishep.log_kv import format_two_column_log


_STOP = "__JARVIS_SAMPLE_ARCHIVE_STOP__"

def _count_jsonl_records(path: str) -> int:
    if not os.path.exists(path):
        return 0
    count = 0
    try:
        with open(path, "r", encoding="utf-8") as f1:
            for line in f1:
                if line.strip():
                    count += 1
    except Exception:
        return 0
    return count


def _resolve_archive_tar_path(sample_root: str, bucket_name: str) -> str:
    target = os.path.join(sample_root, f"{bucket_name}.tar.gz")
    if not os.path.exists(target):
        return target
    for idx in range(1, 100000):
        candidate = os.path.join(sample_root, f"{bucket_name}__{idx:04d}.tar.gz")
        if not os.path.exists(candidate):
            return candidate
    raise RuntimeError(f"Cannot allocate archive tar path for bucket={bucket_name}")


def _archive_worker(
    task_q,
    status_q,
    sample_root: str,
    manifest_jsonl_path: str,
):
    os.makedirs(sample_root, exist_ok=True)
    os.makedirs(os.path.dirname(manifest_jsonl_path), exist_ok=True)

    summary = {
        "enqueued": 0,
        "packed": 0,
        "skipped": 0,
        "failed": 0,
    }
    records = _count_jsonl_records(manifest_jsonl_path)

    while True:
        item = task_q.get()
        if item == _STOP:
            break

        summary["enqueued"] += 1
        bucket_dir = str(item)
        if not os.path.isdir(bucket_dir):
            summary["skipped"] += 1
            continue

        bucket_name = os.path.basename(bucket_dir.rstrip(os.sep))
        primary_tar = os.path.join(sample_root, f"{bucket_name}.tar.gz")
        if os.path.exists(primary_tar):
            shutil.rmtree(bucket_dir, ignore_errors=True)
            summary["skipped"] += 1
            continue
        out_path = _resolve_archive_tar_path(sample_root, bucket_name)
        tmp_path = out_path + ".tmp"
        mode = "tar_gz_prune"

        try:
            with tarfile.open(tmp_path, "w:gz") as tar:
                tar.add(bucket_dir, arcname=bucket_name)
            os.replace(tmp_path, out_path)
            shutil.rmtree(bucket_dir, ignore_errors=True)

            payload = {
                "bucket": bucket_name,
                "source_path": bucket_dir,
                "archive_path": out_path,
                "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "mode": mode,
            }
            with open(manifest_jsonl_path, "a", encoding="utf-8") as f1:
                f1.write(json.dumps(payload, ensure_ascii=True) + "\n")
            records += 1

            summary["packed"] += 1
        except Exception:
            summary["failed"] += 1
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except Exception:
                pass

    summary["manifest_path"] = manifest_jsonl_path
    summary["records"] = int(records)
    status_q.put(summary)


class SampleArchiveManager:
    def __init__(
        self,
        sample_root: str,
        logger=None,
        enabled: bool = True,
        queue_size: int = 256,
    ):
        self.sample_root = os.path.abspath(sample_root)
        self.database_root = os.path.join(os.path.dirname(self.sample_root), "DATABASE")
        self.manifest_jsonl_path = os.path.join(self.database_root, "archive_manifest.jsonl")
        self.logger = logger
        self.enabled = bool(enabled)
        self.queue_size = int(queue_size)

        self._ctx = mp.get_context("spawn")
        self._task_q = None
        self._status_q = None
        self._proc = None
        self._scheduled = set()

    def start(self):
        if not self.enabled:
            return
        if self._proc is not None and self._proc.is_alive():
            return

        os.makedirs(self.sample_root, exist_ok=True)
        os.makedirs(self.database_root, exist_ok=True)
        self._task_q = self._ctx.Queue(maxsize=self.queue_size)
        self._status_q = self._ctx.Queue(maxsize=8)
        self._proc = self._ctx.Process(
            target=_archive_worker,
            args=(
                self._task_q,
                self._status_q,
                self.sample_root,
                self.manifest_jsonl_path,
            ),
            daemon=True,
        )
        self._proc.start()
        self._log(
            "warning",
            format_two_column_log(
                "Sample archive worker started",
                [("pid", self._proc.pid)],
            ),
        )

    def enqueue_bucket_dir(
        self,
        bucket_dir: str,
        blocking: bool = False,
        timeout: float = 5.0,
        force: bool = False,
    ) -> bool:
        if not self.enabled:
            return False
        self.start()
        if self._task_q is None:
            return False

        bucket_dir = os.path.abspath(bucket_dir)
        if not bucket_dir.startswith(self.sample_root):
            return False
        if not os.path.isdir(bucket_dir):
            return False
        if (not force) and bucket_dir in self._scheduled:
            return False

        try:
            if blocking:
                self._task_q.put(bucket_dir, block=True, timeout=float(timeout))
            else:
                self._task_q.put_nowait(bucket_dir)
            self._scheduled.add(bucket_dir)
            self._log(
                "info",
                format_two_column_log(
                    "Sample archive enqueue",
                    [("bucket_dir", bucket_dir)],
                ),
            )
            return True
        except queue.Full:
            try:
                self._task_q.put(bucket_dir, block=True, timeout=float(timeout))
                self._scheduled.add(bucket_dir)
                self._log(
                    "warning",
                    format_two_column_log(
                        "Sample archive queue full; blocking enqueue succeeded",
                        [("bucket_dir", bucket_dir)],
                    ),
                )
                return True
            except Exception:
                self._log(
                    "warning",
                    format_two_column_log(
                        "Sample archive queue full; skip enqueue",
                        [("bucket_dir", bucket_dir)],
                    ),
                )
                return False
        except Exception as exc:
            self._log(
                "error",
                format_two_column_log(
                    "Sample archive enqueue failed",
                    [("bucket_dir", bucket_dir), ("error", exc)],
                ),
            )
            return False

    def _list_existing_bucket_dirs(self):
        try:
            entries = os.listdir(self.sample_root)
        except Exception:
            return []
        bucket_dirs = []
        for name in sorted(entries):
            if not str(name).isdigit():
                continue
            path = os.path.join(self.sample_root, name)
            if os.path.isdir(path):
                bucket_dirs.append(path)
        return bucket_dirs

    def enqueue_all_existing_buckets(self):
        if not self.enabled:
            return
        for bucket_dir in self._list_existing_bucket_dirs():
            self.enqueue_bucket_dir(bucket_dir, blocking=True, timeout=30.0, force=True)

    def shutdown(self, wait: bool = True, timeout: float = 120.0):
        if not self.enabled:
            return
        if self._proc is None:
            return
        if self._task_q is not None:
            deadline = time.time() + 30.0
            while True:
                try:
                    self._task_q.put(_STOP, block=True, timeout=1.0)
                    break
                except queue.Full:
                    if time.time() >= deadline:
                        self._log("warning", "Sample archive worker stop signal timeout")
                        break
                except Exception:
                    break
        if wait:
            self._proc.join(timeout=timeout)
        if self._proc.is_alive():
            self._log("warning", "Sample archive worker did not exit in timeout")

        summary = None
        if self._status_q is not None:
            try:
                summary = self._status_q.get_nowait()
            except Exception:
                summary = None
        if summary is not None:
            self._log(
                "warning",
                format_two_column_log(
                    "Sample archive summary",
                    [
                        ("enqueued", summary.get("enqueued", 0)),
                        ("packed", summary.get("packed", 0)),
                        ("skipped", summary.get("skipped", 0)),
                        ("failed", summary.get("failed", 0)),
                        ("records", summary.get("records", 0)),
                        ("manifest", summary.get("manifest_path", self.manifest_jsonl_path)),
                    ],
                ),
            )

        self._proc = None
        self._task_q = None
        self._status_q = None
        self._scheduled.clear()

    def _log(self, level: str, message: str):
        if self.logger is None:
            return
        fn = getattr(self.logger, level, None)
        if callable(fn):
            fn(message)
