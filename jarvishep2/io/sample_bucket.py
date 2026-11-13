#!/usr/bin/env python3
"""SAMPLE bucket layout + tar packing helpers (V1 sample_directory parity).

Layout::

    SAMPLE/
      000001/
        <uuid>/...
        <uuid>/...
      000001.tar.gz   # after Archiver packs a sealed, idle bucket
"""

from __future__ import annotations

import gzip
import json
import os
import re
import shutil
import tarfile
import threading
import time
from collections.abc import Callable, Mapping
from typing import Any

from jarvishep2.yaml_types import require_bool


_MANIFEST_LOCK = threading.Lock()


SAMPLE_DIRECTORY_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "limit": 200,
    "width": 6,
    "pack": True,
    "start_bucket": 1,
}


def normalize_sample_directory(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize the EnvReqs.V2.sample_directory settings."""
    cfg = dict(SAMPLE_DIRECTORY_DEFAULTS)
    if not isinstance(raw, Mapping):
        return cfg
    if "enabled" in raw:
        cfg["enabled"] = require_bool(
            raw.get("enabled"), field="EnvReqs.V2.sample_directory.enabled"
        )
    if "pack" in raw:
        cfg["pack"] = require_bool(
            raw.get("pack"), field="EnvReqs.V2.sample_directory.pack"
        )
    for key, default in (("limit", 200), ("width", 6), ("start_bucket", 1)):
        if key not in raw:
            continue
        try:
            value = int(raw.get(key))
        except (TypeError, ValueError):
            value = int(default)
        cfg[key] = max(1, value)
    return cfg


def format_bucket_name(bucket_id: int, *, width: int = 6) -> str:
    return f"{int(bucket_id):0{max(1, int(width))}d}"


def bucket_dir_path(sample_root: str, bucket_id: int, *, width: int = 6) -> str:
    return os.path.join(os.path.abspath(sample_root), format_bucket_name(bucket_id, width=width))


def parse_bucket_id(bucket_dir: str | int | None) -> int | None:
    if bucket_dir is None:
        return None
    if isinstance(bucket_dir, int):
        return int(bucket_dir)
    name = os.path.basename(str(bucket_dir).rstrip(os.sep))
    if not name.isdigit():
        return None
    try:
        return int(name)
    except ValueError:
        return None


def archive_manifest_path(sample_root: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(sample_root)), "DATABASE", "archive_manifest.jsonl")


def resume_bucket_start(sample_root: str, *, start_bucket: int) -> int:
    """Reserve bucket numbers already used in this scan's output."""
    highest = 0
    try:
        with os.scandir(sample_root) as entries:
            for entry in entries:
                archive = re.fullmatch(r"([0-9]+)(?:__[0-9]+)?\.tar\.gz(?:\.tmp)?", entry.name)
                bucket = archive.group(1) if archive else entry.name
                highest = max(highest, parse_bucket_id(bucket) or 0)
    except FileNotFoundError:
        pass

    try:
        with open(archive_manifest_path(sample_root), "rb") as manifest:
            for line in manifest:
                try:
                    record = json.loads(line)
                except ValueError:
                    # An interrupted append may leave an incomplete record.
                    continue
                if isinstance(record, dict):
                    highest = max(highest, parse_bucket_id(record.get("bucket")) or 0)
    except FileNotFoundError:
        pass

    return max(start_bucket, highest + 1)


def resolve_archive_tar_path(sample_root: str, bucket_name: str) -> str:
    """Return ``SAMPLE/<bucket>.tar.gz`` or a collision-safe ``__NNNN`` suffix."""
    target = os.path.join(sample_root, f"{bucket_name}.tar.gz")
    if not os.path.exists(target):
        return target
    for idx in range(1, 100_000):
        candidate = os.path.join(sample_root, f"{bucket_name}__{idx:04d}.tar.gz")
        if not os.path.exists(candidate):
            return candidate
    raise RuntimeError(f"Cannot allocate archive tar path for bucket={bucket_name}")


class _BucketOutput:
    """Check for shutdown as tar writes each chunk, including large files."""

    def __init__(self, handle: Any, cancelled: Callable[[], bool] | None) -> None:
        self.handle = handle
        self.cancelled = cancelled

    def write(self, data: bytes) -> int:
        if self.cancelled is not None and self.cancelled():
            raise InterruptedError("SAMPLE bucket packing cancelled")
        return self.handle.write(data)


def pack_bucket_dir(
    bucket_dir: str,
    *,
    sample_root: str | None = None,
    manifest_jsonl_path: str | None = None,
    prune: bool = True,
    cancelled: Callable[[], bool] | None = None,
) -> str | None:
    """Tar a sealed SAMPLE bucket directory. Returns archive path or None if skipped."""
    source = os.path.abspath(str(bucket_dir))
    if not os.path.isdir(source):
        return None
    root = os.path.abspath(sample_root or os.path.dirname(source))
    bucket_name = os.path.basename(source.rstrip(os.sep))
    out_path = resolve_archive_tar_path(root, bucket_name)
    tmp_path = out_path + ".tmp"
    try:
        with gzip.open(tmp_path, "wb") as compressed:
            with tarfile.open(fileobj=_BucketOutput(compressed, cancelled), mode="w|") as tar:
                tar.add(source, arcname=bucket_name)
        if cancelled is not None and cancelled():
            raise InterruptedError("SAMPLE bucket packing cancelled")
        os.replace(tmp_path, out_path)
    except Exception:
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except OSError:
            pass
        raise

    if prune:
        shutil.rmtree(source, ignore_errors=True)

    if manifest_jsonl_path:
        parent = os.path.dirname(manifest_jsonl_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        payload = {
            "bucket": bucket_name,
            "source_path": source,
            "archive_path": out_path,
            "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "mode": "tar_gz_prune" if prune else "tar_gz",
        }
        with _MANIFEST_LOCK, open(manifest_jsonl_path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, ensure_ascii=True) + "\n")
    return out_path


__all__ = [
    "SAMPLE_DIRECTORY_DEFAULTS",
    "archive_manifest_path",
    "bucket_dir_path",
    "format_bucket_name",
    "normalize_sample_directory",
    "pack_bucket_dir",
    "parse_bucket_id",
    "resolve_archive_tar_path",
    "resume_bucket_start",
]
