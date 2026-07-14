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

import json
import os
import shutil
import tarfile
import time
from collections.abc import Mapping
from typing import Any


SAMPLE_DIRECTORY_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "limit": 200,
    "width": 6,
    "pack": True,
    "start_bucket": 1,
}


def normalize_sample_directory(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize Scan/EnvReqs sample_directory settings."""
    cfg = dict(SAMPLE_DIRECTORY_DEFAULTS)
    if not isinstance(raw, Mapping):
        return cfg
    if "enabled" in raw:
        cfg["enabled"] = bool(raw.get("enabled"))
    if "pack" in raw:
        cfg["pack"] = bool(raw.get("pack"))
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


def pack_bucket_dir(
    bucket_dir: str,
    *,
    sample_root: str | None = None,
    manifest_jsonl_path: str | None = None,
    prune: bool = True,
) -> str | None:
    """Tar a sealed SAMPLE bucket directory. Returns archive path or None if skipped."""
    source = os.path.abspath(str(bucket_dir))
    if not os.path.isdir(source):
        return None
    root = os.path.abspath(sample_root or os.path.dirname(source))
    bucket_name = os.path.basename(source.rstrip(os.sep))
    primary = os.path.join(root, f"{bucket_name}.tar.gz")
    if os.path.exists(primary):
        if prune:
            shutil.rmtree(source, ignore_errors=True)
        return primary

    out_path = resolve_archive_tar_path(root, bucket_name)
    tmp_path = out_path + ".tmp"
    try:
        with tarfile.open(tmp_path, "w:gz") as tar:
            tar.add(source, arcname=bucket_name)
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
        with open(manifest_jsonl_path, "a", encoding="utf-8") as handle:
            handle.write(json.dumps(payload, ensure_ascii=True) + "\n")
    return out_path


__all__ = [
    "SAMPLE_DIRECTORY_DEFAULTS",
    "bucket_dir_path",
    "format_bucket_name",
    "normalize_sample_directory",
    "pack_bucket_dir",
    "parse_bucket_id",
    "resolve_archive_tar_path",
]
