#!/usr/bin/env python3
"""Staging handoff helpers for Worker → Archiver Layer-2 I/O (WP-D4)."""

from __future__ import annotations

import os
import shutil
from typing import Any

_VALID_MOVE_STRATEGIES = frozenset({"move", "copy"})


def normalize_move_strategy(value: object | None) -> str:
    strategy = str(value or "move").strip().lower()
    return strategy if strategy in _VALID_MOVE_STRATEGIES else "move"


def resolve_staging_dir(
    task_result_dir: str,
    *,
    configured: str | None = None,
) -> str:
    """Return the absolute staging root for finished Samples."""
    if configured and str(configured).strip():
        path = os.path.abspath(os.path.expanduser(str(configured)))
    else:
        path = os.path.join(os.path.abspath(task_result_dir), "staging")
    os.makedirs(path, exist_ok=True)
    return path


def move_tree(src: str, dst: str, *, strategy: str = "move") -> str:
    """Move or copy *src* to *dst*, preferring same-volume rename."""
    source = os.path.abspath(str(src))
    target = os.path.abspath(str(dst))
    if not os.path.lexists(source):
        raise FileNotFoundError(f"handoff source does not exist: {source}")
    if source in {"", "/", os.path.sep}:
        raise ValueError(f"refusing to move unsafe source path: {src!r}")

    normalized = normalize_move_strategy(strategy)
    parent = os.path.dirname(target)
    if parent:
        os.makedirs(parent, exist_ok=True)

    if os.path.lexists(target):
        if os.path.isdir(target) and not os.path.islink(target):
            shutil.rmtree(target)
        else:
            os.remove(target)

    if normalized == "copy":
        if os.path.isdir(source) and not os.path.islink(source):
            shutil.copytree(source, target)
            shutil.rmtree(source)
        else:
            shutil.copy2(source, target)
            os.remove(source)
        return target

    try:
        os.rename(source, target)
    except OSError:
        if os.path.isdir(source) and not os.path.islink(source):
            shutil.move(source, target)
        else:
            shutil.move(source, target)
    return target


def stage_sample_dir(work_dir: str, staging_root: str, sample_uuid: str) -> str:
    """Move a materialized Sample work dir into ``staging/<uuid>``."""
    staging_path = os.path.join(os.path.abspath(staging_root), str(sample_uuid))
    return move_tree(work_dir, staging_path, strategy="move")


def archive_staging_to_sample(
    staging_path: str,
    sample_root: str,
    sample_uuid: str,
    *,
    strategy: str = "move",
) -> str:
    """Persist staged products under ``SAMPLE/<uuid>/`` (idempotent)."""
    destination = os.path.join(os.path.abspath(sample_root), str(sample_uuid))
    if os.path.isdir(destination):
        return destination
    return move_tree(staging_path, destination, strategy=strategy)


def list_product_names(sample_dir: str) -> list[str]:
    """Return sorted product filenames for archive metadata."""
    root = os.path.abspath(str(sample_dir))
    if not os.path.isdir(root):
        return []
    return sorted(
        name
        for name in os.listdir(root)
        if not name.startswith(".")
    )


def same_volume_move_preserves_inode(src: str, dst: str) -> bool:
    """Return True when *dst* reuses the inode from a same-volume rename."""
    if not os.path.lexists(dst):
        return False
    try:
        return os.stat(src, follow_symlinks=False).st_ino == os.stat(
            dst, follow_symlinks=False
        ).st_ino
    except FileNotFoundError:
        return os.stat(dst, follow_symlinks=False).st_ino > 0


__all__ = [
    "archive_staging_to_sample",
    "list_product_names",
    "move_tree",
    "normalize_move_strategy",
    "resolve_staging_dir",
    "same_volume_move_preserves_inode",
    "stage_sample_dir",
]