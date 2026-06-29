#!/usr/bin/env python3
"""Configurable file deletion helpers (WP-D3.3)."""

from __future__ import annotations

import os
import shutil
import subprocess

DEFAULT_DELETE_METHOD = "shutil"
_VALID_DELETE_METHODS = frozenset({"shutil", "rm"})


def normalize_delete_method(value: object | None) -> str:
    """Normalize ``Runtime.FileOperation.delete_method`` with v1-equivalent default."""
    method = str(value or DEFAULT_DELETE_METHOD).strip().lower()
    if method not in _VALID_DELETE_METHODS:
        return DEFAULT_DELETE_METHOD
    return method


def delete_path(path: str, *, method: str = DEFAULT_DELETE_METHOD, missing_ok: bool = False) -> None:
    """Delete a file or directory using the configured backend."""
    normalized = str(method or DEFAULT_DELETE_METHOD).strip().lower()
    if normalized not in _VALID_DELETE_METHODS:
        allowed = ", ".join(sorted(_VALID_DELETE_METHODS))
        raise ValueError(f"invalid delete_method '{method}'; allowed: {allowed}")

    target = os.path.abspath(os.path.expanduser(str(path)))
    if target in {"", "/", os.path.sep}:
        raise ValueError(f"refusing to delete unsafe path: {path!r}")

    if not os.path.lexists(target):
        if missing_ok:
            return
        raise FileNotFoundError(f"delete_path target does not exist: {target}")

    if normalized == "shutil":
        _delete_with_shutil(target)
        return

    completed = subprocess.run(
        ["rm", "-rf", target],
        capture_output=True,
        text=True,
        check=False,
    )
    if int(completed.returncode) != 0:
        stderr = (completed.stderr or "").strip()
        raise RuntimeError(
            f"rm delete failed for '{target}' (exit {completed.returncode}): {stderr}"
        )


def delete_paths(
    paths: list[str] | tuple[str, ...],
    *,
    method: str = DEFAULT_DELETE_METHOD,
    missing_ok: bool = True,
) -> None:
    """Delete multiple paths with the same backend."""
    for path in paths:
        text = str(path).strip()
        if text:
            delete_path(text, method=method, missing_ok=missing_ok)


def _delete_with_shutil(target: str) -> None:
    if os.path.islink(target) or os.path.isfile(target):
        os.remove(target)
    elif os.path.isdir(target):
        shutil.rmtree(target)
    else:
        os.remove(target)


__all__ = [
    "DEFAULT_DELETE_METHOD",
    "delete_path",
    "delete_paths",
    "normalize_delete_method",
]