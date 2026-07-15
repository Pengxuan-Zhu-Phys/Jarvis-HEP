#!/usr/bin/env python3
"""HEP-owned file operations (delete + SAMPLE save/copy).

Per DESIGN_PORTAL_IO_2.0: Portal owns format serialization only; **save/copy/delete
policy lives in V2**. YAML ``save: true`` is unchanged — execution is here (and
optionally on a dedicated FileOperation process — see ``file_operation_service``).
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

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


def sample_artifact_filename(source_path: str, *, module: str | None = None) -> str:
    """V1 naming: ``basename`` or ``basename@Module``."""
    base = os.path.basename(str(source_path))
    mod = str(module or "").strip()
    return f"{base}@{mod}" if mod else base


def save_io_copy(
    source_path: str,
    sample_save_dir: str,
    *,
    module: str | None = None,
    save: bool = False,
    as_temp: bool = False,
    content: str | bytes | None = None,
) -> str | None:
    """Copy a calculator IO file into SAMPLE (V1 ``save`` / ``.temp`` policy).

    Parameters
    ----------
    source_path:
        Runtime file just written/read by Portal.
    sample_save_dir:
        ``SAMPLE/<bucket>/<uuid>`` (or equivalent).
    save:
        YAML ``save: true`` — land under sample_save_dir.
    as_temp:
        V1 output default: when not ``save``, still copy under ``.temp/``.
    content:
        Optional pre-read text/bytes; if None, read from ``source_path``.

    Returns
    -------
    Absolute path of the SAMPLE copy, or None if neither save nor as_temp.
    """
    do_save = bool(save)
    do_temp = bool(as_temp) and not do_save
    if not do_save and not do_temp:
        return None
    root = str(sample_save_dir or "").strip()
    if not root:
        return None

    base_dir = Path(root)
    if do_temp:
        base_dir = base_dir / ".temp"
    base_dir.mkdir(parents=True, exist_ok=True)

    filename = sample_artifact_filename(source_path, module=module)
    target = base_dir / filename
    src = os.path.abspath(str(source_path))

    if content is not None:
        if isinstance(content, bytes):
            target.write_bytes(content)
        else:
            target.write_text(str(content), encoding="utf-8")
    else:
        if not os.path.isfile(src):
            return None
        shutil.copy2(src, target)
    return str(target.resolve())


def apply_io_save_policy(
    *,
    source_path: str,
    sample_save_dir: str | None,
    module: str | None,
    spec: dict[str, Any],
    direction: str,
    content: str | bytes | None = None,
) -> str | None:
    """Apply YAML ``save`` policy after Portal format R/W.

    - input: copy only when ``save: true``
    - output: copy when ``save: true`` else under ``.temp`` (V1)
    """
    save = bool(spec.get("save", False))
    direction_l = str(direction or "").strip().lower()
    as_temp = direction_l == "output" and not save
    if not save and not as_temp:
        return None
    if not sample_save_dir:
        return None
    return save_io_copy(
        source_path,
        str(sample_save_dir),
        module=module,
        save=save,
        as_temp=as_temp,
        content=content,
    )


__all__ = [
    "DEFAULT_DELETE_METHOD",
    "apply_io_save_policy",
    "delete_path",
    "delete_paths",
    "normalize_delete_method",
    "sample_artifact_filename",
    "save_io_copy",
]
