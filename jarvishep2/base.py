#!/usr/bin/env python3
"""Path and token helpers for Jarvis-HEP V2 (WP-D3.1)."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

_PROJECT_MARKERS = (".jarvis-project.json", "jarvis.project.yaml")
_TASK_ROOT_ENV_VARS = ("JARVIS_HEP_TASK_ROOT", "JHEP_TASK_ROOT")


def env_task_root() -> str | None:
    for env_name in _TASK_ROOT_ENV_VARS:
        value = os.getenv(env_name, "").strip()
        if value:
            return os.path.abspath(os.path.expanduser(value))
    return None


def infer_project_root(start: str | None = None) -> str:
    """Walk up from *start* to find a Jarvis project anchor."""
    path = Path(start or os.getcwd()).expanduser().resolve()
    for candidate in [path, *path.parents]:
        for marker in _PROJECT_MARKERS:
            if (candidate / marker).exists():
                return str(candidate)
        if candidate.name.lower() == "bin":
            return str(candidate.parent)
    return str(path)


def expand_j(text: str, *, project_root: str) -> str:
    """Replace ``&J/`` / ``&J`` with the absolute project root."""
    if text is None:
        return ""
    raw = str(text)
    normalized = raw.replace("\\", "/")
    if normalized.startswith("&J/") and "/src/card/" in normalized:
        raise ValueError(
            f"Legacy card path prefix is no longer supported: {raw}. "
            "Use project-local packaged copies such as '&J/deps/...'."
        )
    if "&J" in raw:
        root = str(project_root).rstrip("/")
        raw = raw.replace("&J/", root + "/").replace("&J", root)
    return raw


def decode_path(
    path: str | None,
    *,
    project_root: str,
    base_dir: str | None = None,
) -> str:
    """Resolve ``&J/``, ``~``, and relative paths to an absolute path."""
    if path is None or not isinstance(path, str):
        return str(path or "")

    resolved = expand_j(path, project_root=project_root)
    if "~" in resolved:
        resolved = os.path.expanduser(resolved)
    if "://" in resolved:
        return resolved
    if resolved.startswith("&"):
        raise ValueError(f"Unable to resolve path: {resolved}")
    if os.path.isabs(resolved):
        return os.path.abspath(resolved)

    anchor = base_dir or project_root
    return os.path.abspath(os.path.join(anchor, resolved))


def scan_output_root(*, project_root: str, scan_name: str) -> str:
    return os.path.join(project_root, "outputs", str(scan_name).strip() or "default")


__all__ = [
    "decode_path",
    "env_task_root",
    "expand_j",
    "infer_project_root",
    "scan_output_root",
]