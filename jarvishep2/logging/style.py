#!/usr/bin/env python3
"""Load process/sample log presentation styles from package card/logging.yaml."""

from __future__ import annotations

import os
from copy import deepcopy
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

# Defaults identical to the historical hard-coded V1 visual contract.
_DEFAULT_STYLE: dict[str, Any] = {
    "version": 1,
    "process": {
        "date_format": "%m-%d %H:%M:%S",
        "millis_digits": 3,
        "bullet": "·•·",
        # DATABASE / samples.hdf5 progress (Archiver flush counter).
        "hdf5_writer_bullet": "Ϡ",
        "hdf5_writer_module": "Jarvis-HEP.DataRecorder",
        "head": "\n{bullet} {module} \n\t-> {timestamp} - [{level}] >>> \n",
        "append_context": True,
    },
    "sample": {
        "date_format": "%m-%d %H:%M:%S.%f",
        "frac_digits": 3,
        "bullet": "·•·",
        "head": "\n{bullet} {module} \n\t-> {timestamp} - [{level}] >>> \n",
        "double_newline_gap": True,
    },
}


def default_logging_style_path() -> str:
    """Path to packaged ``card/logging.yaml``."""
    here = Path(__file__).resolve().parent.parent  # jarvishep2/
    return str(here / "card" / "logging.yaml")


def _merge_section(base: dict[str, Any], overlay: Mapping[str, Any] | None) -> dict[str, Any]:
    out = dict(base)
    if not isinstance(overlay, Mapping):
        return out
    for key, value in overlay.items():
        out[str(key)] = value
    return out


def normalize_logging_style(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Merge user/package YAML onto built-in defaults."""
    style = deepcopy(_DEFAULT_STYLE)
    if not isinstance(raw, Mapping):
        return style
    if "version" in raw:
        style["version"] = raw["version"]
    for section in ("process", "sample"):
        section_raw = raw.get(section)
        if isinstance(section_raw, Mapping):
            style[section] = _merge_section(dict(style[section]), section_raw)
    return style


def load_logging_style(path: str | None = None) -> dict[str, Any]:
    """Load style from YAML path (or package default). Missing file → defaults."""
    resolved = os.path.abspath(os.path.expanduser(str(path or default_logging_style_path())))
    if not os.path.isfile(resolved):
        return deepcopy(_DEFAULT_STYLE)
    try:
        import yaml
    except Exception:
        return deepcopy(_DEFAULT_STYLE)
    try:
        with open(resolved, encoding="utf-8") as handle:
            payload = yaml.safe_load(handle) or {}
    except Exception:
        return deepcopy(_DEFAULT_STYLE)
    if not isinstance(payload, Mapping):
        return deepcopy(_DEFAULT_STYLE)
    return normalize_logging_style(payload)


@lru_cache(maxsize=4)
def get_logging_style(path: str | None = None) -> dict[str, Any]:
    """Cached style loader. Pass ``path`` or None for package default."""
    return load_logging_style(path)


def clear_logging_style_cache() -> None:
    get_logging_style.cache_clear()


def process_style(path: str | None = None) -> dict[str, Any]:
    return dict(get_logging_style(path).get("process") or {})


def sample_style(path: str | None = None) -> dict[str, Any]:
    return dict(get_logging_style(path).get("sample") or {})


__all__ = [
    "clear_logging_style_cache",
    "default_logging_style_path",
    "get_logging_style",
    "load_logging_style",
    "normalize_logging_style",
    "process_style",
    "sample_style",
]
