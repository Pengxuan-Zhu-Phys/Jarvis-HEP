#!/usr/bin/env python3
"""Shared helpers for task-YAML contracts."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any

# Re-export issue type only at type-check time to avoid cycles; callers pass issues lists.


def is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def as_mapping(value: Any) -> Mapping[str, Any] | None:
    return value if isinstance(value, Mapping) else None


def sampling_block(config: Mapping[str, Any]) -> Mapping[str, Any] | None:
    raw = config.get("Sampling")
    return raw if isinstance(raw, Mapping) else None


def is_check_modules_mode(config: Mapping[str, Any]) -> bool:
    sampling = sampling_block(config) or {}
    mode = str(sampling.get("mode") or "").strip().lower()
    return mode in {"check_modules", "check-modules"}


def unknown_keys(
    block: Mapping[str, Any],
    allowed: Iterable[str],
) -> list[str]:
    allowed_set = frozenset(str(k) for k in allowed)
    return sorted(str(k) for k in block.keys() if str(k) not in allowed_set)


def try_int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def try_float(value: Any) -> float | None:
    if isinstance(value, bool):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def try_bool(value: Any) -> bool | None:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"1", "true", "yes", "on"}:
            return True
        if text in {"0", "false", "no", "off"}:
            return False
    return None
