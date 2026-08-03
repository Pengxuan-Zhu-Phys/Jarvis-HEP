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


def sampling_bounds(sampling: Mapping[str, Any] | None) -> Mapping[str, Any]:
    """Return ``Sampling.Bounds`` (empty mapping when absent).

    All method-specific sampler knobs live under Bounds.  Top-level
    ``Sampling`` only holds Method, Variables, selection, LogLikelihood,
    Nuisance, FeedbackReturn, and Bounds itself.
    """
    if not isinstance(sampling, Mapping):
        return {}
    raw = sampling.get("Bounds")
    return raw if isinstance(raw, Mapping) else {}


def bounds_seed(bounds: Mapping[str, Any] | None, default: int = 0) -> int:
    """Read Seed / seed / rseed from Bounds (canonical V2 location)."""
    if not isinstance(bounds, Mapping):
        return int(default)
    for key in ("Seed", "seed", "rseed"):
        if key in bounds and bounds.get(key) is not None:
            try:
                return int(bounds.get(key) or 0)
            except (TypeError, ValueError):
                return int(default)
    return int(default)


def is_check_modules_mode(config: Mapping[str, Any]) -> bool:
    sampling = sampling_block(config) or {}
    mode = str(sampling.get("mode") or "").strip().lower()
    return mode in {"check_modules", "check-modules"}


# Retired Sampling top-level keys (pre-Bounds placement).  Completely rejected.
RETIRED_SAMPLING_TOP_LEVEL: dict[str, str] = {
    "Radius": "Sampling.Bounds.Radius",
    "MaxAttempt": "Sampling.Bounds.MaxAttempt",
    "MaxWorker": "Sampling.Bounds.MaxWorker",
    "Point number": "Sampling.Bounds['Point number']",
    "point_number": "Sampling.Bounds.point_number",
    "Seed": "Sampling.Bounds.Seed",
    "seed": "Sampling.Bounds.seed",
    "CSV": "Sampling.Bounds (path / uuid_column / variables)",
    "AdaptiveBridson": "Sampling.Bounds",
    "adaptive_bridson": "Sampling.Bounds",
}


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
