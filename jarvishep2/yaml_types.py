#!/usr/bin/env python3
"""Strict primitive conversions for user-authored YAML values."""

from __future__ import annotations

from typing import Any


def require_bool(value: Any, *, field: str) -> bool:
    """Require a native YAML boolean for a user-facing boolean field."""
    if isinstance(value, bool):
        return value
    raise ValueError(
        f"{field} must be a YAML boolean (true or false), got {value!r}"
    )
