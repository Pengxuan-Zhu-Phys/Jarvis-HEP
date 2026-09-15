"""Shared formatting rules for Monitor metrics."""

from __future__ import annotations

import math


def format_sample_rate(samples_per_sec: float | None) -> str:
    """Format a sample rate using the Monitor-wide adaptive unit policy."""
    if samples_per_sec is None or not math.isfinite(samples_per_sec):
        return "—"
    per_sec = max(0.0, float(samples_per_sec))
    if per_sec >= 1.0:
        return f"{per_sec:0.1f} / sec"
    per_min = per_sec * 60.0
    if per_min >= 1.0:
        return f"{per_min:0.1f} / min"
    return f"{per_min * 60.0:0.1f} / hour"


__all__ = ["format_sample_rate"]
