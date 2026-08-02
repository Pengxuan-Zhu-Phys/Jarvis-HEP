"""Small deterministic operator used by crash/resume integration tests."""

from __future__ import annotations

import time
from typing import Any


def slow_sum(x: Any, y: Any, *, delay: float = 0.002) -> dict[str, float]:
    """Return a stable observable after a short, intentional delay."""
    time.sleep(max(0.0, float(delay)))
    return {"z": float(x) + float(y)}


def slow_circle_r2(x: Any, y: Any, *, delay: float = 0.01) -> dict[str, float]:
    """Return the deterministic circle fixture slowly enough to interrupt."""
    time.sleep(max(0.0, float(delay)))
    dx = float(x) - 0.5
    dy = float(y) - 0.5
    return {"r2": dx * dx + dy * dy}


__all__ = ["slow_circle_r2", "slow_sum"]
