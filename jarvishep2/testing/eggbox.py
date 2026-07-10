#!/usr/bin/env python3
"""V1-compatible eggbox operators (spawn-safe package path)."""

from __future__ import annotations

import numpy as np


def eggbox_numpy(x: float, y: float, logger=None, **_params: object) -> float:
    if logger is not None:
        logger.info("EggBox input loaded: x -> %s, y -> %s", x, y)
    z = (np.sin(x) * np.cos(y) + 2.0) ** 5
    return float(z)


def eggbox2d_numpy(x: float, y: float, logger=None, **_params: object) -> dict[str, float]:
    return {"z": eggbox_numpy(x, y, logger=logger)}


def circle_r2(x: float, y: float, logger=None, **_params: object) -> dict[str, float]:
    """Level-set fixture: f = (x - 0.5)^2 + (y - 0.5)^2 on the unit box."""
    return {"r2": float((float(x) - 0.5) ** 2 + (float(y) - 0.5) ** 2)}