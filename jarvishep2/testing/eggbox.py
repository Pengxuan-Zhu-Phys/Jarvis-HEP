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


def ellipse_s(
    x: float,
    y: float,
    logger=None,
    *,
    a: float = 0.25,
    b: float = 0.15,
    **_params: object,
) -> dict[str, float]:
    """Level-set fixture: s = ((x-0.5)/a)^2 + ((y-0.5)/b)^2; contour at s=1."""
    xx = (float(x) - 0.5) / float(a)
    yy = (float(y) - 0.5) / float(b)
    return {"s": float(xx * xx + yy * yy)}


def circle_r2_region_fail(
    x: float,
    y: float,
    logger=None,
    **_params: object,
) -> dict[str, float]:
    """Circle fixture that fails for x < 0.15 (exercises failed_regions)."""
    if float(x) < 0.15:
        raise RuntimeError("synthetic failure region x < 0.15")
    return circle_r2(x, y, logger=logger)