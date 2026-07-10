#!/usr/bin/env python3
"""Spawn-safe operators for AdaptiveLevelSet tests."""

from __future__ import annotations


def circle_r2(x: float, y: float, logger=None, **_params: object) -> dict[str, float]:
    """f = (x - 0.5)^2 + (y - 0.5)^2 for a circle level-set around the unit box centre."""
    return {"r2": float((float(x) - 0.5) ** 2 + (float(y) - 0.5) ** 2)}


__all__ = ["circle_r2"]
