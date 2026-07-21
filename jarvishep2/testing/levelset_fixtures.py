#!/usr/bin/env python3
"""Spawn-safe operators for AdaptiveBridson tests.

Single source: re-exports from ``eggbox`` so Workers and tests share one path.
"""

from __future__ import annotations

from jarvishep2.testing.eggbox import (
    circle_r2,
    circle_r2_region_fail,
    ellipse_s,
    hypersphere_r,
    sphere_r,
)

__all__ = [
    "circle_r2",
    "circle_r2_region_fail",
    "ellipse_s",
    "hypersphere_r",
    "sphere_r",
]
