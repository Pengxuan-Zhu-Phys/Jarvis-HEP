#!/usr/bin/env python3
"""Built-in functions for likelihood/expression evaluation in Jarvis-HEP V2."""

from __future__ import annotations

import math
from typing import Any

import numpy as np


def LogGauss(xx: Any, mean: Any, err: Any) -> float:
    """Log of a Gaussian (V1-compatible)."""
    x = float(np.asarray(xx).item() if isinstance(xx, np.generic) else xx)
    m = float(np.asarray(mean).item() if isinstance(mean, np.generic) else mean)
    e = float(np.asarray(err).item() if isinstance(err, np.generic) else err)
    return float(-0.5 * ((x - m) / e) ** 2)


NUMERIC_MODULES: dict[str, Any] = {
    "sin": np.sin,
    "cos": np.cos,
    "exp": np.exp,
    "log": np.log,
    "sqrt": math.sqrt,
    "LogGauss": LogGauss,
}


__all__ = ["LogGauss", "NUMERIC_MODULES"]