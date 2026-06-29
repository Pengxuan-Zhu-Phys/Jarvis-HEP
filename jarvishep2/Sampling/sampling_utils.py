#!/usr/bin/env python3
"""Shared helpers for V2 samplers."""

from __future__ import annotations

from typing import Any, Mapping

import sympy as sp

from jarvishep2.Sampling.variables import Variable


class BoolConversionError(ValueError):
    pass


def evaluate_selection(expression: str | None, variables: Mapping[str, Any]) -> bool:
    if expression is None:
        return True
    if not isinstance(variables, Mapping):
        raise BoolConversionError("Selection variables must be a mapping.")

    symbols = {name: sp.symbols(name) for name in variables.keys()}
    try:
        expr = sp.sympify(expression, locals=symbols)
        evaluated = expr.subs(variables)
        return bool(evaluated)
    except Exception as exc:
        raise BoolConversionError(
            f"Cannot evaluate selection expression '{expression}' as boolean."
        ) from exc


def map_u_to_physical(u_coords, variables: list[Variable]) -> dict[str, float]:
    import numpy as np

    coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
    mapped: dict[str, float] = {}
    for index, var in enumerate(variables):
        mapped[var.name] = var.map_standard_random_to_distribution(float(coords[index]))
    return mapped


def map_row_to_physical(row, variables: list[Variable]) -> dict[str, float]:
    mapped: dict[str, float] = {}
    for index, var in enumerate(variables):
        length = float(var.parameters.get("length", 1.0) or 1.0)
        std_rand = float(row[index]) / length
        mapped[var.name] = var.map_standard_random_to_distribution(std_rand)
    return mapped


def row_to_u_coords(row, variables: list[Variable]):
    import numpy as np

    coords = []
    for index, var in enumerate(variables):
        length = float(var.parameters.get("length", 1.0) or 1.0)
        coords.append(float(row[index]) / length)
    return np.asarray(coords, dtype=np.float64)


__all__ = [
    "BoolConversionError",
    "evaluate_selection",
    "map_row_to_physical",
    "map_u_to_physical",
    "row_to_u_coords",
]