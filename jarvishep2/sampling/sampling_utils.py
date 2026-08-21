#!/usr/bin/env python3
"""Shared helpers for V2 samplers."""

from __future__ import annotations

from typing import Any, Mapping

from jarvishep2.sampling.variables import Variable
from jarvishep2.expression import ExpressionContext


class BoolConversionError(ValueError):
    pass


_SELECTION_EXPRESSIONS = ExpressionContext()
_SELECTION_OPERAS_LOADED = False


def evaluate_selection(
    expression: str | None,
    variables: Mapping[str, Any],
    *,
    context: ExpressionContext | None = None,
) -> bool:
    if expression is None:
        return True
    if not isinstance(variables, Mapping):
        raise BoolConversionError("Selection variables must be a mapping.")

    try:
        normalized_variables = {str(name): value for name, value in variables.items()}
        available_names = tuple(sorted(normalized_variables))
        global _SELECTION_EXPRESSIONS, _SELECTION_OPERAS_LOADED
        if context is None:
            from jarvishep2.operas_functions import (
                build_operas_expression_context,
                expression_uses_operas_function,
            )

            if (
                not _SELECTION_OPERAS_LOADED
                and expression_uses_operas_function(expression)
            ):
                _SELECTION_EXPRESSIONS = build_operas_expression_context(required=True)
                _SELECTION_OPERAS_LOADED = True
            context = _SELECTION_EXPRESSIONS
        compiled = context.compile(
            str(expression),
            symbols=available_names,
        )
        return bool(compiled.evaluate(normalized_variables))
    except Exception as exc:
        raise BoolConversionError(
            f"Cannot evaluate selection expression '{expression}' as boolean."
        ) from exc


def map_u_to_physical(u_coords, variables: list[Variable]) -> dict[str, float]:
    """Distribution-only mapping (no derive). Prefer :class:`MapperPipeline.map`.

    Kept as the internal fast path used when a sampler has not yet attached a
    pipeline, and for unit tests that exercise Variables without Mapper.
    """
    import numpy as np

    coords = np.asarray(u_coords, dtype=np.float64).reshape(-1)
    if len(coords) != len(variables):
        raise ValueError(
            f"u_coords length {len(coords)} must equal variable count {len(variables)}"
        )
    mapped: dict[str, float] = {}
    for index, var in enumerate(variables):
        mapped[var.name] = var.map_standard_random_to_distribution(float(coords[index]))
    return mapped


def map_row_to_physical(row, variables: list[Variable]) -> dict[str, float]:
    """Distribution-only row mapping. Prefer pipeline.map(row_to_u_coords(...))."""
    if len(row) != len(variables):
        raise ValueError(
            f"row length {len(row)} must equal variable count {len(variables)}"
        )
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


def physical_from_u(u_coords, variables: list[Variable], pipeline=None) -> dict[str, float]:
    """Map unit-cube coords via pipeline when present, else distribution-only."""
    if pipeline is not None:
        return pipeline.map(u_coords)
    return map_u_to_physical(u_coords, variables)


__all__ = [
    "BoolConversionError",
    "evaluate_selection",
    "map_row_to_physical",
    "map_u_to_physical",
    "physical_from_u",
    "row_to_u_coords",
]
