#!/usr/bin/env python3
"""Minimal JSON IO helpers for CalculatorModule (parity-project scope)."""

from __future__ import annotations

import json
import math
import os
from collections.abc import Mapping
from typing import Any

import numpy as np
import sympy as sp


def _to_json_compatible(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Mapping):
        return {str(k): _to_json_compatible(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_json_compatible(v) for v in value]
    return value


def _coerce_numeric_param(value: Any) -> float | None:
    if isinstance(value, bool):
        return float(int(value))
    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(np.asarray(value).item() if isinstance(value, np.generic) else value)
    if isinstance(value, str):
        try:
            return float(value)
        except ValueError:
            return None
    return None


def _evaluate_dump_value(var: Mapping[str, Any], param_values: Mapping[str, Any]) -> Any:
    expression = str(var.get("expression", "")).strip()
    if expression:
        expr = sp.sympify(expression, locals={"Pi": sp.pi, "pi": sp.pi})
        substitutions: dict[Any, float] = {}
        missing: list[str] = []
        for symbol in expr.free_symbols:
            name = str(symbol)
            if name in {"Pi", "pi"}:
                substitutions[symbol] = math.pi
                continue
            numeric = _coerce_numeric_param(param_values.get(name))
            if numeric is None:
                missing.append(name)
                continue
            substitutions[symbol] = numeric
        if missing:
            raise ValueError(f"Dump expression '{expression}' misses parameters: {sorted(missing)}")
        return float(expr.evalf(subs=substitutions))
    name = str(var.get("name", ""))
    return param_values.get(name, "MISSING_VALUE")


def write_json_input(
    path: str,
    *,
    actions: list[Mapping[str, Any]],
    param_values: Mapping[str, Any],
) -> dict[str, Any]:
    """Write a JSON input file using Dump actions (V1-compatible expressions)."""
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    else:
        data: dict[str, Any] = {}

    for action in actions:
        if str(action.get("type", "")).strip() != "Dump":
            continue
        for var in action.get("variables", []) or []:
            if not isinstance(var, Mapping):
                continue
            name = str(var.get("name", ""))
            value = _evaluate_dump_value(var, param_values)
            data[name] = _to_json_compatible(value)

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=4)
    return {}


def read_json_output(
    path: str,
    *,
    variables: list[Mapping[str, Any]],
) -> dict[str, Any]:
    """Read observables from a JSON output file."""
    with open(path, "r", encoding="utf-8") as handle:
        content = json.load(handle)
    observables: dict[str, Any] = {}
    for var in variables:
        if not isinstance(var, Mapping) or "name" not in var:
            continue
        name = str(var["name"])
        entry = str(var.get("entry", name))
        current: Any = content
        for part in entry.split("."):
            if not isinstance(current, Mapping) or part not in current:
                current = None
                break
            current = current[part]
        observables[name] = current
    return observables


__all__ = ["read_json_output", "write_json_input"]