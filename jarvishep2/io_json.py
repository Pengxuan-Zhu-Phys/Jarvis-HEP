#!/usr/bin/env python3
"""Minimal JSON IO helpers for CalculatorModule (parity-project scope)."""

from __future__ import annotations

import json
import os
from collections.abc import Mapping
from typing import Any

import numpy as np


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


def write_json_input(
    path: str,
    *,
    actions: list[Mapping[str, Any]],
    param_values: Mapping[str, Any],
) -> dict[str, Any]:
    """Write a JSON input file using Dump actions (EchoCalc parity shape)."""
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
            value = param_values.get(name, "MISSING_VALUE")
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