#!/usr/bin/env python3
"""Execution-plan helpers for Worker-side workflow dispatch."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.sample import ExecutionStep


def build_execution_plan(
    *,
    calculator_modules: Sequence[Mapping[str, Any]] | None = None,
    opera_modules: Sequence[Mapping[str, Any]] | None = None,
    include_likelihood: bool = True,
) -> list[ExecutionStep]:
    """Build a layered plan: calculators → operas → likelihood."""
    steps: list[ExecutionStep] = []
    layer = 0

    for index, module in enumerate(calculator_modules or []):
        name = str(module.get("name", f"Calculator{index}"))
        steps.append(ExecutionStep(type="calculator", name=name, layer=layer))

    calc_count = len(list(calculator_modules or []))
    opera_count = len(list(opera_modules or []))
    if calc_count:
        layer += 1

    for index, module in enumerate(opera_modules or []):
        name = str(module.get("name", f"Operas{index}"))
        steps.append(ExecutionStep(type="opera", name=name, layer=layer))

    if opera_count and (calc_count or opera_count):
        layer += 1
    elif calc_count and include_likelihood:
        layer += 1

    if include_likelihood and (calc_count or opera_count):
        steps.append(ExecutionStep(type="likelihood", name="LogLikelihood", layer=layer))
    return steps


def build_opera_execution_plan(
    opera_modules: Sequence[Mapping[str, Any]],
    *,
    include_likelihood: bool = True,
) -> list[ExecutionStep]:
    """Build a simple layer-0 opera plan followed by an optional likelihood step."""
    return build_execution_plan(
        opera_modules=opera_modules,
        include_likelihood=include_likelihood,
    )


def group_by_layer(steps: Sequence[ExecutionStep]) -> list[list[ExecutionStep]]:
    """Group execution steps by ascending layer index."""
    buckets: dict[int, list[ExecutionStep]] = defaultdict(list)
    for step in steps:
        buckets[int(step.layer)].append(step)
    return [buckets[layer] for layer in sorted(buckets)]


def execution_plan_template(
    opera_modules: Sequence[Mapping[str, Any]] | None = None,
    *,
    calculator_modules: Sequence[Mapping[str, Any]] | None = None,
    include_likelihood: bool = True,
) -> list[dict[str, Any]]:
    """JSON-serializable execution plan for Redis transport."""
    return [
        step.to_dict()
        for step in build_execution_plan(
            calculator_modules=calculator_modules,
            opera_modules=opera_modules,
            include_likelihood=include_likelihood,
        )
    ]


__all__ = [
    "build_execution_plan",
    "build_opera_execution_plan",
    "execution_plan_template",
    "group_by_layer",
]