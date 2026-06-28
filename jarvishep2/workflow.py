#!/usr/bin/env python3
"""Execution-plan helpers for Worker-side workflow dispatch."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.sample import ExecutionStep


def build_opera_execution_plan(
    opera_modules: Sequence[Mapping[str, Any]],
    *,
    include_likelihood: bool = True,
) -> list[ExecutionStep]:
    """Build a simple layer-0 opera plan followed by an optional likelihood step."""
    steps: list[ExecutionStep] = []
    for index, module in enumerate(opera_modules):
        name = str(module.get("name", f"Operas{index}"))
        steps.append(ExecutionStep(type="opera", name=name, layer=0))
    if include_likelihood:
        steps.append(ExecutionStep(type="likelihood", name="LogLikelihood", layer=1))
    return steps


def group_by_layer(steps: Sequence[ExecutionStep]) -> list[list[ExecutionStep]]:
    """Group execution steps by ascending layer index."""
    buckets: dict[int, list[ExecutionStep]] = defaultdict(list)
    for step in steps:
        buckets[int(step.layer)].append(step)
    return [buckets[layer] for layer in sorted(buckets)]


def execution_plan_template(
    opera_modules: Sequence[Mapping[str, Any]],
    *,
    include_likelihood: bool = True,
) -> list[dict[str, Any]]:
    """JSON-serializable execution plan for Redis transport."""
    return [
        step.to_dict()
        for step in build_opera_execution_plan(
            opera_modules,
            include_likelihood=include_likelihood,
        )
    ]


__all__ = [
    "build_opera_execution_plan",
    "execution_plan_template",
    "group_by_layer",
]