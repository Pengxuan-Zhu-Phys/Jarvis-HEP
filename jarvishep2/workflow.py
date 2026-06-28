#!/usr/bin/env python3
"""Execution-plan helpers for Worker-side workflow dispatch."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.sample import ExecutionStep


def _normalize_required_modules(value: Any) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value.strip()] if value.strip() else []
    if isinstance(value, Mapping):
        return []
    return [str(item).strip() for item in value if str(item).strip()]


def resolve_module_layers(modules: Sequence[Mapping[str, Any]] | None) -> dict[str, int]:
    """Assign layer indices from ``required_modules`` (same layer = concurrent-eligible)."""
    names: list[str] = []
    deps: dict[str, set[str]] = {}
    for index, module in enumerate(modules or []):
        if not isinstance(module, Mapping):
            continue
        name = str(module.get("name", f"Module{index}")).strip()
        if not name:
            continue
        names.append(name)
        deps[name] = set(_normalize_required_modules(module.get("required_modules")))

    layers: dict[str, int] = {}
    resolved: set[str] = set()
    remaining = set(names)
    layer_idx = 0
    while remaining:
        ready = sorted(
            name for name in remaining if deps.get(name, set()).issubset(resolved)
        )
        if not ready:
            ready = sorted(remaining)
        for name in ready:
            layers[name] = layer_idx
            resolved.add(name)
            remaining.discard(name)
        layer_idx += 1
    return layers


def max_layer_width(layers: Mapping[str, int]) -> int:
    """Return the widest layer in a name→layer map."""
    if not layers:
        return 1
    counts: dict[int, int] = defaultdict(int)
    for layer in layers.values():
        counts[int(layer)] += 1
    return max(counts.values())


def concurrency_groups(steps: Sequence[ExecutionStep]) -> list[list[str]]:
    """Module names grouped by execution layer."""
    grouped = group_by_layer(steps)
    return [[step.name for step in layer] for layer in grouped]


def build_execution_plan(
    *,
    calculator_modules: Sequence[Mapping[str, Any]] | None = None,
    opera_modules: Sequence[Mapping[str, Any]] | None = None,
    include_likelihood: bool = True,
) -> list[ExecutionStep]:
    """Build a layered plan: calculators → operas → likelihood."""
    steps: list[ExecutionStep] = []
    calc_layers = resolve_module_layers(calculator_modules)
    opera_layers = resolve_module_layers(opera_modules)

    calc_count = len(calc_layers)
    opera_count = len(opera_layers)
    max_calc_layer = max(calc_layers.values()) if calc_layers else -1
    opera_base = max_calc_layer + 1 if calc_count else 0

    for index, module in enumerate(calculator_modules or []):
        name = str(module.get("name", f"Calculator{index}"))
        steps.append(
            ExecutionStep(type="calculator", name=name, layer=calc_layers.get(name, 0))
        )

    for index, module in enumerate(opera_modules or []):
        name = str(module.get("name", f"Operas{index}"))
        steps.append(
            ExecutionStep(
                type="opera",
                name=name,
                layer=opera_base + opera_layers.get(name, 0),
            )
        )

    likelihood_layer = opera_base
    if opera_count:
        likelihood_layer = opera_base + max(opera_layers.values()) + 1
    elif calc_count:
        likelihood_layer = max_calc_layer + 1
        if include_likelihood:
            likelihood_layer += 1

    if include_likelihood and (calc_count or opera_count):
        steps.append(
            ExecutionStep(type="likelihood", name="LogLikelihood", layer=likelihood_layer)
        )
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
    "concurrency_groups",
    "execution_plan_template",
    "group_by_layer",
    "max_layer_width",
    "resolve_module_layers",
]