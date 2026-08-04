#!/usr/bin/env python3
"""Sampling.Mapper load-time contracts (D22 / JV2-MAP-*)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.mapper import (
    STATISTICAL_METHODS,
    MapperError,
    _is_affine_in_sample_vars,
    build_mapper_spec_from_config,
)

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

_MAPPER_EXAMPLE = """Sampling:
  Mapper:
    x: "cos(t)"
    y: "sin(t)" """


def _opera_output_names(config: Mapping[str, Any]) -> set[str]:
    names: set[str] = set()
    operas = config.get("Operas") if isinstance(config.get("Operas"), Mapping) else {}
    modules = operas.get("Modules") if isinstance(operas, Mapping) else None
    if not isinstance(modules, list):
        return names
    for mod in modules:
        if not isinstance(mod, Mapping):
            continue
        for out in mod.get("output") or []:
            if isinstance(out, Mapping):
                name = str(out.get("name") or "").strip()
                if name:
                    names.add(name)
            elif isinstance(out, str) and out.strip():
                names.add(out.strip())
    return names


def _calculator_dump_names(config: Mapping[str, Any]) -> set[str]:
    """Best-effort dump variable names from calculator execution.input actions."""
    names: set[str] = set()
    calcs = config.get("Calculators") if isinstance(config.get("Calculators"), Mapping) else {}
    modules = calcs.get("Modules") if isinstance(calcs, Mapping) else None
    if not isinstance(modules, list):
        return names
    for mod in modules:
        if not isinstance(mod, Mapping):
            continue
        execution = mod.get("execution") if isinstance(mod.get("execution"), Mapping) else {}
        for io in list(execution.get("input") or []) + list(execution.get("output") or []):
            if not isinstance(io, Mapping):
                continue
            for action in io.get("actions") or []:
                if not isinstance(action, Mapping):
                    continue
                for var in action.get("variables") or []:
                    if isinstance(var, Mapping):
                        name = str(var.get("name") or "").strip()
                        if name:
                            names.add(name)
            for var in io.get("variables") or []:
                if isinstance(var, Mapping):
                    name = str(var.get("name") or "").strip()
                    if name:
                        names.add(name)
    return names


def validate_mapper(
    config: Mapping[str, Any],
    *,
    method: str | None = None,
) -> list[ValidationIssue]:
    """Validate optional ``Sampling.Mapper`` (pure; no Redis / processes)."""
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    sampling = config.get("Sampling")
    if not isinstance(sampling, Mapping):
        return issues
    if "Mapper" not in sampling:
        return issues

    method_name = (
        str(method).strip()
        if method is not None
        else str(sampling.get("Method") or sampling.get("method") or "").strip()
    )

    try:
        spec = build_mapper_spec_from_config(config)
    except MapperError as exc:
        issues.append(
            issue(
                "error",
                exc.code,
                exc.path,
                exc.message,
                suggestion=_suggestion_for(exc.code),
                example=_MAPPER_EXAMPLE if exc.code.startswith("JV2-MAP-00") else None,
            )
        )
        return issues

    n_vars = len(spec.variables)
    n_mapped = len(spec.derive_order)
    if n_mapped == 0:
        return issues

    sample_names = [v.name for v in spec.variables]
    is_stat = method_name in STATISTICAL_METHODS

    if is_stat and n_mapped > n_vars:
        issues.append(
            issue(
                "warning",
                "JV2-MAP-050",
                "Sampling.Mapper",
                (
                    f"dimension expansion: {n_mapped} Mapper parameter(s) from "
                    f"{n_vars} sampling variable(s) with Method {method_name!r}; "
                    "samples lie on a lower-dimensional manifold and the prior is "
                    "defined on Sampling.Variables — supply a Jacobian in "
                    "LogLikelihood if inferring in the mapped space"
                ),
                suggestion=(
                    "Use coverage Methods (Bridson/Random/Grid) for manifold scans, "
                    "or add the Jacobian correction to LogLikelihood yourself."
                ),
            )
        )
    elif is_stat and n_mapped > 0:
        nonlinear = [
            name
            for name in spec.derive_order
            if not _is_affine_in_sample_vars(spec.derive_exprs[name], sample_names)
        ]
        if nonlinear:
            issues.append(
                issue(
                    "warning",
                    "JV2-MAP-051",
                    "Sampling.Mapper",
                    (
                        f"nonlinear reparameterization via Mapper {nonlinear} with "
                        f"Method {method_name!r}; prior remains on Sampling.Variables"
                    ),
                    suggestion=(
                        "If the mapped parameters are the inference coordinates, "
                        "include the absolute Jacobian of the transform in LogLikelihood."
                    ),
                )
            )

    collide = set(spec.derive_order) & (
        _opera_output_names(config) | _calculator_dump_names(config)
    )
    if collide:
        issues.append(
            issue(
                "warning",
                "JV2-MAP-052",
                "Sampling.Mapper",
                (
                    f"Mapper name(s) {sorted(collide)} also appear as Opera output / "
                    "calculator dump names; merge_observables will overwrite the "
                    "observable while Sample.params keeps the mapper value"
                ),
                suggestion="Rename the Mapper parameter or the Opera/calculator output.",
            )
        )

    return issues


def _suggestion_for(code: str) -> str:
    return {
        "JV2-MAP-001": (
            "Use Sampling.Mapper as a mapping of name → expression string, "
            'for example: Mapper: {x: "cos(t)", y: "sin(t)"}.'
        ),
        "JV2-MAP-002": (
            "Reference only Sampling.Variables names and earlier Mapper names "
            "(observables / LogL / uuid are not visible to Mapper)."
        ),
        "JV2-MAP-003": "Rename the Mapper parameter so it does not shadow a Variable.",
        "JV2-MAP-004": "Break the cycle among Mapper expressions.",
        "JV2-MAP-005": "Avoid reserved constants Pi/pi/PI/E/Inf as Mapper names.",
        "JV2-MAP-006": "Avoid reserved DATABASE columns uuid/sample_index/status/LogL.",
        "JV2-MAP-007": "Fix the sympy expression syntax for this Mapper entry.",
        "JV2-MAP-010": (
            "Remove Sampling.Mapper for Method CSV (v1); parameters come from CSV columns."
        ),
    }.get(code, "Fix Sampling.Mapper according to the error message.")


__all__ = ["validate_mapper"]
