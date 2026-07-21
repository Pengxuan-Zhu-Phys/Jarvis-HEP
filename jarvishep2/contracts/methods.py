#!/usr/bin/env python3
"""Per-Method Sampling contracts (required keys + Bounds dispatch)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.contracts.common import try_float, try_int
from jarvishep2.contracts.nested import validate_nested_bounds
from jarvishep2.contracts.variables import validate_variables

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

# Methods that need Sampling.Variables (CSV is fixed-set from file).
_METHODS_NEED_VARIABLES: frozenset[str] = frozenset(
    {
        "Bridson",
        "Random",
        "Grid",
        "AdaptiveBridson",
        "MCMC",
        "AMMCMC",
        "AM",
        "DRAM",
        "EnsembleMCMC",
        "Ensemble",
        "DEMCMC",
        "PTMCMC",
        "PT",
        "PTEnsemble",
        "Dynesty",
        "MultiNest",
    }
)

_NESTED_METHODS: frozenset[str] = frozenset({"Dynesty", "MultiNest"})


def method_needs_variables(method: str) -> bool:
    return str(method).strip() in _METHODS_NEED_VARIABLES


def validate_method_sampling(
    sampling: Mapping[str, Any],
    *,
    method: str,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    method = str(method).strip()

    require_length = method == "Bridson"
    require_num = method == "Grid"
    require_vars = method_needs_variables(method)

    if require_vars or "Variables" in sampling:
        issues.extend(
            validate_variables(
                sampling,
                method=method,
                require_nonempty=require_vars,
                require_length=require_length,
                require_num=require_num,
            )
        )

    if method == "Bridson":
        if "Radius" not in sampling:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-010",
                    "Sampling.Radius",
                    "Bridson requires Sampling.Radius (positive number)",
                )
            )
        else:
            radius = try_float(sampling.get("Radius"))
            if radius is None or radius <= 0:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-011",
                        "Sampling.Radius",
                        f"expected positive number, got {sampling.get('Radius')!r}",
                    )
                )
        if "MaxAttempt" not in sampling:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-012",
                    "Sampling.MaxAttempt",
                    "Bridson requires Sampling.MaxAttempt (integer ≥ 1)",
                )
            )
        else:
            k = try_int(sampling.get("MaxAttempt"))
            if k is None or k < 1:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-013",
                        "Sampling.MaxAttempt",
                        f"expected integer ≥ 1, got {sampling.get('MaxAttempt')!r}",
                    )
                )

    if method == "Random":
        point_number = sampling.get("Point number", sampling.get("point_number"))
        if point_number is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-020",
                    "Sampling['Point number']",
                    "Random requires Sampling['Point number'] (or point_number)",
                )
            )
        else:
            n = try_int(point_number)
            if n is None or n < 1:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-021",
                        "Sampling['Point number']",
                        f"expected integer ≥ 1, got {point_number!r}",
                    )
                )

    if method == "CSV":
        csv_cfg = sampling.get("CSV")
        if not isinstance(csv_cfg, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-030",
                    "Sampling.CSV",
                    "CSV sampler requires Sampling.CSV mapping with path",
                )
            )
        else:
            raw_path = str(csv_cfg.get("path") or "").strip()
            if not raw_path:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-031",
                        "Sampling.CSV.path",
                        "Sampling.CSV.path is required",
                    )
                )

    if method in _NESTED_METHODS:
        bounds = sampling.get("Bounds")
        if bounds is None:
            # Bounds optional (defaults apply) — no error.
            pass
        elif not isinstance(bounds, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-BND-000",
                    "Sampling.Bounds",
                    f"expected a mapping, got {type(bounds).__name__}",
                )
            )
        else:
            issues.extend(
                validate_nested_bounds(
                    bounds,
                    method=method,
                )
            )

    return issues


__all__ = [
    "method_needs_variables",
    "validate_method_sampling",
]
