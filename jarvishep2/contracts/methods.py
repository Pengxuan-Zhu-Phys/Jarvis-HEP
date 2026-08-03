#!/usr/bin/env python3
"""Per-Method Sampling contracts (required keys + Bounds dispatch)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.contracts.common import (
    RETIRED_SAMPLING_TOP_LEVEL,
    sampling_bounds,
    try_float,
    try_int,
)
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

_MCMC_METHODS: frozenset[str] = frozenset(
    {
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
    }
)


def method_needs_variables(method: str) -> bool:
    return str(method).strip() in _METHODS_NEED_VARIABLES


def _reject_retired_sampling_top_level(
    sampling: Mapping[str, Any],
) -> list[ValidationIssue]:
    """Old design put sampler knobs under Sampling; V2 requires Bounds only."""
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    for key, destination in RETIRED_SAMPLING_TOP_LEVEL.items():
        if key not in sampling:
            continue
        issues.append(
            issue(
                "error",
                "JV2-MTH-001",
                f"Sampling.{key}",
                f"Sampling.{key} is not a V2 setting; place it under {destination}",
                hint=(
                    "Method-specific sampler settings belong exclusively under "
                    "Sampling.Bounds.  Top-level Sampling holds Method, Variables, "
                    "selection, LogLikelihood, Nuisance, FeedbackReturn, and Bounds."
                ),
                suggestion=f"Move this value to {destination} and remove Sampling.{key}.",
            )
        )
    return issues


def validate_method_sampling(
    sampling: Mapping[str, Any],
    *,
    method: str,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    method = str(method).strip()

    issues.extend(_reject_retired_sampling_top_level(sampling))

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

    bounds = sampling.get("Bounds")
    bounds_map = sampling_bounds(sampling) if bounds is not None else None

    if bounds is not None and not isinstance(bounds, Mapping):
        issues.append(
            issue(
                "error",
                "JV2-BND-000",
                "Sampling.Bounds",
                f"expected a mapping, got {type(bounds).__name__}",
            )
        )
        return issues

    if method == "Bridson":
        if bounds is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-010",
                    "Sampling.Bounds",
                    "Bridson requires Sampling.Bounds with Radius and MaxAttempt",
                )
            )
        else:
            assert bounds_map is not None
            if "Radius" not in bounds_map:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-010",
                        "Sampling.Bounds.Radius",
                        "Bridson requires Sampling.Bounds.Radius (positive number)",
                    )
                )
            else:
                radius = try_float(bounds_map.get("Radius"))
                if radius is None or radius <= 0:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-011",
                            "Sampling.Bounds.Radius",
                            f"expected positive number, got {bounds_map.get('Radius')!r}",
                        )
                    )
            if "MaxAttempt" not in bounds_map:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-012",
                        "Sampling.Bounds.MaxAttempt",
                        "Bridson requires Sampling.Bounds.MaxAttempt (integer ≥ 1)",
                    )
                )
            else:
                k = try_int(bounds_map.get("MaxAttempt"))
                if k is None or k < 1:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-013",
                            "Sampling.Bounds.MaxAttempt",
                            f"expected integer ≥ 1, got {bounds_map.get('MaxAttempt')!r}",
                        )
                    )

    if method == "Random":
        if bounds is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-020",
                    "Sampling.Bounds",
                    "Random requires Sampling.Bounds with 'Point number' (or point_number)",
                )
            )
        else:
            assert bounds_map is not None
            point_number = bounds_map.get(
                "Point number", bounds_map.get("point_number")
            )
            if point_number is None:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-020",
                        "Sampling.Bounds['Point number']",
                        "Random requires Sampling.Bounds['Point number'] (or point_number)",
                    )
                )
            else:
                n = try_int(point_number)
                if n is None or n < 1:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-021",
                            "Sampling.Bounds['Point number']",
                            f"expected integer ≥ 1, got {point_number!r}",
                        )
                    )

    if method == "CSV":
        if bounds is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-030",
                    "Sampling.Bounds",
                    "CSV sampler requires Sampling.Bounds mapping with path",
                )
            )
        else:
            assert bounds_map is not None
            raw_path = str(bounds_map.get("path") or "").strip()
            if not raw_path:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-031",
                        "Sampling.Bounds.path",
                        "Sampling.Bounds.path is required",
                    )
                )

    if method == "AdaptiveBridson":
        if bounds is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-040",
                    "Sampling.Bounds",
                    "AdaptiveBridson requires Sampling.Bounds with "
                    "target_expression and target_value",
                )
            )
        else:
            assert bounds_map is not None
            expr = str(bounds_map.get("target_expression") or "").strip()
            if not expr:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-041",
                        "Sampling.Bounds.target_expression",
                        "AdaptiveBridson requires Sampling.Bounds.target_expression",
                    )
                )
            if "target_value" not in bounds_map:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-042",
                        "Sampling.Bounds.target_value",
                        "AdaptiveBridson requires Sampling.Bounds.target_value",
                    )
                )

    if method in _NESTED_METHODS:
        if bounds is None:
            # Bounds optional (defaults apply) — no error.
            pass
        else:
            assert bounds_map is not None
            issues.extend(
                validate_nested_bounds(
                    bounds_map,
                    method=method,
                )
            )

    if method in _MCMC_METHODS and bounds is not None and not isinstance(bounds, Mapping):
        # Already reported above.
        pass

    return issues


__all__ = [
    "method_needs_variables",
    "validate_method_sampling",
]
