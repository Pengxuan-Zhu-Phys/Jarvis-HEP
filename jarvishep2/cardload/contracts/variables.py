#!/usr/bin/env python3
"""Sampling.Variables distribution contracts (D14 L1)."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.cardload.contracts.common import try_float, unknown_keys

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

DISTRIBUTION_TYPES: frozenset[str] = frozenset(
    {
        "Flat", "Log", "Logit", "Normal", "Log-Normal", "Binomial", "Poisson",
        "Beta", "Exponential", "Gamma",
    }
)

# ``num`` and ``length`` describe the sampling method's unit-cube geometry,
# not the variable probability law.  They are therefore accepted for every
# distribution and required by Grid / Bridson respectively.
_METHOD_PARAMETER_KEYS = frozenset({"length", "num"})

# Per-type allowed parameter keys (closed).
# Public aliases (D24): man / tools import PARAMS_* without touching private names.
_PARAMS_ALLOWED: dict[str, frozenset[str]] = {
    "Flat": frozenset({"min", "max"}) | _METHOD_PARAMETER_KEYS,
    "Log": frozenset({"min", "max"}) | _METHOD_PARAMETER_KEYS,
    "Logit": frozenset({"location", "scale"}) | _METHOD_PARAMETER_KEYS,
    "Normal": frozenset({"mean", "stddev"}) | _METHOD_PARAMETER_KEYS,
    "Log-Normal": frozenset({"mean", "stddev"}) | _METHOD_PARAMETER_KEYS,
    "Binomial": frozenset({"n", "p"}) | _METHOD_PARAMETER_KEYS,
    "Poisson": frozenset({"lambda"}) | _METHOD_PARAMETER_KEYS,
    "Beta": frozenset({"alpha", "beta"}) | _METHOD_PARAMETER_KEYS,
    "Exponential": frozenset({"rate"}) | _METHOD_PARAMETER_KEYS,
    "Gamma": frozenset({"shape", "scale"}) | _METHOD_PARAMETER_KEYS,
}

_PARAMS_REQUIRED: dict[str, frozenset[str]] = {
    "Flat": frozenset({"min", "max"}),
    "Log": frozenset({"min", "max"}),
    "Logit": frozenset({"location", "scale"}),
    "Normal": frozenset({"mean", "stddev"}),
    "Log-Normal": frozenset({"mean", "stddev"}),
    "Binomial": frozenset({"n", "p"}),
    "Poisson": frozenset({"lambda"}),
    "Beta": frozenset({"alpha", "beta"}),
    "Exponential": frozenset({"rate"}),
    "Gamma": frozenset({"shape", "scale"}),
}

# D24 public surface for ``Jarvis man`` (same objects; do not mutate).
PARAMS_ALLOWED: dict[str, frozenset[str]] = _PARAMS_ALLOWED
PARAMS_REQUIRED: dict[str, frozenset[str]] = _PARAMS_REQUIRED

_VARIABLE_ITEM_KEYS = frozenset({"name", "description", "distribution"})
_DISTRIBUTION_KEYS = frozenset({"type", "parameters"})


def validate_variables(
    sampling: Mapping[str, Any],
    *,
    method: str | None,
    require_nonempty: bool,
    require_length: bool = False,
    require_num: bool = False,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    raw_vars = sampling.get("Variables")

    if raw_vars is None:
        if require_nonempty:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-001",
                    "Sampling.Variables",
                    "Sampling.Variables is required for this Method and must be a non-empty list",
                )
            )
        return issues

    if not isinstance(raw_vars, list):
        issues.append(
            issue(
                "error",
                "JV2-VAR-001",
                "Sampling.Variables",
                f"expected a list, got {type(raw_vars).__name__}",
            )
        )
        return issues

    if require_nonempty and len(raw_vars) == 0:
        issues.append(
            issue(
                "error",
                "JV2-VAR-001",
                "Sampling.Variables",
                "must define at least one variable",
            )
        )
        return issues

    valid_count = 0
    for index, item in enumerate(raw_vars):
        path = f"Sampling.Variables[{index}]"
        if not isinstance(item, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-010",
                    path,
                    f"expected a mapping, got {type(item).__name__}",
                )
            )
            continue

        extra = unknown_keys(item, _VARIABLE_ITEM_KEYS)
        if extra:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-011",
                    path,
                    f"unknown key(s): {', '.join(extra)}; "
                    f"allowed: {', '.join(sorted(_VARIABLE_ITEM_KEYS))}",
                )
            )

        name = item.get("name")
        if not isinstance(name, str) or not name.strip():
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-012",
                    f"{path}.name",
                    "name is required and must be a non-empty string "
                    "(nameless entries are not silently dropped)",
                )
            )
            continue

        dist = item.get("distribution")
        if dist is None:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-020",
                    f"{path}.distribution",
                    "distribution is required",
                )
            )
            continue
        if not isinstance(dist, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-020",
                    f"{path}.distribution",
                    f"expected a mapping, got {type(dist).__name__}",
                )
            )
            continue

        dist_extra = unknown_keys(dist, _DISTRIBUTION_KEYS)
        if dist_extra:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-021",
                    f"{path}.distribution",
                    f"unknown key(s): {', '.join(dist_extra)}; "
                    f"allowed: {', '.join(sorted(_DISTRIBUTION_KEYS))}",
                )
            )

        dist_type_raw = dist.get("type")
        if dist_type_raw is None:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-002",
                    f"{path}.distribution.type",
                    "type is required; "
                    f"expected one of: {', '.join(sorted(DISTRIBUTION_TYPES))}",
                )
            )
            continue
        dist_type = str(dist_type_raw).strip()
        if dist_type not in DISTRIBUTION_TYPES:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-002",
                    f"{path}.distribution.type",
                    f"{dist_type!r} is not supported; "
                    f"expected one of: {', '.join(sorted(DISTRIBUTION_TYPES))} "
                    f"(names are case-sensitive)",
                )
            )
            continue

        params = dist.get("parameters")
        if params is None:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-030",
                    f"{path}.distribution.parameters",
                    "parameters is required",
                )
            )
            continue
        if not isinstance(params, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-030",
                    f"{path}.distribution.parameters",
                    f"expected a mapping, got {type(params).__name__}",
                )
            )
            continue

        allowed = _PARAMS_ALLOWED[dist_type]
        required = set(_PARAMS_REQUIRED[dist_type])
        if require_length:
            required.add("length")
        if require_num:
            required.add("num")

        missing = sorted(required - set(params.keys()))
        if missing:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-031",
                    f"{path}.distribution.parameters",
                    f"missing required key(s) for type {dist_type}: {', '.join(missing)}",
                )
            )

        extra_params = unknown_keys(params, allowed)
        if extra_params:
            issues.append(
                issue(
                    "error",
                    "JV2-VAR-032",
                    f"{path}.distribution.parameters",
                    f"unknown key(s) for type {dist_type}: {', '.join(extra_params)}; "
                    f"allowed: {', '.join(sorted(allowed))}",
                )
            )

        # Numeric sanity for range distributions.
        if dist_type in {"Flat", "Log"}:
            lo = try_float(params.get("min")) if "min" in params else None
            hi = try_float(params.get("max")) if "max" in params else None
            if "min" in params and lo is None:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-033",
                        f"{path}.distribution.parameters.min",
                        f"expected a number, got {params.get('min')!r}",
                    )
                )
            if "max" in params and hi is None:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-033",
                        f"{path}.distribution.parameters.max",
                        f"expected a number, got {params.get('max')!r}",
                    )
                )
            if lo is not None and hi is not None and not (lo < hi):
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-034",
                        f"{path}.distribution.parameters",
                        f"require min < max, got min={lo}, max={hi}",
                    )
                )

        if require_length and "length" in params:
            length = try_float(params.get("length"))
            if length is None or not math.isfinite(length) or length <= 0:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-035",
                        f"{path}.distribution.parameters.length",
                        f"Bridson requires positive length, got {params.get('length')!r}",
                    )
                )
        if require_num and "num" in params:
            num = params.get("num")
            num_value = try_float(num)
            if (
                num_value is None
                or not math.isfinite(num_value)
                or not num_value.is_integer()
                or num_value < 1
            ):
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-036",
                        f"{path}.distribution.parameters.num",
                        f"Grid requires integer num ≥ 1, got {num!r}",
                    )
                )

        if dist_type in {"Normal", "Log-Normal"} and "stddev" in params:
            std = try_float(params.get("stddev"))
            if std is None or std <= 0:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-037",
                        f"{path}.distribution.parameters.stddev",
                        f"expected positive number, got {params.get('stddev')!r}",
                    )
                )

        if dist_type == "Binomial":
            n_raw = params.get("n")
            n_value = try_float(n_raw) if "n" in params else None
            if (
                n_value is None
                or not math.isfinite(n_value)
                or not n_value.is_integer()
                or n_value < 0
            ):
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-038",
                        f"{path}.distribution.parameters.n",
                        f"expected a non-negative integer, got {n_raw!r}",
                    )
                )
            p_raw = params.get("p")
            p_value = try_float(p_raw) if "p" in params else None
            if p_value is None or not math.isfinite(p_value) or not 0.0 <= p_value <= 1.0:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-039",
                        f"{path}.distribution.parameters.p",
                        f"expected a probability in [0, 1], got {p_raw!r}",
                    )
                )

        if dist_type == "Poisson":
            value = try_float(params.get("lambda")) if "lambda" in params else None
            if value is None or not math.isfinite(value) or value < 0:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-040",
                        f"{path}.distribution.parameters.lambda",
                        f"expected a non-negative number, got {params.get('lambda')!r}",
                    )
                )

        positive_pairs = {
            "Beta": ("alpha", "beta"),
            "Exponential": ("rate",),
            "Gamma": ("shape", "scale"),
        }
        for parameter in positive_pairs.get(dist_type, ()):
            value = try_float(params.get(parameter)) if parameter in params else None
            if value is None or not math.isfinite(value) or value <= 0:
                issues.append(
                    issue(
                        "error",
                        "JV2-VAR-041",
                        f"{path}.distribution.parameters.{parameter}",
                        f"expected a positive number, got {params.get(parameter)!r}",
                    )
                )

        valid_count += 1

    if require_nonempty and valid_count == 0 and raw_vars:
        # All entries invalid — already reported per entry; keep a summary only if none.
        pass

    _ = method  # reserved for method-specific variable rules
    return issues


__all__ = ["DISTRIBUTION_TYPES", "validate_variables"]
