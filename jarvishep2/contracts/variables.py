#!/usr/bin/env python3
"""Sampling.Variables distribution contracts (D14 L1)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.contracts.common import try_float, unknown_keys

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

DISTRIBUTION_TYPES: frozenset[str] = frozenset(
    {"Flat", "Log", "Logit", "Normal", "Log-Normal"}
)

# Per-type allowed parameter keys (closed).
_PARAMS_ALLOWED: dict[str, frozenset[str]] = {
    "Flat": frozenset({"min", "max", "length", "num"}),
    "Log": frozenset({"min", "max", "length", "num"}),
    "Logit": frozenset({"location", "scale"}),
    "Normal": frozenset({"mean", "stddev"}),
    "Log-Normal": frozenset({"mean", "stddev"}),
}

_PARAMS_REQUIRED: dict[str, frozenset[str]] = {
    "Flat": frozenset({"min", "max"}),
    "Log": frozenset({"min", "max"}),
    "Logit": frozenset({"location", "scale"}),
    "Normal": frozenset({"mean", "stddev"}),
    "Log-Normal": frozenset({"mean", "stddev"}),
}

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
        if require_length and dist_type in {"Flat", "Log"}:
            required.add("length")
        if require_num and dist_type in {"Flat", "Log"}:
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
                if length is None or length <= 0:
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
                try:
                    num_i = int(num)  # type: ignore[arg-type]
                except (TypeError, ValueError):
                    num_i = -1
                if isinstance(num, bool) or num_i < 1:
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

        valid_count += 1

    if require_nonempty and valid_count == 0 and raw_vars:
        # All entries invalid — already reported per entry; keep a summary only if none.
        pass

    _ = method  # reserved for method-specific variable rules
    return issues


__all__ = ["DISTRIBUTION_TYPES", "validate_variables"]
