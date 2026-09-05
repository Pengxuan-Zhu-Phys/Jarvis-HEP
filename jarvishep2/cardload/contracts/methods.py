#!/usr/bin/env python3
"""Per-Method Sampling contracts (required keys + Bounds dispatch)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.cardload.contracts.common import (
    RETIRED_SAMPLING_TOP_LEVEL,
    sampling_bounds,
    try_float,
    try_int,
)
from jarvishep2.cardload.contracts.nested import validate_nested_bounds
from jarvishep2.cardload.contracts.variables import validate_variables
from jarvishep2.sampler_catalog import (
    mcmc_names,
    methods_need_variables,
    nested_names,
)

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

# Views of ``sampler_catalog`` (D25.1). CSV is the only method without Variables.
_METHODS_NEED_VARIABLES: frozenset[str] = methods_need_variables()
_NESTED_METHODS: frozenset[str] = nested_names()
_MCMC_METHODS: frozenset[str] = mcmc_names()


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
                    "Bridson requires Sampling.Bounds with radius and max_attempt",
                )
            )
        else:
            assert bounds_map is not None
            if "radius" not in bounds_map:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-010",
                        "Sampling.Bounds.radius",
                        "Bridson requires Sampling.Bounds.radius (positive number)",
                    )
                )
            else:
                radius = try_float(bounds_map.get("radius"))
                if radius is None or radius <= 0:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-011",
                            "Sampling.Bounds.radius",
                            f"expected positive number, got {bounds_map.get('radius')!r}",
                        )
                    )
            if "max_attempt" not in bounds_map:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-012",
                        "Sampling.Bounds.max_attempt",
                        "Bridson requires Sampling.Bounds.max_attempt (integer ≥ 1)",
                    )
                )
            else:
                k = try_int(bounds_map.get("max_attempt"))
                if k is None or k < 1:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-013",
                            "Sampling.Bounds.max_attempt",
                            f"expected integer ≥ 1, got {bounds_map.get('max_attempt')!r}",
                        )
                    )

    if method == "Random":
        if bounds is None:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-020",
                    "Sampling.Bounds",
                    "Random requires Sampling.Bounds.point_number",
                )
            )
        else:
            assert bounds_map is not None
            point_number = bounds_map.get("point_number")
            if point_number is None:
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-020",
                        "Sampling.Bounds.point_number",
                        "Random requires Sampling.Bounds.point_number",
                    )
                )
            else:
                n = try_int(point_number)
                if n is None or n < 1:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-021",
                            "Sampling.Bounds.point_number",
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

    if method in _MCMC_METHODS:
        issues.extend(_validate_mcmc_bounds(method, bounds, bounds_map))

    return issues


_ENSEMBLE_MCMC = frozenset({"EnsembleMCMC"})
_DE_MCMC = frozenset({"DEMCMC"})
_PT_MCMC = frozenset({"PTMCMC", "PTEnsemble"})
_ADAPTIVE_MCMC = frozenset({"AMMCMC", "DRAM"})

# Per-method Bounds codes. Each YAML Method name has exactly one engine.
_MCMC_CORE = {
    "ToyMCMC": {
        "missing": "JV2-MTH-050",
        "num_chains": "JV2-MTH-051",
        "num_iters": "JV2-MTH-052",
        "proposal_scale": "JV2-MTH-053",
        "proposal_scale_len": "JV2-MTH-055",
        "seed": "JV2-MTH-056",
        "min_chains": 1,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "MCMC": {
        "missing": "JV2-MTH-080",
        "num_chains": "JV2-MTH-081",
        "num_iters": "JV2-MTH-082",
        "proposal_scale": "JV2-MTH-083",
        "proposal_scale_len": "JV2-MTH-085",
        "seed": "JV2-MTH-086",
        "min_chains": 1,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "AMMCMC": {
        "missing": "JV2-MTH-090",
        "num_chains": "JV2-MTH-091",
        "num_iters": "JV2-MTH-092",
        "proposal_scale": "JV2-MTH-093",
        "proposal_scale_len": "JV2-MTH-094",
        "seed": "JV2-MTH-095",
        "min_chains": 1,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "DRAM": {
        "missing": "JV2-MTH-100",
        "num_chains": "JV2-MTH-101",
        "num_iters": "JV2-MTH-102",
        "proposal_scale": "JV2-MTH-103",
        "proposal_scale_len": "JV2-MTH-104",
        "seed": "JV2-MTH-105",
        "min_chains": 1,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "EnsembleMCMC": {
        "missing": "JV2-MTH-110",
        "num_chains": "JV2-MTH-111",
        "num_iters": "JV2-MTH-112",
        "proposal_scale": "JV2-MTH-113",
        "proposal_scale_len": "JV2-MTH-114",
        "seed": "JV2-MTH-115",
        "min_chains": 2,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "DEMCMC": {
        "missing": "JV2-MTH-120",
        "num_chains": "JV2-MTH-121",
        "num_iters": "JV2-MTH-122",
        "proposal_scale": "JV2-MTH-123",
        "proposal_scale_len": "JV2-MTH-124",
        "seed": "JV2-MTH-125",
        "min_chains": 2,
        "required": ("num_chains", "num_iters", "proposal_scale"),
    },
    "PTMCMC": {
        "missing": "JV2-MTH-060",
        "num_chains": "JV2-MTH-061",
        "num_iters": "JV2-MTH-062",
        "proposal_scale": "JV2-MTH-063",
        "proposal_scale_len": "JV2-MTH-066",
        "seed": "JV2-MTH-071",
        "min_chains": 2,
        "required": (
            "num_chains",
            "num_iters",
            "proposal_scale",
            "temperature_ladder",
        ),
    },
}
_MCMC_CORE["PTEnsemble"] = _MCMC_CORE["PTMCMC"]

_REQUIRED_KEY_HELP = {
    "num_chains": "a positive integer",
    "num_iters": "a positive integer",
    "proposal_scale": "a positive number or list",
    "temperature_ladder": "a strictly increasing list starting at 1.0",
}


def _validate_mcmc_bounds(
    method: str,
    bounds: Any,
    bounds_map: Mapping[str, Any] | None,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    spec = _MCMC_CORE.get(method)
    if spec is None:
        return []
    label = method
    issues: list[ValidationIssue] = []
    if bounds is None:
        required = ", ".join(spec["required"])
        issues.append(
            issue(
                "error",
                spec["missing"],
                "Sampling.Bounds",
                f"{label} requires Sampling.Bounds with {required}",
            )
        )
        return issues
    if bounds_map is None:
        return issues

    min_chains = int(spec["min_chains"])
    help_text = dict(_REQUIRED_KEY_HELP)
    if min_chains > 1:
        help_text["num_chains"] = f"an integer >= {min_chains}"
    code_for = {
        "num_chains": spec["num_chains"],
        "num_iters": spec["num_iters"],
        "proposal_scale": spec["proposal_scale"],
        "temperature_ladder": "JV2-MTH-064",
    }
    for key in spec["required"]:
        if key not in bounds_map:
            issues.append(
                issue(
                    "error",
                    code_for[key],
                    f"Sampling.Bounds.{key}",
                    f"{label} requires Sampling.Bounds.{key} ({help_text[key]})",
                )
            )

    nchains = try_int(bounds_map.get("num_chains"))
    if nchains is not None and nchains < min_chains:
        issues.append(
            issue(
                "error",
                spec["num_chains"],
                "Sampling.Bounds.num_chains",
                f"{label} requires num_chains >= {min_chains}, "
                f"got {bounds_map.get('num_chains')!r}",
            )
        )

    niters = try_int(bounds_map.get("num_iters"))
    if niters is not None and niters < 1:
        issues.append(
            issue(
                "error",
                spec["num_iters"],
                "Sampling.Bounds.num_iters",
                f"expected a positive integer, got {bounds_map.get('num_iters')!r}",
            )
        )

    raw_scale = bounds_map.get("proposal_scale")
    scales = raw_scale if isinstance(raw_scale, (list, tuple)) else [raw_scale]
    parsed_scales = [try_float(value) for value in scales]
    if raw_scale is not None and (
        not parsed_scales or any(value is None or value <= 0 for value in parsed_scales)
    ):
        issues.append(
            issue(
                "error",
                spec["proposal_scale"],
                "Sampling.Bounds.proposal_scale",
                f"{label} proposal_scale values must all be positive numbers",
            )
        )
    elif (
        nchains is not None
        and isinstance(raw_scale, (list, tuple))
        and len(raw_scale) not in {1, nchains}
    ):
        issues.append(
            issue(
                "error",
                spec["proposal_scale_len"],
                "Sampling.Bounds.proposal_scale",
                f"{label} proposal_scale must have length 1 or num_chains ({nchains}), "
                f"got {len(raw_scale)}",
            )
        )

    raw_seed = bounds_map.get("seed")
    if raw_seed is not None:
        seed = try_int(raw_seed)
        if seed is None or seed < 0:
            issues.append(
                issue(
                    "error",
                    spec["seed"],
                    "Sampling.Bounds.seed",
                    f"{label} seed must be a non-negative integer, got {raw_seed!r}",
                )
            )

    if method in _PT_MCMC:
        issues.extend(_validate_pt_ladder(label, bounds_map, nchains))
    if method in _ADAPTIVE_MCMC:
        issues.extend(_validate_adapt_knobs(label, bounds_map, method))
    if method == "DRAM":
        issues.extend(_validate_dram_knobs(label, bounds_map))
    if method in _ENSEMBLE_MCMC or method == "PTEnsemble":
        issues.extend(_validate_stretch_a(label, bounds_map, method))
    if method in _DE_MCMC:
        issues.extend(_validate_de_knobs(label, bounds_map))
    return issues


def _validate_pt_ladder(
    label: str,
    bounds_map: Mapping[str, Any],
    nchains: int | None,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    raw_ladder = bounds_map.get("temperature_ladder")
    if raw_ladder is not None:
        if not isinstance(raw_ladder, (list, tuple)):
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-064",
                    "Sampling.Bounds.temperature_ladder",
                    f"{label} temperature_ladder must be a list of positive numbers",
                )
            )
        else:
            temps = [try_float(value) for value in raw_ladder]
            if not temps or any(value is None or value <= 0 for value in temps):
                issues.append(
                    issue(
                        "error",
                        "JV2-MTH-064",
                        "Sampling.Bounds.temperature_ladder",
                        f"{label} temperature_ladder values must all be positive",
                    )
                )
            else:
                assert all(value is not None for value in temps)
                if nchains is not None and len(temps) != nchains:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-067",
                            "Sampling.Bounds.temperature_ladder",
                            f"{label} temperature_ladder length must equal "
                            f"num_chains ({nchains}), got {len(temps)}",
                        )
                    )
                if abs(float(temps[0]) - 1.0) > 1e-12:
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-068",
                            "Sampling.Bounds.temperature_ladder",
                            f"{label} temperature_ladder[0] must be 1.0 (cold), "
                            f"got {temps[0]!r}",
                        )
                    )
                if any(not (float(hi) > float(lo)) for lo, hi in zip(temps, temps[1:])):
                    issues.append(
                        issue(
                            "error",
                            "JV2-MTH-069",
                            "Sampling.Bounds.temperature_ladder",
                            f"{label} temperature_ladder must be strictly increasing",
                        )
                    )

    raw_exchange = bounds_map.get("exchange_interval")
    if raw_exchange is not None:
        exchange = try_int(raw_exchange)
        if exchange is None or exchange < 1:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-070",
                    "Sampling.Bounds.exchange_interval",
                    f"{label} exchange_interval must be an integer >= 1, "
                    f"got {raw_exchange!r}",
                )
            )
    return issues


def _validate_adapt_knobs(
    label: str,
    bounds_map: Mapping[str, Any],
    method: str,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    start_code = "JV2-MTH-106" if method == "DRAM" else "JV2-MTH-096"
    eps_code = "JV2-MTH-106" if method == "DRAM" else "JV2-MTH-097"
    scale_code = "JV2-MTH-106" if method == "DRAM" else "JV2-MTH-098"
    issues: list[ValidationIssue] = []
    for key, code in (
        ("adapt_start_iter", start_code),
        ("adapt_window", start_code),
    ):
        raw = bounds_map.get(key)
        if raw is None:
            continue
        value = try_int(raw)
        if value is None or value < 1:
            issues.append(
                issue(
                    "error",
                    code,
                    f"Sampling.Bounds.{key}",
                    f"{label} {key} must be an integer >= 1, got {raw!r}",
                )
            )
    raw_eps = bounds_map.get("adapt_eps")
    if raw_eps is not None:
        eps = try_float(raw_eps)
        if eps is None or eps < 0:
            issues.append(
                issue(
                    "error",
                    eps_code,
                    "Sampling.Bounds.adapt_eps",
                    f"{label} adapt_eps must be >= 0, got {raw_eps!r}",
                )
            )
    raw_scale = bounds_map.get("adapt_scale")
    if raw_scale is not None:
        scale = try_float(raw_scale)
        if scale is None or scale <= 0:
            issues.append(
                issue(
                    "error",
                    scale_code,
                    "Sampling.Bounds.adapt_scale",
                    f"{label} adapt_scale must be a positive number, got {raw_scale!r}",
                )
            )
    return issues


def _validate_dram_knobs(
    label: str, bounds_map: Mapping[str, Any]
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    raw_steps = bounds_map.get("dr_steps")
    if raw_steps is not None:
        steps = try_int(raw_steps)
        if steps is None or steps < 1:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-107",
                    "Sampling.Bounds.dr_steps",
                    f"{label} dr_steps must be an integer >= 1, got {raw_steps!r}",
                )
            )
    raw_factors = bounds_map.get("dr_scale_factors")
    if raw_factors is not None:
        values = (
            list(raw_factors)
            if isinstance(raw_factors, (list, tuple))
            else [raw_factors]
        )
        parsed = [try_float(item) for item in values]
        if not parsed or any(item is None or item <= 0 for item in parsed):
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-108",
                    "Sampling.Bounds.dr_scale_factors",
                    f"{label} dr_scale_factors must be positive numbers, got {raw_factors!r}",
                )
            )
    return issues


def _validate_stretch_a(
    label: str, bounds_map: Mapping[str, Any], method: str
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    code = "JV2-MTH-072" if method == "PTEnsemble" else "JV2-MTH-116"
    raw = bounds_map.get("stretch_a")
    if raw is None:
        return []
    value = try_float(raw)
    if value is None or value <= 1.0:
        return [
            issue(
                "error",
                code,
                "Sampling.Bounds.stretch_a",
                f"{label} stretch_a must be > 1, got {raw!r}",
            )
        ]
    return []


def _validate_de_knobs(
    label: str, bounds_map: Mapping[str, Any]
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    for key in ("de_gamma", "de_noise"):
        raw = bounds_map.get(key)
        if raw is None:
            continue
        value = try_float(raw)
        if value is None or value < 0:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-126",
                    f"Sampling.Bounds.{key}",
                    f"{label} {key} must be >= 0, got {raw!r}",
                )
            )
    raw_cr = bounds_map.get("de_crossover")
    if raw_cr is not None:
        value = try_float(raw_cr)
        if value is None or value < 0 or value > 1:
            issues.append(
                issue(
                    "error",
                    "JV2-MTH-126",
                    "Sampling.Bounds.de_crossover",
                    f"{label} de_crossover must be in [0, 1], got {raw_cr!r}",
                )
            )
    return issues


__all__ = [
    "method_needs_variables",
    "validate_method_sampling",
]
