#!/usr/bin/env python3
"""Dynesty / MultiNest Bounds contracts (D14 L2)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.contracts.common import try_float, try_int, unknown_keys

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

# Import Nested allow-lists from the sampler module (single source of truth).
from jarvishep2.Sampling.dynesty_sampler import (
    DYNAMIC_RUN_NESTED_KEYS,
    NESTED_CONSTRUCTOR_USER_KEYS,
    STATIC_RUN_NESTED_KEYS,
)

# Canonical + documented aliases for top-level Bounds keys.
_BOUNDS_META_ALLOWED: frozenset[str] = frozenset(
    {
        "nlive",
        "n_live",
        "rseed",
        "seed",
        "Seed",
        # Bounds.dynamic / Dynamic removed: Method selects the engine
        # (Dynesty → always DynamicNestedSampler, MultiNest → always static).
        "dlogz",
        "dlogz_init",
        "run_nested",
        "sampler",
        "Sampler",
        "constructor",
        "Constructor",
    }
)

# Flat constructor keys may also appear at Bounds top level (official dynesty style).
_BOUNDS_TOP_ALLOWED: frozenset[str] = _BOUNDS_META_ALLOWED | NESTED_CONSTRUCTOR_USER_KEYS

_RUN_NESTED_ALLOWED: frozenset[str] = STATIC_RUN_NESTED_KEYS | DYNAMIC_RUN_NESTED_KEYS


def validate_nested_bounds(
    bounds: Mapping[str, Any],
    *,
    method: str,
) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    path = "Sampling.Bounds"

    # Retired YAML switch: engine is selected only by Sampling.Method.
    for dyn_key in ("dynamic", "Dynamic"):
        if dyn_key in bounds:
            issues.append(
                issue(
                    "error",
                    "JV2-BND-012",
                    f"{path}.{dyn_key}",
                    f"Bounds.{dyn_key} is not a V2 setting; the engine is fixed by "
                    f"Sampling.Method (Dynesty → always DynamicNestedSampler, "
                    f"MultiNest → always static NestedSampler). Remove Bounds.{dyn_key}.",
                    hint=(
                        "Use Method: Dynesty for dynamic nested sampling, or "
                        "Method: MultiNest for static NestedSampler — not a dynamic flag."
                    ),
                )
            )

    extra = unknown_keys(bounds, _BOUNDS_TOP_ALLOWED)
    # dynamic already reported with a specific code — drop from generic unknown list.
    extra = [k for k in extra if k not in {"dynamic", "Dynamic"}]
    if extra:
        issues.append(
            issue(
                "error",
                "JV2-BND-001",
                path,
                f"unknown key(s): {', '.join(extra)}; "
                f"see Dynesty/MultiNest Bounds allow-list "
                f"(meta + official constructor keys)",
                hint="Remove typos or move pass-through kwargs under Bounds.sampler / run_nested",
            )
        )

    # nlive
    if "nlive" in bounds or "n_live" in bounds:
        raw_nlive = bounds.get("nlive", bounds.get("n_live"))
        nlive = try_int(raw_nlive)
        if nlive is None or nlive < 2:
            issues.append(
                issue(
                    "error",
                    "JV2-BND-010",
                    f"{path}.nlive",
                    f"expected integer ≥ 2, got {raw_nlive!r}",
                )
            )

    # dlogz optional numbers
    for key in ("dlogz", "dlogz_init"):
        if key not in bounds:
            continue
        val = try_float(bounds.get(key))
        if val is None:
            issues.append(
                issue(
                    "error",
                    "JV2-BND-013",
                    f"{path}.{key}",
                    f"expected a number, got {bounds.get(key)!r}",
                )
            )

    # run_nested block
    run_raw = bounds.get("run_nested")
    if run_raw is not None:
        if not isinstance(run_raw, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-BND-020",
                    f"{path}.run_nested",
                    f"expected a mapping, got {type(run_raw).__name__}",
                )
            )
        else:
            run_extra = unknown_keys(run_raw, _RUN_NESTED_ALLOWED)
            if run_extra:
                issues.append(
                    issue(
                        "error",
                        "JV2-BND-021",
                        f"{path}.run_nested",
                        f"unknown key(s): {', '.join(run_extra)}",
                    )
                )

    # sampler / constructor blocks
    for block_key in ("sampler", "Sampler", "constructor", "Constructor"):
        block = bounds.get(block_key)
        if block is None:
            continue
        if not isinstance(block, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-BND-030",
                    f"{path}.{block_key}",
                    f"expected a mapping, got {type(block).__name__}",
                )
            )
            continue
        # HEP-injected keys are stripped at runtime; reject them here for clarity.
        injected = {"loglikelihood", "prior_transform", "ndim", "pool", "rstate"}
        bad_injected = sorted(str(k) for k in block.keys() if str(k) in injected)
        if bad_injected:
            issues.append(
                issue(
                    "error",
                    "JV2-BND-031",
                    f"{path}.{block_key}",
                    f"key(s) are injected by HEP and must not appear in YAML: "
                    f"{', '.join(bad_injected)}",
                )
            )
        block_extra = unknown_keys(block, NESTED_CONSTRUCTOR_USER_KEYS)
        # Filter out injected already reported.
        block_extra = [k for k in block_extra if k not in injected]
        if block_extra:
            issues.append(
                issue(
                    "error",
                    "JV2-BND-032",
                    f"{path}.{block_key}",
                    f"unknown constructor key(s): {', '.join(block_extra)}",
                )
            )

    return issues


__all__ = ["validate_nested_bounds"]
