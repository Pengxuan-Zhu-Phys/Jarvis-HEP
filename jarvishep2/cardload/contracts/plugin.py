#!/usr/bin/env python3
"""Sampling contract for plug-in methods (external sampler S1, design §5.3–§5.4).

Never imports or runs plug-in code: everything comes from the registry entry
(or, for a shipped plug-in, its file loaded once by the catalog).

* JV2-PLG-001 (info): the plug-in has no ``OPTIONS_SCHEMA``; its options are
  not checked. Leaving the schema out is the author's choice, not a problem.
* JV2-PLG-002 (error): the loaded plug-in cannot be used (snapshot missing or
  changed, or the registry is unreadable).
* JV2-PLG-003 (error): a Jarvis-HEP-reserved Bounds key is invalid.
* JV2-PLG-004 (error): an option breaks ``OPTIONS_SCHEMA`` or is not JSON data.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue


def validate_plugin_sampling(
    sampling: Mapping[str, Any],
    *,
    method: str,
) -> list[ValidationIssue]:
    from jarvis_hep_sampler.errors import PluginConfigError
    from jarvis_hep_sampler.options import (
        freeze_options,
        parse_reserved,
        schema_issues,
        split_bounds,
    )

    from jarvishep2.sampler_catalog import SOURCE_USER, plugin_record
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    man = f"Run: Jarvis man sampler.{method}"
    try:
        record = plugin_record(method)
    except Exception as exc:  # a shipped plug-in that fails to load is a bug
        return [issue("error", "JV2-PLG-002", "Sampling.Method", f"{method}: {exc}")]
    if record is None:
        return issues

    if record.source == SOURCE_USER:
        from jarvis_hep_sampler.registry import SamplerRegistry

        registry = SamplerRegistry()
        entry = registry.get(method)
        problems = registry.problems(entry) if entry is not None else ["registry entry vanished"]
        for problem in problems:
            issues.append(
                issue(
                    "error",
                    "JV2-PLG-002",
                    "Sampling.Method",
                    f"{method}: {problem}",
                    suggestion=f"Load the sampler again: Jarvis sampler load {record.source_path}",
                )
            )

    bounds = sampling.get("Bounds")
    if bounds is None:
        bounds = {}
    if not isinstance(bounds, Mapping):
        return issues  # JV2-BND-000 is reported by the caller
    reserved, options = split_bounds(bounds, record.profile)
    try:
        parse_reserved(reserved, record.profile)
    except PluginConfigError as exc:
        issues.append(
            issue(
                "error",
                "JV2-PLG-003",
                f"Sampling.{exc.yaml_path or 'Bounds'}",
                exc.message,
                hint=man,
                suggestion="seed and point_number belong to Jarvis-HEP; use integers "
                "(seed >= 0, point_number >= 1).",
            )
        )
    try:
        frozen = freeze_options(options)
    except PluginConfigError as exc:
        issues.append(
            issue(
                "error",
                "JV2-PLG-004",
                f"Sampling.{exc.yaml_path or 'Bounds'}",
                exc.message,
                hint=man,
                suggestion="Sampler options must be plain YAML data: strings, numbers, "
                "true/false, null, lists and mappings.",
            )
        )
        return issues

    if record.options_schema is None:
        if options:
            issues.append(
                issue(
                    "info",
                    "JV2-PLG-001",
                    "Sampling.Bounds",
                    f"{method} has no OPTIONS_SCHEMA; its options "
                    f"({', '.join(sorted(map(str, options)))}) are passed through unchecked",
                    hint=man,
                    suggestion="Nothing to fix. The sampler's author may add an "
                    "OPTIONS_SCHEMA to have these options checked.",
                )
            )
        return issues
    for problem in schema_issues(record.options_schema, frozen):
        issues.append(
            issue(
                "error",
                "JV2-PLG-004",
                f"Sampling.{problem.location()}",
                problem.message,
                hint=man,
                suggestion=f"Use the options {method} accepts; Jarvis man sampler.{method} "
                "lists them.",
            )
        )
    return issues


__all__ = ["validate_plugin_sampling"]
