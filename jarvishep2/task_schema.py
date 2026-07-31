"""File-composed JSON Schema validation for the stable Jarvis2 card surface."""

from __future__ import annotations

import json
import re
from difflib import get_close_matches
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from jsonschema import Draft202012Validator
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT202012


SCHEMA_DIR = Path(__file__).with_name("schema")
SCHEMA_PATH = SCHEMA_DIR / "task-card-v2.schema.json"
MANIFEST_PATH = SCHEMA_DIR / "manifest.json"
_NUMERIC_STRING = re.compile(r"[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")


def schema_catalog_lint_errors() -> list[str]:
    """Return authoring errors in the bundled schema catalog.

    ``x-jarvis-zone`` documents who owns an object surface.  Keeping this
    check next to the loader makes a missing zone a CI failure, rather than a
    silent weakening of the strict-card boundary.
    """
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    errors: list[str] = []

    def walk(node: Any, location: str) -> None:
        if isinstance(node, Mapping):
            if node.get("type") == "object" and "x-jarvis-zone" not in node:
                errors.append(f"{location}: object schema has no x-jarvis-zone")
            zone = node.get("x-jarvis-zone")
            if zone is not None and zone not in {"closed", "delegated", "open"}:
                errors.append(f"{location}: invalid x-jarvis-zone {zone!r}")
            if zone == "open" and not isinstance(node.get("x-jarvis-open-reason"), str):
                errors.append(f"{location}: open zone has no x-jarvis-open-reason")
            for key, child in node.items():
                walk(child, f"{location}/{key}")
        elif isinstance(node, list):
            for index, child in enumerate(node):
                walk(child, f"{location}/{index}")

    for relative_name in manifest["schema_files"]:
        path = SCHEMA_DIR / relative_name
        with path.open(encoding="utf-8") as handle:
            walk(json.load(handle), relative_name)
    return errors


@lru_cache(maxsize=1)
def _schema_catalog() -> tuple[dict[str, Any], Registry]:
    """Load the manifest-selected, bundled schemas into a local-only registry."""
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    lint_errors = schema_catalog_lint_errors()
    if lint_errors:
        raise RuntimeError("Invalid Jarvis2 schema catalog:\n" + "\n".join(lint_errors))

    schemas: list[dict[str, Any]] = []
    for relative_name in manifest["schema_files"]:
        path = SCHEMA_DIR / relative_name
        with path.open(encoding="utf-8") as handle:
            schema = json.load(handle)
        Draft202012Validator.check_schema(schema)
        schemas.append(schema)

    registry = Registry().with_resources(
        (schema["$id"], Resource.from_contents(schema, default_specification=DRAFT202012))
        for schema in schemas
    )
    return manifest, registry


@lru_cache(maxsize=1)
def task_card_validator() -> Draft202012Validator:
    """Return the composed root-card validator."""
    manifest, registry = _schema_catalog()
    root = registry.contents(manifest["root"])
    return Draft202012Validator(root, registry=registry)


def _validator_for(schema_uri: str) -> Draft202012Validator:
    _, registry = _schema_catalog()
    return Draft202012Validator(registry.contents(schema_uri), registry=registry)


def _path(parts: Any, prefix: str = "$") -> str:
    rendered = prefix
    for part in parts:
        rendered += f"[{part}]" if isinstance(part, int) else f".{part}"
    return rendered


def _schema_error_guidance(error: Any, prefix: str) -> tuple[str, str | None]:
    """Translate jsonschema's technical validators into an editable YAML fix."""
    schema = error.schema if isinstance(error.schema, Mapping) else {}
    example = schema.get("x-jarvis-example")
    if not isinstance(example, str):
        example = None

    if error.validator == "required":
        missing = re.search(r"'([^']+)' is a required property", error.message)
        field = missing.group(1) if missing else "the required field"
        if example is None:
            if ".execution.input" in prefix:
                example = "- name: input\n  path: input.json\n  type: JSON"
            elif ".execution.output" in prefix:
                example = "- name: output\n  path: output.json\n  type: JSON"
        return f"Add {field!r} at this location using the expected YAML type.", example
    if error.validator == "additionalProperties":
        unknown = re.search(r"\('([^']+)' was unexpected\)", error.message)
        key = unknown.group(1) if unknown else "the unexpected key"
        allowed = schema.get("properties", {})
        names = ", ".join(sorted(allowed)) if isinstance(allowed, Mapping) else ""
        close = get_close_matches(key, allowed, n=1, cutoff=0.6) if isinstance(allowed, Mapping) else []
        rename = f" Did you mean {close[0]!r}?" if close else ""
        suffix = f" Allowed keys: {names}." if names else ""
        return f"Remove or rename {key!r}.{rename}{suffix}", example
    if error.validator in {"enum", "const"}:
        allowed = schema.get("enum")
        if error.validator == "const":
            allowed = [schema.get("const")]
        values = ", ".join(repr(value) for value in allowed) if isinstance(allowed, list) else "the required value"
        return f"Replace this value with {values}.", example
    if error.validator == "type":
        expected = schema.get("type", "the required type")
        if expected == "string" and isinstance(error.instance, bool):
            return (
                'Quote YAML boolean-like values when a string is required, for example: value: "on".',
                example,
            )
        return f"Replace this value with a YAML {expected}.", example
    if error.validator in {"minimum", "exclusiveMinimum", "maximum", "minLength", "minItems"}:
        return "Adjust the value so it satisfies the bound stated in the error message.", example
    if error.validator in {"anyOf", "oneOf"}:
        required_fields = sorted({
            match.group(1)
            for child in error.context
            if child.validator == "required"
            for match in [re.search(r"'([^']+)' is a required property", child.message)]
            if match
        })
        if required_fields:
            return (
                "Use one accepted form; add one of the alternative required fields: "
                + ", ".join(required_fields) + ".",
                example,
            )
        return "Use one complete accepted form for this object; inspect the listed schema alternatives.", example
    return "Correct this value to satisfy the stated schema constraint.", example


def _issues_for(validator: Draft202012Validator, value: Any, prefix: str) -> list[Any]:
    from jarvishep2.task_validation import issue

    errors = sorted(validator.iter_errors(value), key=lambda err: list(err.absolute_path))
    issues = []
    for error in errors:
        suggestion, example = _schema_error_guidance(error, prefix)
        issues.append(issue(
            "error", "JV2-SCH-001", _path(error.absolute_path, prefix), error.message,
            hint="See docs/task-card-schema.md for the strict card interface.",
            suggestion=suggestion, example=example,
        ))
    return issues


def _validate_selected_io(config: Mapping[str, Any]) -> list[Any]:
    """Validate Portal I/O vocabulary, then enrich it with bundled schemas.

    Portal owns supported format names.  The manifest is intentionally only an
    optional local refinement layer, so a newer Portal format cannot be
    rejected merely because Jarvis2 has not yet bundled a detailed schema.
    """
    from jarvishep2.task_validation import issue
    from jarvishep2.io_portal import available_io_formats

    manifest, _ = _schema_catalog()
    issues: list[Any] = []
    calculators = config.get("Calculators")
    if not isinstance(calculators, Mapping):
        return issues
    modules = calculators.get("Modules")
    if not isinstance(modules, list):
        return issues

    for module_index, module in enumerate(modules):
        if not isinstance(module, Mapping):
            continue
        execution = module.get("execution")
        if not isinstance(execution, Mapping):
            continue
        for direction in ("input", "output"):
            entries = execution.get(direction)
            if not isinstance(entries, list):
                continue
            formats = manifest["io"].get(direction, {})
            portal_formats = set(available_io_formats(direction))
            for entry_index, entry in enumerate(entries):
                prefix = f"$.Calculators.Modules[{module_index}].execution.{direction}[{entry_index}]"
                if not isinstance(entry, Mapping):
                    continue
                raw_kind = entry.get("type")
                kind = raw_kind.strip() if isinstance(raw_kind, str) else raw_kind
                schema_uri = formats.get(kind) if isinstance(kind, str) else None
                if not isinstance(kind, str) or kind not in portal_formats:
                    supported = ", ".join(sorted(portal_formats))
                    issues.append(issue(
                        "error", "JV2-SCH-002", f"{prefix}.type",
                        f"unsupported {direction} format {raw_kind!r}; Portal supports: {supported}",
                        hint="Use a format registered by Jarvis-HEP-Portal for this direction.",
                        suggestion="Replace the type with one of the listed Portal formats.",
                    ))
                    continue
                if not isinstance(schema_uri, str):
                    # Portal accepts it; a detailed Jarvis2 schema is optional.
                    continue
                normalized_entry = dict(entry)
                normalized_entry["type"] = kind
                issues.extend(_issues_for(_validator_for(schema_uri), normalized_entry, prefix))
    return issues


def _numeric_string_issues(config: Mapping[str, Any]) -> list[Any]:
    """Warn when an accepted YAML numeric string can be written canonically."""
    from jarvishep2.task_validation import issue

    warnings: list[Any] = []
    numeric_paths = (
        ("Sampling", "Radius"),
        ("Sampling", "MaxAttempt"),
        ("Sampling", "MaxWorker"),
        ("Sampling", "Point number"),
        ("Sampling", "point_number"),
        ("Sampling", "Seed"),
        ("Sampling", "seed"),
    )
    for parts in numeric_paths:
        parent: Any = config
        for part in parts[:-1]:
            parent = parent.get(part) if isinstance(parent, Mapping) else None
        value = parent.get(parts[-1]) if isinstance(parent, Mapping) else None
        if isinstance(value, str) and _NUMERIC_STRING.fullmatch(value):
            path = "$" + "".join(f".{part}" for part in parts)
            warnings.append(issue(
                "warning", "JV2-SCH-003", path,
                f"numeric string {value!r} is accepted for compatibility",
                hint="YAML may parse scientific notation such as 1e-5 as a string.",
                suggestion=f"Use the canonical numeric form: {value.replace('E', 'e')}",
            ))

    sampling = config.get("Sampling")
    adaptive = sampling.get("AdaptiveBridson") if isinstance(sampling, Mapping) else None
    if not isinstance(adaptive, Mapping) and isinstance(sampling, Mapping):
        adaptive = sampling.get("adaptive_bridson")
    if isinstance(adaptive, Mapping):
        numeric_keys = {
            "target_value", "outer_half_width", "min_radius", "initial_radius",
            "refinement_factor", "bridge_span_factor", "max_generations",
            "max_points", "max_new_per_generation", "k_ref", "knn_k",
            "core_half_width", "threshold", "function_tolerance", "final_half_width",
            "outer_shrink_factor", "core_spacing_factor", "min_cores_for_coverage",
            "quiet_fill_passes", "max_fill_passes", "endpoint_stall_passes",
            "endpoint_omni_probes",
        }
        for key, value in adaptive.items():
            if key not in numeric_keys:
                continue
            if isinstance(value, str) and _NUMERIC_STRING.fullmatch(value):
                warnings.append(issue(
                    "warning", "JV2-SCH-003", f"$.Sampling.AdaptiveBridson.{key}",
                    f"numeric string {value!r} is accepted for compatibility",
                    hint="YAML may parse scientific notation such as 1e-5 as a string.",
                    suggestion=f"Use the canonical numeric form: {value.replace('E', 'e')}",
                ))
    return warnings


def _open_zone_issues(config: Mapping[str, Any]) -> list[Any]:
    """Make intentionally unfinished surfaces visible without rejecting cards."""
    from jarvishep2.task_validation import issue

    sampling = config.get("Sampling")
    if not isinstance(sampling, Mapping):
        return []
    reasons = {
        "Control": "RLTPMCMC controller has not been migrated to the V2 contract.",
        "Diagnostics": "RLTPMCMC diagnostics has not been migrated to the V2 contract.",
        "PPO": "RLTPMCMC PPO settings have not been migrated to the V2 contract.",
        "Reward": "RLTPMCMC reward settings have not been migrated to the V2 contract.",
    }
    present = [key for key in reasons if key in sampling]
    if not present:
        return []
    return [issue(
        "warning", "JV2-SCH-004", "$.Sampling." + ", ".join(present),
        "open validation zone: " + " ".join(reasons[key] for key in present),
        hint="Their nested keys are passed through without Jarvis2 schema checking.",
        suggestion="Use the documented sampler settings, or migrate these blocks before relying on strict validation.",
    )]


def validate_task_card_schema(config: Mapping[str, Any]) -> list[Any]:
    """Return normal V2 diagnostics from the root and selected local schemas."""
    issues = _issues_for(task_card_validator(), dict(config), "$")
    manifest, _ = _schema_catalog()

    sampling = config.get("Sampling")
    if isinstance(sampling, Mapping):
        raw_method = sampling.get("Method")
        method = raw_method.strip() if isinstance(raw_method, str) else raw_method
        schema_uri = manifest["sampling_methods"].get(method) if isinstance(method, str) else None
        if schema_uri is not None:
            normalized_sampling = dict(sampling)
            normalized_sampling["Method"] = method
            issues.extend(_issues_for(_validator_for(schema_uri), normalized_sampling, "$.Sampling"))

    issues.extend(_validate_selected_io(config))
    issues.extend(_numeric_string_issues(config))
    issues.extend(_open_zone_issues(config))
    return issues


__all__ = [
    "MANIFEST_PATH", "SCHEMA_PATH", "schema_catalog_lint_errors",
    "task_card_validator", "validate_task_card_schema",
]
