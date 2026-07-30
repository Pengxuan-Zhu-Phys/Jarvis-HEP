"""File-composed JSON Schema validation for the stable Jarvis2 card surface."""

from __future__ import annotations

import json
import re
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping

from jsonschema import Draft202012Validator
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT202012


SCHEMA_DIR = Path(__file__).with_name("schema")
SCHEMA_PATH = SCHEMA_DIR / "task-card-v2.schema.json"
MANIFEST_PATH = SCHEMA_DIR / "manifest.json"


@lru_cache(maxsize=1)
def _schema_catalog() -> tuple[dict[str, Any], Registry]:
    """Load the manifest-selected, bundled schemas into a local-only registry."""
    with MANIFEST_PATH.open(encoding="utf-8") as handle:
        manifest = json.load(handle)

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
        suffix = f" Allowed keys: {names}." if names else ""
        return f"Remove or rename {key!r}.{suffix}", example
    if error.validator in {"enum", "const"}:
        allowed = schema.get("enum")
        if error.validator == "const":
            allowed = [schema.get("const")]
        values = ", ".join(repr(value) for value in allowed) if isinstance(allowed, list) else "the required value"
        return f"Replace this value with {values}.", example
    if error.validator == "type":
        expected = schema.get("type", "the required type")
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
    """Dispatch calculator I/O items through manifest-selected format schemas."""
    from jarvishep2.task_validation import issue

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
            formats = manifest["io"][direction]
            for entry_index, entry in enumerate(entries):
                prefix = f"$.Calculators.Modules[{module_index}].execution.{direction}[{entry_index}]"
                if not isinstance(entry, Mapping):
                    continue
                raw_kind = entry.get("type")
                kind = raw_kind.strip() if isinstance(raw_kind, str) else raw_kind
                schema_uri = formats.get(kind) if isinstance(kind, str) else None
                if schema_uri is None:
                    issues.append(issue(
                        "error", "JV2-SCH-002", f"{prefix}.type",
                        f"unsupported {direction} format {raw_kind!r}; register a local schema in schema/manifest.json",
                        hint="Add the format schema file and its manifest entry together with its Portal adapter.",
                    ))
                    continue
                normalized_entry = dict(entry)
                normalized_entry["type"] = kind
                issues.extend(_issues_for(_validator_for(schema_uri), normalized_entry, prefix))
    return issues


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
    return issues


__all__ = ["MANIFEST_PATH", "SCHEMA_PATH", "task_card_validator", "validate_task_card_schema"]
