"""Actionable, stable guidance for task-card validation diagnostics.

This module intentionally has no dependency on the validation-report classes so
every contract can use it without creating an import cycle.
"""

from __future__ import annotations

import re


_VARIABLE_EXAMPLE = """Variables:
  - name: x
    distribution:
      type: Flat
      parameters: {min: 0.0, max: 1.0}"""

_METHOD_EXAMPLE = """Sampling:
  Method: Random
  Bounds:
    Point number: 100
    Seed: 0
  Variables:
    - name: x
      distribution:
        type: Flat
        parameters: {min: 0.0, max: 1.0}"""

_GUIDANCE_BY_PREFIX: tuple[tuple[str, str, str | None], ...] = (
    ("Sampling.Method", "Set Method to one of the methods listed in the error, then add that method's required Bounds fields.", _METHOD_EXAMPLE),
    ("Sampling.Variables", "Provide a non-empty YAML list of variable mappings; each variable needs name and distribution.", _VARIABLE_EXAMPLE),
    ("Sampling.Bounds.path", "Provide the path to the fixed-point CSV file under Bounds.", "Bounds:\n  path: points.csv"),
    ("Sampling.Bounds", "Put method-specific sampler knobs under Bounds (Radius, Point number, Seed, CSV path, dynesty kwargs, …).", "Bounds:\n  Point number: 100\n  Seed: 0"),
    ("Sampling.CSV", "CSV settings moved to Sampling.Bounds; use Bounds.path (and optional variables).", "Bounds:\n  path: points.csv"),
    ("Sampling.Radius", "Radius moved to Sampling.Bounds.Radius.", "Bounds:\n  Radius: 0.1\n  MaxAttempt: 30"),
    ("Sampling.MaxAttempt", "MaxAttempt moved to Sampling.Bounds.MaxAttempt.", "Bounds:\n  Radius: 0.1\n  MaxAttempt: 30"),
    ("Sampling.Point number", "Point number moved to Sampling.Bounds['Point number'].", "Bounds:\n  Point number: 100"),
    ("Sampling.Seed", "Seed moved to Sampling.Bounds.Seed.", "Bounds:\n  Seed: 0"),
    ("Sampling.AdaptiveBridson", "AdaptiveBridson knobs moved under Sampling.Bounds (no nested AdaptiveBridson block).", "Bounds:\n  target_expression: f\n  target_value: 0.12"),
    ("Calculators.Modules", "Use a module mapping with name; add execution when this module runs an external calculator.", "Modules:\n  - name: MyCalculator\n    execution:\n      path: .\n      commands: [\"./run.sh\"]"),
    ("Operas.Modules", "Use an Opera module mapping with name and operator.", "Modules:\n  - name: likelihood\n    operator: package.module.function"),
    ("EnvReqs.V2.redis", "Use a Redis mapping with host, port, and db only when overriding defaults.", "redis:\n  host: 127.0.0.1\n  port: 6379\n  db: 0"),
    ("EnvReqs.V2", "Use a mapping and keep only documented V2 runtime settings.", "V2:\n  workers: 4\n  batch_size: 1"),
    ("Scan.sample_directory", "Use a mapping; numeric bucket settings must be positive integers.", "sample_directory:\n  enabled: true\n  limit: 200\n  width: 6"),
    ("Calculators.Archiver", "Use an Archiver mapping and choose a supported mode.", "Archiver:\n  mode: rolling"),
    ("Calculators.Cleanup", "Use a Cleanup mapping and choose a supported strategy.", "Cleanup:\n  strategy: never"),
)


def guidance_for(code: str, path: str, message: str) -> tuple[str, str | None]:
    """Return a concrete correction direction for every V2 diagnostic."""
    if code == "JV2-SCH-002":
        return (
            "Correct type to a registered format, or add its local schema and Portal adapter together.",
            "type: JSON",
        )
    if code in {"JV2-MTH-002", "JV2-MTH-003", "JV2-MTH-040"}:
        return "Choose an implemented Method, then provide its required Bounds configuration.", _METHOD_EXAMPLE
    if code == "JV2-MTH-001":
        return (
            "Move this sampler setting under Sampling.Bounds; top-level Sampling no longer accepts method knobs.",
            "Bounds:\n  Seed: 0",
        )
    if code in {"JV2-MTH-010", "JV2-MTH-011"}:
        return "Set Sampling.Bounds.Radius to a positive number for Bridson sampling.", "Bounds:\n  Radius: 0.1\n  MaxAttempt: 30"
    if code in {"JV2-MTH-012", "JV2-MTH-013"}:
        return "Set Sampling.Bounds.MaxAttempt to an integer of at least 1 for Bridson sampling.", "Bounds:\n  Radius: 0.1\n  MaxAttempt: 30"
    if code in {"JV2-MTH-020", "JV2-MTH-021"}:
        return "Set Sampling.Bounds['Point number'] (or point_number) to an integer of at least 1.", "Bounds:\n  Point number: 100"
    if code in {"JV2-MTH-030", "JV2-MTH-031"}:
        return "Provide the path to the fixed-point CSV file under Sampling.Bounds.", "Bounds:\n  path: points.csv"
    if code in {"JV2-MTH-041", "JV2-MTH-042"}:
        return "Set target_expression and target_value under Sampling.Bounds for AdaptiveBridson.", "Bounds:\n  target_expression: f\n  target_value: 0.12"
    # Parameter-level codes must win over the Sampling.Variables list prefix
    # (D21.14: JV2-VAR-031 used to suggest "provide a non-empty variable list").
    if code == "JV2-VAR-031":
        missing = re.findall(r"missing required key\(s\) for type \w+:\s*(.+)$", message)
        keys = missing[0].strip() if missing else "the keys named in the error"
        return (
            f"Add the required distribution parameter(s) under parameters: {keys}.",
            "parameters: {min: 0.0, max: 1.0}",
        )
    if code == "JV2-VAR-032":
        return (
            "Rename or remove the unknown distribution parameter; only the allowed "
            "keys listed in the error are valid for this type.",
            None,
        )
    for prefix, suggestion, example in _GUIDANCE_BY_PREFIX:
        if path.startswith(prefix):
            return suggestion, example

    if "unknown key" in message or "unsupported setting" in message:
        return "Remove the misspelled key or rename it to one of the allowed keys listed in the error.", None
    if "expected one of:" in message or "invalid; expected one of:" in message:
        return "Replace the value with one of the allowed values listed in the error; spelling and case matter.", None
    if "expected a mapping" in message:
        return "Replace this scalar/list with an indented YAML mapping (key: value pairs).", None
    if "expected a list" in message:
        return "Replace this value with an indented YAML list whose items start with '- '.", None
    if "expected integer" in message or "expected non-negative integer" in message:
        return "Replace the value with an integer that satisfies the stated range.", None
    if "expected number" in message or "positive number" in message:
        return "Replace the value with a numeric YAML scalar that satisfies the stated range.", None
    if "required" in message or "missing" in message:
        fields = re.findall(r"(?:'([^']+)'|\b([A-Za-z][A-Za-z0-9_ ]*)\b)", message)
        names = [left or right for left, right in fields]
        suffix = f" Add: {names[0]}" if names else ""
        return f"Add the required field at this location.{suffix}", None
    return "Correct the value so it satisfies the constraint stated in the error message.", None


__all__ = ["guidance_for"]
