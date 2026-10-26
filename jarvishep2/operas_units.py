"""Host bindings for JO units. See Jarvis-Books/Docs/Jarvis-HEP/design/operas-units.md.

Unit arithmetic and numeric validation belong exclusively to Jarvis-Operas.
This module handles paths, host diagnostics and serializable execution identity.
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

DATABASE_KEY = "_jarvis_units"
RESULT_KEY = "unit_execution"
BRIDGE_VERSION = 1


class UnitBridgeError(ValueError):
    def __init__(
        self,
        message: str,
        *,
        code: str = "JV2-UNIT-002",
        path: str = "$.Operas.Modules",
        jo_code: str | None = None,
    ) -> None:
        self.code, self.path, self.jo_code = code, path, jo_code
        self.suggestion = (
            "Check the binding and registered unit contract: Jarvis man operas"
        )
        super().__init__(
            f"{code} {path}: {message}"
            + (f" [{jo_code}]" if jo_code else "")
            + f". {self.suggestion}"
        )

    def diagnostic(self) -> dict[str, Any]:
        return dict(
            code=self.code, path=self.path, jo_code=self.jo_code, message=str(self)
        )


def enabled(config: Mapping[str, Any]) -> bool:
    return "unit_policy" in config or any(
        isinstance(item, Mapping) and "unit" in item
        for side in ("input", "output")
        for item in specs(config, side)
    )


def specs(config, side):
    value = config.get(side)
    return value if isinstance(value, list) else []


def modules(config: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    block = config.get("Operas") or {}
    return (
        [m for m in specs(block, "Modules") if isinstance(m, Mapping)]
        if isinstance(block, Mapping)
        else []
    )


def pointer_token(value: str) -> str:
    return str(value).replace("~", "~0").replace("/", "~1")


def entry_pointer(entry: str) -> str:
    if not isinstance(entry, str):
        raise ValueError("entry must be a string")
    if entry.startswith("/"):
        if re.search(r"~(?![01])", entry):
            raise ValueError("invalid JSON Pointer escape; use ~0 for ~ and ~1 for /")
        return entry
    return "/" + "/".join(pointer_token(p) for p in entry.split(".")) if entry else ""


def read_entry(data: Any, entry: str) -> Any:
    current = data
    for token in entry_pointer(entry).split("/")[1:]:
        key = token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, Mapping) and key in current:
            current = current[key]
        elif (
            isinstance(current, (list, tuple))
            and re.fullmatch(r"0|[1-9][0-9]*", key)
            and int(key) < len(current)
        ):
            current = current[int(key)]
        else:
            raise KeyError(f"missing entry {entry!r}")
    return current


def static_issues(config: Mapping[str, Any]) -> list[Any]:
    """No discovery or user imports: metadata checks belong to run preflight."""
    from jarvishep2.task_validation import issue

    issues = []
    active = any(enabled(m) for m in modules(config))
    for index, module in enumerate(modules(config)):
        base = f"$.Operas.Modules[{index}]"
        if not enabled(module):
            for i, spec in enumerate(specs(module, "output")):
                if isinstance(spec, Mapping) and spec.get("entry") == "":
                    issues.append(
                        issue(
                            "error",
                            "JV2-UNIT-001",
                            f"{base}.output[{i}].entry",
                            "empty entry requires an enabled unit module",
                            suggestion="Set unit_policy: canonical or use a nonempty entry.",
                        )
                    )
            continue
        issues.append(
            issue(
                "info",
                "JV2-UNIT-001",
                base,
                "Unit structure checked; operator metadata and availability require run/check preflight.",
                suggestion="Inspect the declaration with Jarvis operas info NAME; run Jarvis check for execution.",
            )
        )
        seen_inputs: set[str] = set()
        for side in ("input", "output"):
            seen: set[str] = set()
            for i, spec in enumerate(specs(module, side)):
                if not isinstance(spec, Mapping):
                    continue
                path = f"{base}.{side}[{i}]"
                name = spec.get("name")
                if not isinstance(name, str):
                    continue  # The closed schema reports malformed names.
                reason = None
                if name in seen:
                    reason = f"duplicate {side} binding {name!r}"
                seen.add(name)
                if side == "input":
                    seen_inputs.add(name)
                    if "entry" in spec and "expression" in spec:
                        reason = "input entry and expression are mutually exclusive"
                elif "expression" in spec:
                    reason = "unit output bindings use entry, not expression"
                try:
                    pointer = entry_pointer(spec.get("entry", name))
                    if not pointer and side == "input":
                        reason = "empty entry is only supported for output"
                except ValueError as exc:
                    reason = str(exc)
                if reason:
                    issues.append(
                        issue(
                            "error",
                            "JV2-UNIT-001",
                            path,
                            reason,
                            suggestion="Use one unique name and a valid entry; see Jarvis man operas.",
                        )
                    )
                if (
                    "unit" in spec
                    and isinstance(spec["unit"], str)
                    and spec["unit"].strip()
                ):
                    try:
                        from jarvis_operas import conversion_factor

                        conversion_factor(spec["unit"], spec["unit"])
                    except Exception as exc:
                        issues.append(
                            issue(
                                "error",
                                "JV2-UNIT-001",
                                path + ".unit",
                                str(exc),
                                jo_code=getattr(exc, "code", None),
                                suggestion="Use a supported JO unit symbol; see Jarvis operas info NAME.",
                            )
                        )
        kwargs = module.get("kwargs")
        overlap = seen_inputs.intersection(
            kwargs if isinstance(kwargs, Mapping) else {}
        )
        if overlap:
            issues.append(
                issue(
                    "error",
                    "JV2-UNIT-001",
                    base,
                    f"input/kwargs overlap: {sorted(overlap)}",
                    suggestion="Declare each parameter in either input or kwargs, once.",
                )
            )
    if active:
        # A collision anywhere in user-authored keys/names cannot be fixed by
        # silently replacing the user's value when Archiver injects metadata.
        def walk(value: Any, path: str) -> None:
            if isinstance(value, Mapping):
                for k, v in value.items():
                    if k == DATABASE_KEY or (k == "name" and v == DATABASE_KEY):
                        issues.append(
                            issue(
                                "error",
                                "JV2-UNIT-005",
                                path,
                                f"{DATABASE_KEY} is reserved for unit metadata",
                                suggestion="Rename the user field.",
                            )
                        )
                    walk(v, f"{path}.{k}")
            elif isinstance(value, list):
                for i, v in enumerate(value):
                    walk(v, f"{path}[{i}]")

        walk(getattr(config, "raw_task_card", config), "$")
    return issues


@dataclass
class OperaResult:
    values: dict[str, Any]
    units: dict[str, str | None] = field(default_factory=dict)
    call: dict[str, Any] | None = None


def prepare_plan(config: Mapping[str, Any]) -> Any:
    from jarvishep2.operas_functions import ensure_operas_registry_discovered

    operator = str(config["operator"])
    try:
        from jarvis_operas import UnitCallPlan
        from importlib.metadata import version
        from packaging.version import Version

        if Version(version("Jarvis-Operas")) < Version("1.5.1"):
            raise ImportError("JO unit integration needs Jarvis-Operas>=1.5.1")
        if not all(hasattr(UnitCallPlan, n) for n in ("run_inputs", "arun_inputs")):
            raise ImportError("JO unit container API is unavailable")
        registry = ensure_operas_registry_discovered()
        declaration = registry.get(operator)
        if not declaration.metadata.get("unit_contract"):
            raise UnitBridgeError(
                f"{operator}: registered unit_contract is required",
                jo_code="contract_missing",
            )
        return registry.prepare(
            operator,
            strict=config.get("unit_policy") == "strict",
            input_units={
                s["name"]: s["unit"] for s in config.get("input", []) if "unit" in s
            },
            output_units={
                s["name"]: s["unit"] for s in config.get("output", []) if "unit" in s
            },
            output_bindings={
                s["name"]: entry_pointer(s.get("entry", s["name"]))
                for s in config.get("output", [])
            },
        )
    except UnitBridgeError:
        raise
    except Exception as exc:
        raise UnitBridgeError(
            f"{operator}: {exc}", jo_code=getattr(exc, "code", "contract_missing")
        ) from exc


def plan_identity(plan: Any) -> dict[str, Any]:
    # registry.info exposes a detached, JSON-compatible normalized contract.
    from jarvishep2.operas_functions import ensure_operas_registry_discovered

    info = ensure_operas_registry_discovered().info(plan.declaration.full_name)
    return {"operator": plan.declaration.full_name, "contract": info["unit_contract"]}


def digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def execution_identity(config: Mapping[str, Any]) -> dict[str, Any] | None:
    active = [m for m in modules(config) if enabled(m)]
    if not active:
        return None
    problems = [i for i in static_issues(config) if i.level == "error"]
    if problems:
        first = problems[0]
        raise UnitBridgeError(first.message, code=first.code, path=first.path)
    from importlib.metadata import version
    from jarvishep2.workflow import resolve_module_layers

    layers = resolve_module_layers(modules(config))
    owners: dict[tuple[int, str], str] = {}
    declarations = []
    calculators = config.get("Calculators") or {}
    calc_modules = (
        specs(calculators, "Modules") if isinstance(calculators, Mapping) else []
    )
    calc_layers = resolve_module_layers(calc_modules)
    calc_owners = {}
    for module in calc_modules:
        execution = module.get("execution") or {}
        for output in execution.get("output", []):
            for variable in output.get("variables", []):
                name = variable.get("name")
                if not name:
                    continue
                key = (calc_layers[module["name"]], name)
                if key in calc_owners and calc_owners[key] != module["name"]:
                    raise UnitBridgeError(
                        f"same-layer calculator output {name!r}", code="JV2-UNIT-001"
                    )
                calc_owners[key] = module["name"]
    for index, module in enumerate(modules(config)):
        name = module.get("name", str(index))
        for spec in module.get("output", []):
            key = (layers[name], spec["name"])
            if key in owners:
                raise UnitBridgeError(
                    f"same-layer output {key[1]!r}: {owners[key]}, {name}",
                    code="JV2-UNIT-001",
                )
            owners[key] = name
        if enabled(module):
            try:
                declarations.append(plan_identity(prepare_plan(module)))
            except UnitBridgeError as exc:
                raise UnitBridgeError(
                    str(exc),
                    code=exc.code,
                    path=f"$.Operas.Modules[{index}]",
                    jo_code=exc.jo_code,
                ) from exc
    payload = dict(
        bridge_version=BRIDGE_VERSION,
        operas=copy.deepcopy(config.get("Operas")),
        declarations=declarations,
        jo_version=version("Jarvis-Operas"),
    )
    try:
        identity = digest(payload)
    except (TypeError, ValueError) as exc:
        raise UnitBridgeError(
            f"unit configuration must contain finite JSON-compatible values: {exc}",
            code="JV2-UNIT-001",
        ) from exc
    return dict(
        schema_version=1,
        mode="units",
        execution_identity=identity,
        identity_payload=payload,
    )
