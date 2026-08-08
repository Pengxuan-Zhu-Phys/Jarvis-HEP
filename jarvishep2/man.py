#!/usr/bin/env python3
"""Jarvis man: runtime YAML writing manual (DESIGN_YAML_MAN_2.0).

Output is always English. Man only renders surfaces the validator knows about (K1).
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Mapping

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.task_schema import (
    resolve_schema_ref,
    schema_by_id,
    schema_manifest,
)

_DOMAINS = frozenset({"yaml", "sampler", "calculator", "operas", "tokens", "example"})
_ROOT_BLOCKS = ("Scan", "Sampling", "Calculators", "Operas", "EnvReqs", "LibDeps")
_CALC_TOPICS = ("module", "execution", "modes", "pools")
_ARRAY_INDEX = re.compile(r"\[(?:\d+|\*)?\]")
_EXAMPLE_CATALOG: dict[str, str] = {
    "bridson": "bin/quickstart_bridson_operas.yaml",
    "csv": "bin/quickstart_csv_operas.yaml",
    "dynesty": "bin/sampling/Sampling_Dynesty_Simple.yaml",
    "dynesty-full": "bin/sampling/Sampling_Dynesty_Full.yaml",
    "multinest": "bin/sampling/Sampling_MultiNest_Simple.yaml",
    "multinest-full": "bin/sampling/Sampling_MultiNest_Full.yaml",
    "calculator": "bin/quickstart_calculator.yaml",
}

# Diagnostic codes most relevant to sampler Bounds (subset; D24.10 expands fully).
_METHOD_DIAGNOSTICS: dict[str, list[str]] = {
    "Bridson": ["JV2-MTH-010", "JV2-MTH-011", "JV2-MTH-012", "JV2-MTH-013"],
    "Random": ["JV2-MTH-001", "JV2-MTH-020", "JV2-MTH-021", "JV2-SCH-001"],
    "Grid": ["JV2-VAR-001", "JV2-SCH-001"],
    "CSV": ["JV2-MTH-030", "JV2-MTH-031", "JV2-SCH-001"],
    "AdaptiveBridson": ["JV2-MTH-041", "JV2-MTH-042", "JV2-SCH-001"],
    "Dynesty": ["JV2-BND-001", "JV2-BND-012", "JV2-BND-030", "JV2-SCH-001"],
    "MultiNest": ["JV2-BND-001", "JV2-BND-012", "JV2-BND-020", "JV2-SCH-001"],
}


def _console(*, force_terminal: bool | None = None) -> Console:
    return Console(
        stderr=False,
        force_terminal=force_terminal,
        highlight=False,
        soft_wrap=True,
    )


def normalize_yaml_path(raw: str) -> str:
    """Normalize user path input to ``$.A.B[].C`` form."""
    text = str(raw or "").strip()
    if text.startswith("$."):
        text = text[2:]
    elif text.startswith("$"):
        text = text[1:].lstrip(".")
    text = _ARRAY_INDEX.sub("[]", text)
    parts = [p for p in text.split(".") if p]
    return "$." + ".".join(parts) if parts else "$"


def _schema_status(schema: Mapping[str, Any]) -> str:
    status = str(schema.get("x-jarvis-status") or "stable").strip().lower()
    return "unstable" if status == "unstable" else "stable"


def _collect_properties(schema: Mapping[str, Any]) -> dict[str, Any]:
    """Merge ``properties`` across ``allOf`` / ``$ref`` expansions (shallow)."""
    props: dict[str, Any] = {}
    node = resolve_schema_ref(schema)
    if isinstance(node.get("properties"), Mapping):
        props.update(dict(node["properties"]))
    for item in node.get("allOf") or []:
        if isinstance(item, Mapping):
            props.update(_collect_properties(item))
    return props


def _required_keys(schema: Mapping[str, Any]) -> set[str]:
    node = resolve_schema_ref(schema)
    required: set[str] = set()
    for name in node.get("required") or []:
        required.add(str(name))
    for item in node.get("allOf") or []:
        if isinstance(item, Mapping):
            required |= _required_keys(item)
    # anyOf required alternatives (e.g. Point number | point_number)
    for item in node.get("anyOf") or []:
        if isinstance(item, Mapping):
            for name in item.get("required") or []:
                required.add(str(name))
    return required


def _type_label(schema: Mapping[str, Any]) -> str:
    ref = str(schema.get("$ref") or "")
    if ref.endswith("positiveNumeric"):
        return "number >0"
    if ref.endswith("integerish"):
        return "integer"
    if ref.endswith("numeric"):
        return "number"
    if ref.endswith("nonemptyString"):
        return "string"
    node = resolve_schema_ref(schema)
    if "const" in node:
        return f"const {node['const']!r}"
    if "enum" in node:
        return "enum"
    t = node.get("type")
    if isinstance(t, list):
        return "|".join(str(x) for x in t)
    if t:
        return str(t)
    if node.get("oneOf") or node.get("anyOf"):
        # Numeric unions (common.json) after full resolve.
        if node.get("x-jarvis-numeric") or any(
            isinstance(item, Mapping) and item.get("x-jarvis-numeric")
            for item in (node.get("oneOf") or [])
        ):
            return "number"
        return "union"
    return "any"


def _key_entries(schema: Mapping[str, Any]) -> list[dict[str, Any]]:
    props = _collect_properties(schema)
    required = _required_keys(schema)
    keys: list[dict[str, Any]] = []
    for name, raw in props.items():
        raw_map = raw if isinstance(raw, Mapping) else {}
        node = resolve_schema_ref(raw_map)
        desc = str(raw_map.get("description") or node.get("description") or "")
        keys.append(
            {
                "name": name,
                "type": _type_label(raw_map),
                "required": name in required,
                "default": node.get("default"),
                "aliases": [],
                "description": desc,
                "minimum": node.get("minimum"),
                "exclusive": bool(node.get("exclusiveMinimum") is not None)
                or "positiveNumeric" in str(raw_map.get("$ref") or ""),
            }
        )
    return keys


def _bounds_schema_for_method(method: str) -> tuple[dict[str, Any], str]:
    manifest = schema_manifest()
    methods = dict(manifest.get("sampling_methods") or {})
    # Case-insensitive method match
    uri = methods.get(method)
    if uri is None:
        for name, target in methods.items():
            if name.casefold() == method.casefold():
                method = name
                uri = target
                break
    if uri is None:
        raise KeyError(method)
    root = schema_by_id(str(uri))
    props = _collect_properties(root)
    bounds = props.get("Bounds")
    if isinstance(bounds, Mapping):
        return resolve_schema_ref(bounds), method
    return {}, method


def _method_page(method: str) -> dict[str, Any]:
    from jarvishep2.distributor import Distributor

    manifest = schema_manifest()
    methods = dict(manifest.get("sampling_methods") or {})
    uri = None
    canonical = method
    for name, target in methods.items():
        if name.casefold() == method.casefold():
            canonical = name
            uri = target
            break
    if uri is None:
        raise KeyError(method)
    root = schema_by_id(str(uri))
    status = _schema_status(root)
    bounds_schema, _ = _bounds_schema_for_method(canonical)
    keys = _key_entries(bounds_schema) if status == "stable" else []
    try:
        from jarvishep2.distributor import STATELESS_METHODS, Distributor as Dist

        resume = Dist.get_resume_status(canonical)
        is_stateless = canonical in STATELESS_METHODS
    except Exception:
        resume = "unknown"
        is_stateless = False

    summary = str(root.get("description") or f"Sampling.Method = {canonical}")
    if status == "unstable":
        summary = (
            f"STATUS: not finalised — {canonical} Bounds knobs are not locked yet. "
            "Jarvis validate does not reject unknown Bounds keys for this method."
        )
    page = {
        "path": "$.Sampling.Bounds",
        "context": {"method": canonical},
        "zone": str(bounds_schema.get("x-jarvis-zone") or root.get("x-jarvis-zone") or "delegated"),
        "status": status,
        "summary": summary,
        "keys": keys,
        "examples": [str(root.get("x-jarvis-example"))] if root.get("x-jarvis-example") else [],
        "diagnostics": list(_METHOD_DIAGNOSTICS.get(canonical, [])),
        "see_also": [
            "Jarvis man sampler",
            "Jarvis man yaml Sampling.Variables",
            f"Jarvis man example {canonical.casefold()}",
        ],
        "further_reading": None,
        "capabilities": {
            "stateless": is_stateless,
            "resume": resume,
        },
    }
    return page


def _sampler_index() -> dict[str, Any]:
    from jarvishep2.distributor import STATELESS_METHODS, Distributor

    manifest = schema_manifest()
    rows = []
    for name in sorted(manifest.get("sampling_methods") or {}):
        root = schema_by_id(str(manifest["sampling_methods"][name]))
        status = _schema_status(root)
        rows.append(
            {
                "method": name,
                "status": status,
                "stateless": name in STATELESS_METHODS,
                "resume": Distributor.get_resume_status(name),
                "summary": str(root.get("description") or ""),
            }
        )
    return {
        "path": "$.Sampling.Method",
        "context": {},
        "zone": "closed",
        "status": "stable",
        "summary": "Sampling methods registered in the schema catalog.",
        "keys": [],
        "methods": rows,
        "examples": [],
        "diagnostics": [],
        "see_also": ["Jarvis man sampler Bridson", "Jarvis man yaml Sampling"],
        "further_reading": None,
    }


def _root_overview() -> dict[str, Any]:
    root = schema_by_id(str(schema_manifest()["root"]))
    props = _collect_properties(root)
    blocks = []
    for name in _ROOT_BLOCKS:
        node = props.get(name)
        if not isinstance(node, Mapping):
            blocks.append({"name": name, "zone": "?", "description": ""})
            continue
        resolved = resolve_schema_ref(node)
        blocks.append(
            {
                "name": name,
                "zone": str(resolved.get("x-jarvis-zone") or node.get("x-jarvis-zone") or "?"),
                "description": str(resolved.get("description") or ""),
            }
        )
    return {
        "path": "$",
        "context": {},
        "zone": str(root.get("x-jarvis-zone") or "closed"),
        "status": "stable",
        "summary": "Jarvis-HEP V2 task card top-level blocks.",
        "keys": [
            {
                "name": b["name"],
                "type": "object",
                "required": False,
                "default": None,
                "aliases": [],
                "description": f"zone={b['zone']}. {b['description']}".strip(),
            }
            for b in blocks
        ],
        "examples": [],
        "diagnostics": ["JV2-SCH-001"],
        "see_also": [
            "Jarvis man sampler",
            "Jarvis man calculator",
            "Jarvis man operas",
            "Jarvis man tokens",
        ],
        "further_reading": None,
        "blocks": blocks,
    }


def _path_page(path: str) -> dict[str, Any]:
    """Resolve a YAML path against the root card schema (best-effort)."""
    norm = normalize_yaml_path(path)
    parts = [p for p in norm[2:].split(".") if p] if norm.startswith("$.") else []
    # Special case Sampling.Bounds without method context
    if len(parts) >= 2 and parts[0].casefold() == "sampling" and parts[1].casefold() == "bounds":
        return {
            "path": "$.Sampling.Bounds",
            "context": {},
            "zone": "delegated",
            "status": "stable",
            "summary": (
                "Method-specific sampler knobs. Shape depends on Sampling.Method; "
                "use Jarvis man sampler <Method> for the closed key list."
            ),
            "keys": [],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man sampler", "Jarvis man sampler Bridson"],
            "further_reading": None,
        }
    if len(parts) >= 2 and parts[0].casefold() == "sampling" and parts[1].casefold() == "variables":
        from jarvishep2.contracts.variables import DISTRIBUTION_TYPES, PARAMS_REQUIRED

        keys = []
        for dtype in sorted(DISTRIBUTION_TYPES):
            req = sorted(PARAMS_REQUIRED.get(dtype, ()))
            keys.append(
                {
                    "name": f"distribution.type={dtype}",
                    "type": "object",
                    "required": True,
                    "default": None,
                    "aliases": [],
                    "description": "Required parameters: " + ", ".join(req),
                }
            )
        return {
            "path": "$.Sampling.Variables",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": "List of sampling variables; each item has name + distribution.",
            "keys": keys,
            "examples": [
                "Variables:\n  - name: x\n    distribution:\n      type: Flat\n      parameters: {min: 0.0, max: 1.0}"
            ],
            "diagnostics": ["JV2-VAR-001"],
            "see_also": ["Jarvis man sampler Bridson"],
            "further_reading": None,
        }

    # Walk root properties for first-level block
    root = schema_by_id(str(schema_manifest()["root"]))
    props = _collect_properties(root)
    if not parts:
        return _root_overview()
    head = parts[0]
    match = None
    for name in props:
        if name.casefold() == head.casefold():
            match = name
            break
    if match is None:
        raise KeyError(path)
    node = resolve_schema_ref(props[match])
    # Descend roughly into properties / items
    cursor = node
    trail = [match]
    for part in parts[1:]:
        bare = part.replace("[]", "")
        cursor = resolve_schema_ref(cursor)
        sub = _collect_properties(cursor)
        items = cursor.get("items")
        if bare.casefold() == "modules" and "Modules" in sub:
            cursor = resolve_schema_ref(sub["Modules"])
            trail.append("Modules")
            if isinstance(cursor.get("items"), Mapping):
                cursor = resolve_schema_ref(cursor["items"])
            continue
        # items schema for arrays
        if isinstance(items, Mapping) and not sub:
            cursor = resolve_schema_ref(items)
            sub = _collect_properties(cursor)
        hit = None
        for name in sub:
            if name.casefold() == bare.casefold() or name.casefold() == part.casefold():
                hit = name
                break
        if hit is None and isinstance(items, Mapping):
            cursor = resolve_schema_ref(items)
            sub = _collect_properties(cursor)
            for name in sub:
                if name.casefold() == bare.casefold():
                    hit = name
                    break
        if hit is None:
            # $defs on calculators
            defs = cursor.get("$defs") if isinstance(cursor.get("$defs"), Mapping) else {}
            for name, child in defs.items():
                if name.casefold() == bare.casefold():
                    hit = name
                    cursor = resolve_schema_ref(child)
                    trail.append(name)
                    break
            else:
                raise KeyError(path)
            continue
        cursor = resolve_schema_ref(sub[hit])
        trail.append(hit)

    keys = _key_entries(cursor)
    return {
        "path": "$." + ".".join(trail),
        "context": {},
        "zone": str(cursor.get("x-jarvis-zone") or "delegated"),
        "status": _schema_status(cursor),
        "summary": str(cursor.get("description") or f"YAML path {'.$'.join([''] + trail)}"),
        "keys": keys,
        "examples": [str(cursor["x-jarvis-example"])] if cursor.get("x-jarvis-example") else [],
        "diagnostics": [],
        "see_also": ["Jarvis man yaml"],
        "further_reading": None,
    }


def _io_format_page(direction: str, fmt: str) -> dict[str, Any]:
    manifest = schema_manifest()
    io_map = dict((manifest.get("io") or {}).get(direction) or {})
    # Runtime format list (A10): Portal may expose formats beyond the catalog.
    runtime_names: set[str] = set()
    try:
        from jarvishep2.io_portal import available_io_formats

        runtime_names = {str(x) for x in (available_io_formats(direction) or [])}
    except Exception:
        runtime_names = set(io_map)

    uri = None
    canonical = fmt
    for name, target in io_map.items():
        if name.casefold() == fmt.casefold():
            canonical = name
            uri = target
            break
    if uri is None:
        for name in runtime_names:
            if name.casefold() == fmt.casefold():
                canonical = name
                break

    keys: list[dict[str, Any]] = []
    example = None
    summary = f"{direction} format {canonical}"
    further = None
    status = "stable"
    if uri is None:
        summary = (
            f"Format {canonical!r} is available from Portal at runtime, but this Jarvis "
            "build has no bundled field schema for it (field manual missing in this build)."
        )
        further = {
            "command": f"Jarvis portal man {canonical}",
            "topic": "runtime adapter behaviour (supplementary; not a substitute for YAML keys)",
        }
    else:
        schema = schema_by_id(str(uri))
        keys = _key_entries(schema)
        # Provide minimal fallback descriptions when schema prose is thin.
        _IO_FALLBACK = {
            "name": "Logical name of this IO slot (referenced by actions).",
            "path": "File path relative to execution.path (tokens allowed).",
            "type": f"Portal format name ({canonical}).",
            "save": "When true, keep the file after the sample finishes.",
            "actions": "Ordered read/write actions for this file.",
        }
        for key in keys:
            if not key.get("description"):
                key["description"] = _IO_FALLBACK.get(str(key.get("name")), "")
        example = schema.get("x-jarvis-example")
        summary = str(schema.get("description") or summary)
    return {
        "path": f"$.Calculators.Modules[].execution.{direction}[]",
        "context": {"type": canonical, "direction": direction},
        "zone": "delegated",
        "status": status,
        "summary": summary,
        "keys": keys,
        "examples": [str(example)] if example else [],
        "diagnostics": ["JV2-SCH-002"],
        "see_also": ["Jarvis man calculator execution", "Jarvis man calculator"],
        "further_reading": further,
        "runtime_formats": sorted(runtime_names),
    }


def _calculator_page(topic: str | None, *, io_type: str | None = None) -> dict[str, Any]:
    if topic is None:
        return {
            "path": "$.Calculators",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "Calculator modules: write inputs (Portal IO) -> run commands -> "
                "read outputs -> optional Operas modules."
            ),
            "keys": [
                {"name": t, "type": "topic", "required": False, "default": None, "aliases": [], "description": ""}
                for t in _CALC_TOPICS
            ],
            "examples": [],
            "diagnostics": [],
            "see_also": [
                "Jarvis man calculator module",
                "Jarvis man calculator execution",
                "Jarvis man calculator modes",
                "Jarvis man calculator pools",
                "Jarvis man operas",
                "Jarvis man tokens",
            ],
            "further_reading": None,
            "chain": [
                "1. Calculators.Modules[] declare name/path/execution",
                "2. execution.input[] writes files via Portal formats (JSON, SLHA, ...)",
                "3. execution.commands[] run the external program",
                "4. execution.output[] reads files back into observables",
                "5. Operas.Modules[] consume observables (optional)",
            ],
        }
    topic_l = topic.casefold()
    if topic_l == "execution":
        if io_type:
            # Prefer input schema; fall back to output
            try:
                return _io_format_page("input", io_type)
            except Exception:
                return _io_format_page("output", io_type)
        # List formats from runtime + catalog
        try:
            from jarvishep2.io_portal import available_io_formats

            input_fmts = sorted(available_io_formats("input") or [])
            output_fmts = sorted(available_io_formats("output") or [])
        except Exception:
            manifest = schema_manifest()
            input_fmts = sorted((manifest.get("io") or {}).get("input") or {})
            output_fmts = sorted((manifest.get("io") or {}).get("output") or {})
        return {
            "path": "$.Calculators.Modules[].execution",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "execution block: path, commands, input[], output[]. "
                "Use --type FMT for a full field table of one Portal format."
            ),
            "keys": [
                {
                    "name": "path",
                    "type": "string",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Working directory for commands (tokens allowed).",
                },
                {
                    "name": "commands",
                    "type": "array",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Shell commands or {cmd,cwd} objects to run.",
                },
                {
                    "name": "input[]",
                    "type": "array",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Portal input file specs. Formats: " + ", ".join(input_fmts),
                },
                {
                    "name": "output[]",
                    "type": "array",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Portal output file specs. Formats: " + ", ".join(output_fmts),
                },
            ],
            "examples": [
                "execution:\n  path: &J/calculators/runtime\n  commands: [python run.py]\n  input:\n    - name: in\n      path: input.json\n      type: JSON\n      actions:\n        - type: Dump\n          variables:\n            - {name: x, expression: x}\n  output:\n    - name: out\n      path: output.json\n      type: JSON\n      actions:\n        - type: Load\n          variables:\n            - {name: y, entry: y}"
            ],
            "diagnostics": [],
            "see_also": [
                "Jarvis man calculator execution --type JSON",
                "Jarvis man operas",
            ],
            "further_reading": None,
            "input_formats": input_fmts,
            "output_formats": output_fmts,
        }
    if topic_l == "module":
        calc = schema_by_id("https://jarvis-hep.org/schema/v2/core/calculators.json")
        fields = resolve_schema_ref((calc.get("$defs") or {}).get("moduleFields") or {})
        return {
            "path": "$.Calculators.Modules[]",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": "One calculator module entry under Calculators.Modules[].",
            "keys": _key_entries(fields),
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man calculator execution", "Jarvis man calculator modes"],
            "further_reading": None,
        }
    if topic_l == "modes":
        return {
            "path": "$.Calculators.Modules[].modes",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "Multimode modules share one physical pack; each mode has its own "
                "execution and is referenced with dot-qualified names."
            ),
            "keys": [
                {
                    "name": "modes[]",
                    "type": "array",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Alternative execution profiles; mutually exclusive with top-level execution.",
                }
            ],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man calculator module"],
            "further_reading": None,
        }
    if topic_l == "pools":
        return {
            "path": "$.Calculators.Pools",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": "Named integer pools for concurrent calculator capacity.",
            "keys": [
                {
                    "name": "<pool-name>",
                    "type": "integer",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Max concurrent slots for modules assigned to this pool.",
                }
            ],
            "examples": ["Pools:\n  susy: 4"],
            "diagnostics": [],
            "see_also": ["Jarvis man calculator"],
            "further_reading": None,
        }
    raise KeyError(topic)


def _operas_page() -> dict[str, Any]:
    operas = schema_by_id("https://jarvis-hep.org/schema/v2/core/operas.json")
    modules = (operas.get("properties") or {}).get("Modules")
    item: Mapping[str, Any] = {}
    if isinstance(modules, Mapping) and isinstance(modules.get("items"), Mapping):
        item = modules["items"]
    return {
        "path": "$.Operas.Modules[]",
        "context": {},
        "zone": str(item.get("x-jarvis-zone") or "closed"),
        "status": "stable",
        "summary": (
            "Operas.Modules[] writes derived observables. operator is a registered "
            "Jarvis-Operas name (see Jarvis operas info for signatures)."
        ),
        "keys": _key_entries(item),
        "examples": [
            "Operas:\n  Modules:\n    - name: chi2\n      operator: math.add\n      input:\n        - {name: a, expression: x}\n        - {name: b, expression: y}\n      output:\n        - {name: s, expression: result}"
        ],
        "diagnostics": ["JV2-OPR-002"],
        "see_also": ["Jarvis man calculator", "Jarvis man tokens"],
        "further_reading": {
            "command": "Jarvis operas info math.add",
            "topic": "what a specific operator does",
        },
    }


def _tokens_page() -> dict[str, Any]:
    keys = [
        {"name": "&J/...", "type": "token", "required": False, "default": None, "aliases": [], "description": "Project root (standalone project or cwd)."},
        {"name": "@Sdir", "type": "token", "required": False, "default": None, "aliases": [], "description": "Per-sample working directory under SAMPLE/."},
        {"name": "@PackID", "type": "token", "required": False, "default": None, "aliases": [], "description": "Physical calculator pack id for multimode affinity."},
        {"name": "@SampleID", "type": "token", "required": False, "default": None, "aliases": [], "description": "Current sample uuid."},
        {"name": "${LibDeps:NAME}", "type": "token", "required": False, "default": None, "aliases": [], "description": "Path from LibDeps-built shared libraries."},
        {"name": "${Scan:KEY}", "type": "token", "required": False, "default": None, "aliases": [], "description": "Value from the Scan block (unknown keys may fall back silently)."},
        {"name": "@{ROOT path}", "type": "token", "required": False, "default": None, "aliases": [], "description": "ROOT file path helper where supported."},
    ]
    return {
        "path": "$",
        "context": {"domain": "tokens"},
        "zone": "open",
        "status": "stable",
        "summary": "Path and expression tokens used in task YAML.",
        "keys": keys,
        "examples": ["path: &J/calculators/runtime\ncommand: ${LibDeps:flexiblesusy}/bin/Spectrum"],
        "diagnostics": [],
        "see_also": ["Jarvis man calculator execution", "Jarvis man yaml LibDeps"],
        "further_reading": None,
    }


def _example_page(name: str | None) -> dict[str, Any]:
    root = Path(__file__).with_name("project_template")
    if name is None:
        return {
            "path": "$",
            "context": {"domain": "example"},
            "zone": "open",
            "status": "stable",
            "summary": "Minimal task cards shipped with the package (same files as project create).",
            "keys": [
                {
                    "name": key,
                    "type": "example",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": rel,
                }
                for key, rel in _EXAMPLE_CATALOG.items()
            ],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man example bridson"],
            "further_reading": None,
        }
    key = name.casefold()
    rel = _EXAMPLE_CATALOG.get(key)
    if rel is None:
        raise KeyError(name)
    path = root / rel
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    return {
        "path": "$",
        "context": {"example": key},
        "zone": "open",
        "status": "stable",
        "summary": f"Shipped template card: {rel}",
        "keys": [],
        "examples": [text],
        "diagnostics": [],
        "see_also": ["Jarvis man yaml", "Jarvis validate"],
        "further_reading": None,
        "body": text,
    }


def _center_page() -> dict[str, Any]:
    return {
        "path": "$",
        "context": {"domain": "center"},
        "zone": "open",
        "status": "stable",
        "summary": "Jarvis man — YAML writing manuals for humans and coding agents.",
        "keys": [
            {"name": "yaml", "type": "domain", "required": False, "default": None, "aliases": [], "description": "Task-card structure and path lookup"},
            {"name": "sampler", "type": "domain", "required": False, "default": None, "aliases": [], "description": "Sampling.Method Bounds knobs and capabilities"},
            {"name": "calculator", "type": "domain", "required": False, "default": None, "aliases": [], "description": "Calculators.Modules / execution / modes / pools"},
            {"name": "operas", "type": "domain", "required": False, "default": None, "aliases": [], "description": "Operas.Modules[] YAML shape"},
            {"name": "tokens", "type": "domain", "required": False, "default": None, "aliases": [], "description": "&J / @Sdir / ${LibDeps:} / ${Scan:}"},
            {"name": "example", "type": "domain", "required": False, "default": None, "aliases": [], "description": "Emit a minimal runnable card"},
        ],
        "examples": [],
        "diagnostics": [],
        "see_also": [
            "Jarvis man yaml",
            "Jarvis man sampler Bridson",
            "Jarvis man calculator execution --type JSON",
            "Jarvis validate TASK.yaml",
        ],
        "further_reading": {
            "command": "Jarvis portal man",
            "topic": "Portal adapter runtime behaviour (not YAML field lists)",
        },
    }


def resolve_man_request(
    tokens: list[str],
    *,
    io_type: str | None = None,
    code: str | None = None,
) -> dict[str, Any]:
    """Resolve argv tokens after ``man`` into a page dict."""
    if code:
        from jarvishep2.man_codes import man_target_for_code, resolve_tokens_for_code

        entry = man_target_for_code(str(code))
        if entry is None:
            raise KeyError(code)
        tokens = resolve_tokens_for_code(str(code))
        # Fall through to normal token resolution.

    args = [str(t) for t in tokens if str(t).strip()]
    # Optional leading yaml prefix (A4)
    if args and args[0].casefold() == "yaml" and len(args) > 1 and args[1].casefold() in _DOMAINS - {"yaml"}:
        args = args[1:]

    if not args:
        return _center_page()

    head = args[0].casefold()
    rest = args[1:]

    if head == "yaml":
        if not rest:
            return _root_overview()
        try:
            return _path_page(".".join(rest))
        except KeyError:
            choices = list(_ROOT_BLOCKS) + ["Bounds", "Variables", "Modules", "Method"]
            suggestions = get_close_matches(rest[0], choices, n=5, cutoff=0.4)
            raise KeyError(
                f"{rest[0]} (did you mean: {', '.join(suggestions) or 'nothing close'})"
            ) from None
    if head == "sampler":
        if not rest:
            return _sampler_index()
        return _method_page(rest[0])
    if head == "calculator":
        topic = rest[0] if rest else None
        return _calculator_page(topic, io_type=io_type)
    if head == "operas":
        return _operas_page()
    if head == "tokens":
        return _tokens_page()
    if head == "example":
        return _example_page(rest[0] if rest else None)

    # Bare path or unknown domain → try YAML path
    try:
        return _path_page(".".join(args))
    except KeyError:
        # did-you-mean for top-level blocks / domains
        choices = list(_DOMAINS) + list(_ROOT_BLOCKS) + list(schema_manifest().get("sampling_methods") or {})
        suggestions = get_close_matches(args[0], choices, n=5, cutoff=0.5)
        raise KeyError(f"{args[0]} (did you mean: {', '.join(suggestions) or 'nothing close'})") from None


def _print_page(page: Mapping[str, Any], *, as_json: bool) -> None:
    if as_json:
        # Stable public keys only (+ optional extras for indexes)
        payload = {
            "path": page.get("path"),
            "context": page.get("context") or {},
            "zone": page.get("zone"),
            "status": page.get("status"),
            "summary": page.get("summary"),
            "keys": page.get("keys") or [],
            "examples": page.get("examples") or [],
            "diagnostics": page.get("diagnostics") or [],
            "see_also": page.get("see_also") or [],
            "further_reading": page.get("further_reading"),
        }
        for extra in (
            "methods",
            "blocks",
            "chain",
            "capabilities",
            "input_formats",
            "output_formats",
            "runtime_formats",
            "body",
        ):
            if extra in page:
                payload[extra] = page[extra]
        json.dump(payload, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return

    console = _console()
    title = str(page.get("path") or "Jarvis man")
    ctx = page.get("context") or {}
    if ctx.get("method"):
        title = f"{title} · {ctx['method']}"
    if ctx.get("type"):
        title = f"{title} · type={ctx['type']}"
    status = str(page.get("status") or "stable")
    summary = str(page.get("summary") or "")
    if status == "unstable":
        summary = f"[bold red]STATUS: not finalised[/]\n{summary}"
    console.print(Panel(summary, title=title, box=box.ROUNDED, border_style="cyan"))

    if page.get("methods"):
        table = Table(title="Methods", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("Method", style="bold")
        table.add_column("Status")
        table.add_column("Stateless")
        table.add_column("Resume")
        table.add_column("Summary", overflow="fold")
        for row in page["methods"]:
            table.add_row(
                str(row.get("method")),
                str(row.get("status")),
                "yes" if row.get("stateless") else "no",
                str(row.get("resume")),
                str(row.get("summary") or ""),
            )
        console.print(table)

    if page.get("chain"):
        console.print(Panel("\n".join(str(x) for x in page["chain"]), title="Calculator chain", box=box.ROUNDED))

    if page.get("capabilities"):
        cap = page["capabilities"]
        console.print(
            Panel(
                f"stateless: {cap.get('stateless')}\nresume: {cap.get('resume')}",
                title="Capabilities",
                box=box.ROUNDED,
            )
        )

    keys = list(page.get("keys") or [])
    if keys:
        table = Table(title="Keys", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("Name", style="bold green")
        table.add_column("Type")
        table.add_column("Required")
        table.add_column("Description", overflow="fold")
        for key in keys:
            table.add_row(
                str(key.get("name")),
                str(key.get("type")),
                "yes" if key.get("required") else "optional",
                str(key.get("description") or "—"),
            )
        console.print(table)
    elif status == "unstable":
        console.print(
            Panel(
                "No Keys table: Bounds vocabulary is not finalised for this method.",
                title="Keys",
                border_style="red",
                box=box.ROUNDED,
            )
        )

    examples = list(page.get("examples") or [])
    if page.get("body"):
        examples = [str(page["body"])]
    if examples:
        body = examples[0]
        console.print(Panel(body, title="YAML", box=box.ROUNDED, border_style="green"))

    diags = list(page.get("diagnostics") or [])
    if diags:
        console.print(
            Panel("\n".join(diags), title="Diagnostics that fire here", box=box.ROUNDED)
        )

    see = list(page.get("see_also") or [])
    further = page.get("further_reading")
    if further:
        see.append(f"Further reading: {further.get('command')} ({further.get('topic')})")
    if see:
        console.print(Panel("\n".join(str(x) for x in see), title="See also", box=box.ROUNDED))


def build_man_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="Jarvis man",
        description="YAML writing manuals (English). Domains: yaml sampler calculator operas tokens example.",
    )
    parser.add_argument("tokens", nargs="*", help="Domain and/or path tokens")
    parser.add_argument("--json", action="store_true", help="Machine-readable page")
    parser.add_argument("--code", metavar="JV2-XXX-NNN", help="Jump from a diagnostic code")
    parser.add_argument(
        "--type",
        dest="io_type",
        metavar="FMT",
        help="Portal format for calculator execution field tables (e.g. JSON)",
    )
    return parser


def dispatch_man(argv: list[str] | None = None) -> int:
    parser = build_man_parser()
    args = parser.parse_args(list(argv or []))
    try:
        page = resolve_man_request(
            list(args.tokens or []),
            io_type=getattr(args, "io_type", None),
            code=getattr(args, "code", None),
        )
    except KeyError as exc:
        print(f"Unknown man topic: {exc}", file=sys.stderr)
        print("Try: Jarvis man   or   Jarvis man sampler Bridson", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - defensive
        print(f"Jarvis man failed: {exc}", file=sys.stderr)
        return 1
    _print_page(page, as_json=bool(args.json))
    return 0


__all__ = [
    "build_man_parser",
    "dispatch_man",
    "normalize_yaml_path",
    "resolve_man_request",
]
