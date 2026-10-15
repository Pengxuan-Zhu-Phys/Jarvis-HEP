"""Read-only schema and shipped YAML exports for coding agents.

These adapters print installed package artifacts. They never create a project,
write a card, execute a template, or retrieve a schema from its URI.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any
from urllib.parse import unquote, urldefrag

import yaml

from jarvishep2 import __version__
from jarvishep2.cardload.task_schema import SCHEMA_DIR, schema_manifest
from jarvishep2.cli.man import resolve_man_request


CONTRACT_VERSION = 1
SCHEMA_SCOPE = (
    "bundled task-card schemas; runtime/delegated/plugin metadata require man and validate"
)
_TEMPLATE_ROOT = Path(__file__).resolve().parents[1] / "project_template"
# Facts about shipped examples, not a classifier for arbitrary dotted names.
_SHIPPED_IMPORT_OPERATORS = {"jarvishep2.testing.eggbox.eggbox2d_numpy"}
_DOCUMENTED_REGISTERED_OPERATORS = {"helper.eggbox2d"}


class _UsageError(ValueError):
    pass


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise _UsageError(message)


def _emit_json(payload: Any) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False))


def _emit_error(
    parser: argparse.ArgumentParser, message: str, *, as_json: bool, exit_code: int
) -> int:
    if as_json:
        _emit_json({
            "ok": False,
            "contract_version": CONTRACT_VERSION,
            "command": parser.prog,
            "status": "error",
            "error": {"code": "JV2-AGENT-001", "message": message},
            "usage": parser.format_usage().strip(),
        })
    else:
        print(f"{parser.prog}: {message}", file=sys.stderr)
        print(parser.format_usage().strip(), file=sys.stderr)
    return exit_code


def schema_bundle() -> dict[str, Any]:
    """Return all manifest-selected documents, keyed by their local schema IDs."""
    manifest = schema_manifest()
    schemas: dict[str, Any] = {}
    for relative_path in manifest["schema_files"]:
        document = json.loads((SCHEMA_DIR / relative_path).read_text(encoding="utf-8"))
        schema_id = document["$id"]
        if schema_id in schemas:
            raise RuntimeError(f"Duplicate bundled schema ID: {schema_id}")
        schemas[schema_id] = document
    if manifest["root"] not in schemas:
        raise RuntimeError("Bundled root schema is missing from the manifest")
    return {
        "contract_version": CONTRACT_VERSION,
        "package_version": __version__,
        "root_schema_id": manifest["root"],
        "manifest": manifest,
        "schemas": schemas,
        "scope": SCHEMA_SCOPE,
    }


def _schema_selection(bundle: Mapping[str, Any], reference: str) -> Any:
    """Select an exact local document or RFC 6901 URI-fragment JSON Pointer."""
    try:
        schema_id, fragment = urldefrag(reference)
    except ValueError:
        raise _UsageError(f"Invalid schema reference: {reference}") from None
    if schema_id not in bundle["schemas"]:
        raise _UsageError(f"Unknown bundled schema ID: {schema_id}. List IDs with Jarvis agent schema --list --json.")
    document = bundle["schemas"][schema_id]
    if not fragment:
        return document
    if re.search(r"%(?![0-9a-fA-F]{2})", fragment):
        raise _UsageError(f"Invalid percent escape in schema reference: {reference}")
    try:
        pointer = unquote(fragment, errors="strict")
    except UnicodeError:
        raise _UsageError(f"Invalid UTF-8 JSON Pointer: {reference}") from None
    if not pointer.startswith("/"):
        raise _UsageError(f"Schema fragment must be a JSON Pointer beginning with /: {reference}")
    node: Any = document
    for encoded_part in pointer[1:].split("/"):
        if re.search(r"~(?:[^01]|$)", encoded_part):
            raise _UsageError(f"Invalid JSON Pointer escape: {reference}")
        part = encoded_part.replace("~1", "/").replace("~0", "~")
        try:
            if isinstance(node, Mapping):
                node = node[part]
            elif isinstance(node, list) and re.fullmatch(r"0|[1-9][0-9]*", part):
                node = node[int(part)]
            else:
                raise KeyError(part)
        except (KeyError, IndexError, ValueError):
            raise _UsageError(f"Unknown bundled schema fragment: {reference}") from None
    return node


def _schema_index(bundle: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "contract_version": CONTRACT_VERSION,
        "package_version": __version__,
        "root_schema_id": bundle["root_schema_id"],
        "scope": SCHEMA_SCOPE,
        "schema_ids": list(bundle["schemas"]),
        "sampling_methods": bundle["manifest"]["sampling_methods"],
        "io": bundle["manifest"]["io"],
        "select_argv": ["Jarvis", "agent", "schema", "--id", "<SchemaID-or-JSONPointerRef>", "--json"],
        "note": "Select an exact local $id or append #/... for a JSON Pointer fragment. All IDs and references resolve from this installed bundle; never fetch their URLs.",
    }


def dispatch_schema(argv: list[str] | None = None) -> int:
    """Print a local schema bundle, small index, or exact document/fragment."""
    parser = _Parser(
        prog="Jarvis agent schema",
        description=(
            "Export installed task-card JSON schemas without network access. "
            "JSON is the default. Structural schemas still require man and validate."
        ),
        add_help=False,
        allow_abbrev=False,
    )
    parser.add_argument("--json", action="store_true", help="JSON output (the default)")
    selector = parser.add_mutually_exclusive_group()
    selector.add_argument("--id", metavar="URI", help="Select an exact bundled $id, optionally followed by #/... JSON Pointer")
    selector.add_argument("--list", action="store_true", help="Print a small local schema ID/method/IO index")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help")
    raw_argv = list(argv or [])
    try:
        args = parser.parse_args(raw_argv)
        if args.help:
            if args.json:
                _emit_json({"contract_version": CONTRACT_VERSION, "help": parser.format_help()})
            else:
                parser.print_help()
            return 0
        bundle = schema_bundle()
        payload: Any = _schema_index(bundle) if args.list else bundle
        if args.id is not None:
            payload = _schema_selection(bundle, args.id)
        _emit_json(payload)
        return 0
    except _UsageError as exc:
        return _emit_error(parser, str(exc), as_json=True, exit_code=2)
    except Exception as exc:
        return _emit_error(parser, f"Schema export failed: {exc}", as_json=True, exit_code=1)


def _walk_values(node: Any, path: str = "$"):
    if isinstance(node, Mapping):
        for key, value in node.items():
            yield from _walk_values(value, f"{path}.{key}")
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from _walk_values(value, f"{path}[{index}]")
    else:
        yield path, node


def _prerequisites(document: Mapping[str, Any], artifact: Mapping[str, Any]) -> list[dict[str, Any]]:
    paths = [(path, value) for path, value in _walk_values(document)
             if isinstance(value, str) and value.startswith("&J/")]
    requirements: list[dict[str, Any]] = [{
        "kind": "project_scaffold",
        "required": bool(paths),
        "command": "Jarvis project create NAME",
        "argv": ["Jarvis", "project", "create", "NAME"],
        "note": (
            "Use an existing scaffold or create a new project only when authorized. "
            "The scaffold copies bin/, data/, and deps/. Save the card inside the project; "
            "&J/ resolves to its project root. This export does not create any files."
        ),
    }]
    missing = [name for name in ("Scan", "Sampling", "EnvReqs") if name not in document]
    if missing:
        requirements.append({
            "kind": "composition",
            "missing_root_blocks": missing,
            "alternative_root_blocks": ["Calculators", "Operas"],
            "note": (
                "This shipped sampling example is a fragment. Merge its Sampling block into "
                "a complete task card with Scan, EnvReqs, and at least one of Calculators or "
                "Operas, preserving the intended physics. "
                "Then run explicit Jarvis validate TASK.yaml --json."
            ),
            "manual_argv": ["Jarvis", "man", "yaml", "--json"],
            "base_template_argv": ["Jarvis", "agent", "template", "calculator", "--json"],
            "steps": [
                "Read a complete base card and its prerequisites. The calculator template is a synthetic demo; preserve the chosen Calculators/Operas and EnvReqs when adapting an existing scientific card.",
                "Replace the base document's entire root Sampling mapping with this fragment's Sampling mapping. Do not append a second Sampling key or concatenate two YAML documents.",
                "Review Scan.name, parameter ranges, likelihood and budget against the authorized intent. Save the composed mapping as <ComposedTask.yaml> inside its project.",
                "Validate the composed task, then follow diagnostics and validate again. Do not validate or run the fragment alone.",
            ],
            "operation": {"action": "replace_mapping", "source_path": "$.Sampling", "target_path": "$.Sampling"},
        })
    for path, value in paths:
        relative_path = value[3:]
        requirements.append({
            "kind": "project_path", "yaml_path": path, "value": value,
            "relative_path": relative_path,
            "provided_by_scaffold": (_TEMPLATE_ROOT / relative_path).exists(),
        })
    sampling = document.get("Sampling") or {}
    method = sampling.get("Method")
    if method:
        requirements.append({
            "kind": "sampler", "name": method,
            "manual_argv": ["Jarvis", "man", f"sampler.{method}", "--json"],
        })
    if str(method).casefold() == "csv":
        bounds = sampling.get("Bounds") or {}
        # These are direct identifier inputs of the actual shipped card. This
        # does not guess free variables of arbitrary expressions or operators.
        input_columns = sorted({
            value for path, value in _walk_values(document)
            if (path.endswith(".expression") or path.endswith(".entry"))
            and ".input[" in path and isinstance(value, str) and value.isidentifier()
        })
        bundled_header = None
        csv_path = bounds.get("path")
        if isinstance(csv_path, str) and csv_path.startswith("&J/"):
            bundled_csv = _TEMPLATE_ROOT / csv_path[3:]
            if bundled_csv.is_file():
                with bundled_csv.open(encoding=bounds.get("encoding") or "utf-8", newline="") as handle:
                    bundled_header = next(csv.reader(handle, delimiter=bounds.get("delimiter") or ","), [])
        requirements.append({
            "kind": "csv", "path": bounds.get("path"),
            "uuid_column": bounds.get("uuid_column"),
            "input_columns": input_columns,
            "input_columns_scope": "Direct identifier references in the shipped card's input bindings; no general expression or callable introspection.",
            "bundled_header": bundled_header,
            "note": "Supply the CSV at this path with the columns referenced by the task. The scaffold supplies data/points.csv.",
        })
    for module in (document.get("Calculators") or {}).get("Modules") or []:
        execution = module.get("execution") or {}
        requirements.append({
            "kind": "calculator", "name": module.get("name"), "path": module.get("path"),
            "execution_path": execution.get("path"), "commands": execution.get("commands") or [],
            "io_paths": [
                {"direction": direction, "path": entry.get("path"), "type": entry.get("type")}
                for direction in ("input", "output") for entry in execution.get(direction) or []
            ],
            "note": (
                "Ensure the declared paths and command executables are available. "
                "These shipped synthetic calculators use the project's bin/ and python3; "
                "the scaffold does not install an external calculator. @Sdir is a per-sample "
                "runtime directory. run/check execute these commands; template/validate do not."
            ),
            "manual_argv": ["Jarvis", "man", "calculator.execution", "--json"],
        })
    for module in (document.get("Operas") or {}).get("Modules") or []:
        operator = module.get("operator")
        row = {
            "kind": "operator", "name": module.get("name"), "operator": operator,
            "input": module.get("input") or [], "output": module.get("output") or [],
            "signature_verified": False,
            "manual_argv": ["Jarvis", "man", "operas", "--json"],
        }
        if operator in _SHIPPED_IMPORT_OPERATORS:
            row.update({
                "reference_kind": "python_import_path",
                "module": operator.rpartition(".")[0], "attribute": operator.rpartition(".")[2],
                "note": "This shipped operator is a Python import path, not a registered Operas catalog name. Jarvis operas info only queries registered names and cannot inspect this callable. Template metadata does not import it. The callable must be available and return the mapping declared by output bindings; inspect its source only within the authorized scope before adapting it.",
            })
        elif operator in _DOCUMENTED_REGISTERED_OPERATORS:
            row.update({
                "reference_kind": "registered_name",
                "query_argv": ["Jarvis", "operas", "info", operator],
                "note": "This name is documented in the Operas registry. Query its installed declaration before adapting it; metadata export does not execute the operator.",
            })
        else:
            row.update({
                "reference_kind": "registered_or_import_path",
                "registry_query_argv": ["Jarvis", "operas", "info", operator],
                "note": "Dotted syntax alone cannot distinguish a registered name from a Python import path. The registry query checks only registered names; a miss does not disprove an importable callable. Template metadata does not import user callables or claim their signatures are verified.",
            })
        requirements.append(row)
    requirements.append({
        "kind": "validation",
        "argv": ["Jarvis", "validate", "<ComposedTask.yaml>" if missing else str(artifact["path"]), "--json"],
        "note": "Validate the saved, composed task from its project. Validation success does not prove calculator execution or physical correctness. run/check require execution authorization.",
    })
    return requirements


def template_document(name: str) -> dict[str, Any]:
    """Read one exact catalog name through the same source used by Jarvis man."""
    catalog = resolve_man_request(["example"])["list_rows"]
    key = name.casefold()
    if key not in {row["name"] for row in catalog}:
        raise _UsageError(f"Unknown template: {name}. List names with Jarvis agent template --json.")
    page = resolve_man_request([f"example.{key}"])
    body = page.get("body")
    if not isinstance(body, str) or not body.strip():
        raise RuntimeError(f"Bundled template is missing or empty: {key}")
    document = yaml.safe_load(body)
    if not isinstance(document, Mapping):
        raise RuntimeError(f"Bundled template is not a YAML mapping: {key}")
    artifact = dict(page["artifact"])
    if not all(block in document for block in ("Scan", "Sampling", "EnvReqs")):
        artifact["kind"] = "sampling_fragment" if set(document) == {"Sampling"} else "task_card_fragment"
    return {
        "contract_version": CONTRACT_VERSION,
        "package_version": __version__, "name": key,
        "status": page.get("status"), "summary": page.get("summary"),
        "body": body, "artifact": artifact,
        "prerequisites": _prerequisites(document, artifact),
    }


def template_catalog() -> dict[str, Any]:
    rows = resolve_man_request(["example"])["list_rows"]
    return {
        "contract_version": CONTRACT_VERSION, "package_version": __version__,
        "templates": [{
            "name": row["name"],
            "artifact": template_document(row["name"])["artifact"],
            "argv": ["Jarvis", "agent", "template", row["name"], "--json"],
        } for row in rows],
        "note": "Read prerequisites before saving or composing a template. Sampling fragments need a complete task card. This command only prints installed artifacts.",
    }


def dispatch_template(argv: list[str] | None = None) -> int:
    """List templates or print exact YAML; JSON adds metadata and prerequisites."""
    parser = _Parser(
        prog="Jarvis agent template", add_help=False, allow_abbrev=False,
        description="List shipped templates, or print exact YAML without writing or executing it. Sampling fragments must be composed into complete task cards.",
    )
    parser.add_argument("name", nargs="?", metavar="NAME", help="Exact shipped template name")
    parser.add_argument("--json", action="store_true", help="Include metadata, artifact and prerequisites")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help")
    raw_argv = list(argv or [])
    as_json = any(token == "--json" or token.startswith("--json=") for token in raw_argv)
    try:
        args = parser.parse_args(raw_argv)
        if args.help:
            if args.json:
                _emit_json({"contract_version": CONTRACT_VERSION, "help": parser.format_help()})
            else:
                parser.print_help()
            return 0
        payload = template_catalog() if args.name is None else template_document(args.name)
        if args.json:
            _emit_json(payload)
        elif args.name is not None:
            sys.stdout.write(payload["body"])
        else:
            print("Shipped templates (Jarvis agent template NAME --json shows prerequisites):")
            for entry in payload["templates"]:
                artifact = entry["artifact"]
                print(f"  {entry['name']:<20} {artifact['kind']}: {artifact['path']}")
        return 0
    except _UsageError as exc:
        return _emit_error(parser, str(exc), as_json=as_json, exit_code=2)
    except Exception as exc:
        return _emit_error(parser, f"Template export failed: {exc}", as_json=as_json, exit_code=1)


__all__ = ["dispatch_schema", "dispatch_template", "schema_bundle", "template_catalog", "template_document"]
