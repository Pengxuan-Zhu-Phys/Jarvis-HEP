#!/usr/bin/env python3
"""SDK-owned sampler CLI and rendering, with HEP catalog data hooks."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.text import Text

from jarvishep2.run_outcome import EXIT_USAGE

_JSON_FLAG = "--json"
SAMPLER_COMMAND_HELP = {
    "list": "List loaded samplers and methods shipped with Jarvis-HEP",
    "info": "Show one loaded sampler or method shipped with Jarvis-HEP",
}


def dispatch_sampler(argv: list[str] | None = None) -> int:
    args = list(argv or [])
    from jarvishep2.sampler_sdk import SamplerSDKError, require_sampler_sdk
    try:
        require_sampler_sdk()
        from jarvis_hep_sampler.cli import CatalogExtension, main
        from jarvis_hep_sampler.errors import PluginError
    except (SamplerSDKError, ImportError) as exc:
        return _command_error(args, str(exc), type(exc).__name__)
    try:
        return int(main(
            args, prog="Jarvis sampler", command_help=SAMPLER_COMMAND_HELP,
            command_description=SAMPLER_COMMAND_HELP,
            catalog=CatalogExtension(_list_result, _info_result, (("jarvis_hep", "Jarvis-HEP"),)),
        ))
    except (PluginError, OSError) as exc:
        return _command_error(args, str(exc), type(exc).__name__)


def _report():
    from jarvis_hep_sampler.errors import RegistryError
    from jarvishep2.sampler_catalog import registry_report
    report = registry_report()
    if report.error is not None:
        raise RegistryError(report.error)
    return report


def _apply_entry(payload: dict[str, Any], entry: Any) -> None:
    if entry.status != "ok":
        payload["status"] = entry.status
    if entry.problems:
        payload["problems"] = list(entry.problems)
    if entry.repair_command:
        payload["repair_command"] = entry.repair_command


def _list_result(payload: dict[str, Any]) -> dict[str, Any]:
    from jarvishep2.sampler_catalog import SOURCE_USER
    report = _report()
    entries = {entry.method: entry for entry in report.entries if entry.source == SOURCE_USER}
    for item in payload["samplers"]:
        entry = entries.get(item["method"])
        if entry is not None:
            _apply_entry(item, entry)
    payload["jarvis_hep"] = _builtin_records()
    return payload


def _info_result(method: str, payload: dict[str, Any] | None) -> dict[str, Any] | None:
    from jarvishep2.sampler_catalog import SOURCE_USER
    report = _report()
    if payload is not None:
        entry = next((entry for entry in report.entries
                      if entry.method == payload["method"] and entry.source == SOURCE_USER), None)
        if entry is not None:
            _apply_entry(payload, entry)
        return payload
    spec = _internal_spec(method)
    if spec is None:
        return None
    payload = _internal_info(spec)
    conflicts = []
    for entry in report.matching(method):
        if entry.source == SOURCE_USER and entry.status == "name-conflict":
            item = {"method": entry.method, "source": entry.source, "status": entry.status}
            _apply_entry(item, entry)
            conflicts.append(item)
    if conflicts:
        payload["conflicts"] = conflicts
    return payload


def _command_error(args: list[str], message: str, kind: str) -> int:
    command, as_json = _command(args)
    path = f"sampler {command[0]}" if command else "sampler"
    hint = "Use Jarvis sampler <command> -h for usage."
    try:
        from jarvis_hep_sampler.cli import emit_error
    except ImportError:
        # The dependency itself is unavailable, so its renderer cannot run.
        if as_json:
            _print_json({
                "command": path, "error": {"type": kind, "message": message, "hint": hint},
                "ok": False,
            })
        else:
            Console(stderr=True, highlight=False).print(Panel(
                Group(Text(message), Text(""), Text(hint, style="dim")),
                title=Text("Jarvis sampler", style="bold"),
                box=box.ROUNDED, border_style="red", padding=(0, 1),
            ))
        return EXIT_USAGE
    return emit_error(
        "Jarvis sampler", message, command=path, error_type=kind,
        hint=hint, as_json=as_json, code=EXIT_USAGE,
    )


def _command(args: list[str]) -> tuple[list[str], bool]:
    """Drop ``--json`` so ``list`` and ``info`` match with the flag in any position."""
    return [arg for arg in args if arg != _JSON_FLAG], _JSON_FLAG in args


def _internal_spec(name: str):
    """A sampler shipped in Jarvis-HEP, or None when only jsampler knows it."""
    from jarvishep2.sampler_catalog import SOURCE_USER, canonical_name, spec_for

    exact = canonical_name(name)
    if exact is None:
        return None
    spec = spec_for(exact)
    if spec is None or spec.source == SOURCE_USER:
        return None
    return spec


def _display_source(path: str) -> str:
    """Package-relative path, so a shipped file does not widen the panel."""
    raw = Path(path)
    package = Path(__file__).resolve().parents[1]
    try:
        return str(raw.resolve().relative_to(package))
    except ValueError:
        return raw.name


def _profile(spec) -> str:
    from jarvishep2.sampler_catalog import SOURCE_BUILTIN_PLUGIN, plugin_record

    if spec.source == SOURCE_BUILTIN_PLUGIN:
        record = plugin_record(spec.name)
        return record.profile if record is not None else spec.family
    if spec.stateless and not spec.statistical:
        return "independent"
    return spec.family


def _builtin_records() -> list[dict[str, str]]:
    """Rows for the Jarvis-HEP panel, in the same order as the text table."""
    from jarvishep2 import __version__
    from jarvishep2.sampler_catalog import (
        SOURCE_BUILTIN_PLUGIN,
        builtin_plugin_specs,
        builtin_specs,
        plugin_record,
    )

    rows: list[dict[str, str]] = []
    for spec in (*builtin_specs(), *builtin_plugin_specs()):
        if spec.source == SOURCE_BUILTIN_PLUGIN:
            record = plugin_record(spec.name)
            version = record.version if record is not None else ""
        else:
            version = __version__
        rows.append({
            "method": spec.name,
            "version": version,
            "profile": _profile(spec),
            "source": "Jarvis-HEP",
            "status": "built-in",
        })
    return rows


def _internal_info(spec) -> dict[str, Any]:
    """One built-in sampler. Keys follow ``jsampler info --json``, plus family and manual."""
    from jarvishep2 import __version__
    from jarvishep2.sampler_catalog import SOURCE_BUILTIN_PLUGIN, plugin_record

    record = plugin_record(spec.name) if spec.source == SOURCE_BUILTIN_PLUGIN else None
    if record is not None:
        version = record.version
        profile = record.profile
        description = record.description or ""
        statistical = bool(record.statistical)
        options_schema = record.options_schema is not None
        functions = list(record.hooks)
        source = _display_source(record.source_path)
    else:
        version = __version__
        profile = _profile(spec)
        description = ""
        statistical = bool(spec.statistical)
        options_schema = False
        functions = []
        source = "Jarvis-HEP"
    return {
        "command": "sampler info",
        "description": description,
        "family": spec.family,
        "functions": functions,
        "manual": f"Jarvis man sampler.{spec.name}",
        "method": spec.name,
        "ok": True,
        "options_schema": options_schema,
        "profile": profile,
        "source": source,
        "statistical": statistical,
        "status": "built-in",
        "version": version,
    }


def _print_json(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False))


__all__ = ["dispatch_sampler"]
