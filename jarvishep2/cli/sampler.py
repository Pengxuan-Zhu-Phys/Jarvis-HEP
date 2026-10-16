#!/usr/bin/env python3
"""``Jarvis sampler`` forwards to the ``jsampler`` command.

The command surface lives in Jarvis-HEP-Sampler. HEP adds its shipped methods
to read-only list and info results.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path
from typing import Any

from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.run_outcome import EXIT_OK, EXIT_USAGE

_JSON_FLAG = "--json"
SAMPLER_COMMAND_HELP = {
    "list": "List loaded samplers and methods shipped with Jarvis-HEP",
    "info": "Show one loaded sampler or method shipped with Jarvis-HEP",
}


def _call_sdk(sampler_main: Any, args: list[str]) -> int:
    return int(sampler_main(
        args, prog="Jarvis sampler",
        command_help=SAMPLER_COMMAND_HELP,
        command_description=SAMPLER_COMMAND_HELP,
    ))


def dispatch_sampler(argv: list[str] | None = None) -> int:
    """Handle ``Jarvis sampler``, which is ``jsampler`` with a Jarvis program name."""
    args = list(argv or [])
    from jarvishep2.sampler_sdk import SDK_DISTRIBUTION, SDK_SPECIFIER, SamplerSDKError, require_sampler_sdk

    try:
        require_sampler_sdk()
    except SamplerSDKError as exc:
        return _command_error(args, str(exc), type(exc).__name__)
    try:
        from jarvis_hep_sampler.cli import main as sampler_main
    except ImportError as exc:
        return _command_error(
            args,
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is required for Jarvis sampler. "
            f"Install it with pip install -U '{SDK_DISTRIBUTION}{SDK_SPECIFIER}'. "
            f"Detail: {exc}",
            type(exc).__name__,
        )
    from jarvis_hep_sampler.errors import PluginError

    try:
        return _dispatch_sdk(args, sampler_main)
    except (PluginError, OSError) as exc:
        return _command_error(args, str(exc), type(exc).__name__)
    except SystemExit as exc:
        code = int(exc.code or 0)
        if code and _JSON_FLAG in args:
            return _command_error(args, "Invalid sampler arguments; use Jarvis sampler <command> --help.", "usage_error")
        return code


def _dispatch_sdk(args: list[str], sampler_main: Any) -> int:
    command, as_json = _command(args)
    from jarvishep2.sampler_catalog import SOURCE_USER, registry_report

    if command == ["list"]:
        report = registry_report()
        if as_json:
            return _print_list_json(sampler_main, args, report)
        return _print_list_text(sampler_main, args, report)
    if len(command) == 2 and command[0] == "info":
        report = registry_report()
        if report.error is not None:
            return _call_sdk(sampler_main, args)
        exact_user = next(
            (item for item in report.entries
             if item.source == SOURCE_USER and item.method == command[1]),
            None,
        )
        spec = None if exact_user is not None else _internal_spec(command[1])
        if spec is not None:
            if as_json:
                _print_json(_internal_info(spec))
            else:
                _print_internal_info(spec)
            return EXIT_OK
        if exact_user is not None and as_json:
            return _print_info_json(sampler_main, args, exact_user)
        if exact_user is not None:
            return _print_info_text(sampler_main, args, exact_user)
        code = _call_sdk(sampler_main, args)
        return code
    return _call_sdk(sampler_main, args)


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


def _print_internal_list() -> None:
    table = Table(box=box.SIMPLE_HEAVY, show_header=True, header_style="bold dim")
    for name in ("Method", "Version", "Profile", "Source", "Status"):
        table.add_column(name.upper(), no_wrap=True)
    for row in _builtin_records():
        table.add_row(
            Text(row["method"], style="bold #c8c8ff"),
            row["version"],
            row["profile"],
            row["source"],
            Text(row["status"], style="bold cyan"),
        )
    _console().print(
        Panel(
            table,
            title=Text("Jarvis-HEP", style="bold"),
            box=box.ROUNDED,
            border_style="dim",
            padding=(0, 1),
        )
    )


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


def _print_internal_info(spec) -> None:
    from jarvishep2.sampler_catalog import SOURCE_BUILTIN_PLUGIN

    payload = _internal_info(spec)
    if spec.source == SOURCE_BUILTIN_PLUGIN:
        rows: list[tuple[str, Any]] = [
            ("version", payload["version"]),
            ("profile", payload["profile"]),
            ("description", payload["description"] or "—"),
            ("statistical", "yes" if payload["statistical"] else "no"),
            ("options schema", "yes" if payload["options_schema"] else "no"),
            ("functions", ", ".join(payload["functions"])),
            ("source", payload["source"]),
            ("status", payload["status"]),
        ]
    else:
        rows = [
            ("version", payload["version"]),
            ("profile", payload["profile"]),
            ("family", payload["family"]),
            ("statistical", "yes" if payload["statistical"] else "no"),
            ("source", payload["source"]),
            ("status", payload["status"]),
        ]
    _panel(f"Sampler · {spec.name}", rows, note=f"Manual: {payload['manual']}")


def _sdk_json(sampler_main: Any, args: list[str]) -> tuple[int, dict[str, Any] | None, str]:
    buffer = io.StringIO()
    sdk_args = args if _JSON_FLAG in args else [*args, _JSON_FLAG]
    with contextlib.redirect_stdout(buffer):
        code = _call_sdk(sampler_main, sdk_args)
    raw = buffer.getvalue()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        payload = None
    return code, payload if isinstance(payload, dict) else None, raw


def _user_statuses(report: Any) -> dict[str, Any]:
    from jarvishep2.sampler_catalog import SOURCE_USER

    return {entry.method: entry for entry in report.entries if entry.source == SOURCE_USER}


def _print_list_text(sampler_main: Any, args: list[str], report: Any) -> int:
    code, payload, _ = _sdk_json(sampler_main, args)
    if code != EXIT_OK or payload is None or payload.get("ok") is not True:
        return _call_sdk(sampler_main, args)
    by_name = _user_statuses(report)
    table = Table(box=box.SIMPLE_HEAVY, show_header=True, header_style="bold dim")
    for name in ("Method", "Version", "Profile", "Source", "Status"):
        table.add_column(name.upper(), overflow="fold")
    for item in payload.get("samplers", []):
        entry = by_name.get(item["method"])
        status = entry.status if entry is not None and entry.status != "ok" else item["status"]
        table.add_row(
            Text(item["method"], style="bold #c8c8ff"),
            item["version"], item["profile"],
            Text(item["source"], style="dim"),
            Text(status, style="bold red" if status != "ok" else "bold green"),
        )
    note = Text(f"Registry: {payload.get('registry', '')}", style="dim")
    body: Any = Group(table, Text(""), note)
    _console().print(Panel(body, title=Text("Loaded samplers", style="bold"),
                           box=box.ROUNDED, border_style="dim", padding=(0, 1)))
    for entry in by_name.values():
        if entry.repair_command:
            _console().print(f"{entry.method}: {entry.repair_command}")
    _print_internal_list()
    return code


def _print_info_text(sampler_main: Any, args: list[str], entry: Any) -> int:
    code, payload, _ = _sdk_json(sampler_main, args)
    if code != EXIT_OK or payload is None or payload.get("ok") is not True:
        return _call_sdk(sampler_main, args)
    status = entry.status if entry.status != "ok" else payload.get("status", "ok")
    rows: list[tuple[str, Any]] = [
        ("version", payload.get("version", "")),
        ("profile", payload.get("profile", "")),
        ("description", payload.get("description") or "—"),
        ("statistical", "yes" if payload.get("statistical") else "no"),
        ("options schema", "yes" if payload.get("options_schema") else "no"),
        ("functions", ", ".join(payload.get("functions", []))),
        ("source", payload.get("source", "")),
        ("snapshot", payload.get("snapshot", "")),
        ("sha256", payload.get("sha256", "")),
        ("loaded at", payload.get("loaded_at", "")),
        ("status", status),
    ]
    if entry.problems:
        rows.append(("problem", "; ".join(entry.problems)))
    if entry.repair_command:
        rows.append(("repair", entry.repair_command))
    _panel(f"Sampler · {entry.method}", rows)
    return code


def _print_info_json(sampler_main: Any, args: list[str], entry: Any) -> int:
    code, payload, raw = _sdk_json(sampler_main, args)
    if payload is None:
        sys.stdout.write(raw)
        return code
    if payload.get("ok") is True:
        payload["status"] = entry.status if entry.status != "ok" else payload.get("status", "ok")
        if entry.problems:
            payload["problems"] = list(entry.problems)
        if entry.repair_command:
            payload["repair_command"] = entry.repair_command
        _print_json(payload)
        return code
    sys.stdout.write(raw)
    return code


def _print_list_json(sampler_main: Any, args: list[str], report: Any) -> int:
    """Keep jsampler's registry object and add the shipped methods beside it."""
    code, payload, raw = _sdk_json(sampler_main, args)
    if payload is None:
        sys.stdout.write(raw)
        return code
    if payload.get("ok") is True:
        by_name = _user_statuses(report)
        for item in payload.get("samplers", []):
            entry = by_name.get(item.get("method"))
            if entry is None:
                continue
            if entry.status != "ok":
                item["status"] = entry.status
            if entry.problems:
                item["problems"] = list(entry.problems)
            if entry.repair_command:
                item["repair_command"] = entry.repair_command
        payload["jarvis_hep"] = _builtin_records()
        _print_json(payload)
        return code
    sys.stdout.write(raw)
    return code


def _print_json(payload: dict[str, Any]) -> None:
    print(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False))


def _console() -> Console:
    return Console()


def _panel(title: str, rows: list[tuple[str, Any]], *, note: str = "") -> None:
    grid = Table.grid(padding=(0, 1))
    grid.add_column(style="bold dim", no_wrap=True)
    grid.add_column(overflow="fold")
    for key, value in rows:
        grid.add_row(key.upper(), Text(str(value)))
    body: Any = grid
    if note:
        body = Group(grid, Text(""), Text(note, style="dim"))
    _console().print(
        Panel(
            body,
            title=Text(title, style="bold"),
            box=box.ROUNDED,
            border_style="cyan",
            padding=(0, 1),
        )
    )


__all__ = ["dispatch_sampler"]
