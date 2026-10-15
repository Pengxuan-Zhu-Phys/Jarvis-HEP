#!/usr/bin/env python3
"""``Jarvis sampler`` forwards to the ``jsampler`` command.

The command surface lives in Jarvis-HEP-Sampler. This module is the same
entrypoint pattern as ``Jarvis plot`` and ``jplot``, except ``list`` and
``info`` also show the samplers that ship inside Jarvis-HEP.
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

_HELP_FLAGS = frozenset({"-h", "--help"})
_JSON_FLAG = "--json"


def dispatch_sampler(argv: list[str] | None = None) -> int:
    """Handle ``Jarvis sampler``, which is ``jsampler`` with a Jarvis program name."""
    args = list(argv or [])
    if not args or args[0] in _HELP_FLAGS or any(arg in _HELP_FLAGS for arg in args):
        from jarvishep2.cli.help import _render_help_for_args, build_sampler_help_parser

        return _render_help_for_args(build_sampler_help_parser(), args or ["-h"])
    from jarvishep2.sampler_sdk import SDK_DISTRIBUTION, SDK_SPECIFIER, SamplerSDKError, require_sampler_sdk

    try:
        require_sampler_sdk()
    except SamplerSDKError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    try:
        from jarvis_hep_sampler.cli import main as sampler_main
    except ImportError as exc:
        print(
            f"{SDK_DISTRIBUTION}{SDK_SPECIFIER} is required for `Jarvis sampler`. "
            f"Install it with `pip install -U \"{SDK_DISTRIBUTION}{SDK_SPECIFIER}\"`. "
            f"Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    command, as_json = _command(args)
    if command == ["list"]:
        if as_json:
            return _print_list_json(sampler_main, args)
        code = int(sampler_main(args, prog="Jarvis sampler"))
        _print_internal_list()
        return code
    if len(command) == 2 and command[0] == "info":
        spec = _internal_spec(command[1])
        if spec is not None:
            if as_json:
                _print_json(_internal_info(spec))
            else:
                _print_internal_info(spec)
            return EXIT_OK
    return int(sampler_main(args, prog="Jarvis sampler"))


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
        "command": "info",
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


def _print_list_json(sampler_main: Any, args: list[str]) -> int:
    """Keep jsampler's registry object and add the shipped methods beside it."""
    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        code = int(sampler_main(args, prog="Jarvis sampler"))
    text = buffer.getvalue()
    try:
        payload = json.loads(text)
    except json.JSONDecodeError:
        sys.stdout.write(text)
        return code
    if isinstance(payload, dict) and payload.get("ok") is True:
        payload["jarvis_hep"] = _builtin_records()
        _print_json(payload)
        return code
    sys.stdout.write(text)
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
