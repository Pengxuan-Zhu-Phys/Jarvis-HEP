#!/usr/bin/env python3
"""``Jarvis sampler`` forwards to the ``jsampler`` command.

The command surface lives in Jarvis-HEP-Sampler. This module is the same
entrypoint pattern as ``Jarvis plot`` and ``jplot``, except ``list`` and
``info`` also show the samplers that ship inside Jarvis-HEP.
"""

from __future__ import annotations

import sys
from typing import Any

from jarvishep2.run_outcome import EXIT_OK, EXIT_USAGE


def dispatch_sampler(argv: list[str] | None = None) -> int:
    """Handle ``Jarvis sampler``, which is ``jsampler`` with a Jarvis program name."""
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
    args = list(argv or [])
    if args == ["list"]:
        code = int(sampler_main(args, prog="Jarvis sampler"))
        _print_internal_list()
        return code
    if len(args) == 2 and args[0] == "info":
        spec = _internal_spec(args[1])
        if spec is not None:
            _print_internal_info(spec)
            return EXIT_OK
    return int(sampler_main(args, prog="Jarvis sampler"))


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


def _profile(spec) -> str:
    from jarvishep2.sampler_catalog import SOURCE_BUILTIN_PLUGIN, plugin_record

    if spec.source == SOURCE_BUILTIN_PLUGIN:
        record = plugin_record(spec.name)
        return record.profile if record is not None else spec.family
    if spec.stateless and not spec.statistical:
        return "independent"
    return spec.family


def _print_internal_list() -> None:
    from jarvishep2 import __version__
    from jarvishep2.sampler_catalog import (
        SOURCE_BUILTIN_PLUGIN,
        builtin_plugin_specs,
        builtin_specs,
        plugin_record,
    )

    rows: list[tuple[str, str, str, str, str]] = []
    for spec in (*builtin_specs(), *builtin_plugin_specs()):
        if spec.source == SOURCE_BUILTIN_PLUGIN:
            record = plugin_record(spec.name)
            version = record.version if record is not None else ""
            source = record.source_path if record is not None else "Jarvis-HEP"
        else:
            version = __version__
            source = "Jarvis-HEP"
        rows.append((spec.name, version, _profile(spec), source, "built-in"))
    print()
    print("Jarvis-HEP:")
    _print_rows(rows)


def _print_internal_info(spec) -> None:
    from jarvishep2 import __version__
    from jarvishep2.sampler_catalog import SOURCE_BUILTIN_PLUGIN, plugin_record

    if spec.source == SOURCE_BUILTIN_PLUGIN:
        record = plugin_record(spec.name)
        rows: list[tuple[str, Any]] = [
            ("version", record.version if record is not None else ""),
            ("profile", _profile(spec)),
            ("description", (record.description if record is not None else "") or "—"),
            ("statistical", "yes" if spec.statistical else "no"),
            ("options schema", "yes" if record is not None and record.options_schema is not None else "no"),
            ("functions", ", ".join(record.hooks) if record is not None else ""),
            ("source", record.source_path if record is not None else "Jarvis-HEP"),
            ("status", "built-in"),
        ]
    else:
        rows = [
            ("version", __version__),
            ("profile", _profile(spec)),
            ("family", spec.family),
            ("statistical", "yes" if spec.statistical else "no"),
            ("source", "Jarvis-HEP"),
            ("status", "built-in"),
            ("manual", f"Jarvis man sampler.{spec.name}"),
        ]
    _table(f"Sampler {spec.name} ->", rows)


def _print_rows(rows: list[tuple[str, str, str, str, str]]) -> None:
    header = ("Method", "Version", "Profile", "Source", "Status")
    widths = [max(len(str(row[i])) for row in [header, *rows]) for i in range(len(header))]
    for row in [header, tuple("-" * width for width in widths), *rows]:
        line = "  ".join(str(cell).ljust(width) for cell, width in zip(row, widths, strict=True))
        print(line.rstrip())


def _table(title: str, rows: list[tuple[str, Any]]) -> None:
    print(title)
    width = max((len(str(key)) for key, _ in rows), default=0)
    for key, value in rows:
        print(f"  {str(key).ljust(width)}  {value}")


__all__ = ["dispatch_sampler"]
