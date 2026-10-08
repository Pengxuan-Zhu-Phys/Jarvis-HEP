#!/usr/bin/env python3
"""``Jarvis sampler``: create, check, load and manage your own samplers.

The plug-in format, the template, the offline check and the registry
(``~/.jarvis/samplers``) belong to Jarvis-HEP-Sampler; this module is the CLI
over them. After ``Jarvis sampler load``, ``Sampling.Method: <Name>`` works
like a built-in method in every new terminal. Design:
``docs/jarvis-hep-sampler/technical-design.md`` §5.1.
"""

from __future__ import annotations

import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from jarvishep2.cli.help import _render_help_for_args, build_sampler_help_parser
from jarvishep2.log_kv import format_two_column_log
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE

_HELP_FLAGS = frozenset({"-h", "--help"})


def _error(message: str) -> int:
    print(f"[Jarvis] {message}", file=sys.stderr)
    return EXIT_USAGE


def _reserved() -> frozenset[str]:
    from jarvishep2.sampler_catalog import reserved_method_names

    return reserved_method_names()


def _table(title: str, rows: list[tuple[str, Any]]) -> None:
    print(format_two_column_log(title, rows))


# ---------------------------------------------------------------- init
def _run_init(method: str, output: str | None) -> int:
    from jarvis_hep_sampler import write_template
    from jarvis_hep_sampler.errors import PluginError

    try:
        path = write_template(method, output, reserved_names=_reserved())
    except (PluginError, FileExistsError) as exc:
        return _error(str(exc))
    print(f"Created {path} from the independent sampler template.")
    print(f"Edit draw(), then run: Jarvis sampler load {path}")
    return EXIT_OK


# --------------------------------------------------------------- check
def _load_yaml_mapping(path: str, what: str) -> Mapping[str, Any]:
    import yaml

    with open(Path(path).expanduser(), encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, Mapping):
        raise ValueError(f"{what} {path} must contain a YAML mapping")
    return data


def _check_target(target: str) -> str | None:
    """A file path, or the file a loaded / shipped sampler runs from."""
    if Path(target).expanduser().is_file():
        return target
    from jarvishep2.sampler_catalog import plugin_record

    record = plugin_record(target)
    return record.snapshot_path if record is not None else None


def _run_check(args: Any) -> int:
    from jarvis_hep_sampler.options import parse_reserved, split_bounds
    from jarvis_hep_sampler.testing import check_plugin

    target = _check_target(str(args.target))
    if target is None:
        return _error(f"{args.target} is neither a file nor a loaded sampler")

    variable_names: list[str] | None = None
    options: dict[str, Any] = {}
    seed = 0
    point_number = 100
    overrides: list[str] = []
    try:
        if args.card:
            from jarvishep2.sampling.variables import load_variables
            from jarvishep2.task_config import load_task_yaml

            config = load_task_yaml(args.card)
            variable_names = [variable.name for variable in load_variables(config)]
            bounds = dict((config.get("Sampling") or {}).get("Bounds") or {})
            reserved, options = split_bounds(bounds)
            settings = parse_reserved(reserved)
            seed = settings.seed
            point_number = int(settings.point_number)
        if args.options:
            options = dict(_load_yaml_mapping(args.options, "options file"))
            if args.card:
                overrides.append("options (replaced by --options)")
    except Exception as exc:
        return _error(str(exc))
    if args.seed is not None:
        seed = int(args.seed)
        if args.card:
            overrides.append("seed")
    if args.points is not None:
        point_number = int(args.points)
        if args.card:
            overrides.append("point_number (changes the total candidate budget)")
    dimension = None if variable_names else int(args.dimension or 2)

    _table("Offline check settings ->", [
        ("variables", ", ".join(variable_names) if variable_names else dimension),
        ("seed", seed), ("point number", point_number),
        ("card overrides", ", ".join(overrides) if overrides else "none"),
    ])

    report = check_plugin(
        target,
        point_number=point_number,
        seed=seed,
        variable_names=variable_names,
        dimension=dimension,
        options=options,
        reserved_names=() if target != args.target else _reserved(),
    )
    _table(f"{report.method or args.target} check ->", report.rows())
    return EXIT_OK if report.ok else EXIT_RUN_FAILED


# ---------------------------------------------------------------- load
def _run_load(file: str) -> int:
    from jarvis_hep_sampler import SamplerRegistry
    from jarvis_hep_sampler.errors import PluginError

    try:
        result = SamplerRegistry().load(file, reserved_names=_reserved())
    except PluginError as exc:
        print(f"[Jarvis] Cannot load {file}: {exc}", file=sys.stderr)
        cause = exc.__cause__
        if cause is not None:
            print(f"[Jarvis] Caused by {type(cause).__name__}: {cause}", file=sys.stderr)
        return EXIT_USAGE
    entry = result.entry
    rows: list[tuple[str, Any]] = [
        ("version", entry.version),
        ("profile", entry.profile),
        ("options schema", "yes" if entry.options_schema is not None else "no"),
        ("source", entry.source_path),
        ("snapshot", SamplerRegistry().snapshot_path(entry)),
    ]
    if result.replaced is not None:
        old = result.replaced
        rows.append(
            (
                "replaced",
                f"{old.version} ({old.short_hash}) -> {entry.version} ({entry.short_hash})"
                if result.changed
                else "unchanged (same file)",
            )
        )
    _table(f"Loaded sampler {entry.method} ->", rows)
    print(f"Use it in a task card as  Sampling.Method: {entry.method}")
    return EXIT_OK


# ---------------------------------------------------------------- list
def _status(registry: Any, entry: Any, reserved: frozenset[str]) -> str:
    if entry.method.casefold() in {name.casefold() for name in reserved}:
        return "shadowed by a built-in method"
    problems = registry.problems(entry)
    if problems:
        return f"broken: {problems[0]}"
    return registry.source_state(entry)


def _run_list() -> int:
    from jarvis_hep_sampler import SamplerRegistry
    from jarvis_hep_sampler.errors import RegistryError

    from jarvishep2.sampler_catalog import builtin_plugin_specs, plugin_record

    registry = SamplerRegistry()
    rows: list[tuple[str, str, str, str, str]] = []
    for spec in builtin_plugin_specs():
        record = plugin_record(spec.name)
        rows.append((spec.name, record.version if record else "", "independent", "shipped", "ok"))
    try:
        entries = registry.entries()
    except RegistryError as exc:
        _print_rows(rows)
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    reserved = _reserved()
    for method, entry in entries.items():
        rows.append(
            (method, entry.version, entry.profile, entry.source_path, _status(registry, entry, reserved))
        )
    _print_rows(rows)
    print(f"Registry: {registry.home}")
    return EXIT_OK


def _print_rows(rows: list[tuple[str, str, str, str, str]]) -> None:
    header = ("Method", "Version", "Profile", "Source", "Status")
    widths = [max(len(str(row[i])) for row in [header, *rows]) for i in range(len(header))]
    for row in [header, tuple("-" * width for width in widths), *rows]:
        print("  ".join(str(cell).ljust(width) for cell, width in zip(row, widths)).rstrip())


# ---------------------------------------------------------------- info
def _run_info(method: str) -> int:
    from jarvis_hep_sampler import SamplerRegistry

    from jarvishep2.sampler_catalog import SOURCE_USER, plugin_record, spec_for

    spec = spec_for(method)
    record = plugin_record(method)
    if spec is None or record is None:
        return _error(f"no plug-in sampler named {method!r}; see: Jarvis sampler list")
    rows: list[tuple[str, Any]] = [
        ("version", record.version),
        ("profile", record.profile),
        ("description", record.description or "—"),
        ("statistical", "yes" if record.statistical else "no"),
        ("options schema", "yes" if record.options_schema is not None else "no"),
        ("functions", ", ".join(record.hooks)),
        ("source", record.source_path),
    ]
    if record.source == SOURCE_USER:
        registry = SamplerRegistry()
        entry = registry.get(method)
        rows += [
            ("snapshot", record.snapshot_path),
            ("sha256", record.source_sha256),
            ("loaded at", record.loaded_at),
            ("status", _status(registry, entry, _reserved()) if entry else "missing"),
        ]
    else:
        rows.append(("status", "shipped with Jarvis-HEP"))
    _table(f"Sampler {spec.name} ->", rows)
    return EXIT_OK


# -------------------------------------------------------------- remove
def _run_remove(method: str, *, yes: bool) -> int:
    from jarvis_hep_sampler import SamplerRegistry
    from jarvis_hep_sampler.errors import RegistryError

    if method.casefold() in {name.casefold() for name in _reserved()}:
        return _error(f"{method} is a built-in method and cannot be removed")
    registry = SamplerRegistry()
    if not yes:
        if not sys.stdin.isatty():
            return _error("pass --yes to remove a sampler non-interactively")
        answer = input(f"Remove sampler {method} from {registry.home}? [y/N] ").strip().lower()
        if answer not in {"y", "yes"}:
            print("Nothing removed.")
            return EXIT_OK
    try:
        removed = registry.remove(method)
    except RegistryError as exc:
        return _error(str(exc))
    print(f"Removed sampler {removed.method}; your file {removed.source_path} is unchanged.")
    return EXIT_OK


# ------------------------------------------------------------ dispatch
def dispatch_sampler(argv: list[str] | None = None) -> int:
    """Handle ``Jarvis sampler init|check|load|list|info|remove``."""
    args = list(argv or [])
    parser = build_sampler_help_parser()
    if not args or args[0] in _HELP_FLAGS:
        parser.print_help()
        return EXIT_OK
    if any(arg in _HELP_FLAGS for arg in args):
        return _render_help_for_args(parser, args)
    try:
        parsed = parser.parse_args(args)
    except SystemExit as exc:
        return int(exc.code or EXIT_USAGE)

    command = parsed.sampler_command
    if command == "init":
        return _run_init(parsed.method, parsed.output)
    if command == "check":
        return _run_check(parsed)
    if command == "load":
        return _run_load(parsed.file)
    if command == "list":
        return _run_list()
    if command == "info":
        return _run_info(parsed.method)
    if command == "remove":
        return _run_remove(parsed.method, yes=bool(parsed.yes))
    parser.print_help()
    return EXIT_USAGE


__all__ = ["dispatch_sampler"]
