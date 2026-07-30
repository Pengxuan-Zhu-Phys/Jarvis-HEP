#!/usr/bin/env python3
"""Jarvis2 CLI surface: one-intent subcommands + legacy flag aliases (D11.1 / D11.2)."""

from __future__ import annotations

import argparse
import contextlib
import io
import os
import shutil
import sys
import warnings
from typing import Any

import click
import typer.rich_utils
from rich import box
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.core import Jarvis2Core
from jarvishep2.dashboard import attach_reader, format_monitor_view
from jarvishep2.factory import TaskFactory
from jarvishep2.plot_bridge import PlotBridgeError, run_plot
from jarvishep2.redis_queue import INTERNAL_REDIS_CONFIG, RedisQueue
from jarvishep2.run_outcome import (
    EXIT_INTERRUPTED,
    EXIT_OK,
    EXIT_RUN_FAILED,
    EXIT_USAGE,
    RunOutcome,
)


_SUBCOMMANDS = frozenset(
    {
        "run",
        "check",
        "validate",
        "convert",
        "monitor",
        "plot",
        "gen-plot-yaml",
        "portal",
        "operas",
        "project",
        "ps",
        "kill",
    }
)
_PROJECT_COMMANDS = frozenset(
    {"create", "pack", "browse", "list", "fetch", "info", "encrypt"}
)
_PACK_MODE_FLAGS = {
    "--share": "share",
    "--repro": "repro",
    "--full": "full",
}
_PACK_MANIFEST_FLAG = "--man"
_PACK_ENCRYPT_FLAGS = frozenset({"--encrypt", "--encrypted"})
_HELP_FLAGS = frozenset({"-h", "--help"})
_LEGACY_OPTION_DESTS = frozenset(
    {
        "monitor",
        "plot",
        "validate",
        "strict",
        "check_modules",
        "convert",
        "resume",
        "skip_draw_flowchart",
    }
)
_COMMAND_HELP_PANELS = {
    "run": "Scan workflow",
    "check": "Scan workflow",
    "validate": "Scan workflow",
    "convert": "Data export",
    "gen-plot-yaml": "Plots",
    "plot": "Plots",
    "monitor": "Runtime control",
    "ps": "Runtime control",
    "kill": "Runtime control",
    "portal": "Extensions & projects",
    "operas": "Extensions & projects",
    "project": "Extensions & projects",
}
_COMMAND_HELP_ORDER = {
    # Scan lifecycle
    "validate": 10,
    "check": 20,
    "run": 30,
    "convert": 40,
    # Data export and plot lifecycle
    "gen-plot-yaml": 50,
    "plot": 60,
    # Local runtime lifecycle
    "monitor": 70,
    "ps": 80,
    "kill": 90,
    # Project / extension entry points
    "project": 100,
    "portal": 110,
    "operas": 120,
    # Nested helper workflows
    "create": 130,
    "pack": 140,
    "encrypt": 150,
    "list": 160,
    "browse": 170,
    "fetch": 180,
    "info": 190,
}


_HELP_PRIMARY_COLUMN_WIDTH = 24
_HELP_ALIAS_COLUMN_WIDTH = 6


class JarvisHelpGroup(click.Group):
    """Click command group with an intentional workflow-oriented help order."""

    def list_commands(self, ctx: click.Context) -> list[str]:
        names = super().list_commands(ctx)
        return sorted(
            names,
            key=lambda name: (_COMMAND_HELP_ORDER.get(name, 1_000), name),
        )


def _fixed_option_label(
    param: click.Option | click.Argument, ctx: click.Context
) -> tuple[Text, Text]:
    """Return the common primary/alias columns used by every help panel."""
    if isinstance(param, click.Argument):
        return (
            typer.rich_utils.highlighter(param.name or ""),
            typer.rich_utils.metavar_highlighter(param.make_metavar(ctx=ctx)),
        )

    long_options = [option for option in param.opts if option.startswith("--")]
    short_options = [option for option in param.opts if not option.startswith("--")]
    label = ",".join(long_options or short_options)
    aliases = ",".join(short_options if long_options else [])
    metavar = param.make_metavar(ctx=ctx)
    if metavar != "BOOLEAN":
        label = f"{label} {metavar}".strip()
    return (
        typer.rich_utils.highlighter(label),
        typer.rich_utils.highlighter(aliases),
    )


def _print_fixed_options_panel(
    *,
    name: str,
    params: list[click.Option] | list[click.Argument],
    ctx: click.Context,
    markup_mode: typer.rich_utils.MarkupMode,
    console: Any,
) -> None:
    """Render options in the same fixed columns as command panels."""
    if not params:
        return
    table = Table(
        highlight=True,
        show_header=False,
        expand=True,
        box=getattr(box, typer.rich_utils.STYLE_OPTIONS_TABLE_BOX, None),
        show_lines=typer.rich_utils.STYLE_OPTIONS_TABLE_SHOW_LINES,
        leading=typer.rich_utils.STYLE_OPTIONS_TABLE_LEADING,
        border_style=typer.rich_utils.STYLE_OPTIONS_TABLE_BORDER_STYLE,
        row_styles=typer.rich_utils.STYLE_OPTIONS_TABLE_ROW_STYLES,
        pad_edge=typer.rich_utils.STYLE_OPTIONS_TABLE_PAD_EDGE,
        padding=typer.rich_utils.STYLE_OPTIONS_TABLE_PADDING,
    )
    table.add_column(width=_HELP_PRIMARY_COLUMN_WIDTH, no_wrap=True)
    table.add_column(width=_HELP_ALIAS_COLUMN_WIDTH, no_wrap=True)
    table.add_column(justify="left", ratio=1)
    for param in params:
        primary, aliases = _fixed_option_label(param, ctx)
        table.add_row(
            primary,
            aliases,
            typer.rich_utils._get_parameter_help(
                param=param, ctx=ctx, markup_mode=markup_mode
            ),
        )
    console.print(
        Panel(
            table,
            border_style=typer.rich_utils.STYLE_OPTIONS_PANEL_BORDER,
            title=name,
            title_align=typer.rich_utils.ALIGN_OPTIONS_PANEL,
        )
    )


def _print_fixed_commands_panel(
    *,
    name: str,
    commands: list[click.Command],
    markup_mode: typer.rich_utils.MarkupMode,
    console: Any,
    cmd_len: int,
) -> None:
    """Render command groups with the option panel's fixed column grid."""
    if not commands:
        return
    table = Table(
        highlight=False,
        show_header=False,
        expand=True,
        box=getattr(box, typer.rich_utils.STYLE_COMMANDS_TABLE_BOX, None),
        show_lines=typer.rich_utils.STYLE_COMMANDS_TABLE_SHOW_LINES,
        leading=typer.rich_utils.STYLE_COMMANDS_TABLE_LEADING,
        border_style=typer.rich_utils.STYLE_COMMANDS_TABLE_BORDER_STYLE,
        row_styles=typer.rich_utils.STYLE_COMMANDS_TABLE_ROW_STYLES,
        pad_edge=typer.rich_utils.STYLE_COMMANDS_TABLE_PAD_EDGE,
        padding=typer.rich_utils.STYLE_COMMANDS_TABLE_PADDING,
    )
    table.add_column(
        style=typer.rich_utils.STYLE_COMMANDS_TABLE_FIRST_COLUMN,
        width=_HELP_PRIMARY_COLUMN_WIDTH,
        no_wrap=True,
    )
    table.add_column(width=_HELP_ALIAS_COLUMN_WIDTH, no_wrap=True)
    table.add_column(justify="left", ratio=1)
    for command in commands:
        helptext = command.short_help or command.help or ""
        table.add_row(
            Text(command.name or ""),
            Text(),
            typer.rich_utils._make_command_help(
                help_text=helptext, markup_mode=markup_mode
            ),
        )
    console.print(
        Panel(
            table,
            border_style=typer.rich_utils.STYLE_COMMANDS_PANEL_BORDER,
            title=name,
            title_align=typer.rich_utils.ALIGN_COMMANDS_PANEL,
        )
    )


class JarvisArgumentParser(argparse.ArgumentParser):
    """Argparse parser rendered by Jarvis-Lit's Typer/Rich help formatter.

    ``argparse`` remains the source of truth for parsing and legacy aliases.
    Help is converted to a lightweight Click command tree and sent through the
    same Typer formatter used by Jarvis-Lit, so its rounded panels and layout
    are shared rather than reimplemented here.
    """

    _RICH_CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}

    @staticmethod
    def _metavar(action: argparse.Action) -> str:
        if action.type is int:
            return "<int>"
        if action.type is float:
            return "<float>"
        return "<str>"

    def _click_params(self, *, legacy: bool | None = None) -> list[click.Parameter]:
        params: list[click.Parameter] = []
        for action in self._actions:
            if action.help == argparse.SUPPRESS or action.dest == "help":
                continue
            if isinstance(action, argparse._SubParsersAction):
                continue
            is_legacy = action.dest in _LEGACY_OPTION_DESTS
            if legacy is not None and is_legacy is not legacy:
                continue
            if action.option_strings:
                params.append(
                    click.Option(
                        action.option_strings,
                        help=action.help,
                        is_flag=action.nargs == 0,
                        required=bool(action.required),
                        metavar=None if action.nargs == 0 else self._metavar(action),
                    )
                )
            else:
                argument = click.Argument(
                    [action.dest],
                    required=bool(action.required),
                    metavar=self._metavar(action),
                )
                # Click arguments have no help constructor parameter, whereas
                # argparse stores it on the action. Preserve it for Rich.
                argument.help = action.help
                params.append(argument)
        return params

    def _click_command(
        self, name: str | None = None, *, include_legacy: bool = True
    ) -> click.Command:
        commands: dict[str, click.Command] = {}
        for action in self._actions:
            if not isinstance(action, argparse._SubParsersAction):
                continue
            short_help = {
                choice.dest: choice.help or ""
                for choice in getattr(action, "_choices_actions", [])
            }
            for command_name, command_parser in action.choices.items():
                command = command_parser._click_command(command_name)
                command.short_help = short_help.get(command_name)
                panel = _COMMAND_HELP_PANELS.get(command_name)
                if panel:
                    setattr(command, typer.rich_utils._RICH_HELP_PANEL_NAME, panel)
                commands[command_name] = command
        command_name = name or self.prog.rsplit(" ", maxsplit=1)[-1]
        if commands:
            return JarvisHelpGroup(
                name=command_name,
                help=self.description,
                params=self._click_params(legacy=None if include_legacy else False),
                commands=commands,
                context_settings=self._RICH_CONTEXT_SETTINGS,
            )
        return click.Command(
            name=command_name,
            help=self.description,
            params=self._click_params(legacy=None if include_legacy else False),
            context_settings=self._RICH_CONTEXT_SETTINGS,
        )

    def format_help(self) -> str:
        root_help = self.prog == "Jarvis2"
        command = self._click_command(include_legacy=not root_help)
        terminal_stdout = sys.stdout
        output = io.StringIO()
        old_force_terminal = typer.rich_utils.FORCE_TERMINAL
        old_max_width = typer.rich_utils.MAX_WIDTH
        old_options_panel = typer.rich_utils._print_options_panel
        old_commands_panel = typer.rich_utils._print_commands_panel
        try:
            typer.rich_utils.FORCE_TERMINAL = terminal_stdout.isatty()
            typer.rich_utils.MAX_WIDTH = shutil.get_terminal_size().columns
            typer.rich_utils._print_options_panel = _print_fixed_options_panel
            typer.rich_utils._print_commands_panel = _print_fixed_commands_panel
            with contextlib.redirect_stdout(output):
                typer.rich_utils.rich_format_help(
                    obj=command,
                    ctx=click.Context(
                        command,
                        info_name=self.prog,
                        **self._RICH_CONTEXT_SETTINGS,
                    ),
                    markup_mode=typer.rich_utils.MARKUP_MODE_RICH,
                )
                if root_help:
                    legacy_params = self._click_params(legacy=True)
                    typer.rich_utils._print_options_panel(
                        name="Legacy options",
                        params=legacy_params,
                        ctx=click.Context(command, info_name=self.prog),
                        markup_mode=typer.rich_utils.MARKUP_MODE_RICH,
                        console=typer.rich_utils._get_rich_console(),
                    )
        finally:
            typer.rich_utils.FORCE_TERMINAL = old_force_terminal
            typer.rich_utils.MAX_WIDTH = old_max_width
            typer.rich_utils._print_options_panel = old_options_panel
            typer.rich_utils._print_commands_panel = old_commands_panel
        return output.getvalue()


def normalize_argv(argv: list[str] | None) -> list[str]:
    """Rewrite legacy bare-YAML invocations into explicit subcommands.

    Argparse subparsers treat the first positional as the command name, so
    ``Jarvis2 task.yaml`` would otherwise fail as an invalid choice.
    """
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        return args
    head = args[0]
    if head in _SUBCOMMANDS or head.startswith("-"):
        return args

    # First token is a path (legacy).
    if "--plot" in args:
        rest = [a for a in args[1:] if a != "--plot"]
        return ["plot", head, *rest]
    if "--check-modules" in args:
        rest = [a for a in args[1:] if a != "--check-modules"]
        return ["check", head, *rest]
    if "--validate" in args:
        rest = [a for a in args[1:] if a != "--validate"]
        return ["validate", head, *rest]
    if "--convert" in args:
        rest = [a for a in args[1:] if a != "--convert"]
        return ["convert", head, *rest]
    # Bare path (+ optional --resume) → run
    return ["run", *args]


def build_parser() -> argparse.ArgumentParser:
    parser = JarvisArgumentParser(
        prog="Jarvis2",
        description=(
            "Run and validate distributed HEP scans, prepare plot YAML, "
            "and manage their local runtime."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Intents (preferred):\n"
            "  Jarvis2 run TASK.yaml [--resume]\n"
            "  Jarvis2 check TASK.yaml\n"
            "  Jarvis2 validate TASK.yaml [--strict] [--json]\n"
            "  Jarvis2 convert TASK.yaml\n"
            "  Jarvis2 monitor\n"
            "  Jarvis2 plot PLOT.yaml\n"
            "  Jarvis2 portal …            # same CLI as jportal (V2 registry)\n"
            "  Jarvis2 operas list|info NAME\n"
            "  Jarvis2 project create|pack|list|browse|fetch|info …\n"
            "  Jarvis2 ps                  # list running Jarvis* processes\n"
            "  Jarvis2 kill [--yes]        # kill them (asks for confirmation)\n"
            "  Jarvis2 -v / --version      # logo + authors + package version\n"
            "\n"
            "Legacy aliases (still accepted):\n"
            "  Jarvis2 TASK.yaml [--resume]\n"
            "  Jarvis2 TASK.yaml --check-modules\n"
            "  Jarvis2 TASK.yaml --validate\n"
            "  Jarvis2 TASK.yaml --convert\n"
            "  Jarvis2 --monitor\n"
            "  Jarvis2 PLOT.yaml --plot   (deprecated)\n"
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Print Jarvis-HEP logo, author information and runtime package version",
    )
    sub = parser.add_subparsers(
        dest="command", required=False, parser_class=JarvisArgumentParser
    )

    def _add_logging_flags(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--console-level",
            default="WARNING",
            metavar="LEVEL",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help=(
                "Minimum screen severity: [cyan]DEBUG[/cyan], [green]INFO[/green], "
                "[yellow]WARNING[/yellow], [red]ERROR[/red], or [bold red]CRITICAL[/bold red] "
                "(default: WARNING). Files always retain DEBUG and above."
            ),
        )
        p.add_argument(
            "--silence",
            "-s",
            action="store_true",
            help="Disable console log output (files still written).",
        )

    run_p = sub.add_parser("run", help="Run a distributed scan task YAML")
    run_p.add_argument("task_yaml", help="Path to scan task YAML")
    run_p.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    run_p.add_argument(
        "--skip-draw-flowchart",
        action="store_true",
        help="Skip workflow flowchart.json / flowchart.png export",
    )
    run_p.add_argument(
        "--strict",
        action="store_true",
        help="Treat config validation warnings as errors",
    )
    _add_logging_flags(run_p)
    run_p.add_argument(
        "--debug",
        "-d",
        action="store_true",
        help="Enable DEBUG logs on the console",
    )

    check_p = sub.add_parser("check", help="Run fixed-point calculator smoke (check-modules)")
    check_p.add_argument("task_yaml", help="Path to check-modules task YAML")
    check_p.add_argument(
        "--strict",
        action="store_true",
        help="Treat config validation warnings as errors",
    )
    _add_logging_flags(check_p)

    validate_p = sub.add_parser(
        "validate",
        help="Validate task YAML (no Redis / workers); D13.9 config gate",
    )
    validate_p.add_argument("task_yaml", help="Path to scan task YAML")
    validate_p.add_argument(
        "--strict",
        action="store_true",
        help="Treat warnings as errors",
    )
    validate_p.add_argument(
        "--json",
        action="store_true",
        dest="validate_json",
        help="Emit a minimal JSON report on stdout",
    )

    convert_p = sub.add_parser(
        "convert",
        help="Refresh CSV snapshots when project DATABASE HDF5 files change",
    )
    convert_p.add_argument("task_yaml", help="Path to scan task YAML")

    sub.add_parser("monitor", help="Print one monitor snapshot and exit")

    plot_p = sub.add_parser("plot", help="Render a JarvisPLOT scene YAML")
    plot_p.add_argument("plot_yaml", help="Path to plot scene YAML")

    scene_p = sub.add_parser(
        "gen-plot-yaml",
        help="Generate images/<scan>/plot.yaml from a finished run task YAML",
    )
    scene_p.add_argument("task_yaml", help="Path to the finished scan task YAML")

    # ``portal`` is handled by argv passthrough in main() so that
    # ``Jarvis2 portal man|file|-h|-v`` matches the standalone jportal CLI.
    sub.add_parser(
        "portal",
        help="Jarvis-Portal CLI (same as jportal; uses V2 format registry)",
        add_help=False,
    )

    operas_p = sub.add_parser("operas", help="Jarvis-Operas discovery helpers")
    operas_sub = operas_p.add_subparsers(
        dest="operas_command", required=False, parser_class=JarvisArgumentParser
    )
    operas_sub.add_parser("list", help="List registered Operas operators")
    operas_info = operas_sub.add_parser("info", help="Show one Operas operator")
    operas_info.add_argument("name", help="Operator name or dotted path")

    # ``project`` uses argv passthrough (pack has many optional flags).
    sub.add_parser(
        "project",
        help="Standalone project tools (create/pack/browse/fetch/info)",
        add_help=False,
    )

    sub.add_parser(
        "ps",
        help="List running Jarvis OS processes (control/workers/archiver/redis titles)",
    )

    kill_p = sub.add_parser(
        "kill",
        help="Kill running Jarvis OS processes (prompts for confirmation)",
    )
    kill_p.add_argument(
        "--yes",
        "-y",
        action="store_true",
        help="Do not ask for confirmation",
    )
    kill_p.add_argument(
        "--no-force",
        action="store_true",
        help="Only SIGTERM; do not escalate to SIGKILL",
    )

    # ---- legacy top-level flags (paths go through normalize_argv → subcommands) ----
    parser.add_argument(
        "--monitor",
        action="store_true",
        help="(legacy) Print one monitor snapshot; prefer `Jarvis2 monitor`",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="(legacy) Treat path as plot scene; prefer `Jarvis2 plot`",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="(legacy) Validate task YAML only; prefer `Jarvis2 validate`",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="(legacy/global) Treat config validation warnings as errors",
    )
    parser.add_argument(
        "--check-modules",
        action="store_true",
        help="(legacy) Check-modules path; prefer `Jarvis2 check`",
    )
    parser.add_argument(
        "--convert",
        action="store_true",
        help="(legacy) Convert DATABASE HDF5 to CSV; prefer `Jarvis2 convert`",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="(legacy with bare TASK.yaml via normalize) Resume from checkpoint",
    )
    parser.add_argument(
        "--skip-draw-flowchart",
        action="store_true",
        help="(legacy) Skip flowchart export",
    )
    parser.add_argument(
        "--pid",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    return parser


def build_project_help_parser() -> JarvisArgumentParser:
    """Document project tools without changing their free-form dispatcher."""
    parser = JarvisArgumentParser(
        prog="Jarvis2 project",
        description="Create, package, browse, fetch, and inspect standalone Jarvis projects.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(
        dest="project_command", required=False, parser_class=JarvisArgumentParser
    )

    create = sub.add_parser("create", help="Scaffold a local project")
    create.add_argument("name", help="New project directory name")

    pack = sub.add_parser("pack", help="Package a project for sharing or release")
    pack.add_argument("path", nargs="?", help="Project path or pack manifest YAML")
    pack.add_argument("--share", action="store_true", help="Build a shareable package")
    pack.add_argument("--repro", action="store_true", help="Build a reproducible package")
    pack.add_argument("--full", action="store_true", help="Include full project contents")
    pack.add_argument("--man", action="store_true", help="Write a pack manifest instead")
    pack.add_argument("--encrypt", action="store_true", help="Encrypt the finished package")
    pack.add_argument("--key", help="Encryption key")

    for command, help_text in (
        ("list", "List projects in the official library"),
        ("browse", "Alias for list"),
    ):
        sub.add_parser(command, help=help_text)

    fetch = sub.add_parser("fetch", help="Download an official project")
    fetch.add_argument("name", help="Official project name")
    fetch.add_argument("--key", help="Key for a restricted project")

    info = sub.add_parser("info", help="Show official project details")
    info.add_argument("name", help="Official project name")

    encrypt = sub.add_parser("encrypt", help="Encrypt an existing project archive")
    encrypt.add_argument("archive", help="Existing .tar.gz archive")
    encrypt.add_argument("--key", required=True, help="Encryption key")
    return parser


def _render_help_for_args(parser: argparse.ArgumentParser, args: list[str]) -> int:
    """Render contextual parser help while keeping dispatchers integer-returning."""
    try:
        parser.parse_args(args)
    except SystemExit as exc:
        return int(exc.code)
    parser.print_help()
    return EXIT_OK


def _mode_conflict_message(modes: list[str]) -> str:
    return (
        "conflicting CLI intents: "
        + ", ".join(modes)
        + "; use a single subcommand "
        "(run|check|validate|convert|monitor|plot|gen-plot-yaml|portal|operas|project)"
    )


def resolve_intent(args: argparse.Namespace) -> tuple[str, argparse.Namespace]:
    """Normalize subcommand + legacy flags into one intent name.

    Raises SystemExit-ready ValueError with usage semantics for conflicts.
    """
    if getattr(args, "pid", None) is not None:
        raise ValueError(
            "--pid is not implemented and has been retired; "
            "use `Jarvis2 monitor` (Redis) or attach by scan outputs"
        )

    command = getattr(args, "command", None)
    legacy_modes: list[str] = []
    if args.monitor:
        legacy_modes.append("monitor")
    if args.plot:
        legacy_modes.append("plot")
    if args.check_modules:
        legacy_modes.append("check")
    if getattr(args, "validate", False):
        legacy_modes.append("validate")
    if getattr(args, "convert", False):
        legacy_modes.append("convert")

    if command:
        if legacy_modes:
            conflict = [command, *legacy_modes]
            if len(set(conflict)) > 1:
                raise ValueError(_mode_conflict_message(sorted(set(conflict))))
        return str(command), args

    if len(legacy_modes) > 1:
        raise ValueError(_mode_conflict_message(sorted(legacy_modes)))
    if len(legacy_modes) == 1:
        return legacy_modes[0], args
    raise ValueError(
        "Provide a task YAML or a subcommand "
        "(run|check|validate|convert|monitor|plot|portal|operas|project). "
        "See Jarvis2 -h."
    )


def run_monitor(
    *,
    factory: TaskFactory | None = None,
    redis: RedisQueue | None = None,
) -> int:
    reader = attach_reader(factory=factory, redis=redis)
    view = reader.read()
    if not view.has_active_scan():
        print("No active scan found.", file=sys.stderr)
        return EXIT_RUN_FAILED
    sys.stdout.write(format_monitor_view(view))
    return EXIT_OK


def dispatch_version() -> int:
    """Print V1-style logo banner with injected runtime version (no scan init)."""
    from jarvishep2.versioning import render_logo_with_version

    print(render_logo_with_version())
    return EXIT_OK


def dispatch_monitor(_args: argparse.Namespace) -> int:
    redis_config = dict(INTERNAL_REDIS_CONFIG)
    factory = TaskFactory(redis_config)
    try:
        factory.init_redis()
        if factory.get_monitor_snapshot():
            return run_monitor(factory=factory)
    except Exception:
        pass

    redis = RedisQueue(redis_config)
    try:
        redis.connect()
    except Exception as exc:
        print(f"Unable to connect to Redis for monitor attach: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    try:
        return run_monitor(redis=redis)
    finally:
        redis.close()


def dispatch_plot(plot_yaml: str, *, legacy: bool = False) -> int:
    if not plot_yaml:
        print("Plot YAML is required (usage: Jarvis2 plot scene.yaml).", file=sys.stderr)
        return EXIT_USAGE
    if legacy:
        print(
            "warning: `Jarvis2 <plot.yaml> --plot` is deprecated; use `Jarvis2 plot <plot.yaml>`",
            file=sys.stderr,
        )
    try:
        return int(run_plot(plot_yaml))
    except ImportError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE
    except PlotBridgeError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED


def dispatch_gen_plot_yaml(task_yaml: str) -> int:
    """Generate (but never render) JarvisPLOT YAML for one finished scan."""
    if not task_yaml:
        print("Task YAML is required (usage: Jarvis2 gen-plot-yaml TASK.yaml).", file=sys.stderr)
        return EXIT_USAGE
    try:
        core = Jarvis2Core()
        core.load_task_yaml(task_yaml, validate=False)
        from jarvishep2.base import infer_project_root_from_task_result_dir
        from jarvishep2.plot_scene import emit_jplot_scan_levelset_yaml

        task_result_dir = str(core.config.get("task_result_dir") or "").strip()
        scan_name = str(core.config.get("scan_name") or "scan").strip() or "scan"
        if not task_result_dir:
            raise ValueError("task YAML does not resolve a task_result_dir")
        project_root = str(
            core.config.get("project_root") or core.config.get("task_root") or ""
        ).strip() or infer_project_root_from_task_result_dir(task_result_dir)
        output = emit_jplot_scan_levelset_yaml(
            task_result_dir,
            scan_name=scan_name,
            project_root=project_root,
            config=core.config,
            yaml_basename="plot.yaml",
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except (ValueError, OSError) as exc:
        print(f"gen-plot-yaml failed: {exc}", file=sys.stderr)
        return EXIT_USAGE
    if not output:
        print("No plot YAML could be generated: run outputs were not found.", file=sys.stderr)
        return EXIT_RUN_FAILED
    print(output)
    return EXIT_OK


def dispatch_portal(portal_argv: list[str] | None = None) -> int:
    """Forward to Jarvis-Portal's CLI with the **V2** registry surface.

    Interface matches standalone ``jportal``::

        Jarvis2 portal man
        Jarvis2 portal man json
        Jarvis2 portal man slha
        Jarvis2 portal file.yaml
        Jarvis2 portal -h
        Jarvis2 portal -v

    Convenience alias: ``Jarvis2 portal formats`` → ``man`` (list formats).
    """
    argv = list(portal_argv or [])
    if argv == ["formats"]:
        argv = ["man"]
    try:
        from jarvis_portal.cli import main as portal_main
        from jarvis_portal.v2 import create_default_registry as v2_registry
    except ImportError as exc:
        print(
            "Jarvis-Portal is required for `Jarvis2 portal`. "
            "Install with `pip install -U Jarvis-HEP-Portal` "
            f"(or `pip install -e ../Jarvis-Portal`). Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    try:
        return int(
            portal_main(
                argv,
                registry_factory=v2_registry,
                prog="Jarvis2 portal",
            )
        )
    except Exception as exc:
        print(f"Jarvis2 portal failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED


def _human_bytes(value: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    amount = float(max(0, int(value)))
    idx = 0
    while amount >= 1024.0 and idx < len(units) - 1:
        amount /= 1024.0
        idx += 1
    if idx == 0:
        return f"{int(amount)} {units[idx]}"
    return f"{amount:.2f} {units[idx]}"


def _print_kv_block(title: str, rows: list[tuple[str, object]]) -> None:
    print(title)
    width = max((len(str(key)) for key, _ in rows), default=8)
    for key, value in rows:
        print(f"  {str(key):<{width}}  {value}")


def _print_project_help() -> None:
    build_project_help_parser().print_help()


def _run_project_create(project_name: str) -> int:
    from jarvishep2.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold

    try:
        project_root = create_project_scaffold(project_name, cwd=os.getcwd())
    except ValueError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except FileExistsError as exc:
        print(f"[Jarvis2] Project directory already exists: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    print(f"[Jarvis2] Project scaffold created at: {project_root}")
    print(f"[Jarvis2] Created folders: {', '.join(PROJECT_SUBDIRS)}")
    print("[Jarvis2] Try: Jarvis2 run bin/quickstart_bridson_operas.yaml")
    return EXIT_OK


def _run_project_pack(
    project_path: str | None,
    profile: str,
    *,
    encrypt: bool = False,
    key: str | None = None,
) -> int:
    from jarvishep2.project_crypto import ProjectCryptoError, encrypt_file, resolve_fetch_key
    from jarvishep2.project_packager import ProjectPackError, create_project_package

    if encrypt and not resolve_fetch_key(key):
        print(
            "[Jarvis2] --encrypt requires --key KEY or JARVIS_PROJECT_FETCH_KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE

    try:
        report = create_project_package(project_root=project_path, profile=profile)
    except ProjectPackError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis2] Failed to package project: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    archive_path = report.archive_path
    encrypted_path = None
    if encrypt:
        try:
            encrypted_path = archive_path + ".jenc"
            if archive_path.endswith(".tar.gz"):
                encrypted_path = archive_path[: -len(".tar.gz")] + ".tar.gz.jenc"
            encrypt_file(archive_path, encrypted_path, key=str(resolve_fetch_key(key)))
        except ProjectCryptoError as exc:
            print(f"[Jarvis2] {exc}", file=sys.stderr)
            return EXIT_RUN_FAILED

    rows = [
        ("Project root", report.project_root),
        ("Archive", report.archive_path),
        ("Mode", report.profile),
        ("Packed files", report.included_files),
        ("Skipped files", report.excluded_files),
        ("Payload size", _human_bytes(report.total_bytes)),
    ]
    if encrypted_path:
        rows.append(("Encrypted", encrypted_path))
        rows.append(("Fetch with", f"Jarvis2 project fetch <name> --key <key>"))
    _print_kv_block("Jarvis2 project package created", rows)
    return EXIT_OK


def _run_project_encrypt(archive_path: str, *, key: str | None) -> int:
    from jarvishep2.project_crypto import ProjectCryptoError, encrypt_file, resolve_fetch_key

    resolved = resolve_fetch_key(key)
    if not resolved:
        print(
            "[Jarvis2] encrypt requires --key KEY or JARVIS_PROJECT_FETCH_KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE
    path = os.path.abspath(os.path.expanduser(archive_path))
    if not os.path.isfile(path):
        print(f"[Jarvis2] Archive not found: {path}", file=sys.stderr)
        return EXIT_USAGE
    if path.endswith(".jenc"):
        print("[Jarvis2] File already looks encrypted (.jenc)", file=sys.stderr)
        return EXIT_USAGE
    out = path + ".jenc"
    if path.endswith(".tar.gz"):
        out = path[: -len(".tar.gz")] + ".tar.gz.jenc"
    try:
        encrypt_file(path, out, key=resolved)
    except ProjectCryptoError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    _print_kv_block(
        "Jarvis2 project archive encrypted",
        [
            ("Input", path),
            ("Encrypted", out),
            ("Scheme", "openssl-aes-256-cbc"),
            ("Users fetch with", "Jarvis2 project fetch NAME --key …"),
        ],
    )
    return EXIT_OK


def _run_project_pack_manifest(project_path: str | None, profile: str) -> int:
    from jarvishep2.project_packager import ProjectPackError, create_project_pack_manifest

    try:
        report = create_project_pack_manifest(project_root=project_path, profile=profile)
    except ProjectPackError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis2] Failed to write project pack manifest: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    _print_kv_block(
        "Jarvis2 project pack manifest created",
        [
            ("Project root", report.project_root),
            ("Manifest", report.manifest_path),
            ("Pack ID", report.pack_id),
            ("Mode", report.profile),
            ("Manifest files", report.included_files),
            ("Manifest excludes", report.excluded_files),
            ("Output", report.output),
        ],
    )
    return EXIT_OK


def _run_project_pack_from_manifest(manifest_path: str) -> int:
    from jarvishep2.project_packager import (
        ProjectPackError,
        create_project_package_from_manifest,
    )

    try:
        report = create_project_package_from_manifest(manifest_path)
    except ProjectPackError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis2] Failed to package project from manifest: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    _print_kv_block(
        "Jarvis2 project package created from manifest",
        [
            ("Project root", report.project_root),
            ("Archive", report.archive_path),
            ("Mode", report.profile),
            ("Packed files", report.included_files),
            ("Skipped files", report.excluded_files),
            ("Payload size", _human_bytes(report.total_bytes)),
        ],
    )
    return EXIT_OK


def _run_project_browse() -> int:
    from jarvishep2.official_project_library import (
        OfficialLibraryError,
        format_project_list_table,
        list_official_projects,
    )

    try:
        projects = list_official_projects()
    except OfficialLibraryError as exc:
        print(f"[Jarvis2] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    if not projects:
        print("[Jarvis2] No verified projects are listed in the official library.")
        return EXIT_OK

    print("[Jarvis2] Official library (catalog: Jarvis-Examples GitHub JSON)")
    print(format_project_list_table(projects))
    return EXIT_OK


def _run_project_info(project_name: str) -> int:
    from jarvishep2.official_project_library import (
        OfficialLibraryError,
        OfficialProjectNotFoundError,
        get_official_project,
    )

    try:
        project = get_official_project(project_name)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except OfficialLibraryError as exc:
        print(f"[Jarvis2] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    needs_key = bool(project.get("requires_key"))
    print(f"[Jarvis2] Official project: {project['name']}")
    print(f"[Jarvis2] Summary: {project.get('summary') or 'N/A'}")
    print(f"[Jarvis2] Category: {project.get('category') or 'N/A'}")
    print(f"[Jarvis2] Access: {project.get('access') or 'public'}")
    print(f"[Jarvis2] Key required: {'yes' if needs_key else 'no'}")
    if needs_key and project.get("encryption_hint"):
        print(f"[Jarvis2] Key hint: {project.get('encryption_hint')}")
    if needs_key:
        print(f"[Jarvis2] Encryption: {project.get('encryption_scheme') or 'openssl-aes-256-cbc'}")
    print(f"[Jarvis2] Entrypoint: {project.get('entrypoint') or 'N/A'}")
    print(f"[Jarvis2] Compatibility notes: {project.get('compatibility_notes') or 'None'}")
    return EXIT_OK


def _run_project_fetch(project_name: str, *, key: str | None = None) -> int:
    from jarvishep2.official_project_library import (
        OfficialLibraryError,
        OfficialProjectFetchError,
        OfficialProjectNotFoundError,
        fetch_official_project,
    )

    try:
        report = fetch_official_project(project_name, key=key)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except OfficialProjectFetchError as exc:
        print(f"[Jarvis2] {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    except OfficialLibraryError as exc:
        print(f"[Jarvis2] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    print(
        f"[Jarvis2] Official project '{report.project_name}' fetched to: {report.target_dir}"
    )
    print(f"[Jarvis2] Access: {report.access}")
    print(f"[Jarvis2] Entrypoint: {report.entrypoint or 'N/A'}")
    return EXIT_OK


def _parse_fetch_arguments(tokens: list[str]) -> tuple[str, str | None] | int:
    """Parse ``fetch <name> [--key KEY]``."""
    if not tokens or (len(tokens) == 1 and tokens[0] in _HELP_FLAGS):
        print(
            "usage: Jarvis2 project fetch <name> [--key KEY]\n"
            "  restricted projects need --key or JARVIS_PROJECT_FETCH_KEY\n"
        )
        return 0 if tokens and tokens[0] in _HELP_FLAGS else EXIT_USAGE

    name: str | None = None
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis2] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis2] Unsupported option for project fetch: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if name is None:
            name = tok
            i += 1
            continue
        print(f"[Jarvis2] Unexpected argument for project fetch: {tok}", file=sys.stderr)
        return EXIT_USAGE

    if not name:
        print("[Jarvis2] Usage: Jarvis2 project fetch <name> [--key KEY]", file=sys.stderr)
        return EXIT_USAGE
    return name, key


def _looks_like_yaml_path(path: str | None) -> bool:
    if path is None:
        return False
    return path.lower().endswith((".yaml", ".yml"))


def _parse_pack_arguments(
    tokens: list[str],
) -> tuple[str | None, str, bool, bool, bool, str | None] | int:
    path: str | None = None
    mode_flag: str | None = None
    manifest_only = False
    encrypt = False
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in _HELP_FLAGS:
            print(
                "usage: Jarvis2 project pack [path] [--share|--repro|--full] [--man]\n"
                "       Jarvis2 project pack [path] --encrypt --key KEY\n"
                "       Jarvis2 project pack <pack_manifest.yaml>\n"
            )
            return 0
        if tok == _PACK_MANIFEST_FLAG:
            manifest_only = True
            i += 1
            continue
        if tok in _PACK_ENCRYPT_FLAGS:
            encrypt = True
            i += 1
            continue
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis2] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok in _PACK_MODE_FLAGS:
            if mode_flag is not None:
                print(
                    "[Jarvis2] project pack modes are mutually exclusive: "
                    "--share, --repro, --full",
                    file=sys.stderr,
                )
                return EXIT_USAGE
            mode_flag = tok
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis2] Unsupported option for project pack: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if path is None:
            path = tok
            i += 1
            continue
        print(f"[Jarvis2] Unexpected argument for project pack: {tok}", file=sys.stderr)
        return EXIT_USAGE

    if manifest_only and _looks_like_yaml_path(path):
        print(
            "[Jarvis2] --man writes a new manifest from a project path, not from a manifest file",
            file=sys.stderr,
        )
        return EXIT_USAGE

    is_manifest_input = _looks_like_yaml_path(path)
    if is_manifest_input and mode_flag is not None:
        print(
            "[Jarvis2] Manifest packing does not accept --share, --repro, or --full",
            file=sys.stderr,
        )
        return EXIT_USAGE
    if is_manifest_input and encrypt:
        print(
            "[Jarvis2] --encrypt is for packing a project path; "
            "encrypt a finished archive with: Jarvis2 project encrypt FILE --key K",
            file=sys.stderr,
        )
        return EXIT_USAGE

    profile = _PACK_MODE_FLAGS.get(mode_flag or "", "share")
    return path, profile, manifest_only, is_manifest_input, encrypt, key


def _parse_encrypt_arguments(tokens: list[str]) -> tuple[str, str | None] | int:
    if not tokens or (len(tokens) == 1 and tokens[0] in _HELP_FLAGS):
        print("usage: Jarvis2 project encrypt <archive.tar.gz> --key KEY\n")
        return 0 if tokens and tokens[0] in _HELP_FLAGS else EXIT_USAGE
    path: str | None = None
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis2] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis2] Unsupported option for project encrypt: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if path is None:
            path = tok
            i += 1
            continue
        print(f"[Jarvis2] Unexpected argument: {tok}", file=sys.stderr)
        return EXIT_USAGE
    if not path:
        print(
            "[Jarvis2] Usage: Jarvis2 project encrypt <archive.tar.gz> --key KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE
    return path, key


def dispatch_project(project_argv: list[str] | None = None) -> int:
    """Handle ``Jarvis2 project create|pack|browse|fetch|info`` (D12.5)."""
    args = list(project_argv or [])
    if not args or args[0] in _HELP_FLAGS:
        _print_project_help()
        return EXIT_OK
    if any(arg in _HELP_FLAGS for arg in args):
        return _render_help_for_args(build_project_help_parser(), args)

    command = args[0]
    rest = args[1:]
    if command not in _PROJECT_COMMANDS:
        print(f"[Jarvis2] Unknown project command: {command}", file=sys.stderr)
        _print_project_help()
        return EXIT_USAGE

    if command == "create":
        if len(rest) == 1 and rest[0] in _HELP_FLAGS:
            print("usage: Jarvis2 project create <name>\n")
            return EXIT_OK
        if len(rest) != 1 or rest[0].startswith("-"):
            print("[Jarvis2] Usage: Jarvis2 project create <name>", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_create(rest[0])

    if command == "pack":
        parsed = _parse_pack_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        path, profile, manifest_only, is_manifest_input, encrypt, key = parsed
        if manifest_only:
            return _run_project_pack_manifest(path, profile)
        if is_manifest_input:
            return _run_project_pack_from_manifest(str(path))
        return _run_project_pack(path, profile, encrypt=encrypt, key=key)

    if command == "encrypt":
        parsed = _parse_encrypt_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        archive, key = parsed
        return _run_project_encrypt(archive, key=key)

    if command in {"browse", "list"}:
        if rest:
            if len(rest) == 1 and rest[0] in _HELP_FLAGS:
                print("usage: Jarvis2 project list|browse\n")
                return EXIT_OK
            print("[Jarvis2] Usage: Jarvis2 project list|browse", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_browse()

    if command == "fetch":
        parsed = _parse_fetch_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        name, key = parsed
        return _run_project_fetch(name, key=key)

    if command == "info":
        if len(rest) == 1 and rest[0] in _HELP_FLAGS:
            print("usage: Jarvis2 project info <name>\n")
            return EXIT_OK
        if len(rest) != 1 or rest[0].startswith("-"):
            print("[Jarvis2] Usage: Jarvis2 project info <name>", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_info(rest[0])

    return EXIT_USAGE


def dispatch_operas(args: argparse.Namespace) -> int:
    cmd = getattr(args, "operas_command", None)
    if cmd is None:
        print("Usage: Jarvis2 operas list | Jarvis2 operas info NAME", file=sys.stderr)
        return EXIT_USAGE
    try:
        from jarvishep2.operas_functions import ensure_operas_registry_discovered

        registry = ensure_operas_registry_discovered()
    except Exception as exc:
        print(f"Unable to load Operas registry: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    names: list[str] = []
    for attr in ("list_functions", "list", "names", "keys"):
        if hasattr(registry, attr):
            try:
                raw = getattr(registry, attr)
                names = list(raw() if callable(raw) else raw)
                break
            except Exception:
                continue
    if not names and isinstance(registry, dict):
        names = sorted(str(k) for k in registry.keys())
    names = [str(n) for n in names]

    if cmd == "list":
        for name in sorted(names):
            print(name)
        return EXIT_OK if names else EXIT_RUN_FAILED

    if cmd == "info":
        name = str(getattr(args, "name", "") or "").strip()
        if not name:
            print("operas info requires NAME", file=sys.stderr)
            return EXIT_USAGE
        target = None
        if hasattr(registry, "get"):
            try:
                target = registry.get(name)
            except Exception:
                target = None
        if target is None and isinstance(registry, dict):
            target = registry.get(name)
        if target is None and name not in names:
            print(f"operator not found: {name}", file=sys.stderr)
            return EXIT_RUN_FAILED
        print(f"name: {name}")
        if target is not None:
            print(f"repr: {target!r}")
        else:
            print("registered: true")
        return EXIT_OK

    print(f"Unknown operas subcommand: {cmd}", file=sys.stderr)
    return EXIT_USAGE


def _print_outcome(outcome: RunOutcome) -> None:
    """Human-readable one-line summary on stderr (stdout stays clean for pipes)."""
    ncall = None
    if isinstance(outcome.extras, dict):
        ncall = outcome.extras.get("ncall")
    ncall_bit = f" ncall={ncall}" if ncall is not None else ""
    print(
        f"RunOutcome status={outcome.status} "
        f"submitted={outcome.submitted} completed={outcome.completed} "
        f"failed={outcome.failed} archived={outcome.archived}"
        f"{ncall_bit} "
        f"exit={outcome.exit_code}",
        file=sys.stderr,
    )
    if outcome.error:
        print(f"  error_type={outcome.error_type} error={outcome.error}", file=sys.stderr)


def dispatch_validate(
    task_yaml: str,
    *,
    strict: bool = False,
    as_json: bool = False,
    check_modules: bool | None = None,
) -> int:
    """Load + validate a task card; never start Redis or workers."""
    import json

    from jarvishep2.task_validation import (
        ConfigValidationError,
        format_report,
        format_report_success,
        validate_task_config,
    )

    if not task_yaml:
        print("Task YAML is required.", file=sys.stderr)
        return EXIT_USAGE

    try:
        core = Jarvis2Core()
        # Load without gate, then validate once with full control over output.
        core.load_task_yaml(task_yaml, validate=False)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    report = validate_task_config(
        core.config, strict=strict, check_modules=check_modules
    )
    if as_json:
        payload = {
            "ok": report.ok,
            "task_yaml": str(core.config.get("task_yaml") or task_yaml),
            "scan_name": str(core.config.get("scan_name") or ""),
            "issues": [
                {
                    "level": i.level,
                    "code": i.code,
                    "path": i.path,
                    "message": i.message,
                    "hint": i.hint,
                }
                for i in report.issues
            ],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return EXIT_OK if report.ok else EXIT_USAGE

    if report.ok:
        # Mirror run-path logging: explicit success line on stderr for humans.
        print(
            "Config validation successful. The task YAML meets the V2 contract.",
            file=sys.stderr,
        )
        if report.warnings():
            print(format_report(report), file=sys.stderr)
        else:
            # Keep format_report_success available for tests / quieter tooling.
            _ = format_report_success(report)
        return EXIT_OK

    print(format_report(report), file=sys.stderr)
    return EXIT_USAGE


def dispatch_convert(task_yaml: str) -> int:
    """Refresh DATABASE HDF5 → CSV snapshots when their source content changes."""
    if not task_yaml:
        print("Task YAML is required for convert.", file=sys.stderr)
        return EXIT_USAGE

    try:
        core = Jarvis2Core()
        # Post-run utility: path resolution only; skip full config gate.
        core.load_task_yaml(task_yaml, validate=False)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    results = core.convert()
    if not results:
        return EXIT_RUN_FAILED

    if any(item.get("status") in {"converted", "skipped_unchanged"} for item in results):
        return EXIT_OK
    return EXIT_RUN_FAILED


def dispatch_run(
    task_yaml: str,
    *,
    resume: bool = False,
    check_modules: bool = False,
    skip_draw_flowchart: bool = False,
    console_level: str = "WARNING",
    silence: bool = False,
    debug: bool = False,
    strict: bool = False,
) -> int:
    if not task_yaml:
        print("Task YAML is required.", file=sys.stderr)
        return EXIT_USAGE

    try:
        core = Jarvis2Core()
        core.load_task_yaml(
            task_yaml,
            validate=True,
            strict=strict,
            check_modules=check_modules if check_modules else None,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except ValueError as exc:
        # Includes ConfigValidationError (subclass of ValueError).
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    if debug:
        console_level = "DEBUG"

    core.set_logging_options(
        console_level=console_level,
        silence=silence,
    )

    if skip_draw_flowchart:
        core.config["skip_draw_flowchart"] = True
        core.info["skip_draw_flowchart"] = True

    try:
        if check_modules:
            outcome = core.check_modules()
        else:
            outcome = core.run(resume=resume)
    except NotImplementedError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE
    except (RuntimeError, TimeoutError, OSError) as exc:
        print(f"Run failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    if not isinstance(outcome, RunOutcome):
        # Backward safety if a test double returns a bare int.
        try:
            count = int(outcome)  # type: ignore[arg-type]
        except Exception:
            count = 0
        return EXIT_OK if count > 0 else EXIT_RUN_FAILED

    _print_outcome(outcome)
    return int(outcome.exit_code)


def dispatch(args: argparse.Namespace) -> int:
    if bool(getattr(args, "version", False)):
        return dispatch_version()

    try:
        intent, args = resolve_intent(args)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    if intent == "monitor":
        return dispatch_monitor(args)
    if intent == "plot":
        plot_yaml = getattr(args, "plot_yaml", None) or args.task_yaml
        legacy = getattr(args, "command", None) is None
        return dispatch_plot(str(plot_yaml or ""), legacy=legacy)
    if intent == "gen-plot-yaml":
        return dispatch_gen_plot_yaml(str(getattr(args, "task_yaml", "") or ""))
    if intent == "portal":
        # Prefer argv passthrough via main(); this path is a thin fallback.
        return dispatch_portal([])
    if intent == "operas":
        return dispatch_operas(args)
    if intent == "project":
        return dispatch_project([])
    if intent == "ps":
        from jarvishep2.process_cleanup import list_running_jarvis_cli

        return int(list_running_jarvis_cli())
    if intent == "kill":
        from jarvishep2.process_cleanup import kill_running_jarvis_cli

        return int(
            kill_running_jarvis_cli(
                yes=bool(getattr(args, "yes", False)),
                force=not bool(getattr(args, "no_force", False)),
            )
        )
    if intent == "validate":
        task = getattr(args, "task_yaml", None)
        return dispatch_validate(
            str(task or ""),
            strict=bool(getattr(args, "strict", False)),
            as_json=bool(getattr(args, "validate_json", False)),
        )
    if intent == "convert":
        task = getattr(args, "task_yaml", None)
        return dispatch_convert(str(task or ""))
    if intent == "check":
        task = getattr(args, "task_yaml", None)
        return dispatch_run(
            str(task or ""),
            check_modules=True,
            skip_draw_flowchart=bool(getattr(args, "skip_draw_flowchart", False)),
            console_level=str(getattr(args, "console_level", "WARNING") or "WARNING"),
            silence=bool(getattr(args, "silence", False)),
            strict=bool(getattr(args, "strict", False)),
        )
    if intent == "run":
        task = getattr(args, "task_yaml", None)
        return dispatch_run(
            str(task or ""),
            resume=bool(args.resume),
            check_modules=False,
            skip_draw_flowchart=bool(getattr(args, "skip_draw_flowchart", False)),
            console_level=str(getattr(args, "console_level", "WARNING") or "WARNING"),
            silence=bool(getattr(args, "silence", False)),
            debug=bool(getattr(args, "debug", False)),
            strict=bool(getattr(args, "strict", False)),
        )

    print(f"Unknown intent: {intent}", file=sys.stderr)
    return EXIT_USAGE


def main(argv: list[str] | None = None) -> int:
    from jarvishep2.proc_title import control_title, set_process_title

    set_process_title(control_title())
    raw = list(sys.argv[1:] if argv is None else argv)
    normalized = normalize_argv(raw)
    # Full jportal-compatible surface: do not let argparse eat portal args.
    if normalized and normalized[0] == "portal":
        return dispatch_portal(normalized[1:])
    # Project pack flags are free-form; pass through like portal.
    if normalized and normalized[0] == "project":
        return dispatch_project(normalized[1:])
    return dispatch(build_parser().parse_args(normalized))


if __name__ == "__main__":
    raise SystemExit(main())
