#!/usr/bin/env python3
"""CLI help chrome and argparse schema (D25.9)."""

from __future__ import annotations

import argparse
import contextlib
import io
import shutil
import sys
from typing import Any

import click
import typer.rich_utils
from rich import box
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.run_outcome import EXIT_OK

_LEGACY_OPTION_DESTS = frozenset(
    {
        "monitor",
        "plot",
        "validate",
        "strict",
        "convert",
        "resume",
        "skip_draw_flowchart",
    }
)

_COMMAND_HELP_PANELS = {
    "run": "Scan workflow",
    "check": "Scan workflow",
    "validate": "Scan workflow",
    "man": "YAML authoring",
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
    "man": 5,
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
    "browse": 160,
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
        # Usage benefits from a semantic metavar (``<task_yaml>``), while the
        # Arguments panel already has that name in its first column. Keep the
        # type compact there so the fixed help grid never truncates it.
        metavar = "<str>" if param.name == "task_yaml" else param.make_metavar(ctx=ctx)
        return (
            typer.rich_utils.highlighter(param.name or ""),
            typer.rich_utils.metavar_highlighter(metavar),
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
        if action.dest == "task_yaml":
            return "<task_yaml>"
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
        root_help = self.prog == "Jarvis"
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

def build_parser() -> argparse.ArgumentParser:
    parser = JarvisArgumentParser(
        prog="Jarvis",
        description=(
            "Run and validate distributed HEP scans, prepare plot YAML, "
            "and manage their local runtime.\n\n"
            "Command help: Jarvis COMMAND -h\n"
            "See each command's arguments and options, e.g. Jarvis run -h"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Intents (preferred):\n"
            "  Jarvis run TASK.yaml [--resume]\n"
            "  Jarvis check TASK.yaml\n"
            "  Jarvis validate TASK.yaml [--strict] [--json]\n"
            "  Jarvis man [yaml|sampler|calculator|…]  # YAML writing manuals\n"
            "  Jarvis convert TASK.yaml\n"
            "  Jarvis monitor\n"
            "  Jarvis plot …              # same CLI as jplot\n"
            "  Jarvis portal …            # same CLI as jportal (V2 registry)\n"
            "  Jarvis operas …            # same CLI as jopera\n"
            "  Jarvis project create|pack|browse|fetch|info …\n"
            "  Jarvis ps                  # list running Jarvis* processes\n"
            "  Jarvis kill [--yes]        # kill them (asks for confirmation)\n"
            "  Jarvis -v / --version      # logo + authors + package version\n"
            "  Jarvis --refs              # framework and sampler references\n"
            "\n"
            "Legacy aliases (still accepted):\n"
            "  Jarvis TASK.yaml [--resume]\n"
            "  Jarvis TASK.yaml --validate\n"
            "  Jarvis TASK.yaml --convert\n"
            "  Jarvis --monitor\n"
            "  Jarvis PLOT.yaml --plot   (deprecated)\n"
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        help="Print Jarvis-HEP logo, author information and runtime package version",
    )
    parser.add_argument(
        "--refs",
        action="store_true",
        help="Print framework and built-in sampler references for citation",
    )
    sub = parser.add_subparsers(
        dest="command", required=False, parser_class=JarvisArgumentParser
    )

    def _add_logging_flags(
        p: argparse.ArgumentParser, *, default_console_level: str = "WARNING"
    ) -> None:
        p.add_argument(
            "--console-level",
            default=default_console_level,
            metavar="LEVEL",
            choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            help=(
                "Minimum screen severity: [cyan]DEBUG[/cyan], [green]INFO[/green], "
                "[yellow]WARNING[/yellow], [red]ERROR[/red], or [bold red]CRITICAL[/bold red] "
                f"(default: {default_console_level}). Files always retain DEBUG and above."
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
        "--skip-library-installation",
        action="store_true",
        help="Do not build LibDeps; require every declared library path to already exist",
    )
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

    check_p = sub.add_parser("check", help="Run fixed-point calculator smoke")
    check_p.add_argument("task_yaml", help="Path to the task YAML to smoke-test")
    check_p.add_argument(
        "--skip-library-installation",
        action="store_true",
        help="Do not build LibDeps; require every declared library path to already exist",
    )
    check_p.add_argument(
        "--strict",
        action="store_true",
        help="Treat config validation warnings as errors",
    )
    check_p.add_argument(
        "--timeout",
        type=float,
        default=None,
        metavar="SEC",
        help=(
            "Max seconds to wait for smoke samples to archive (deadline; default "
            "from EnvReqs.V2.check_modules.timeout_sec, else 120)"
        ),
    )
    _add_logging_flags(check_p, default_console_level="DEBUG")

    sub.add_parser(
        "man",
        help=(
            "YAML writing manuals; coding agents: use --json for structured keys, "
            "requirements, examples, and actions"
        ),
        add_help=False,
    )

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

    monitor_p = sub.add_parser(
        "monitor", help="List running scans, or print one selected scan monitor snapshot"
    )
    monitor_p.add_argument(
        "scan_ref",
        nargs="?",
        help="Running scan reference (sticky R1/R2, control PID, or Scan.name)",
    )

    # ``plot`` is handled by argv passthrough in main(), matching jplot.
    sub.add_parser("plot", help="Jarvis-PLOT CLI (same as jplot)", add_help=False)

    scene_p = sub.add_parser(
        "gen-plot-yaml",
        help="Generate images/<scan>/plot.yaml from a finished run task YAML",
    )
    scene_p.add_argument("task_yaml", help="Path to the finished scan task YAML")

    # ``portal`` is handled by argv passthrough in main() so that
    # ``Jarvis portal man|file|-h|-v`` matches the standalone jportal CLI.
    sub.add_parser(
        "portal",
        help="Jarvis-Portal CLI (same as jportal; uses V2 format registry)",
        add_help=False,
    )

    # ``operas`` is handled by argv passthrough in main(), matching jopera.
    sub.add_parser("operas", help="Jarvis-Operas CLI (same as jopera)", add_help=False)

    # ``project`` uses argv passthrough (pack has many optional flags).
    sub.add_parser(
        "project",
        help="Standalone project tools (create/pack/browse/fetch/info)",
        add_help=False,
    )

    ps_p = sub.add_parser(
        "ps",
        help="List process groups, or show one scan / ZP group",
    )
    ps_p.add_argument(
        "scan_ref",
        nargs="?",
        help="Process-group reference (sticky R1/R2, Scan.name, control PID, or ZP)",
    )

    kill_p = sub.add_parser(
        "kill",
        help="List process groups, or terminate one selected scan / ZP group",
    )
    kill_p.add_argument(
        "scan_ref",
        nargs="?",
        help="Process-group reference (sticky R1/R2, Scan.name, control PID, or ZP)",
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
        help="(legacy) List running scans; prefer `Jarvis monitor`",
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="(legacy) Treat path as plot scene; prefer `Jarvis plot`",
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="(legacy) Validate task YAML only; prefer `Jarvis validate`",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="(legacy/global) Treat config validation warnings as errors",
    )
    parser.add_argument(
        "--convert",
        action="store_true",
        help="(legacy) Convert DATABASE HDF5 to CSV; prefer `Jarvis convert`",
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
        prog="Jarvis project",
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

    sub.add_parser("browse", help="Browse projects in the official library")

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

