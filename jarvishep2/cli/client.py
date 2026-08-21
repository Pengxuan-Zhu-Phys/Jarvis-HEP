#!/usr/bin/env python3
"""Jarvis CLI surface: argv router (D25.9).

Help chrome lives in ``jarvishep2.cli.help``, project tools in
``jarvishep2.cli.project``, and validate/run/check/convert in
``jarvishep2.cli.run``. ``main()`` remains the entry point.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from jarvishep2.cli.help import (
    JarvisArgumentParser,
    build_parser,
    build_project_help_parser,
)
from jarvishep2.cli.project import dispatch_project
from jarvishep2.cli.run import (
    _emit_load_diagnostic,
    _log_outcome,
    dispatch_convert,
    dispatch_run,
    dispatch_validate,
)
from jarvishep2.references import render_references
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE

_SUBCOMMANDS = frozenset(
    {
        "run",
        "check",
        "validate",
        "man",
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

def normalize_argv(argv: list[str] | None) -> list[str]:
    """Rewrite legacy bare-YAML invocations into explicit subcommands.

    Argparse subparsers treat the first positional as the command name, so
    ``Jarvis task.yaml`` would otherwise fail as an invalid choice.
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
    if "--validate" in args:
        rest = [a for a in args[1:] if a != "--validate"]
        return ["validate", head, *rest]
    if "--convert" in args:
        rest = [a for a in args[1:] if a != "--convert"]
        return ["convert", head, *rest]
    # Bare path (+ optional --resume) → run
    return ["run", *args]


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
            "use `Jarvis monitor` (Redis) or attach by scan outputs"
        )

    command = getattr(args, "command", None)
    legacy_modes: list[str] = []
    if args.monitor:
        legacy_modes.append("monitor")
    if args.plot:
        legacy_modes.append("plot")
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
        "See Jarvis -h."
    )


def run_monitor(
    *,
    factory: Any | None = None,
    redis: Any | None = None,
) -> int:
    from jarvishep2.dashboard import attach_reader, format_monitor_view

    reader = attach_reader(factory=factory, redis=redis)
    view = reader.read()
    if not view.has_active_scan():
        print("No active scan found.", file=sys.stderr)
        return EXIT_RUN_FAILED
    sys.stdout.write(format_monitor_view(view))
    return EXIT_OK


def dispatch_refs() -> int:
    """Print the framework and built-in sampler citations without a task card."""
    from jarvishep2.versioning import render_logo_with_version

    sys.stdout.write(f"{render_logo_with_version()}\n\n{render_references()}")
    return EXIT_OK


def dispatch_version() -> int:
    """Print V1-style logo banner with injected runtime version (no scan init)."""
    from jarvishep2.versioning import render_logo_with_version

    print(render_logo_with_version())
    return EXIT_OK


def dispatch_monitor(args: argparse.Namespace) -> int:
    from jarvishep2.process_cleanup import (
        format_scan_table,
        list_active_scans,
        print_scan_table,
        resolve_scan_reference,
        runtime_metadata_for_scan,
    )

    scans = list_active_scans()
    scan_ref = str(getattr(args, "scan_ref", "") or "").strip()
    if not scan_ref:
        print_scan_table(scans)
        return EXIT_OK
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print_scan_table(scans)
        return EXIT_USAGE
    print(f"{scan.reference}: {scan.name}")
    metadata = runtime_metadata_for_scan(scan)
    if metadata is None:
        print("Unable to verify this scan's Redis runtime metadata.", file=sys.stderr)
        return EXIT_RUN_FAILED
    if int(metadata.get("control_pid", -1)) not in {proc.pid for proc in scan.processes}:
        print("Redis runtime metadata does not match the selected control process.", file=sys.stderr)
        return EXIT_RUN_FAILED
    redis_config = dict(metadata["redis"])
    from jarvishep2.factory import TaskFactory
    from jarvishep2.redis_queue import RedisQueue

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


def dispatch_plot(plot_argv: list[str] | None = None) -> int:
    """Forward the complete Jarvis-PLOT CLI without duplicating its surface."""
    try:
        from jarvisplot.client import main as plot_main
    except ImportError as exc:
        print(
            "JarvisPLOT is required for `Jarvis plot`. "
            "Install it with `pip install -U JarvisPLOT` "
            f"(or `pip install -e ../Jarvis-PLOT`). Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    try:
        return int(plot_main(list(plot_argv or [])))
    except Exception as exc:
        print(f"Jarvis plot failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED


def dispatch_gen_plot_yaml(task_yaml: str) -> int:
    """Generate (but never render) JarvisPLOT YAML for one finished scan."""
    if not task_yaml:
        print("Task YAML is required (usage: Jarvis gen-plot-yaml TASK.yaml).", file=sys.stderr)
        return EXIT_USAGE
    from jarvishep2.core import Jarvis2Core
    from jarvishep2.task_config import TaskCardLoadError

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
    except TaskCardLoadError as exc:
        _emit_load_diagnostic(
            task_yaml=task_yaml, code=exc.code, message=str(exc),
            suggestion=exc.suggestion, example=exc.example, as_json=False,
        )
        return EXIT_USAGE
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

        Jarvis portal man
        Jarvis portal man json
        Jarvis portal man slha
        Jarvis portal file.yaml
        Jarvis portal -h
        Jarvis portal -v

    Convenience alias: ``Jarvis portal formats`` → ``man`` (list formats).
    """
    argv = list(portal_argv or [])
    if argv == ["formats"]:
        argv = ["man"]
    try:
        from jarvis_portal.cli import main as portal_main
        from jarvis_portal.v2 import create_default_registry as v2_registry
    except ImportError as exc:
        print(
            "Jarvis-Portal is required for `Jarvis portal`. "
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
                prog="Jarvis portal",
            )
        )
    except Exception as exc:
        print(f"Jarvis portal failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED


def dispatch_operas(operas_argv: list[str] | None = None) -> int:
    """Forward the complete Jarvis-Operas CLI without duplicating its surface."""
    try:
        from jarvis_operas.cli import main as operas_main
    except ImportError as exc:
        print(
            "Jarvis-Operas is required for `Jarvis operas`. "
            f"Detail: {exc}",
            file=sys.stderr,
        )
        return EXIT_USAGE
    try:
        return int(operas_main(list(operas_argv or [])))
    except Exception as exc:
        print(f"Jarvis operas failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED


def dispatch(args: argparse.Namespace) -> int:
    if bool(getattr(args, "version", False)):
        return dispatch_version()
    if bool(getattr(args, "refs", False)):
        return dispatch_refs()

    try:
        intent, args = resolve_intent(args)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    if intent == "monitor":
        return dispatch_monitor(args)
    if intent == "plot":
        return dispatch_plot([])
    if intent == "gen-plot-yaml":
        return dispatch_gen_plot_yaml(str(getattr(args, "task_yaml", "") or ""))
    if intent == "portal":
        # Prefer argv passthrough via main(); this path is a thin fallback.
        return dispatch_portal([])
    if intent == "operas":
        return dispatch_operas([])
    if intent == "project":
        return dispatch_project([])
    if intent == "ps":
        from jarvishep2.process_cleanup import list_running_jarvis_cli

        return int(list_running_jarvis_cli(getattr(args, "scan_ref", None)))
    if intent == "kill":
        from jarvishep2.process_cleanup import kill_running_jarvis_cli

        return int(
            kill_running_jarvis_cli(
                scan_ref=getattr(args, "scan_ref", None),
                yes=bool(getattr(args, "yes", False)),
                force=not bool(getattr(args, "no_force", False)),
            )
        )
    if intent == "man":
        # Handled in main() via argv passthrough (same pattern as portal/operas).
        from jarvishep2.man import dispatch_man

        return dispatch_man([])
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
        raw_timeout = getattr(args, "timeout", None)
        check_timeout: float | None
        try:
            check_timeout = float(raw_timeout) if raw_timeout is not None else None
        except (TypeError, ValueError):
            check_timeout = None
        if check_timeout is not None and check_timeout <= 0:
            print("--timeout must be a positive number of seconds.", file=sys.stderr)
            return EXIT_USAGE
        return dispatch_run(
            str(task or ""),
            check_modules=True,
            skip_draw_flowchart=bool(getattr(args, "skip_draw_flowchart", False)),
            console_level=str(getattr(args, "console_level", "DEBUG") or "DEBUG"),
            silence=bool(getattr(args, "silence", False)),
            strict=bool(getattr(args, "strict", False)),
            skip_library_installation=bool(getattr(args, "skip_library_installation", False)),
            check_timeout=check_timeout,
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
            skip_library_installation=bool(getattr(args, "skip_library_installation", False)),
        )

    print(f"Unknown intent: {intent}", file=sys.stderr)
    return EXIT_USAGE


def main(argv: list[str] | None = None) -> int:
    from jarvishep2.proc_title import control_title, set_process_title

    set_process_title(control_title())
    try:
        raw = list(sys.argv[1:] if argv is None else argv)
        if raw[:1] == ["refs"]:
            print("`Jarvis refs` has moved to `Jarvis --refs`.", file=sys.stderr)
            return EXIT_USAGE
        normalized = normalize_argv(raw)
        # Manual center: own argparse (domains / --json / --code / --type).
        if normalized and normalized[0] == "man":
            from jarvishep2.man import dispatch_man

            return dispatch_man(normalized[1:])
        # Full jportal-compatible surface: do not let argparse eat portal args.
        if normalized and normalized[0] == "portal":
            return dispatch_portal(normalized[1:])
        # Full jplot-compatible surface: do not let argparse consume its args.
        if normalized and normalized[0] == "plot":
            return dispatch_plot(normalized[1:])
        # Full jopera-compatible surface: do not let argparse consume its args.
        if normalized and normalized[0] == "operas":
            return dispatch_operas(normalized[1:])
        # Project pack flags are free-form; pass through like portal.
        if normalized and normalized[0] == "project":
            return dispatch_project(normalized[1:])
        return dispatch(build_parser().parse_args(normalized))
    finally:
        # Keep the onboarding hint after all normal command output. It is a
        # stderr diagnostic so commands such as ``Jarvis validate --json`` keep
        # stdout valid for scripts and other tools.
        try:
            from jarvishep2.redis_notice import emit_first_run_redis_notice

            emit_first_run_redis_notice()
        except Exception:
            # A best-effort environment hint must never change Jarvis's exit
            # status or mask the real command error.
            pass


if __name__ == "__main__":
    raise SystemExit(main())
