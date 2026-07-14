#!/usr/bin/env python3
"""Jarvis2 CLI surface: one-intent subcommands + legacy flag aliases (D11.1 / D11.2)."""

from __future__ import annotations

import argparse
import sys
import warnings
from typing import Any

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


_SUBCOMMANDS = frozenset({"run", "check", "monitor", "plot", "portal", "operas"})


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
    # Bare path (+ optional --resume) → run
    return ["run", *args]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="Jarvis2",
        description="Jarvis-HEP V2 CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Intents (preferred):\n"
            "  Jarvis2 run TASK.yaml [--resume]\n"
            "  Jarvis2 check TASK.yaml\n"
            "  Jarvis2 monitor\n"
            "  Jarvis2 plot PLOT.yaml\n"
            "  Jarvis2 portal …            # same CLI as jportal (V2 registry)\n"
            "  Jarvis2 operas list|info NAME\n"
            "\n"
            "Legacy aliases (still accepted):\n"
            "  Jarvis2 TASK.yaml [--resume]\n"
            "  Jarvis2 TASK.yaml --check-modules\n"
            "  Jarvis2 --monitor\n"
            "  Jarvis2 PLOT.yaml --plot   (deprecated)\n"
        ),
    )
    sub = parser.add_subparsers(dest="command", required=False)

    run_p = sub.add_parser("run", help="Run a distributed scan task YAML")
    run_p.add_argument("task_yaml", help="Path to scan task YAML")
    run_p.add_argument("--resume", action="store_true", help="Resume from checkpoint")
    run_p.add_argument(
        "--skip-draw-flowchart",
        action="store_true",
        help="Skip workflow flowchart.json / flowchart.png export",
    )

    check_p = sub.add_parser("check", help="Run fixed-point calculator smoke (check-modules)")
    check_p.add_argument("task_yaml", help="Path to check-modules task YAML")

    sub.add_parser("monitor", help="Print one monitor snapshot and exit")

    plot_p = sub.add_parser("plot", help="Render a JarvisPLOT scene YAML")
    plot_p.add_argument("plot_yaml", help="Path to plot scene YAML")

    # ``portal`` is handled by argv passthrough in main() so that
    # ``Jarvis2 portal man|file|-h|-v`` matches the standalone jportal CLI.
    sub.add_parser(
        "portal",
        help="Jarvis-Portal CLI (same as jportal; uses V2 format registry)",
        add_help=False,
    )

    operas_p = sub.add_parser("operas", help="Jarvis-Operas discovery helpers")
    operas_sub = operas_p.add_subparsers(dest="operas_command", required=False)
    operas_sub.add_parser("list", help="List registered Operas operators")
    operas_info = operas_sub.add_parser("info", help="Show one Operas operator")
    operas_info.add_argument("name", help="Operator name or dotted path")

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
        "--check-modules",
        action="store_true",
        help="(legacy) Check-modules path; prefer `Jarvis2 check`",
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


def _mode_conflict_message(modes: list[str]) -> str:
    return (
        "conflicting CLI intents: "
        + ", ".join(modes)
        + "; use a single subcommand (run|check|monitor|plot|portal|operas)"
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
        "(run|check|monitor|plot|portal|operas). See Jarvis2 -h."
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
    print(
        f"RunOutcome status={outcome.status} "
        f"submitted={outcome.submitted} completed={outcome.completed} "
        f"failed={outcome.failed} archived={outcome.archived} "
        f"exit={outcome.exit_code}",
        file=sys.stderr,
    )
    if outcome.error:
        print(f"  error_type={outcome.error_type} error={outcome.error}", file=sys.stderr)


def dispatch_run(
    task_yaml: str,
    *,
    resume: bool = False,
    check_modules: bool = False,
    skip_draw_flowchart: bool = False,
) -> int:
    if not task_yaml:
        print("Task YAML is required.", file=sys.stderr)
        return EXIT_USAGE

    try:
        core = Jarvis2Core()
        core.load_task_yaml(task_yaml)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED

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
    if intent == "portal":
        # Prefer argv passthrough via main(); this path is a thin fallback.
        return dispatch_portal([])
    if intent == "operas":
        return dispatch_operas(args)
    if intent == "check":
        task = getattr(args, "task_yaml", None)
        return dispatch_run(
            str(task or ""),
            check_modules=True,
            skip_draw_flowchart=bool(getattr(args, "skip_draw_flowchart", False)),
        )
    if intent == "run":
        task = getattr(args, "task_yaml", None)
        return dispatch_run(
            str(task or ""),
            resume=bool(args.resume),
            check_modules=False,
            skip_draw_flowchart=bool(getattr(args, "skip_draw_flowchart", False)),
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
    return dispatch(build_parser().parse_args(normalized))


if __name__ == "__main__":
    raise SystemExit(main())
