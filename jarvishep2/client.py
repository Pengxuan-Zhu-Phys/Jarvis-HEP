#!/usr/bin/env python3
"""Jarvis2 CLI surface: run dispatch, monitor attach, plot bridge, check-modules smoke."""

from __future__ import annotations

import argparse
import sys

from jarvishep2.core import Jarvis2Core
from jarvishep2.dashboard import attach_reader, format_monitor_view
from jarvishep2.factory import TaskFactory
from jarvishep2.plot_bridge import PlotBridgeError, run_plot
from jarvishep2.redis_queue import INTERNAL_REDIS_CONFIG, RedisQueue


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Jarvis2", description="Jarvis-HEP V2 CLI")
    parser.add_argument(
        "task_yaml",
        nargs="?",
        help="Task YAML for distributed runs, or plot YAML when --plot is set",
    )
    parser.add_argument("--monitor", action="store_true", help="Print one monitor snapshot and exit")
    parser.add_argument(
        "--plot",
        action="store_true",
        help="Run Jarvis-PLOT on the given YAML (delegates to the JarvisPLOT package)",
    )
    parser.add_argument("--check-modules", action="store_true", help="Run fixed-point calculator smoke")
    parser.add_argument("--resume", action="store_true", help="Resume from an existing checkpoint")
    parser.add_argument("--pid", type=int, default=None, help="Attach to a running scan by control PID")
    return parser


def run_monitor(
    *,
    factory: TaskFactory | None = None,
    redis: RedisQueue | None = None,
) -> int:
    reader = attach_reader(factory=factory, redis=redis)
    view = reader.read()
    if not view.has_active_scan():
        print("No active scan found.", file=sys.stderr)
        return 1
    sys.stdout.write(format_monitor_view(view))
    return 0


def dispatch_monitor(args: argparse.Namespace) -> int:
    # The V2 broker is internal and always uses the local service.
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
        return 1
    try:
        return run_monitor(redis=redis)
    finally:
        redis.close()


def dispatch_plot(args: argparse.Namespace) -> int:
    if not args.task_yaml:
        print("Plot YAML is required with --plot (usage: Jarvis2 --plot scene.yaml).", file=sys.stderr)
        return 2
    try:
        return int(run_plot(args.task_yaml))
    except ImportError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except PlotBridgeError as exc:
        print(str(exc), file=sys.stderr)
        return 1


def dispatch_run(args: argparse.Namespace) -> int:
    if not args.task_yaml:
        print("Task YAML is required for distributed runs.", file=sys.stderr)
        return 2

    try:
        core = Jarvis2Core()
        core.load_task_yaml(args.task_yaml)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 1
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return 1

    try:
        if args.check_modules:
            count = core.check_modules()
        else:
            count = core.run(resume=args.resume)
    except NotImplementedError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    except (RuntimeError, TimeoutError, OSError) as exc:
        print(f"Run failed: {exc}", file=sys.stderr)
        return 1

    return 0 if count > 0 else 1


def dispatch(args: argparse.Namespace) -> int:
    if args.monitor:
        return dispatch_monitor(args)
    if args.plot:
        return dispatch_plot(args)
    if args.check_modules or args.resume or args.task_yaml:
        return dispatch_run(args)
    print("Provide a task YAML, or use --monitor / --plot.", file=sys.stderr)
    return 2


def main(argv: list[str] | None = None) -> int:
    return dispatch(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
