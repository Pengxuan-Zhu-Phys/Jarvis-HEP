#!/usr/bin/env python3
"""Minimal Jarvis2 CLI surface for monitor attachment (WP-D5.2)."""

from __future__ import annotations

import argparse
import sys
from typing import Any

from jarvishep2.dashboard import SnapshotReader, attach_reader, format_monitor_view
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import RedisQueue


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Jarvis2", description="Jarvis-HEP V2 CLI")
    parser.add_argument("task_yaml", nargs="?", help="Task YAML (reserved for future run dispatch)")
    parser.add_argument("--monitor", action="store_true", help="Print one monitor snapshot and exit")
    parser.add_argument("--pid", type=int, default=None, help="Attach to a running scan by control PID")
    parser.add_argument("--redis-host", default="127.0.0.1")
    parser.add_argument("--redis-port", type=int, default=6379)
    parser.add_argument("--redis-db", type=int, default=0)
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


def dispatch(args: argparse.Namespace) -> int:
    if not args.monitor:
        print("Only --monitor is implemented in this milestone.", file=sys.stderr)
        return 2

    factory = TaskFactory.get_instance()
    if factory.redis is not None and factory.get_monitor_snapshot():
        return run_monitor(factory=factory)

    redis = RedisQueue(
        {
            "host": args.redis_host,
            "port": args.redis_port,
            "db": args.redis_db,
        }
    )
    try:
        redis.connect()
    except Exception as exc:
        print(f"Unable to connect to Redis for monitor attach: {exc}", file=sys.stderr)
        return 1
    try:
        return run_monitor(redis=redis)
    finally:
        redis.close()


def main(argv: list[str] | None = None) -> int:
    return dispatch(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())