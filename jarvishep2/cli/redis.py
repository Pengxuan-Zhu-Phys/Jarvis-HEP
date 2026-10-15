"""Redis environment inspection and copyable setup instructions."""

from __future__ import annotations

import json
import sys

from jarvishep2.redis_notice import build_redis_environment_report
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED


def dispatch_redis(*, topic: str = "info", as_json: bool = False) -> int:
    """Inspect the local executable; install/start topics only print guidance."""

    try:
        report = build_redis_environment_report(probe_version=topic in {"info", "check"})
    except Exception as exc:
        error = {"ok": False, "kind": "redis-environment", "error": {
            "code": "environment_error", "message": str(exc),
        }}
        if as_json:
            print(json.dumps(error, indent=2))
        else:
            print(f"Cannot inspect the Redis environment: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    report["topic"] = topic
    if topic == "check":
        report["ok"] = bool(report["server"]["usable"]) and report["runtime_supported"] is not False
    if as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_report(report, topic=topic)
    return EXIT_RUN_FAILED if not report["ok"] else EXIT_OK


def _print_report(report: dict, *, topic: str) -> None:
    system = report["system"]
    label = system["distro_name"] or system["label"]
    if system["wsl"]:
        label += " (WSL)"
    print(f"Redis environment — {label}")
    server = report["server"]
    if server["found"]:
        print(f"Server executable: {server['path']}")
    else:
        print("Redis-compatible server executable was not found on PATH.")
    if topic in {"info", "check"}:
        if server["version"]:
            print(f"Server version: {server['version']}")
        if server["version_error"]:
            print(f"Version check: {server['version_error']}")
        if server["supports_blmove"] is False:
            print("Version compatibility: this Redis version is too old; Jarvis requires Redis >= 6.2.")
        if topic == "check":
            print("Environment check: " + ("passed" if report["ok"] else "failed"))

    if topic in {"info", "check", "install"}:
        print("\nInstallation instructions (run these yourself):")
        for command in report["install"]["commands"]:
            print(f"  {command}")
        for note in report["install"]["notes"]:
            print(note)
        for url in report["install"]["documentation"]:
            print(f"Documentation: {url}")
    if topic in {"info", "check", "start"}:
        print("\nStartup:")
        for note in report["start"]["notes"]:
            print(note)
        for command in report["start"]["commands"]:
            print(f"  {command}")
    for note in report["notes"]:
        print(note)
    print("\nMore: Jarvis redis check | Jarvis redis install | Jarvis redis start | Jarvis redis --json")
