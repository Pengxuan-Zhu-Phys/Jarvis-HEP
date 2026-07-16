#!/usr/bin/env python3
"""Detect / terminate leftover Jarvis OS processes (orphan scan cleanup).

Process titles are set via ``setproctitle`` (see ``proc_title.py``)::

    Jarvis2:<scan>
    Jarvis2-Worker-N:<scan>
    Jarvis2-Archiver:<scan>
    Jarvis2-FileOperation
    Jarvis-Redis:<scan>

CLI::

    Jarvis2 ps              # list running Jarvis processes
    Jarvis2 kill            # kill them after interactive confirmation
    Jarvis2 kill --yes      # non-interactive kill
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Iterable


# Match OS titles we set. "Jarvis2*" covers control/worker/archiver/file-ops;
# "Jarvis-Redis" is the managed redis-server title.
_TITLE_PREFIXES: tuple[str, ...] = (
    "Jarvis2",
    "Jarvis-Redis",
)


@dataclass(frozen=True)
class JarvisProcess:
    pid: int
    command: str


def _command_is_jarvis(command: str) -> bool:
    text = str(command or "").strip()
    if not text:
        return False
    # Full command string may be the title alone, or "title args…".
    first = text.split(None, 1)[0]
    for prefix in _TITLE_PREFIXES:
        if text.startswith(prefix) or first.startswith(prefix):
            return True
    return False


def list_jarvis_processes(*, exclude_pids: Iterable[int] | None = None) -> list[JarvisProcess]:
    """Return running processes whose title/command looks like a Jarvis component."""
    skip = {int(os.getpid())}
    if exclude_pids:
        skip.update(int(p) for p in exclude_pids)

    try:
        # Portable: PID + full command/args (title after setproctitle).
        completed = subprocess.run(
            ["ps", "-ax", "-o", "pid=,command="],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if int(completed.returncode) != 0:
        return []

    found: list[JarvisProcess] = []
    for line in (completed.stdout or "").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) < 2:
            continue
        try:
            pid = int(parts[0])
        except ValueError:
            continue
        if pid in skip or pid <= 1:
            continue
        command = parts[1].strip()
        if not _command_is_jarvis(command):
            continue
        found.append(JarvisProcess(pid=pid, command=command))
    found.sort(key=lambda item: item.pid)
    return found


def format_process_table(procs: list[JarvisProcess]) -> str:
    if not procs:
        return "No running Jarvis processes.\n"
    width = max(len(str(p.pid)) for p in procs)
    lines = [
        f"Running Jarvis processes ({len(procs)}):",
        f"{'PID':>{width}}  COMMAND",
        f"{'-' * width}  -------",
    ]
    for proc in procs:
        lines.append(f"{proc.pid:>{width}}  {proc.command}")
    return "\n".join(lines) + "\n"


def confirm_kill(count: int, *, yes: bool = False) -> bool:
    """Ask the user before killing. ``--yes`` skips the prompt."""
    if count <= 0:
        return False
    if yes:
        return True
    if not sys.stdin.isatty():
        print(
            "Refusing to kill without confirmation in non-interactive mode; "
            "re-run with: Jarvis2 kill --yes",
            file=sys.stderr,
        )
        return False
    try:
        answer = input(f"Kill {count} Jarvis process(es)? [y/N] ").strip().lower()
    except EOFError:
        return False
    return answer in {"y", "yes"}


def _kill_order_key(proc: JarvisProcess) -> tuple[int, int]:
    """Prefer Workers/FileOp → Archiver → control → Redis (drop deps first)."""
    cmd = proc.command
    if "Worker" in cmd:
        rank = 0
    elif "FileOperation" in cmd:
        rank = 1
    elif "Archiver" in cmd:
        rank = 2
    elif cmd.startswith("Jarvis-Redis"):
        rank = 4
    elif cmd.startswith("Jarvis2:") or cmd == "Jarvis2" or cmd.startswith("Jarvis2 "):
        rank = 3
    else:
        rank = 5
    return (rank, proc.pid)


def _safe_kill(pid: int, sig: int) -> str:
    """Send ``sig`` to ``pid``. Return ok|missing|failed."""
    try:
        os.kill(pid, sig)
        return "ok"
    except ProcessLookupError:
        return "missing"
    except PermissionError:
        return "failed"
    except OSError:
        return "failed"


def kill_jarvis_processes(
    procs: list[JarvisProcess] | None = None,
    *,
    sigterm_grace_sec: float = 2.0,
    force: bool = True,
) -> dict[str, list[int]]:
    """SIGTERM then optional SIGKILL for listed (or current) Jarvis processes.

    Returns ``{"signaled": [...], "killed": [...], "missing": [...], "failed": [...]}``.
    """
    targets = list(procs) if procs is not None else list_jarvis_processes()
    targets = sorted(targets, key=_kill_order_key)
    signaled: list[int] = []
    missing: list[int] = []
    failed: list[int] = []

    for proc in targets:
        # Alive check
        status = _safe_kill(proc.pid, 0)
        if status == "missing":
            missing.append(proc.pid)
            continue
        if status == "failed":
            failed.append(proc.pid)
            continue
        status = _safe_kill(proc.pid, signal.SIGTERM)
        if status == "ok":
            signaled.append(proc.pid)
        elif status == "missing":
            missing.append(proc.pid)
        else:
            failed.append(proc.pid)

    if signaled and sigterm_grace_sec > 0:
        time.sleep(max(0.0, float(sigterm_grace_sec)))

    killed: list[int] = []
    if force:
        for pid in list(signaled):
            status = _safe_kill(pid, 0)
            if status == "missing":
                continue
            if status == "failed":
                failed.append(pid)
                continue
            status = _safe_kill(pid, signal.SIGKILL)
            if status == "ok":
                killed.append(pid)
            elif status == "failed":
                failed.append(pid)

    return {
        "signaled": signaled,
        "killed": killed,
        "missing": missing,
        "failed": sorted(set(failed)),
    }


def list_running_jarvis_cli() -> int:
    """``Jarvis2 ps`` — show running Jarvis processes (scan-time friendly)."""
    procs = list_jarvis_processes()
    print(format_process_table(procs), end="")
    if procs:
        print("To terminate: Jarvis2 kill")
    return 0


def kill_running_jarvis_cli(
    *,
    yes: bool = False,
    force: bool = True,
) -> int:
    """``Jarvis2 kill`` — confirm (unless ``--yes``), then SIGTERM/SIGKILL."""
    procs = list_jarvis_processes()
    print(format_process_table(procs), end="")
    if not procs:
        return 0
    if not confirm_kill(len(procs), yes=yes):
        print("Aborted.")
        return 0
    result = kill_jarvis_processes(procs, force=force)
    print(
        "kill: "
        f"SIGTERM={len(result['signaled'])} "
        f"SIGKILL={len(result['killed'])} "
        f"gone={len(result['missing'])} "
        f"failed={len(result['failed'])}"
    )
    if result["failed"]:
        print(f"  failed pids: {result['failed']}")
        return 1
    left = list_jarvis_processes()
    if left:
        print("Still running after kill:")
        print(format_process_table(left), end="")
        return 1
    print("All matched Jarvis processes terminated.")
    return 0


def cleanup_jarvis_processes(*, kill: bool = False, force: bool = True) -> int:
    """Backward-compatible helper (prefer ``ps`` / ``kill`` CLI)."""
    if kill:
        return kill_running_jarvis_cli(yes=True, force=force)
    return list_running_jarvis_cli()


__all__ = [
    "JarvisProcess",
    "cleanup_jarvis_processes",
    "confirm_kill",
    "format_process_table",
    "kill_jarvis_processes",
    "kill_running_jarvis_cli",
    "list_jarvis_processes",
    "list_running_jarvis_cli",
]
