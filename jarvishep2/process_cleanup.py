#!/usr/bin/env python3
"""Detect / terminate leftover Jarvis OS processes (orphan scan cleanup).

Process titles are set via ``setproctitle`` (see ``proc_title.py``)::

    Jarvis2:<scan>
    Jarvis2-Worker-N:<scan>
    Jarvis2-Archiver:<scan>
    Jarvis2-FileOperation
    Jarvis-Redis:<scan>

CLI: ``Jarvis2 cleanup`` (list) / ``Jarvis2 cleanup --kill``.
"""

from __future__ import annotations

import os
import signal
import subprocess
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
        return "No leftover Jarvis processes found.\n"
    width = max(len(str(p.pid)) for p in procs)
    lines = [f"{'PID':>{width}}  COMMAND", f"{'-' * width}  -------"]
    for proc in procs:
        lines.append(f"{proc.pid:>{width}}  {proc.command}")
    lines.append(f"\n{len(procs)} process(es)")
    return "\n".join(lines) + "\n"


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
    signaled: list[int] = []
    missing: list[int] = []
    failed: list[int] = []

    for proc in targets:
        try:
            os.kill(proc.pid, 0)
        except ProcessLookupError:
            missing.append(proc.pid)
            continue
        except PermissionError:
            failed.append(proc.pid)
            continue
        try:
            os.kill(proc.pid, signal.SIGTERM)
            signaled.append(proc.pid)
        except ProcessLookupError:
            missing.append(proc.pid)
        except PermissionError:
            failed.append(proc.pid)

    if signaled and sigterm_grace_sec > 0:
        time.sleep(max(0.0, float(sigterm_grace_sec)))

    killed: list[int] = []
    if force:
        for pid in list(signaled):
            try:
                os.kill(pid, 0)
            except ProcessLookupError:
                continue
            except PermissionError:
                failed.append(pid)
                continue
            try:
                os.kill(pid, signal.SIGKILL)
                killed.append(pid)
            except ProcessLookupError:
                pass
            except PermissionError:
                failed.append(pid)

    return {
        "signaled": signaled,
        "killed": killed,
        "missing": missing,
        "failed": sorted(set(failed)),
    }


def cleanup_jarvis_processes(*, kill: bool = False, force: bool = True) -> int:
    """CLI helper: list (default) or kill leftovers. Returns process exit code."""
    procs = list_jarvis_processes()
    print(format_process_table(procs), end="")
    if not kill:
        if procs:
            print("Re-run with: Jarvis2 cleanup --kill")
        return 0
    if not procs:
        return 0
    result = kill_jarvis_processes(procs, force=force)
    print(
        "cleanup: "
        f"SIGTERM={len(result['signaled'])} "
        f"SIGKILL={len(result['killed'])} "
        f"gone={len(result['missing'])} "
        f"failed={len(result['failed'])}"
    )
    if result["failed"]:
        print(f"  failed pids: {result['failed']}")
        return 1
    # Verify
    left = list_jarvis_processes()
    if left:
        print("Still running after kill:")
        print(format_process_table(left), end="")
        return 1
    print("All matched Jarvis processes terminated.")
    return 0


__all__ = [
    "JarvisProcess",
    "cleanup_jarvis_processes",
    "format_process_table",
    "kill_jarvis_processes",
    "list_jarvis_processes",
]
