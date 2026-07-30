#!/usr/bin/env python3
"""List / terminate running Jarvis OS processes (``Jarvis2 ps`` / ``Jarvis2 kill``).

Process titles are set via ``setproctitle`` (see ``proc_title.py``)::

    Jarvis2:<scan>
    Jarvis2-Worker-N:<scan>
    Jarvis2-Archiver:<scan>
    Jarvis2-FileOperation
    Jarvis-Redis:<scan>
"""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from typing import Iterable

from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_metadata import read_scan_metadata


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


@dataclass(frozen=True)
class JarvisScan:
    """One running scan, resolved from the titles of its component processes."""

    reference: str
    name: str
    processes: tuple[JarvisProcess, ...]


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


def _scan_name_from_command(command: str) -> str | None:
    """Extract ``<scan>`` from a Jarvis process title ending in ``:<scan>``."""
    title = str(command or "").strip().split(None, 1)[0]
    head, marker, name = title.partition(":")
    if marker and (
        head == "Jarvis2"
        or head.startswith("Jarvis2-Worker-")
        or head == "Jarvis2-Archiver"
        or head == "Jarvis-Redis"
    ):
        return name.split("@", 1)[0].strip() or None
    return None


def runtime_metadata_for_scan(scan: JarvisScan) -> dict[str, object] | None:
    """Read a scan identity record through its advertised Redis endpoint."""
    redis_process = next(
        (proc for proc in scan.processes if proc.command.startswith("Jarvis-Redis:")),
        None,
    )
    if redis_process is None:
        return None
    title = redis_process.command.split(None, 1)[0]
    marker = title.rsplit("@", 1)
    if len(marker) != 2 or "/" not in marker[1]:
        return None
    try:
        port_text, db_text = marker[1].split("/", 1)
        redis_config = {"host": "127.0.0.1", "port": int(port_text), "db": int(db_text)}
    except ValueError:
        return None
    redis = RedisQueue(redis_config)
    try:
        redis.connect()
        path = redis.get_runtime_metadata_path()
        if not path:
            return None
        return read_scan_metadata(path, redis=redis_config, expected_scan=scan.name)
    except Exception:
        return None
    finally:
        try:
            redis.close()
        except Exception:
            pass


def _scan_has_advertised_redis(scan: JarvisScan) -> bool:
    return any(
        "@" in proc.command.split(None, 1)[0]
        for proc in scan.processes
        if proc.command.startswith("Jarvis-Redis:")
    )


def _verified_runtime_metadata(scan: JarvisScan) -> dict[str, object] | None:
    """Require metadata for modern runs, while retaining old orphan recovery."""
    metadata = runtime_metadata_for_scan(scan)
    if metadata is None:
        if _scan_has_advertised_redis(scan):
            print("Refusing unverified runtime: Redis metadata is unavailable.", file=sys.stderr)
        return None
    control_pid = int(metadata.get("control_pid", -1))
    if control_pid not in {proc.pid for proc in scan.processes}:
        print("Refusing unverified runtime: control PID does not match metadata.", file=sys.stderr)
        return None
    return metadata


def list_running_scans(
    procs: Iterable[JarvisProcess] | None = None,
) -> list[JarvisScan]:
    """Group titled control/worker/archiver/Redis processes into scan rows."""
    grouped: dict[str, list[JarvisProcess]] = {}
    for proc in procs if procs is not None else list_jarvis_processes():
        name = _scan_name_from_command(proc.command)
        if name:
            grouped.setdefault(name, []).append(proc)
    scans: list[JarvisScan] = []
    for index, name in enumerate(sorted(grouped, key=str.casefold), start=1):
        scans.append(
            JarvisScan(
                reference=f"R{index}",
                name=name,
                processes=tuple(sorted(grouped[name], key=lambda proc: proc.pid)),
            )
        )
    return scans


def resolve_scan_reference(selector: str, scans: Iterable[JarvisScan]) -> JarvisScan:
    """Resolve an ``R1`` table reference or an exact scan name."""
    text = str(selector or "").strip()
    choices = list(scans)
    for scan in choices:
        if text.upper() == scan.reference.upper() or text == scan.name:
            return scan
    available = ", ".join(f"{scan.reference} ({scan.name})" for scan in choices) or "none"
    raise ValueError(f"Unknown running scan {text!r}; available: {available}")


def format_scan_table(scans: list[JarvisScan]) -> str:
    """Render the compact task chooser shared by monitor / ps / kill."""
    if not scans:
        return "No running Jarvis scan tasks.\n"
    ref_w = max(len(scan.reference) for scan in scans)
    name_w = max(len(scan.name) for scan in scans)
    lines = [
        f"Running Jarvis scan tasks ({len(scans)}):",
        f"{'REF':<{ref_w}}  {'SCAN':<{name_w}}  PROCESSES  PIDS",
        f"{'-' * ref_w}  {'-' * name_w}  ---------  ----",
    ]
    for scan in scans:
        pids = ",".join(str(proc.pid) for proc in scan.processes)
        lines.append(
            f"{scan.reference:<{ref_w}}  {scan.name:<{name_w}}  "
            f"{len(scan.processes):>9}  {pids}"
        )
    lines.append("Select one: Jarvis2 monitor R1 | Jarvis2 ps R1 | Jarvis2 kill R1")
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


def list_running_jarvis_cli(scan_ref: str | None = None) -> int:
    """``Jarvis2 ps [R#]`` — choose a scan or list its component processes."""
    procs = list_jarvis_processes()
    scans = list_running_scans(procs)
    if not scan_ref:
        print(format_scan_table(scans), end="")
        return 0
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print(format_scan_table(scans), end="")
        return 2
    print(f"{scan.reference}: {scan.name}")
    metadata = _verified_runtime_metadata(scan)
    if _scan_has_advertised_redis(scan) and metadata is None:
        return 1
    if metadata is not None:
        print(f"Runtime metadata: {metadata['task_result_dir']}")
    print(format_process_table(list(scan.processes)), end="")
    return 0


def kill_running_jarvis_cli(
    *,
    scan_ref: str | None = None,
    yes: bool = False,
    force: bool = True,
) -> int:
    """``Jarvis2 kill [R#]`` — choose a scan, then confirm termination."""
    procs = list_jarvis_processes()
    scans = list_running_scans(procs)
    if not scan_ref:
        print(format_scan_table(scans), end="")
        return 0
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print(format_scan_table(scans), end="")
        return 2
    targets = list(scan.processes)
    print(f"{scan.reference}: {scan.name}")
    metadata = _verified_runtime_metadata(scan)
    if _scan_has_advertised_redis(scan) and metadata is None:
        return 1
    if metadata is not None:
        print(f"Runtime metadata: {metadata['task_result_dir']}")
    print(format_process_table(targets), end="")
    if not targets:
        return 0
    if not confirm_kill(len(targets), yes=yes):
        print("Aborted.")
        return 0
    result = kill_jarvis_processes(targets, force=force)
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
    left = [proc for proc in list_jarvis_processes() if _scan_name_from_command(proc.command) == scan.name]
    if left:
        print("Still running after kill:")
        print(format_process_table(left), end="")
        return 1
    print("All matched Jarvis processes terminated.")
    return 0


__all__ = [
    "JarvisProcess",
    "JarvisScan",
    "confirm_kill",
    "format_process_table",
    "format_scan_table",
    "kill_jarvis_processes",
    "kill_running_jarvis_cli",
    "list_jarvis_processes",
    "list_running_scans",
    "list_running_jarvis_cli",
    "resolve_scan_reference",
    "runtime_metadata_for_scan",
]
