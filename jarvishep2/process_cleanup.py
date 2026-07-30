#!/usr/bin/env python3
"""List / terminate running Jarvis OS processes (``Jarvis2 ps`` / ``Jarvis2 kill``).

Process titles are set via ``setproctitle`` (see ``proc_title.py``)::

    Jarvis2:<scan>
    Jarvis2-Worker-N:<scan>
    Jarvis2-Archiver:<scan>
    Jarvis2-FileOperation:<scan>
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

from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

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
        or head == "Jarvis2-FileOperation"
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
    """Group titled control/worker/file-operation/archiver/Redis processes into scan rows."""
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


def ensure_scan_name_available(scan_name: str) -> None:
    """Reject a new run when its exact scan name already has a live controller."""
    name = str(scan_name or "").strip()
    if not name:
        raise ValueError("scan name is required")
    for scan in list_running_scans():
        if scan.name != name:
            continue
        if any(proc.command.split(None, 1)[0].startswith("Jarvis2:") for proc in scan.processes):
            raise RuntimeError(
                f"a Jarvis scan named {name!r} is already running; "
                f"inspect it with `Jarvis2 monitor {scan.reference}` or "
                f"`Jarvis2 ps {scan.reference}` before starting another run"
            )


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


def print_scan_table(scans: list[JarvisScan]) -> None:
    """Render the running-task chooser with explicit, safe next actions."""
    if not scans:
        Console().print("[dim]No running Jarvis scan tasks.[/]")
        return
    table = Table(box=box.SIMPLE_HEAVY, show_header=True, header_style="bold dim")
    table.add_column("REF", style="bold green", no_wrap=True)
    table.add_column("SCAN", style="bold #c8c8ff")
    table.add_column("PROCESSES", justify="right")
    table.add_column("PIDS", style="dim")
    for scan in scans:
        table.add_row(
            scan.reference,
            scan.name,
            str(len(scan.processes)),
            ", ".join(str(proc.pid) for proc in scan.processes),
        )
    guidance = Text()
    guidance.append("Choose a task by its REF (for example R1):\n", style="bold")
    guidance.append("  Jarvis2 monitor R1", style="bold cyan")
    guidance.append("  view this task's live scan status\n", style="dim")
    guidance.append("  Jarvis2 ps R1", style="bold #c8c8ff")
    guidance.append(
        "       list this task's control, workers, file operations, archiver, and Redis\n",
        style="dim",
    )
    guidance.append("  Jarvis2 kill R1", style="bold red")
    guidance.append("     terminate ALL processes belonging to this task (asks for confirmation)", style="bold red")
    Console().print(
        Panel(
            Group(table, Text(""), guidance),
            title=f"[bold]Running Jarvis scan tasks[/] [dim]({len(scans)})[/]",
            box=box.ROUNDED,
            border_style="dim",
            padding=(0, 1),
        )
    )


def _process_role_and_detail(proc: JarvisProcess, scan_name: str) -> tuple[str, str]:
    """Turn a long OS title into the concise role shown in the runtime card."""
    title = proc.command.split(None, 1)[0]
    if title.startswith("Jarvis2:"):
        return "Control", "scan coordinator"
    if title.startswith("Jarvis2-Worker-"):
        worker = title.removeprefix("Jarvis2-Worker-").split(":", 1)[0]
        return f"Worker {worker}", "calculator / sampling worker"
    if title.startswith("Jarvis2-FileOperation:"):
        return "FileOperation", "Worker-owned SAMPLE save / copy / delete"
    if title.startswith("Jarvis2-Archiver:"):
        return "Archiver", "writes DATABASE and exports"
    if title.startswith("Jarvis-Redis:"):
        endpoint = proc.command.split(None, 1)[1] if " " in proc.command else ""
        advertised = title.rsplit("@", 1)[-1] if "@" in title else ""
        return "Redis", endpoint or advertised or "runtime broker"
    return "Jarvis process", title.replace(scan_name, "<scan>")


def print_scan_processes(
    scan: JarvisScan,
    metadata: dict[str, object] | None = None,
    *,
    kill_warning: bool = False,
) -> None:
    """Render a selected scan as a concise, readable process card."""
    details = Table.grid(padding=(0, 1))
    details.add_column(style="bold dim", no_wrap=True)
    details.add_column()
    details.add_row("REF", f"[bold green]{scan.reference}[/]")
    details.add_row("SCAN", f"[bold #c8c8ff]{scan.name}[/]")
    if metadata is not None:
        redis = metadata.get("redis")
        if isinstance(redis, dict):
            details.add_row("REDIS", f"{redis.get('host')}:{redis.get('port')}  ·  db {redis.get('db')}")
        output = str(metadata.get("task_result_dir") or "")
        if output:
            details.add_row("OUTPUT", f"[dim]{output}[/]")
    processes = Table(box=box.SIMPLE_HEAVY, header_style="bold dim")
    processes.add_column("PID", justify="right", style="bold cyan", no_wrap=True)
    processes.add_column("ROLE", style="bold #c8c8ff", no_wrap=True)
    processes.add_column("DETAIL", style="dim")
    for proc in scan.processes:
        role, detail = _process_role_and_detail(proc, scan.name)
        processes.add_row(str(proc.pid), role, detail)
    body: list[object] = [details, Text(""), processes]
    if kill_warning:
        body.extend(
            [
                Text(""),
                Text(
                    f"DANGER: Jarvis2 kill {scan.reference} will terminate all {len(scan.processes)} processes above.",
                    style="bold red",
                ),
            ]
        )
    Console().print(
        Panel(
            Group(*body),
            title=f"[bold]Runtime task[/] [dim]·[/] [bold green]{scan.reference}[/]",
            box=box.ROUNDED,
            border_style="red" if kill_warning else "dim",
            padding=(0, 1),
        )
    )


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
        print_scan_table(scans)
        return 0
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print_scan_table(scans)
        return 2
    metadata = _verified_runtime_metadata(scan)
    if _scan_has_advertised_redis(scan) and metadata is None:
        return 1
    print_scan_processes(scan, metadata)
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
        print_scan_table(scans)
        return 0
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print_scan_table(scans)
        return 2
    targets = list(scan.processes)
    metadata = _verified_runtime_metadata(scan)
    if _scan_has_advertised_redis(scan) and metadata is None:
        return 1
    print_scan_processes(scan, metadata, kill_warning=True)
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
    "ensure_scan_name_available",
    "format_process_table",
    "format_scan_table",
    "print_scan_table",
    "print_scan_processes",
    "kill_jarvis_processes",
    "kill_running_jarvis_cli",
    "list_jarvis_processes",
    "list_running_scans",
    "list_running_jarvis_cli",
    "resolve_scan_reference",
    "runtime_metadata_for_scan",
]
