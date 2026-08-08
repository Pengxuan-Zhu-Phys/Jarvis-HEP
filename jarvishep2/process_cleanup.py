"""List / terminate running Jarvis OS processes (``Jarvis ps`` / ``Jarvis kill``).

Process titles are set via ``setproctitle`` (see ``proc_title.py``)::

    Jarvis:<scan>
    Jarvis-Worker-N:<scan>
    Jarvis-Archiver:<scan>
    Jarvis-FileOperation:<scan>
    Jarvis-Redis:<scan>
"""

from __future__ import annotations

import json
import os
import re
import signal
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.redis_queue import RedisQueue
from jarvishep2.runtime_metadata import read_scan_metadata

# Sticky R1/R2 registry path. Override with JARVIS_SCAN_REF_REGISTRY for tests.
_ENV_REF_REGISTRY = "JARVIS_SCAN_REF_REGISTRY"


# Jarvis-Lit is a separate application.  Its executable may appear either as
# a setproctitle-style name or as an absolute path to its Python/app bundle.
_JARVIS_LIT_MARKERS: tuple[str, ...] = (
    "Jarvis-Lit",
    "JarvisLit",
)
# Retired V1 used Jarvis-HEP@scan; never treat it (or sibling CLIs) as V2 runtime.
_FOREIGN_JARVIS_PREFIXES: tuple[str, ...] = (
    "Jarvis-HEP",
    "Jarvis-Lit",
    "JarvisLit",
    "JarvisPLOT",
    "Jarvis-Agent",
    "jarvis-agent",
    "jlit",
    "jplot",
    "jportal",
    "jopera",
)
ZP_REFERENCE = "ZP"
ZP_LABEL = "Zambie Process"

# Control process title is exactly ``Jarvis`` or ``Jarvis:<scan>``.
_CONTROL_HEADS: frozenset[str] = frozenset({"Jarvis"})
_ROLE_HEAD_PREFIXES: tuple[str, ...] = (
    "Jarvis-Worker-",
    "Jarvis-Archiver",
    "Jarvis-FileOperation",
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


def _title_head(command: str) -> str:
    return str(command or "").strip().split(None, 1)[0]


def _is_control_title(title: str) -> bool:
    head = title.split(":", 1)[0]
    return head in _CONTROL_HEADS


def _is_control_command(command: str) -> bool:
    return _is_control_title(_title_head(command))


def _command_is_jarvis(command: str) -> bool:
    text = str(command or "").strip()
    if not text:
        return False
    if _command_is_jarvis_lit(text):
        return False
    first = text.split(None, 1)[0]
    # Managed Redis is part of the runtime even though it shares the Jarvis- prefix.
    if first == "Jarvis-Redis" or first.startswith("Jarvis-Redis:") or first.startswith("Jarvis-Redis@"):
        return True
    if any(
        first == p or first.startswith(p + ":") or first.startswith(p + "@") or first.startswith(p + "-")
        for p in _FOREIGN_JARVIS_PREFIXES
    ):
        return False
    if _is_control_title(first):
        return True
    for prefix in _ROLE_HEAD_PREFIXES:
        bare = prefix.rstrip("-")
        if first == bare or first.startswith(prefix):
            return True
    return False


def _command_is_jarvis_lit(command: str) -> bool:
    """Return whether a command belongs to the separate Jarvis-Lit app."""
    text = str(command or "").strip()
    if not text:
        return False
    first = text.split(None, 1)[0]
    if first in _JARVIS_LIT_MARKERS or first.startswith(("Jarvis-Lit:", "JarvisLit:")):
        return True
    return any(
        f"/{marker}/" in text or f"/{marker}.app/" in text
        for marker in _JARVIS_LIT_MARKERS
    )

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
    if not marker:
        return None
    if (
        head in _CONTROL_HEADS
        or head.startswith("Jarvis-Worker-")
        or head == "Jarvis-FileOperation"
        or head == "Jarvis-Archiver"
        or head == "Jarvis-Redis"
    ):
        return name.split("@", 1)[0].strip() or None
    return None


def _is_unscoped_jarvis_title(command: str) -> bool:
    """Return whether ``command`` is a known Jarvis title without a scan id.

    ZP cleanup is narrow: only titles emitted by Jarvis itself are eligible,
    preventing a third-party command such as ``JarvisHelper`` from being
    killed merely because it shares a prefix.
    """
    title = str(command or "").strip().split(None, 1)[0]
    if title in _CONTROL_HEADS or title in {
        "Jarvis-Archiver",
        "Jarvis-FileOperation",
        "Jarvis-Redis",
    }:
        return True
    if re.fullmatch(r"Jarvis-Worker-\d+", title):
        return True
    return re.fullmatch(r"Jarvis-Redis@\d+/\d+", title) is not None


def list_zombie_processes(
    procs: Iterable[JarvisProcess] | None = None,
) -> list[JarvisProcess]:
    """Return known Jarvis components not owned by a live scan controller.

    This includes both bare legacy titles and components which retain a stale
    ``:<scan>`` suffix after their matching ``Jarvis:<scan>`` control process
    has disappeared.
    """
    candidates = list(procs if procs is not None else list_jarvis_processes())
    controlled_scans = {
        name
        for proc in candidates
        if _is_control_command(proc.command)
        if (name := _scan_name_from_command(proc.command)) is not None
    }
    zombies: list[JarvisProcess] = []
    for proc in candidates:
        name = _scan_name_from_command(proc.command)
        if name is None:
            if _is_unscoped_jarvis_title(proc.command):
                zombies.append(proc)
            continue
        if name not in controlled_scans:
            zombies.append(proc)
    return sorted(zombies, key=lambda proc: proc.pid)


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
    control_alive = any(
        _is_control_command(proc.command)
        for proc in scan.processes
    )
    metadata = runtime_metadata_for_scan(scan)
    if metadata is None:
        if control_alive and _scan_has_advertised_redis(scan):
            print("Refusing unverified runtime: Redis metadata is unavailable.", file=sys.stderr)
            return None
        if not control_alive:
            print(
                "Warning: stale Jarvis runtime has no live control process; "
                "orphan cleanup is allowed.",
                file=sys.stderr,
            )
            return {}
        return None
    control_pid = int(metadata.get("control_pid", -1))
    if control_pid not in {proc.pid for proc in scan.processes}:
        if control_alive:
            print("Refusing unverified runtime: control PID does not match metadata.", file=sys.stderr)
            return None
        print(
            "Warning: control process is gone; treating the verified runtime as stale.",
            file=sys.stderr,
        )
    return metadata


def _control_pid_for_processes(processes: Iterable[JarvisProcess]) -> int | None:
    """Return the live control PID for a scan group, if any."""
    control_pids = [proc.pid for proc in processes if _is_control_command(proc.command)]
    if not control_pids:
        return None
    return min(control_pids)


def _ref_registry_path() -> Path:
    override = str(os.environ.get(_ENV_REF_REGISTRY) or "").strip()
    if override:
        return Path(override).expanduser()
    return Path.home() / ".jarvis" / "scan_refs.json"


def _load_ref_registry() -> dict[str, object]:
    path = _ref_registry_path()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return {"slots": {}}
    if not isinstance(raw, dict):
        return {"slots": {}}
    slots = raw.get("slots")
    if not isinstance(slots, dict):
        raw["slots"] = {}
    return raw


def _save_ref_registry(data: dict[str, object]) -> None:
    path = _ref_registry_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        tmp.replace(path)
    except OSError:
        # Registry is a convenience; listing must still work if the home dir is RO.
        pass


def _sync_sticky_slots(live: dict[int, str]) -> dict[int, int]:
    """Map live control PIDs to sticky slot numbers (1 → R1, 2 → R2, …).

    Rules:
    - A control PID keeps the same slot for its whole lifetime.
    - Dead PIDs are dropped from the registry.
    - New PIDs take the lowest free slot (so you still get short R1/R2 labels).
    - Adding a new scan does **not** renumber existing ones (the original kill-wrong bug).
    """
    data = _load_ref_registry()
    old_slots = data.get("slots")
    if not isinstance(old_slots, dict):
        old_slots = {}

    live_pids = set(live)
    kept: dict[str, dict[str, object]] = {}
    pid_to_slot: dict[int, int] = {}
    for slot_key, meta in old_slots.items():
        if not isinstance(meta, dict):
            continue
        try:
            slot = int(slot_key)
            control_pid = int(meta.get("control_pid"))
        except (TypeError, ValueError):
            continue
        if slot < 1 or control_pid not in live_pids:
            continue
        if control_pid in pid_to_slot:
            continue  # duplicate registry row; keep first
        kept[str(slot)] = {
            "control_pid": control_pid,
            "name": str(live.get(control_pid) or meta.get("name") or ""),
        }
        pid_to_slot[control_pid] = slot

    used = set(pid_to_slot.values())

    def _next_free_slot() -> int:
        slot = 1
        while slot in used:
            slot += 1
        used.add(slot)
        return slot

    # Deterministic assignment for first-seen controls: lower control PID first.
    for control_pid, name in sorted(live.items(), key=lambda item: item[0]):
        if control_pid in pid_to_slot:
            continue
        slot = _next_free_slot()
        kept[str(slot)] = {"control_pid": control_pid, "name": name}
        pid_to_slot[control_pid] = slot

    _save_ref_registry({"slots": kept})
    return pid_to_slot


def _apply_sticky_references(scans: list[JarvisScan]) -> list[JarvisScan]:
    """Attach sticky ``R1``/``R2`` refs to controller-owned groups; orphans keep scan name."""
    live: dict[int, str] = {}
    for scan in scans:
        control_pid = _control_pid_for_processes(scan.processes)
        if control_pid is not None:
            live[control_pid] = scan.name
    pid_to_slot = _sync_sticky_slots(live) if live else {}

    updated: list[JarvisScan] = []
    for scan in scans:
        control_pid = _control_pid_for_processes(scan.processes)
        if control_pid is None:
            updated.append(
                JarvisScan(reference=scan.name, name=scan.name, processes=scan.processes)
            )
            continue
        slot = pid_to_slot[control_pid]
        updated.append(
            JarvisScan(
                reference=f"R{slot}",
                name=scan.name,
                processes=scan.processes,
            )
        )

    def _sort_key(scan: JarvisScan) -> tuple[int, str]:
        if re.fullmatch(r"R\d+", scan.reference):
            return (int(scan.reference[1:]), scan.name.casefold())
        return (10**9, scan.name.casefold())

    return sorted(updated, key=_sort_key)


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
    for name in sorted(grouped, key=str.casefold):
        processes = tuple(sorted(grouped[name], key=lambda proc: proc.pid))
        scans.append(
            JarvisScan(
                reference=name,  # replaced by sticky R# when a control is live
                name=name,
                processes=processes,
            )
        )
    return _apply_sticky_references(scans)


def list_active_scans(
    procs: Iterable[JarvisProcess] | None = None,
) -> list[JarvisScan]:
    """Return only controller-owned scan groups with sticky ``R1``/``R2``/… refs."""
    return [
        scan
        for scan in list_running_scans(procs)
        if _control_pid_for_processes(scan.processes) is not None
    ]


def ensure_scan_name_available(scan_name: str, *, cleanup_stale: bool = False) -> None:
    """Reject a new run when its exact scan name already has a live controller."""
    name = str(scan_name or "").strip()
    if not name:
        raise ValueError("scan name is required")
    for scan in list_running_scans():
        if scan.name != name:
            continue
        control_alive = any(
            _is_control_command(proc.command)
            for proc in scan.processes
        )
        if control_alive:
            raise RuntimeError(
                f"a Jarvis scan named {name!r} is already running; "
                f"inspect it with `Jarvis monitor {scan.reference}` or "
                f"`Jarvis ps {scan.reference}` before starting another run"
            )
        if not cleanup_stale:
            raise RuntimeError(
                f"a stale Jarvis runtime named {name!r} still has orphan processes; "
                f"clean it with `Jarvis kill {scan.reference} --yes` or rerun with --resume"
            )
        result = kill_jarvis_processes(list(scan.processes), force=True)
        if result["failed"]:
            raise RuntimeError(
                f"resume could not clean stale runtime {name!r}; failed pids: {result['failed']}"
            )
        remaining = [
            proc
            for proc in list_jarvis_processes()
            if _scan_name_from_command(proc.command) == name
        ]
        if remaining:
            raise RuntimeError(
                f"resume could not clean all stale processes for {name!r}: "
                + ", ".join(str(proc.pid) for proc in remaining)
            )
        return


def resolve_scan_reference(selector: str, scans: Iterable[JarvisScan]) -> JarvisScan:
    """Resolve sticky ``R1``/``R2``, bare control PID, or exact scan name."""
    text = str(selector or "").strip()
    choices = list(scans)
    bare_pid: int | None = None
    if re.fullmatch(r"[Rr]\d+", text):
        # Prefer sticky slot match (R1 → slot 1). Only fall back to control PID
        # when no sticky row claims that number (rare for large PIDs).
        slot_hit = next(
            (scan for scan in choices if scan.reference.upper() == text.upper()),
            None,
        )
        if slot_hit is not None:
            return slot_hit
        bare_pid = int(text[1:])
    elif text.isdigit():
        bare_pid = int(text)

    for scan in choices:
        if text == scan.name:
            return scan
        if bare_pid is not None and _control_pid_for_processes(scan.processes) == bare_pid:
            return scan
    available = ", ".join(f"{scan.reference} ({scan.name})" for scan in choices) or "none"
    raise ValueError(f"Unknown running scan {text!r}; available: {available}")


def format_scan_table(
    scans: list[JarvisScan],
    zombie_processes: Iterable[JarvisProcess] = (),
) -> str:
    """Render the compact task chooser shared by monitor / ps / kill."""
    zombies = tuple(zombie_processes)
    if not scans and not zombies:
        return "No running Jarvis scan tasks.\n"
    rows: list[tuple[str, str, str, tuple[JarvisProcess, ...]]] = []
    for scan in scans:
        control = _control_pid_for_processes(scan.processes)
        rows.append(
            (
                scan.reference,
                scan.name,
                str(control) if control is not None else "-",
                scan.processes,
            )
        )
    if zombies:
        rows.append((ZP_REFERENCE, ZP_LABEL, "-", zombies))
    ref_w = max(len(reference) for reference, _, _, _ in rows)
    name_w = max(len(name) for _, name, _, _ in rows)
    ctrl_w = max(len(control) for _, _, control, _ in rows)
    lines = [
        f"Running Jarvis process groups ({len(rows)}):",
        f"{'REF':<{ref_w}}  {'SCAN':<{name_w}}  {'CONTROL':<{ctrl_w}}  PROCESSES  PIDS",
        f"{'-' * ref_w}  {'-' * name_w}  {'-' * ctrl_w}  ---------  ----",
    ]
    for reference, name, control, processes in rows:
        pids = ",".join(str(proc.pid) for proc in processes)
        lines.append(
            f"{reference:<{ref_w}}  {name:<{name_w}}  {control:<{ctrl_w}}  "
            f"{len(processes):>9}  {pids}"
        )
    example = scans[0].reference if scans else "R1"
    lines.append(
        f"Select one (sticky REF: same R# while that control lives): "
        f"Jarvis monitor {example} | Jarvis ps {example} | Jarvis kill {example}"
    )
    lines.append("Also: Jarvis kill <Scan.name>  or  Jarvis kill <control_pid>")
    if zombies:
        lines.append("Clean unscoped processes: Jarvis ps ZP | Jarvis kill ZP -y")
    return "\n".join(lines) + "\n"


def print_scan_table(
    scans: list[JarvisScan],
    zombie_processes: Iterable[JarvisProcess] = (),
) -> None:
    """Render the running-task chooser with explicit, safe next actions."""
    zombies = tuple(zombie_processes)
    if not scans and not zombies:
        Console().print("[dim]No running Jarvis scan tasks.[/]")
        return
    table = Table(box=box.SIMPLE_HEAVY, show_header=True, header_style="bold dim")
    table.add_column("REF", style="bold green", no_wrap=True)
    table.add_column("SCAN", style="bold #c8c8ff")
    table.add_column("CONTROL", style="bold cyan", no_wrap=True)
    table.add_column("PROCESSES", justify="right")
    table.add_column("PIDS", style="dim")
    for scan in scans:
        control = _control_pid_for_processes(scan.processes)
        table.add_row(
            scan.reference,
            scan.name,
            str(control) if control is not None else "-",
            str(len(scan.processes)),
            ", ".join(str(proc.pid) for proc in scan.processes),
        )
    if zombies:
        table.add_row(
            ZP_REFERENCE,
            f"{ZP_LABEL} [dim](no live control)[/]",
            "-",
            str(len(zombies)),
            ", ".join(str(proc.pid) for proc in zombies),
        )
    example = scans[0].reference if scans else "R1"
    example_name = scans[0].name if scans else "<Scan.name>"
    guidance = Text()
    guidance.append(
        "Choose by sticky REF (R1/R2 stay put while that control lives), "
        "scan name, or control PID:\n",
        style="bold",
    )
    guidance.append(f"  Jarvis monitor {example}", style="bold cyan")
    guidance.append("  view this task's live scan status\n", style="dim")
    guidance.append(f"  Jarvis ps {example}", style="bold #c8c8ff")
    guidance.append(
        "       list this task's control, workers, file operations, archiver, and Redis\n",
        style="dim",
    )
    guidance.append(f"  Jarvis kill {example}", style="bold red")
    guidance.append(
        "     terminate ALL processes belonging to this task (asks for confirmation)\n",
        style="bold red",
    )
    guidance.append(f"  Jarvis kill {example_name}", style="bold red")
    guidance.append("  same task, selected by Scan.name", style="dim")
    if zombies:
        guidance.append("\n  Jarvis ps ZP", style="bold yellow")
        guidance.append("       inspect Jarvis processes without a live control\n", style="dim")
        guidance.append("  Jarvis kill ZP -y", style="bold red")
        guidance.append("  terminate ALL orphan Jarvis processes", style="bold red")
    Console().print(
        Panel(
            Group(table, Text(""), guidance),
            title=(
                "[bold]Running Jarvis process groups[/] "
                f"[dim]({len(scans) + int(bool(zombies))})[/]"
            ),
            box=box.ROUNDED,
            border_style="dim",
            padding=(0, 1),
        )
    )


def _process_role_and_detail(proc: JarvisProcess, scan_name: str) -> tuple[str, str]:
    """Turn a long OS title into the concise role shown in the runtime card."""
    title = proc.command.split(None, 1)[0]
    if _is_control_title(title):
        return "Control", "scan coordinator"
    if title.startswith("Jarvis-Worker-"):
        worker = title.removeprefix("Jarvis-Worker-").split(":", 1)[0]
        return f"Worker {worker}", "calculator / sampling worker"
    if title.startswith("Jarvis-FileOperation"):
        return "FileOperation", "Worker-owned SAMPLE save / copy / delete"
    if title.startswith("Jarvis-Archiver"):
        return "Archiver", "writes DATABASE and exports"
    if (
        title.startswith("Jarvis-Redis:")
        or title.startswith("Jarvis-Redis@")
        or title == "Jarvis-Redis"
    ):
        endpoint = proc.command.split(None, 1)[1] if " " in proc.command else ""
        advertised = title.rsplit("@", 1)[-1] if "@" in title else ""
        return "Redis", endpoint or advertised or "runtime broker"
    return "Jarvis process", title.replace(scan_name, "<scan>")


def print_zombie_processes(
    processes: Iterable[JarvisProcess],
    *,
    kill_warning: bool = False,
) -> None:
    """Render the fixed ZP group of Jarvis processes without a run id."""
    print_scan_processes(
        JarvisScan(
            reference=ZP_REFERENCE,
            name=ZP_LABEL,
            processes=tuple(processes),
        ),
        metadata=None,
        kill_warning=kill_warning,
    )


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
                    f"DANGER: Jarvis kill {scan.reference} will terminate all {len(scan.processes)} processes above.",
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
            "re-run with: Jarvis kill --yes",
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
    elif _is_control_command(cmd):
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
    """``Jarvis ps [R#|ZP]`` — inspect one scan or unscoped processes."""
    procs = list_jarvis_processes()
    scans = list_active_scans(procs)
    zombies = list_zombie_processes(procs)
    if not scan_ref:
        print_scan_table(scans, zombies)
        return 0
    if str(scan_ref).strip().upper() == ZP_REFERENCE:
        if not zombies:
            print("No ZP (unscoped Jarvis processes) found.")
            return 0
        print_zombie_processes(zombies)
        return 0
    try:
        scan = resolve_scan_reference(scan_ref, scans)
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        print_scan_table(scans, zombies)
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
    """``Jarvis kill [R#|ZP]`` — terminate one scan or unscoped processes."""
    procs = list_jarvis_processes()
    scans = list_active_scans(procs)
    zombies = list_zombie_processes(procs)
    if not scan_ref:
        print_scan_table(scans, zombies)
        return 0
    is_zombie_group = str(scan_ref).strip().upper() == ZP_REFERENCE
    if is_zombie_group:
        if not zombies:
            print("No ZP (unscoped Jarvis processes) found.")
            return 0
        targets = list(zombies)
        scan = None
        print_zombie_processes(targets, kill_warning=True)
    else:
        try:
            scan = resolve_scan_reference(scan_ref, scans)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            print_scan_table(scans, zombies)
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
    if is_zombie_group:
        left = list_zombie_processes()
    else:
        assert scan is not None
        left = [
            proc
            for proc in list_jarvis_processes()
            if _scan_name_from_command(proc.command) == scan.name
        ]
    if left:
        print("Still running after kill:")
        print(format_process_table(left), end="")
        return 1
    print("All matched Jarvis processes terminated.")
    return 0


__all__ = [
    "JarvisProcess",
    "JarvisScan",
    "ZP_LABEL",
    "ZP_REFERENCE",
    "confirm_kill",
    "ensure_scan_name_available",
    "format_process_table",
    "format_scan_table",
    "print_scan_table",
    "print_scan_processes",
    "kill_jarvis_processes",
    "kill_running_jarvis_cli",
    "list_jarvis_processes",
    "list_active_scans",
    "list_running_scans",
    "list_running_jarvis_cli",
    "list_zombie_processes",
    "print_zombie_processes",
    "resolve_scan_reference",
    "runtime_metadata_for_scan",
]
