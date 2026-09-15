"""Project ``list_active_scans`` rows for the splash chooser."""

from __future__ import annotations

from dataclasses import dataclass

from jarvishep2.process_cleanup import (
    JarvisScan,
    _control_pid_for_processes,
    list_active_scans,
    resolve_scan_reference,
)


@dataclass(frozen=True)
class ScanChoice:
    reference: str
    name: str
    control_pid: int | None
    process_count: int
    pids: tuple[int, ...]
    simulated: bool = False


class ScanExitWatch:
    """Observe the attached Core's process identity, including PID reuse."""

    def __init__(self, choice: ScanChoice) -> None:
        self._pid = None if choice.simulated else choice.control_pid
        self._process = None
        self._exited = False
        # Capture identity on attach, before the first periodic check.
        self.has_exited()

    def has_exited(self) -> bool:
        if self._pid is None:
            return False
        if self._exited:
            return True
        try:
            import psutil
        except ImportError:
            return False
        try:
            if self._process is None:
                self._process = psutil.Process(self._pid)
            self._exited = (
                not self._process.is_running()
                or self._process.status() in {psutil.STATUS_ZOMBIE, psutil.STATUS_DEAD}
            )
        except psutil.NoSuchProcess:
            self._exited = True
        except psutil.AccessDenied:
            # Lack of permission is not evidence that the scan exited.
            pass
        return self._exited


def choice_from_scan(scan: JarvisScan) -> ScanChoice:
    control = _control_pid_for_processes(scan.processes)
    return ScanChoice(
        reference=str(scan.reference),
        name=str(scan.name),
        control_pid=control,
        process_count=len(scan.processes),
        pids=tuple(proc.pid for proc in scan.processes),
        simulated=False,
    )


def simulated_choice() -> ScanChoice:
    """In-process fixture for monitor layout work. Not a live scan."""
    control = 44001
    nproc = 52
    return ScanChoice(
        reference="SIM",
        name="simu-iDM_Vector_V1",
        control_pid=control,
        process_count=nproc,
        pids=tuple(range(control, control + nproc)),
        simulated=True,
    )


def list_scan_choices(
    lister: object | None = None,
) -> list[ScanChoice]:
    """OS inventory of live controller-owned scans (same as ``Jarvis ps``)."""
    fetch = list_active_scans if lister is None else lister
    scans = list(fetch())
    return [choice_from_scan(scan) for scan in scans]


def resolve_choice(
    selector: str,
    choices: list[ScanChoice],
    *,
    scans: list[JarvisScan] | None = None,
    lister: object | None = None,
) -> ScanChoice:
    """Resolve sticky REF / name / control PID using the ps resolver."""
    fetch = list_active_scans if lister is None else lister
    live = list(scans if scans is not None else fetch())
    scan = resolve_scan_reference(selector, live)
    return choice_from_scan(scan)


__all__ = [
    "ScanChoice",
    "ScanExitWatch",
    "choice_from_scan",
    "list_scan_choices",
    "resolve_choice",
    "simulated_choice",
]
