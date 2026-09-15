"""Read-only Workspace pages backed by the Collector's standard snapshot."""

from __future__ import annotations

from typing import Any, Mapping

from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import Static

from jarvishep2.monitor.collector import CollectorFrame


def _text(value: Any, fallback: str = "—") -> str:
    text = str(value or "").strip()
    return text or fallback


def _age(value: Any) -> str:
    try:
        return f"{max(0.0, float(value)):.1f}s"
    except (TypeError, ValueError):
        return "—"


def _board_lines(title: str, board: Mapping[str, Any]) -> list[str]:
    if not board:
        return [title, "", "No live board is available."]
    lines = [title, ""]
    for key, value in sorted(board.items()):
        lines.append(f"{str(key):<20} {_text(value)}")
    return lines


class WorkersPane(Horizontal):
    """Known Worker boards only; discovery remains the OS inventory's job."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "workers")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="workers")

    def compose(self) -> ComposeResult:
        yield Static(id="workers-list")
        yield Static(id="workers-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#workers-list", Static)
            right = self.query_one("#workers-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        rows = self.frame.workers
        lines = [f"WORKERS{stale}  known process boards", "", "ID   PID     STATUS    HB       CALC  UUID"]
        if rows:
            for row in rows:
                lines.append(
                    f"{_text(row.get('worker_id')):>3}  "
                    f"{_text(row.get('pid')):<6.6}  "
                    f"{_text(row.get('status'), 'unknown'):<8.8}  "
                    f"{_age(row.get('heartbeat_age_s')):<7}  "
                    f"{_text(row.get('held_calc_n')):>4}  "
                    f"{_text(row.get('current_uuid')):<24.24}"
                )
        else:
            lines.append("—    no Worker board is visible for this scan")
        left.update("\n".join(lines))

        selected = rows[0] if rows else None
        if selected is None:
            right.update("WORKER DETAIL\n\nNo Worker is selected.")
            return
        detail = [
            "WORKER DETAIL",
            "",
            f"worker         {_text(selected.get('worker_id'))}",
            f"pid            {_text(selected.get('pid'))}",
            f"status         {_text(selected.get('status'), 'unknown')}",
            f"uuid           {_text(selected.get('current_uuid'))}",
            f"heartbeat      {_age(selected.get('heartbeat_age_s'))}",
            f"held_calc      {_text(selected.get('held_calc_n'))}",
            f"file_operator  {_text(selected.get('file_operation_pid'))}",
            f"fo_pgid         {_text(selected.get('file_operation_pgid'))}",
            f"fo_mode         {_text(selected.get('file_operation_mode'), 'unknown')}",
            f"alive          {_text(selected.get('alive'), 'unknown')}",
            "",
            "Logs and result payloads are not fetched.",
        ]
        right.update("\n".join(detail))


class FactoryPane(Horizontal):
    """Control-plane boards: Core, Archiver, and managed Redis."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "factory")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="factory")

    def compose(self) -> ComposeResult:
        yield Static(id="factory-list")
        yield Static(id="factory-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#factory-list", Static)
            right = self.query_one("#factory-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        boards = (
            ("Core", self.frame.proc_core),
            ("Archiver", self.frame.proc_archiver),
            ("Redis", self.frame.proc_redis),
        )
        lines = [f"FACTORY{stale}  control plane", "", "ROLE       PID     STATUS       HOST"]
        for role, board in boards:
            lines.append(
                f"{role:<10} {_text(board.get('pid')):<6.6}  "
                f"{_text(board.get('status'), 'unknown'):<11.11}  "
                f"{_text(board.get('host')):<20.20}"
            )
        left.update("\n".join(lines))

        role, selected = next(((role, board) for role, board in boards if board), ("Core", {}))
        right.update("\n".join(_board_lines(f"{role.upper()} DETAIL", selected)))


class SamplerPane(Horizontal):
    """Counters and queue lengths only; queue payloads stay private work items."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "sampler")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="sampler")

    def compose(self) -> ComposeResult:
        yield Static(id="sampler-list")
        yield Static(id="sampler-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#sampler-list", Static)
            right = self.query_one("#sampler-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        queues = self.frame.queues
        stats = self.frame.sample_stats
        lines = [
            f"SAMPLER{stale}  counters + queue depths",
            "",
            f"completed  {_text(stats.get('completed'), '0')}    running  {_text(stats.get('running'), '0')}    failed  {_text(stats.get('failed'), '0')}",
            "",
            "QUEUE              LENGTH",
            f"task               {_text(queues.get('task_queue_length'), '0')}",
            f"archive            {_text(queues.get('archive_queue_length'), '0')}",
            "",
            "Queue payloads are intentionally not read.",
        ]
        left.update("\n".join(lines))

        detail = ["SAMPLER DETAIL", "", "operation counters"]
        if self.frame.op_counts:
            for kind, count in sorted(self.frame.op_counts.items()):
                detail.append(f"{kind:<16} {_text(count, '0')}")
        else:
            detail.append("—    no counters published")
        detail.extend(["", "Metadata-backed method details are staged separately."])
        right.update("\n".join(detail))


class HostPane(Horizontal):
    """Host aggregates and the process inventory already scoped to this scan."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "host")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="host")

    def compose(self) -> ComposeResult:
        yield Static(id="host-list")
        yield Static(id="host-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#host-list", Static)
            right = self.query_one("#host-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        host = self.frame.host
        if not host.get("available"):
            left.update(
                f"HOST{stale}\n\nHost metrics are unavailable: {_text(host.get('reason'), 'unknown error')}"
            )
            right.update("PROCESS DETAIL\n\nNo process is selected.")
            return
        mem_used = float(host.get("memory_used", 0) or 0) / (1024**3)
        mem_total = float(host.get("memory_total", 0) or 0) / (1024**3)
        swap_used = float(host.get("swap_used", 0) or 0) / (1024**3)
        swap_total = float(host.get("swap_total", 0) or 0) / (1024**3)
        load = host.get("load") or ()
        load_text = " / ".join(f"{float(value):.2f}" for value in load) or "—"
        rows = list(host.get("processes") or [])
        lines = [
            f"HOST{stale}  scan process inventory",
            "",
            f"CPU  {float(host.get('cpu_percent', 0) or 0):.1f}%    LOAD  {load_text}",
            f"MEM  {mem_used:.1f}/{mem_total:.1f}G    SWAP  {swap_used:.1f}/{swap_total:.1f}G",
            "",
            "ROLE       PID     CPU%   RSS      FDS  ALIVE",
        ]
        if rows:
            for row in rows:
                rss = float(row.get("rss", 0) or 0) / (1024**2)
                lines.append(
                    f"{_text(row.get('role'), 'process'):<10.10} "
                    f"{_text(row.get('pid')):<6.6}  "
                    f"{float(row.get('cpu_percent', 0) or 0):>4.1f}  "
                    f"{rss:>6.1f}M  {_text(row.get('fds')):>3}  {_text(row.get('alive'))}"
                )
        else:
            lines.append("—    no process from this scan is available")
        left.update("\n".join(lines))

        selected = rows[0] if rows else None
        if selected is None:
            right.update("PROCESS DETAIL\n\nNo process is selected.")
            return
        detail = [
            "PROCESS DETAIL",
            "",
            f"role           {_text(selected.get('role'), 'process')}",
            f"pid            {_text(selected.get('pid'))}",
            f"parent         {_text(selected.get('ppid'))}",
            f"threads        {_text(selected.get('threads'))}",
            f"fds            {_text(selected.get('fds'))}",
            f"alive          {_text(selected.get('alive'))}",
            "",
            f"cmdline        {_text(selected.get('cmdline'))}",
        ]
        right.update("\n".join(detail))


__all__ = ["FactoryPane", "HostPane", "SamplerPane", "WorkersPane"]
