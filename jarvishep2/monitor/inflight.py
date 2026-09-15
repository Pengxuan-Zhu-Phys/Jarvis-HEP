"""Read-only live panes for calculator and sample inflight telemetry."""

from __future__ import annotations

from typing import Any

from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import Static

from jarvishep2.monitor.collector import CollectorFrame


def format_worker_owner(owner: Any) -> str:
    """Render Redis's unpadded Worker id without mistaking pool state for one."""
    try:
        identifier = int(str(owner).strip())
    except (TypeError, ValueError):
        text = str(owner or "").strip()
        return f"worker-{text}" if text else "worker-?"
    return f"worker-{identifier:02d}"


def _age_text(value: Any) -> str:
    try:
        return f"{max(0.0, float(value)):.1f}s"
    except (TypeError, ValueError):
        return "—"


class CalculatorsPane(Horizontal):
    """Known pool occupancy plus owner sidecar for one selected pool only."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "calculators")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="calculators")

    def compose(self) -> ComposeResult:
        yield Static(id="calculators-list")
        yield Static(id="calculators-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#calculators-list", Static)
            right = self.query_one("#calculators-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        rows = self.frame.occupancy
        selected = self.frame.selected_calc
        if not rows:
            left.update(
                f"CALCULATORS{stale}  pack pools\n\nNo pool metadata is available for this scan."
            )
        else:
            lines = [f"CALCULATORS{stale}  pack pools", "", "NAME           BUSY  FREE  TOTAL  UTIL"]
            for name, row in rows.items():
                busy = int(row.get("busy", 0) or 0)
                free = int(row.get("free", 0) or 0)
                slots = int(row.get("slots", 0) or 0)
                util = 100 * busy / max(1, slots)
                marker = ">" if name == selected else " "
                lines.append(
                    f"{marker}{name:<14.14} {busy:>4}  {free:>4}  {slots:>5}  {util:>3.0f}%"
                )
            left.update("\n".join(lines))

        if not selected:
            right.update("SELECTED CALCULATOR\n\nNo calculator pool is selected.")
            return
        owners = self.frame.calc_busy
        detail = [selected, "", "pack   owner"]
        if owners:
            for pack, owner in sorted(owners.items()):
                detail.append(f"{str(pack):<6} {format_worker_owner(owner)}")
        else:
            detail.extend(["", "waiting for next acquire", "", "Owners are intentionally shown only", "after a monitored acquire event."])
        detail.extend(
            [
                "",
                f"hash   hep:monitor:calc:busy:{selected}",
                f"busy   calc:busy:{selected}",
            ]
        )
        right.update("\n".join(detail))


class SamplesPane(Horizontal):
    """Worker-board rows with optional sample overlay fields on the right."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "samples")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="samples")

    def compose(self) -> ComposeResult:
        yield Static(id="samples-list")
        yield Static(id="samples-detail")

    def on_mount(self) -> None:
        self.paint()

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            left = self.query_one("#samples-list", Static)
            right = self.query_one("#samples-detail", Static)
        except Exception:
            return
        stale = "  [STALE]" if self.frame.stale else ""
        stats = self.frame.sample_stats
        active = [row for row in self.frame.workers if str(row.get("current_uuid") or "")]
        lines = [
            f"SAMPLES{stale}  counters + currently running",
            "",
            "Redis has no cheap UUID catalogue. This is not a browser.",
            "",
            f"completed  {int(stats.get('completed', 0) or 0)}    running  {int(stats.get('running', 0) or 0)}    failed  {int(stats.get('failed', 0) or 0)}",
            "",
            "WRK  UUID                              STAT      PID",
        ]
        if active:
            for row in active:
                worker = str(row.get("worker_id") or "?")
                uuid = str(row.get("current_uuid") or "")
                status = str(row.get("status") or "unknown")
                pid = str(row.get("pid") or "—")
                lines.append(f"{worker:>3}  {uuid:<32.32}  {status:<8.8}  {pid}")
        else:
            lines.append("—    no Worker currently advertises a sample")
        left.update("\n".join(lines))

        selected = active[0] if active else None
        if selected is None:
            right.update("RUNNING SAMPLE\n\nNo running Worker is visible on the process boards.")
            return
        worker_id = str(selected.get("worker_id") or "")
        overlay = self.frame.sample_running.get(worker_id) or {}
        step = str(overlay.get("step") or "").strip()
        detail = [
            "RUNNING SAMPLE",
            "",
            f"uuid           {selected.get('current_uuid') or '—'}",
            f"worker         {format_worker_owner(worker_id)}",
            f"pid            {selected.get('pid') or '—'}",
            f"status         {selected.get('status') or 'unknown'}",
            f"hb             {_age_text(selected.get('heartbeat_age_s'))}",
            f"held_calc      {selected.get('held_calc_n') if selected.get('held_calc_n') is not None else '—'}",
        ]
        if step:
            detail.append(f"step           {step}")
        else:
            detail.append("step           unknown until next work event")
        if overlay.get("t0") not in (None, ""):
            detail.append(f"started        {overlay['t0']}")
        detail.extend(["", "result         not fetched (v1)", "logs           not tailed (v1)"])
        right.update("\n".join(detail))


__all__ = ["CalculatorsPane", "SamplesPane", "format_worker_owner"]
