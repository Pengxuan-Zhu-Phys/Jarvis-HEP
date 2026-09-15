"""Fleet-level WORKERS block for Overview; individual evidence stays on page 2."""

from __future__ import annotations

import math
from dataclasses import dataclass

from rich.cells import cell_len

from jarvishep2.monitor.styles import paint


_BACKGROUND_RGB = (16, 18, 22)
_LED_RGB = {
    "green": (53, 201, 138),
    "amber": (246, 211, 63),
    "red": (239, 107, 115),
    "idle": (74, 81, 96),
}


@dataclass(frozen=True)
class WorkerBlockData:
    expected: int | None = None
    present: int | None = None
    process_unknown: int = 0
    busy: int = 0
    idle: int = 0
    assigned: int = 0
    file_ops_expected: int = 0
    file_ops_alive: int = 0
    file_ops_missing: int = 0
    file_ops_unknown: int = 0
    file_ops_inline: int = 0
    heartbeat_recent: int = 0
    heartbeat_stale: int = 0
    heartbeat_unknown: int = 0
    source_stale: bool = False


@dataclass(frozen=True)
class WorkerBlockRender:
    border_title: str
    body: str
    border_subtitle: str


def _clip(text: str, width: int) -> str:
    if width <= 0:
        return ""
    if cell_len(text) <= width:
        return text
    if width == 1:
        return "…"
    result: list[str] = []
    used = 0
    for char in text:
        size = cell_len(char)
        if used + size > width - 1:
            break
        result.append(char)
        used += size
    return "".join(result) + "…"


def _dot(kind: str, phase: float) -> str:
    rgb = _LED_RGB[kind]
    amount = 1.0 if kind == "idle" else 0.58 + 0.42 * (1.0 - math.cos(phase)) / 2.0
    mixed = tuple(
        round(_BACKGROUND_RGB[index] + amount * (rgb[index] - _BACKGROUND_RGB[index]))
        for index in range(3)
    )
    return f"[#{mixed[0]:02x}{mixed[1]:02x}{mixed[2]:02x}]●[/]"


def _row(label: str, metric: str, kind: str, width: int, phase: float) -> str:
    metric = _clip(metric, width)
    label = _clip(f" {label}", max(0, width - cell_len(metric) - 1))
    gap = " " * max(0, width - 1 - cell_len(label) - cell_len(metric))
    return _dot(kind, phase) + paint("panel-title", label) + gap + f"[bold #e6e8eb]{metric}[/]"


def _note(text: str, width: int) -> str:
    text = _clip(text, width)
    return paint("dim", text + " " * max(0, width - cell_len(text)))


def _title(outer_width: int) -> str:
    available = max(len("WORKERS"), int(outer_width) - 6)
    gap = max(1, available - len("WORKERS") - len("FLEET") - 2)
    return f"WORKERS {'─' * gap} FLEET"


def render_workers_block(
    width: int,
    data: WorkerBlockData,
    *,
    outer_width: int | None = None,
    breathe_phase: float = 0.0,
) -> WorkerBlockRender:
    """Render the eight-row fleet view with FileOperator as a first-class path."""
    width = max(1, int(width))
    expected = data.expected
    present = data.present
    process_metric = (
        "—  waiting"
        if expected is None or present is None
        else f"{present} / {expected} present"
    )
    process_kind = "idle" if expected in (None, 0) else "green"
    if data.source_stale or data.process_unknown:
        process_kind = "red" if data.source_stale else "amber"
    elif expected and present is not None and present < expected:
        process_kind = "red"

    execution_metric = f"{data.busy} busy  ·  {data.idle} idle"
    execution_kind = "red" if data.source_stale else "green"
    if not data.busy and not data.idle and expected:
        execution_kind = "amber"

    if data.file_ops_expected:
        file_ops_metric = (
            f"{data.file_ops_alive} / {data.file_ops_expected} alive"
            f"  ·  {data.file_ops_missing} miss"
        )
        file_ops_kind = "green"
        if data.source_stale or data.file_ops_missing:
            file_ops_kind = "red"
        elif data.file_ops_unknown:
            file_ops_kind = "amber"
    elif data.file_ops_inline:
        file_ops_metric = f"—  {data.file_ops_inline} inline"
        file_ops_kind = "idle"
    else:
        file_ops_metric = "—  waiting"
        file_ops_kind = "amber"

    heartbeat_metric = f"{data.heartbeat_stale} stale  ·  {data.heartbeat_recent} recent"
    heartbeat_kind = "green"
    if data.source_stale or data.heartbeat_stale:
        heartbeat_kind = "red"
    elif data.heartbeat_unknown:
        heartbeat_kind = "amber"

    rows = (
        _row("PROCESSES", process_metric, process_kind, width, breathe_phase),
        _note("  OS-admitted Workers · Core target", width),
        _row("EXECUTING", execution_metric, execution_kind, width, breathe_phase + 1.2),
        _note(f"  Worker boards · {data.assigned} samples assigned", width),
        _row("FILE OPS", file_ops_metric, file_ops_kind, width, breathe_phase + 2.4),
        _row("HEARTBEAT", heartbeat_metric, heartbeat_kind, width, breathe_phase + 3.6),
    )
    return WorkerBlockRender(
        border_title=_title(outer_width or width + 4),
        body="\n".join(rows),
        border_subtitle="OS inventory + known Worker boards",
    )


__all__ = ["WorkerBlockData", "WorkerBlockRender", "render_workers_block"]
