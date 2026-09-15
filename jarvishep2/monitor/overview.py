"""Overview pane: rounded fieldset blocks (总 / 分) inside the page frame."""

from __future__ import annotations

from collections import deque
from collections.abc import Callable
from dataclasses import dataclass
import math
import random

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.events import Click, Leave, MouseMove, Resize
from textual.screen import ModalScreen
from textual.widgets import Button, Static

from jarvishep2.monitor.bars import bar_rail, visible_width
from jarvishep2.monitor.calculators import render_calculators_block
from jarvishep2.monitor.collector import CollectorFrame
from jarvishep2.monitor.metrics import format_sample_rate
from jarvishep2.monitor.queues import render_queue_block
from jarvishep2.monitor.samples import (
    SamplesBlockData,
    render_samples_block,
    sample_border_title,
)
from jarvishep2.monitor.workers import render_workers_block
from jarvishep2.monitor.simu import OverviewFrame, simu_overview
from jarvishep2.monitor.styles import paint

SPARK_CHARS = "▁▂▃▄▅▆▇█"
SPARK_ROWS = 4
SPARK_SEPARATOR = "  ·  "
SPARK_METRICS = 2
# Fixed four-row ramp: top is Jarvis light blue, bottom is Jarvis dark blue.
SPARK_ROW_COLORS = ("#73b8f4", "#5393d1", "#336eaf", "#134a8d")
STATUS_WHITE_RGB = (230, 232, 235)
STATUS_YELLOW_RGB = (246, 211, 63)
HEALTH_COMPONENT_WIDTH = 10
# Keep the help trigger at the leading edge. The status text and the right-hand
# fact may change on every refresh, but the click target must never move.
HEALTH_INFO_START = 0
HEALTH_INFO_WIDTH = visible_width("💡")


@dataclass(frozen=True)
class HealthBlockRender:
    """Rendered half-width HEALTH component summary for the simulated Overview."""

    border_title: str
    body: str
    border_subtitle: str
    items: tuple["HealthItem", ...]


@dataclass(frozen=True)
class HealthItem:
    """One interactive HEALTH row and the context used by its help modal."""

    component: str
    state: str
    fact: str
    compact_fact: str


def _spark_column_widths(width: int, count: int) -> list[int]:
    """Return widths that consume the full plot, giving remainder to the last metric."""
    if width <= 0 or count <= 0:
        return []
    usable = max(count, width - len(SPARK_SEPARATOR) * (count - 1))
    base, remainder = divmod(usable, count)
    return [base] * (count - 1) + [base + remainder]


def _fit_spark_series(
    values: tuple[float | None, ...], width: int
) -> list[float | None]:
    """Validate a one-to-one ring-to-column mapping and reverse it for display."""
    if width <= 0 or not values:
        return []
    if len(values) != width:
        raise ValueError(
            f"SPARKS history has {len(values)} columns, expected display width {width}"
        )
    # The ring mutates as [new, old...] and is reversed only for display:
    # left is oldest, right is newest, so every refresh scrolls left.
    return list(reversed(values))


def _balanced_spark_width(width: int, count: int) -> int:
    """Use the full available plot width; the last metric owns any odd remainder."""
    if width <= 0 or count <= 0:
        return 0
    return width


def _spark_area_rows(
    values: tuple[float | None, ...],
    width: int,
    *,
    height: int = SPARK_ROWS,
    scale: float | None = None,
) -> tuple[str, ...]:
    """Render one history ring as a stacked, bottom-filled sparkline."""
    if width <= 0 or height <= 0:
        return tuple("" for _ in range(max(0, height)))
    if not values:
        return tuple(" " * width for _ in range(height))
    sequence = _fit_spark_series(values, width)

    observed = [value for value in sequence if value is not None]
    if not observed:
        return tuple(" " * width for _ in range(height))
    if scale is not None and scale > 0:
        levels = [
            None
            if value is None
            else max(0.0, min(float(height), value / scale * height))
            for value in sequence
        ]
    else:
        lo, hi = min(observed), max(observed)
        span = max(1, hi - lo)
        levels = [
            None
            if value is None
            else (
                float(height)
                if hi == lo and hi != 0
                else max(0.0, min(float(height), (value - lo) / span * height))
            )
            for value in sequence
        ]
    rows: list[str] = []
    for row in range(height - 1, -1, -1):
        chars: list[str] = []
        for level in levels:
            if level is None or level <= row:
                chars.append(" ")
            elif level >= row + 1:
                chars.append("█")
            else:
                fraction = level - row
                index = max(
                    0,
                    min(
                        len(SPARK_CHARS) - 1,
                        int(round(fraction * len(SPARK_CHARS))) - 1,
                    ),
                )
                chars.append(SPARK_CHARS[index])
        rows.append("".join(chars))
    return tuple(rows)


def _edge_align(left: str, right: str, width: int) -> str:
    """Place a metric label at the left edge and its value at the right edge."""
    if width <= 0:
        return ""
    left = left[:width]
    right = right[:width]
    if len(left) + len(right) > width:
        right = right[-max(0, width - len(left)) :]
    return left + right.rjust(max(0, width - len(left)))


def _edge_align_markup(left: str, right: str, width: int) -> str:
    """Edge-align a row whose left side may contain Rich markup."""
    gap = max(0, width - visible_width(left) - visible_width(right))
    return left + " " * gap + right


def _sample_border_title(width: int, method: str) -> str:
    """Compatibility alias for the shared SAMPLES renderer."""
    return sample_border_title(width, method)


def render_adaptive_bridson_samples(width: int, frame: OverviewFrame) -> str:
    """Compatibility wrapper around the shared sampler-specific renderer."""
    return render_samples_block(
        width,
        SamplesBlockData(
            completed=frame.done,
            running=frame.running,
            failed=frame.failed,
            rate=frame.rate,
            sampler=frame.sampler,
        ),
    ).body


def _status_breathe_color(phase: float) -> str:
    """Blend Jarvis white to yellow and back over one breathing cycle."""
    amount = _status_breathe_amount(phase)
    rgb = tuple(
        round(STATUS_WHITE_RGB[index] + amount * (STATUS_YELLOW_RGB[index] - STATUS_WHITE_RGB[index]))
        for index in range(3)
    )
    return f"#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}"


def _status_breathe_amount(phase: float) -> float:
    """Return zero at the white end and one at the yellow end."""
    return (1.0 - math.cos(phase)) / 2.0


def _status_breathe_dots(phase: float) -> str:
    """Render three dots with their white-to-yellow breathing phases offset."""
    offsets = (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)
    return "".join(
        f"[{_status_breathe_color(phase + offset)}]●[/]" for offset in offsets
    )


def render_status_block(
    width: int,
    scan: str,
    mode: str,
    elapsed: str,
    ref: str,
    redis: str,
    *,
    breathe_phase: float = 0.0,
) -> str:
    """Render STATUS as two edge-aligned rows with no internal blank row."""
    first = _edge_align(
        scan,
        f"●●● {mode.lower()}  ·  {elapsed}",
        width,
    )
    if mode.lower() == "running":
        dots = _status_breathe_dots(breathe_phase)
        dots_at = first.index("●●●")
        first = first[:dots_at] + dots + first[dots_at + 3 :]
    second = _edge_align(f"REF  {ref}", f"Redis  {redis}", width)
    return f"{first}\n{second}"


def _health_border_title(
    outer_width: int, verdict: str, healthy: int, applicable: int
) -> str:
    """Keep HEALTH left and the global verdict on the top-right edge."""
    base = "HEALTH"
    right = f"{verdict} · {healthy}/{applicable}"
    available = max(len(base), int(outer_width) - 6)
    if len(base) + len(right) + 2 > available:
        right = verdict[: max(1, available - len(base) - 2)]
    gap = max(1, available - len(base) - len(right) - 2)
    return f"{base} {'─' * gap} {right}"


def _health_row(item: HealthItem, width: int) -> str:
    """Render one HEALTH row with a fixed help column and right-aligned fact."""
    normalized = str(item.state or "UNKNOWN").upper()
    state_role = {
        "HEALTHY": "health-good",
        "DEGRADED": "health-warn",
        "CRITICAL": "health-critical",
        "UNKNOWN": "health-idle",
        "PAUSED": "health-idle",
        "DRAINING": "health-idle",
        "COMPLETE": "health-idle",
    }.get(normalized, "health-idle")
    dot = "●" if normalized in {"HEALTHY", "DEGRADED", "CRITICAL"} else "◌"
    led_role = {
        "HEALTHY": "health-led-good",
        "DEGRADED": "health-led-warn",
        "CRITICAL": "health-led-critical",
    }.get(normalized, "health-idle")
    left = (
        paint("health-info", "💡")
        + " "
        + f"{item.component:<{HEALTH_COMPONENT_WIDTH}}"
        + paint(led_role, dot)
        + " "
        + paint(state_role, normalized)
    )
    # The help marker is deliberately not part of the flowing information
    # string. It remains at the leading edge while all live text refreshes.
    fact = str(item.fact)
    if visible_width(left) + visible_width(fact) > width:
        fact = str(item.compact_fact)
    return _edge_align_markup(left, fact, max(1, width))


def render_health_block(width: int, frame: OverviewFrame) -> HealthBlockRender:
    """Render the four Jarvis-component HEALTH rows for the simulated Overview.

    This uses the lightweight simulation frame. The live collector can adopt the same
    layout after its detailed health evaluator is implemented.
    """
    if frame.health_items:
        items = tuple(HealthItem(*item) for item in frame.health_items)
    else:
        mode = str(frame.mode or "unknown").lower()
        outstanding = frame.task_q > 0 or frame.running > 0 or frame.done < frame.target

        core_state, core_fact = "HEALTHY", "hb 0.3s"
        if mode == "stopping":
            core_state, core_fact = "DRAINING", "stopping"

        factory_state, factory_fact = "HEALTHY", "dispatching"
        if mode == "paused":
            factory_state, factory_fact = "CRITICAL", "watchdog paused"
        elif mode == "degraded":
            factory_state, factory_fact = "DEGRADED", "watchdog degraded"
        elif mode in {"stopping", "draining"} or (
            not outstanding and frame.archive_q > 0
        ):
            factory_state, factory_fact = "DRAINING", "draining"
        elif not outstanding:
            factory_state, factory_fact = "COMPLETE", "complete"

        worker_state, worker_fact = (
            "HEALTHY",
            f"{frame.workers_alive}/{frame.workers_total} alive",
        )
        if mode == "paused":
            worker_state, worker_fact = "PAUSED", "watchdog"
        elif outstanding and frame.workers_total > 0 and frame.workers_alive <= 0:
            worker_state, worker_fact = "CRITICAL", "0 workers alive"
        elif frame.stale > 0:
            worker_state, worker_fact = "DEGRADED", f"{frame.stale} stale"

        archiver_state = "HEALTHY"
        # Keep the queue fact short enough to retain right alignment at 80 columns.
        # "draining" is implicit in HEALTHY and is explained by the block design.
        archiver_fact = "queue empty" if frame.archive_q <= 0 else f"queue {frame.archive_q}"
        items = (
            HealthItem("CORE", core_state, core_fact, "0.3s"),
            HealthItem("FACTORY", factory_state, factory_fact, "run"),
            HealthItem("WORKERS", worker_state, worker_fact, f"{frame.workers_alive}/{frame.workers_total}"),
            HealthItem("ARCHIVER", archiver_state, archiver_fact, f"q {frame.archive_q}"),
        )
    applicable_states = {"HEALTHY", "DEGRADED", "CRITICAL", "UNKNOWN"}
    applicable = [item for item in items if item.state in applicable_states]
    healthy = sum(1 for item in applicable if item.state == "HEALTHY")
    states = {item.state for item in applicable}
    if "CRITICAL" in states:
        verdict = "CRITICAL"
    elif "DEGRADED" in states:
        verdict = "DEGRADED"
    elif "UNKNOWN" in states or not applicable:
        verdict = "UNKNOWN"
    else:
        verdict = "HEALTHY"

    issue = "No active warnings"
    priority = ("CORE", "FACTORY", "WORKERS", "ARCHIVER")
    rank = {"CRITICAL": 0, "DEGRADED": 1, "UNKNOWN": 2}
    active = [item for item in items if item.state in rank]
    if active:
        item = min(
            active,
            key=lambda row: (rank[row.state], priority.index(row.component)),
        )
        if item.state != "HEALTHY":
            prefix = (
                "CRITICAL"
                if item.state == "CRITICAL"
                else "WARN"
                if item.state == "DEGRADED"
                else "UNKNOWN"
            )
            issue = f"{prefix} · {item.component.lower()} {item.fact}"

    return HealthBlockRender(
        border_title=_health_border_title(width, verdict, healthy, len(applicable)),
        body="\n".join(_health_row(item, width) for item in items),
        border_subtitle=issue,
        items=items,
    )


_HEALTH_EXPLANATIONS: dict[str, tuple[str, str, str, str]] = {
    "CORE": (
        "Jarvis2Core, Redis reachability, and the control-lock lease.",
        "Core is coordinating the scan and still owns its control path.",
        "Redis timeout, missing/mismatched lock, or an expired Core heartbeat.",
        "hep:proc:core · hep:proc:redis · PING · control-lock TTL",
    ),
    "FACTORY": (
        "TaskFactory inside Core: dispatch, watchdog, worker replacement, and scan progress.",
        "Tasks are being dispatched or the sampler is making observable progress.",
        "Watchdog fuse paused the scan, workers cannot make progress, or work stalls.",
        "Core Factory projection · task queue · sample counters · sampler status",
    ),
    "WORKERS": (
        "The spawned Worker pool and calculator ownership needed to execute samples.",
        "Usable workers are alive; calculator ownership and resource pressure are sane.",
        "Stale/dead workers, no workers with work outstanding, lost calculator owner, or pressure.",
        "Worker proc/status boards · OS process group · known calculator pools",
    ),
    "ARCHIVER": (
        "The Archiver process and its progress writing completed records to the archive.",
        "The archive queue is empty or draining while written-record count advances.",
        "Archiver dies, its heartbeat expires, or archive backlog stops making write progress.",
        "hep:proc:archiver · archive queue length · Archiver OS process",
    ),
}


def _health_modal_body(item: HealthItem) -> str:
    """Explain one component without leaking raw payloads or implementation detail."""
    watches, healthy, attention, source = _HEALTH_EXPLANATIONS[item.component]
    state = str(item.state).upper()
    return "\n\n".join(
        (
            paint("panel-title", "WHAT IT WATCHES") + f"\n{watches}",
            paint("panel-title", "CURRENT")
            + f"\n{item.component} · {state} · {item.fact}",
            paint("panel-title", "HEALTHY MEANS") + f"\n{healthy}",
            paint("panel-title", "NEEDS ATTENTION WHEN") + f"\n{attention}",
            paint("dim", "SOURCE") + f"\n{source}",
        )
    )


class HealthDetailModal(ModalScreen[None]):
    """Centered explanatory overlay opened from a HEALTH row's 💡 trigger."""

    BINDINGS = [
        Binding("escape", "close", "Close", show=False),
        Binding("q", "close", "Close", show=False),
    ]

    def __init__(self, item: HealthItem) -> None:
        super().__init__()
        self._item = item

    def compose(self) -> ComposeResult:
        with Vertical(id="health-detail-dialog"):
            yield Static(
                paint("panel-title", f"💡 HEALTH · {self._item.component}"),
                id="health-detail-title",
            )
            yield Static(_health_modal_body(self._item), id="health-detail-body")
            yield Button("Close", id="health-detail-close")

    def action_close(self) -> None:
        self.dismiss(None)

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "health-detail-close":
            self.dismiss(None)

    def on_click(self, event: Click) -> None:
        if event.widget is self:
            self.dismiss(None)


class HealthInfoRow(Static):
    """One click-targeted HEALTH row; only its 💡 marker opens the explanation."""

    def __init__(self, component: str) -> None:
        super().__init__(id=f"ov-health-{component.lower()}", classes="health-row")
        self._component = component
        self._item = HealthItem(component, "UNKNOWN", "waiting", "—")
        self._width = 1
        self._info_start = HEALTH_INFO_START
        self._info_width = HEALTH_INFO_WIDTH
        self.tooltip = None

    def set_item(self, item: HealthItem, width: int) -> None:
        self._item = item
        self._width = max(1, width)
        self._refresh()

    def _refresh(self) -> None:
        width = self.content_region.width or self.size.width or self._width
        self._width = max(1, int(width))
        self.update(_health_row(self._item, self._width))

    def on_resize(self, _event: Resize) -> None:
        # A terminal resize can leave a row narrower than the width used by
        # the parent paint pass.  Reflow against the row's own content width
        # so the right edge is never clipped.
        self._refresh()

    def _over_info(self, x: int) -> bool:
        return self._info_start <= x < self._info_start + self._info_width

    def on_mouse_move(self, event: MouseMove) -> None:
        over_info = self._over_info(int(event.x))
        self.set_class(over_info, "health-info-hover")
        self.tooltip = (
            f"💡 Click to explain {self._component}" if over_info else None
        )

    def on_leave(self, _event: Leave) -> None:
        self.set_class(False, "health-info-hover")
        self.tooltip = None

    def on_click(self, event: Click) -> None:
        if not self._over_info(int(event.x)):
            return
        event.stop()
        self.app.push_screen(HealthDetailModal(self._item))


class HealthPanel(Vertical):
    """Fieldset containing four interactive rows while retaining Overview styling."""

    def __init__(self, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.add_class("ov-block")
        self.border_title = "HEALTH"

    def compose(self) -> ComposeResult:
        for component in ("CORE", "FACTORY", "WORKERS", "ARCHIVER"):
            yield HealthInfoRow(component)

    def set_render(self, rendered: HealthBlockRender, width: int) -> None:
        self.border_title = rendered.border_title
        self.border_subtitle = rendered.border_subtitle
        rows = list(self.query(HealthInfoRow))
        for row, item in zip(rows, rendered.items):
            row_width = row.content_region.width or row.size.width or width
            row.set_item(item, int(row_width))


def _colorize_spark_row(text: str, row: int) -> str:
    """Color one complete plot row; both metrics use the same row color."""
    color = SPARK_ROW_COLORS[row]
    return f"[{color}]{text}[/]"


def _peak_value(values: tuple[float | None, ...]) -> float | None:
    """Return the peak observed value, ignoring warm-up empty slots."""
    observed = [value for value in values if value is not None]
    return max(observed) if observed else None


def _spark_peak_subtitle(frame: OverviewFrame) -> str:
    """Build the compact lower-right SPARKS peak summary."""
    sample_peak = format_sample_rate(_peak_value(frame.spark_samples))
    queue_peak = _peak_value(frame.spark_queue)
    queue_text = "—" if queue_peak is None else f"{queue_peak:0.0f}"
    return f"Peak : samples {sample_peak} · queue {queue_text}"


def render_spark_block(
    width: int,
    metrics: tuple[tuple[str, str, tuple[float | None, ...], float], ...],
) -> str:
    """Render the simu SPARKS block as balanced four-row plots."""
    if width <= 0 or not metrics:
        return ""
    separator = SPARK_SEPARATOR
    widths = _spark_column_widths(width, len(metrics))
    plots = [
        _spark_area_rows(values, column_width, scale=scale)
        for (_label, _current, values, scale), column_width in zip(metrics, widths)
    ]
    lines = [
        separator.join(
            _edge_align(label, current, column_width)
            for (label, current, _values, _scale), column_width in zip(metrics, widths)
        ).ljust(width)
    ]
    for row in range(SPARK_ROWS):
        plot_line = separator.join(
            _colorize_spark_row(plot[row].ljust(column_width), row)
            for plot, column_width in zip(plots, widths)
        )
        visible_plot_width = sum(widths) + len(separator) * (len(widths) - 1)
        lines.append(plot_line + " " * max(0, width - visible_plot_width))
    return "\n".join(lines)


def _layers(total: list[str], detail: list[str]) -> str:
    return "\n".join([*total, "", *detail])


def render_resource_row(
    inner_w: int,
    cpu: float | None,
    mem_g: float | None,
    mem_total_g: float | None,
) -> str:
    """LOCKED RESOURCES inner row. See docs/TUI/STYLES.txt."""
    cpu_sfx = " —" if cpu is None else f" {cpu:0.0f}%"
    mem_known = mem_g is not None and mem_total_g is not None and mem_total_g > 0
    mem_sfx = " —" if not mem_known else f" {mem_g:0.1f}/{mem_total_g:0.0f}G"
    inner_w = max(24, inner_w)
    dot_at = inner_w // 2
    left_w = max(8, dot_at - 2)
    right_w = max(8, inner_w - dot_at - 3)
    cpu_rail = max(4, left_w - len("CPU ") - len(cpu_sfx))
    mem_rail = max(4, right_w - len("MEM ") - len(mem_sfx))
    return (
        "CPU "
        + bar_rail(0.0 if cpu is None else cpu / 100.0, cpu_rail)
        + cpu_sfx
        + "  ·  "
        + "MEM "
        + bar_rail(0.0 if not mem_known else mem_g / mem_total_g, mem_rail)
        + mem_sfx
    )


class FieldBlock(Static):
    """Rounded fieldset; ``border_title`` is the label on the top-left edge."""

    def __init__(self, label: str, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.add_class("ov-block")
        self.border_title = label


@dataclass
class _PacmanGhost:
    row: int
    column: int
    direction: int
    role: str
    speed: int


class PacmanGame(Static):
    """A tiny borderless Pac-Man board with independent ghosts and pellet loops."""

    _PACMAN_OPEN = "ᗧ"
    _PACMAN_CLOSED = "◯"
    _GHOST_ROLES = (
        "pacman-ghost-red",
        "pacman-ghost-yellow",
        "pacman-ghost-blue",
    )
    _GHOST_SPEEDS = (1, 2, 1)
    _CAUGHT_TICKS = 7  # 7 × 0.14 s ≈ 1 second.

    def __init__(self, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self._step = 0
        self._shape = (0, 0)
        self._pellets: list[list[bool]] = []
        self._blocks: list[set[int]] = []
        self._restore_queue: deque[tuple[int, int]] = deque()
        self._path: list[tuple[int, int]] = []
        self._path_index = 0
        self._ghosts: list[_PacmanGhost] = []
        self._giant_ticks = 0
        self._caught_ticks = 0
        self._giant: tuple[int, int] | None = None
        self._giant_spawn_ticks = 0
        self._rng = random.Random(0x4A4152564953)

    def on_mount(self) -> None:
        self._paint_lane()
        self.set_interval(0.14, self._tick)

    def on_resize(self, _event: Resize) -> None:
        self._paint_lane()

    def _tick(self) -> None:
        if not self.is_on_screen or self.content_region.width <= 0 or self.content_region.height <= 0:
            return
        self._step += 1
        self._advance_game()
        self._paint_board()

    def _paint_lane(self) -> None:
        # Kept as a small compatibility shim for Textual's first resize pass.
        self._paint_board()

    def _ensure_board(self, width: int, height: int) -> None:
        shape = (max(1, width), max(1, height))
        if shape == self._shape:
            return
        track_width, board_rows = shape
        self._shape = shape
        self._pellets = [[True] * track_width for _ in range(board_rows)]
        self._blocks = [
            set(
                self._rng.sample(
                    range(track_width),
                    k=min(track_width, self._rng.randint(1, 2)),
                )
            )
            for _ in range(board_rows)
        ]
        self._restore_queue.clear()
        self._giant_ticks = 0
        self._caught_ticks = 0
        self._giant = None
        self._giant_spawn_ticks = self._next_giant_delay(track_width)
        # Sweep down, then back up. Every visited row is eaten left-to-right so
        # ᗧ always faces the direction of travel.
        row_order = list(range(board_rows)) + list(range(board_rows - 2, 0, -1))
        self._path = [
            (row, column)
            for row in row_order
            for column in range(track_width)
        ]
        self._path_index = 0
        max_ghost_column = max(0, track_width - 2)
        start_columns = (
            track_width // 4,
            track_width // 2,
            (track_width * 3) // 4,
        )
        self._ghosts = []
        for index, (role, start_column) in enumerate(
            zip(self._GHOST_ROLES, start_columns)
        ):
            row = index % board_rows
            self._ghosts.append(
                _PacmanGhost(
                    row=row,
                    column=self._first_open_ghost_column(
                        row,
                        min(max_ghost_column, start_column),
                        max_ghost_column,
                    ),
                    direction=1 if index % 2 == 0 else -1,
                    role=role,
                    speed=self._GHOST_SPEEDS[index],
                )
            )
        # Even a narrow board must not start a two-cell ghost inside a wall.
        for ghost in self._ghosts:
            self._blocks[ghost.row].difference_update((ghost.column, ghost.column + 1))
        self._blocks[0].discard(0)
        for row in range(board_rows):
            occupied = {
                column
                for ghost in self._ghosts if ghost.row == row
                for column in (ghost.column, ghost.column + 1)
            }
            if row == 0:
                occupied.add(0)
            available = list(set(range(track_width)) - occupied - self._blocks[row])
            if not self._blocks[row] and available:
                self._blocks[row].update(self._rng.sample(available, min(len(available), self._rng.randint(1, 2))))

    def _first_open_ghost_column(self, row: int, start: int, maximum: int) -> int:
        for offset in range(maximum + 1):
            column = (start + offset) % max(1, maximum + 1)
            if column not in self._blocks[row] and column + 1 not in self._blocks[row]:
                return column
        return start

    def _advance_game(self) -> None:
        if not self._path:
            return
        width, _height = self._shape
        row, column = self._path[self._path_index]
        if self._caught_ticks:
            self._move_ghosts(row, column, width)
            self._caught_ticks -= 1
            self._giant_ticks = max(0, self._giant_ticks - 1)
            self._restore_one()
            self._advance_giant_spawn((row, column))
            return

        self._consume(row, column)
        caught = self._handle_collisions(row, column)
        self._move_ghosts(row, column, width)
        caught = caught or self._handle_collisions(row, column)
        next_index = (self._path_index + 1) % len(self._path)
        if not caught:
            self._path_index = next_index
            next_row, next_column = self._path[next_index]
            self._consume(next_row, next_column)
            self._handle_collisions(next_row, next_column)
            # Only an actual row departure starts recovery, including a
            # single-row board wrapping to its first cell.
            if column == width - 1:
                self._queue_row_restore(row)
        self._restore_one()
        self._giant_ticks = max(0, self._giant_ticks - 1)
        self._advance_giant_spawn(self._path[self._path_index])

    def _consume(self, row: int, column: int) -> None:
        self._blocks[row].discard(column)
        self._pellets[row][column] = False
        if self._giant == (row, column):
            self._giant = None
            self._giant_ticks = max(12, self._shape[0] // 2)
            self._giant_spawn_ticks = self._next_giant_delay(self._shape[0])

    def _next_giant_delay(self, width: int) -> int:
        return self._rng.randint(max(6, width // 6), max(12, width // 2))

    def _advance_giant_spawn(self, occupied: tuple[int, int]) -> None:
        if self._giant is not None:
            return
        if self._giant_spawn_ticks > 0:
            self._giant_spawn_ticks -= 1
            return
        width, height = self._shape
        blocked = {occupied}
        for ghost in self._ghosts:
            blocked.add((ghost.row, ghost.column))
            blocked.add((ghost.row, ghost.column + 1))
        for row, blocks in enumerate(self._blocks):
            blocked.update((row, column) for column in blocks)
        choices = [
            (row, column)
            for row in range(height)
            for column in range(width)
            if (row, column) not in blocked
        ]
        if choices:
            self._giant = self._rng.choice(choices)
        self._giant_spawn_ticks = self._next_giant_delay(width)

    def _queue_row_restore(self, row: int) -> None:
        queued = set(self._restore_queue)
        added = False
        for column, present in enumerate(self._pellets[row]):
            if not present and (row, column) not in queued:
                self._restore_queue.append((row, column))
                added = True
        if added:
            self._spawn_breakable_blocks()

    def _spawn_breakable_blocks(self) -> None:
        """Replenish one random row during recovery, keeping at most two walls."""
        width, height = self._shape
        actor = self._path[self._path_index] if self._path else None
        blocked = {actor} if actor is not None else set()
        if self._giant is not None:
            blocked.add(self._giant)
        for ghost in self._ghosts:
            blocked.add((ghost.row, ghost.column))
            blocked.add((ghost.row, ghost.column + 1))

        candidates: list[tuple[int, list[int]]] = []
        for row in range(height):
            available = [
                column
                for column in range(width)
                if column not in self._blocks[row] and (row, column) not in blocked
            ]
            if available and len(self._blocks[row]) < 2:
                candidates.append((row, available))
        if not candidates:
            return
        row, available = self._rng.choice(candidates)
        count = min(len(available), 2 - len(self._blocks[row]), self._rng.randint(1, 2))
        self._blocks[row].update(self._rng.sample(available, k=count))

    def _restore_one(self) -> None:
        if self._restore_queue:
            row, column = self._restore_queue.popleft()
            self._pellets[row][column] = True

    def _move_ghosts(self, pacman_row: int, pacman_column: int, width: int) -> None:
        max_column = max(0, width - 2)
        height = self._shape[1]
        for ghost in self._ghosts:
            if self._step % ghost.speed:
                continue
            if width < 2:
                continue
            if self._giant_ticks and ghost.row == pacman_row:
                away = -1 if ghost.column <= pacman_column else 1
                # A flee decision must not undo a wall bounce every frame.
                lane = ghost.row * (max_column + 1) + ghost.column
                candidate = (lane + away) % (height * (max_column + 1))
                target_row, target_column = divmod(candidate, max_column + 1)
                if not self._blocks[target_row].intersection((target_column, target_column + 1)):
                    ghost.direction = away
            next_row = ghost.row
            next_column = ghost.column + ghost.direction
            if next_column < 0:
                next_row = (ghost.row - 1) % height
                next_column = max_column
            elif next_column > max_column:
                next_row = (ghost.row + 1) % height
                next_column = 0
            if (
                next_column in self._blocks[next_row]
                or next_column + 1 in self._blocks[next_row]
            ):
                ghost.direction *= -1
                continue
            ghost.row = next_row
            ghost.column = next_column

    def _handle_collisions(self, pacman_row: int, pacman_column: int) -> bool:
        if self._caught_ticks:
            return False
        for ghost in self._ghosts:
            if self._shape[0] < 2 or ghost.row != pacman_row or not (
                ghost.column <= pacman_column <= ghost.column + 1
            ):
                continue
            if self._giant_ticks:
                # Respawn only in a free pair of cells; never teleport into a wall.
                free = [
                    (row, column)
                    for row in range(self._shape[1])
                    for column in range(self._shape[0] - 1)
                    if column not in self._blocks[row]
                    and column + 1 not in self._blocks[row]
                    and not (row == pacman_row and column <= pacman_column <= column + 1)
                ]
                if free:
                    ghost.row, ghost.column = self._rng.choice(free)
            else:
                self._caught_ticks = self._CAUGHT_TICKS
                return True
        return False

    def _paint_board(self) -> None:
        width = self.content_region.width or self.size.width
        height = self.content_region.height or self.size.height
        if width <= 0 or height <= 0:
            self.update("")
            return
        self._ensure_board(width, height)
        self.update(self._render_board())

    def _render_board(self) -> str:
        track_width, board_rows = self._shape
        actor_row, pacman_position = self._path[self._path_index]
        cells: list[list[tuple[str, str] | None]] = [
            [
                ("·" if present else " ", "pacman-pellet")
                for present in pellets
            ]
            for pellets in self._pellets
        ]
        for row, blocks in enumerate(self._blocks):
            for column in blocks:
                cells[row][column] = ("▣", "pacman-block")
        if self._giant is not None:
            giant_row, giant_column = self._giant
            cells[giant_row][giant_column] = ("●", "pacman-giant")
        pacman = self._PACMAN_CLOSED if self._step % 2 else self._PACMAN_OPEN
        pacman_role = "pacman-giant" if self._giant_ticks else "pacman"
        actor_display_position = pacman_position
        actor_width = 1
        if self._caught_ticks:
            actor_width = min(2, track_width)
            pacman = ("👾" if actor_width == 2 else "◯") if self._step % 2 else " " * actor_width
            pacman_role = "pacman-caught"
            actor_display_position = min(pacman_position, track_width - actor_width)
        cells[actor_row][actor_display_position] = (pacman, pacman_role)
        if actor_width == 2:
            cells[actor_row][actor_display_position + 1] = None
        occupied = {(actor_row, column) for column in range(actor_display_position, actor_display_position + actor_width)}
        for ghost in self._ghosts:
            position = ghost.column
            if track_width < 2 or position + 1 >= track_width:
                continue
            footprint = {(ghost.row, position), (ghost.row, position + 1)}
            if footprint & occupied or self._blocks[ghost.row].intersection((position, position + 1)):
                continue
            if self._giant is not None:
                giant_row, giant_column = self._giant
                if ghost.row == giant_row and position <= giant_column <= position + 1:
                    continue
            role = "pacman-ghost-flee" if self._giant_ticks else ghost.role
            cells[ghost.row][position] = ("👻", role)
            cells[ghost.row][position + 1] = None
            occupied.update(footprint)

        rows: list[str] = []
        for row in range(board_rows):
            rows.append(
                "".join(
                    "" if cell is None else paint(cell[1], cell[0])
                    for cell in cells[row]
                )
            )
        return "\n".join(rows)


class OverviewPane(Vertical):
    def __init__(
        self,
        frame: OverviewFrame | None = None,
        *,
        on_spark_width_change: Callable[[tuple[int, int]], None] | None = None,
        **kwargs: object,
    ) -> None:
        kwargs.setdefault("id", "overview")
        super().__init__(**kwargs)
        self.frame = frame or simu_overview()
        self._on_spark_width_change = on_spark_width_change
        self._status_breathe_phase = 0.0
        self._queue_breathe_phase = 0.0
        self._worker_breathe_phase = 0.0
        self._calculator_pulse_phase = 0.0

    def compose(self) -> ComposeResult:
        yield FieldBlock("RESOURCES", id="ov-resources")
        yield FieldBlock("STATUS", id="ov-status")
        yield FieldBlock("SPARKS", id="ov-sparks")
        with Horizontal(id="ov-lower-grid"):
            with Vertical(id="ov-lower-left"):
                yield FieldBlock("SAMPLES", id="ov-samples")
                yield FieldBlock("CALCULATORS", id="ov-calcs")
                yield PacmanGame(id="ov-pacman-left")
            with Vertical(id="ov-lower-right"):
                yield HealthPanel(id="ov-health")
                yield FieldBlock("QUEUES", id="ov-queues")
                yield FieldBlock("WORKERS", id="ov-workers")
                yield PacmanGame(id="ov-pacman")

    def on_mount(self) -> None:
        self.paint()
        self.call_after_refresh(self._sync_pacman_height)
        self.set_interval(0.08, self._tick_status_breathe)

    def _sync_pacman_height(self) -> None:
        """Fill the shorter column, measuring blocks without either filler."""
        try:
            left = self.query_one("#ov-lower-left")
            right = self.query_one("#ov-lower-right")
            calculators = self.query_one("#ov-calcs")
            workers = self.query_one("#ov-workers")
            left_game = self.query_one("#ov-pacman-left", PacmanGame)
            right_game = self.query_one("#ov-pacman", PacmanGame)
        except Exception:
            return
        left_blocks_height = calculators.region.bottom - left.region.y
        right_blocks_height = workers.region.y + workers.region.height - right.region.y
        difference = left_blocks_height - right_blocks_height
        for game, desired in ((left_game, max(0, -difference)), (right_game, max(0, difference))):
            if game.region.height != desired:
                game.styles.height = desired

    def _tick_status_breathe(self) -> None:
        if self.frame.mode.lower() != "running":
            return
        self._status_breathe_phase = (self._status_breathe_phase + 0.12) % (2 * math.pi)
        self._queue_breathe_phase = (self._queue_breathe_phase + 0.09) % (2 * math.pi)
        self._worker_breathe_phase = (self._worker_breathe_phase + 0.10) % (2 * math.pi)
        self._calculator_pulse_phase = (self._calculator_pulse_phase + 0.65) % (2 * math.pi)
        try:
            status = self.query_one("#ov-status", FieldBlock)
        except Exception:
            return
        status_width = status.content_region.width or status.size.width or self.size.width or 80
        status.update(
            render_status_block(
                status_width,
                self.frame.scan,
                self.frame.mode,
                self.frame.elapsed,
                self.frame.ref,
                self.frame.redis,
                breathe_phase=self._status_breathe_phase,
            )
        )
        self._paint_queues()
        self._paint_workers()
        self._paint_calculators()

    def _paint_queues(self) -> None:
        """Refresh fixed-position queue LEDs without moving its live values."""
        try:
            queues_widget = self.query_one("#ov-queues", FieldBlock)
        except Exception:
            return
        queues_width = (
            queues_widget.content_region.width or queues_widget.size.width or 52
        )
        queues_render = render_queue_block(
            queues_width,
            self.frame.queue_block,
            outer_width=queues_widget.size.width or queues_width + 4,
            breathe_phase=self._queue_breathe_phase,
        )
        queues_widget.border_title = queues_render.border_title
        queues_widget.border_subtitle = queues_render.border_subtitle
        queues_widget.update(queues_render.body)

    def _paint_workers(self) -> None:
        """Refresh WORKERS LEDs while keeping their click-free rows stationary."""
        try:
            workers_widget = self.query_one("#ov-workers", FieldBlock)
        except Exception:
            return
        workers_width = workers_widget.content_region.width or workers_widget.size.width or 52
        workers_render = render_workers_block(
            workers_width,
            self.frame.worker_block,
            outer_width=workers_widget.size.width or workers_width + 4,
            breathe_phase=self._worker_breathe_phase,
        )
        workers_widget.border_title = workers_render.border_title
        workers_widget.border_subtitle = workers_render.border_subtitle
        workers_widget.update(workers_render.body)

    def _paint_calculators(self) -> None:
        """Fast-pulse only the fixed PackID positions currently occupied."""
        try:
            calculators_widget = self.query_one("#ov-calcs", FieldBlock)
        except Exception:
            return
        width = calculators_widget.content_region.width or calculators_widget.size.width or 52
        rendered = render_calculators_block(
            width,
            self.frame.calculator_block,
            outer_width=calculators_widget.size.width or width + 4,
            pulse_phase=self._calculator_pulse_phase,
        )
        calculators_widget.border_title = rendered.border_title
        calculators_widget.border_subtitle = rendered.border_subtitle
        calculators_widget.update(rendered.body)

    def on_resize(self) -> None:
        self.call_after_refresh(self._sync_pacman_height)
        if self._on_spark_width_change is not None:
            self._on_spark_width_change(self.spark_history_widths())

    def spark_history_widths(self) -> tuple[int, int]:
        """Return the exact sample and queue history widths for SPARKS."""
        try:
            sparks_widget = self.query_one("#ov-sparks", FieldBlock)
        except Exception:
            return (1, 1)
        block_width = sparks_widget.content_region.width
        if block_width <= 0:
            block_width = sparks_widget.size.width or self.size.width or 80
        balanced_width = _balanced_spark_width(block_width, SPARK_METRICS)
        widths = _spark_column_widths(balanced_width, SPARK_METRICS)
        if len(widths) != SPARK_METRICS:
            return (1, 1)
        return tuple(max(1, width) for width in widths)

    def spark_history_width(self) -> int:
        """Return the sample history width for compatibility with older callers."""
        return self.spark_history_widths()[0]

    def set_frame(self, frame: OverviewFrame) -> None:
        self.frame = frame
        self.paint()
        self.call_after_refresh(self._sync_pacman_height)

    def paint(self) -> None:
        try:
            self.query_one("#ov-status", FieldBlock)
        except Exception:
            return
        frame = self.frame
        status_widget = self.query_one("#ov-status", FieldBlock)
        status_width = status_widget.content_region.width or status_widget.size.width or self.size.width or 80
        status_widget.update(
            render_status_block(
                status_width,
                frame.scan,
                frame.mode,
                frame.elapsed,
                frame.ref,
                frame.redis,
                breathe_phase=self._status_breathe_phase,
            )
        )
        samples_widget = self.query_one("#ov-samples", FieldBlock)
        samples_width = (
            samples_widget.content_region.width
            or samples_widget.size.width
            or 52
        )
        samples_render = render_samples_block(
            samples_width,
            SamplesBlockData(
                completed=frame.done,
                running=frame.running,
                failed=frame.failed,
                rate=frame.rate,
                sampler=frame.sampler,
            ),
        )
        samples_widget.border_title = samples_render.border_title
        samples_widget.border_subtitle = samples_render.border_subtitle
        samples_widget.update(samples_render.body)
        self._paint_queues()
        self._paint_workers()
        self._paint_calculators()
        sparks_widget = self.query_one("#ov-sparks", FieldBlock)
        block_width = sparks_widget.content_region.width
        if block_width <= 0:
            block_width = sparks_widget.size.width or self.size.width or 80
        spark_width = _balanced_spark_width(
            max(1, block_width),
            SPARK_METRICS,
        )
        has_samples = any(value is not None for value in frame.spark_samples)
        has_queue = any(value is not None for value in frame.spark_queue)
        sparks_widget.update(
            render_spark_block(
                spark_width,
                (
                    (
                        "SAMPLES",
                        frame.rate if has_samples else "—",
                        frame.spark_samples,
                        20.0 / 60.0,
                    ),
                    (
                        "TASK QUEUE",
                        str(frame.task_q) if has_queue else "—",
                        frame.spark_queue,
                        32.0,
                    ),
                ),
            )
        )
        sparks_widget.border_subtitle = _spark_peak_subtitle(frame)
        res_widget = self.query_one("#ov-resources", FieldBlock)
        inner_w = max(24, res_widget.size.width or 80)
        if frame.resources_available:
            res_widget.update(
                render_resource_row(
                    inner_w,
                    frame.cpu,
                    frame.mem_g,
                    frame.mem_total_g,
                )
            )
            fds = "—" if frame.fds_limit <= 0 else f"{frame.fds}/{frame.fds_limit}"
            res_widget.border_subtitle = (
                f"{frame.cpu:0.0f}% · {frame.mem_g:0.1f}G · {fds}"
            )
        else:
            res_widget.update(render_resource_row(inner_w, None, None, None))
            res_widget.border_subtitle = "host metrics unavailable · fds —"
        health_widget = self.query_one("#ov-health", HealthPanel)
        health_width = (
            health_widget.content_region.width or health_widget.size.width or 52
        )
        health_render = render_health_block(health_width, frame)
        health_widget.set_render(health_render, health_width)


class LiveOverviewPane(Vertical):
    """Honest live overview from the Collector, never the simulated fixture."""

    def __init__(self, **kwargs: object) -> None:
        kwargs.setdefault("id", "overview")
        super().__init__(**kwargs)
        self.frame = CollectorFrame(page="overview")
        self._queue_breathe_phase = 0.0

    def compose(self) -> ComposeResult:
        with Horizontal(id="ov-row1"):
            yield FieldBlock("SAMPLES", id="ov-samples")
            yield FieldBlock("QUEUES", id="ov-queues")
        yield Static(id="live-overview-body")

    def on_mount(self) -> None:
        self.paint()
        self.set_interval(0.08, self._tick_queue_breathe)

    def _tick_queue_breathe(self) -> None:
        """Animate only the fixed queue signal lamps between collector ticks."""
        self._queue_breathe_phase = (self._queue_breathe_phase + 0.09) % (2 * math.pi)
        self._paint_queues()

    def _paint_queues(self) -> None:
        try:
            queues_widget = self.query_one("#ov-queues", FieldBlock)
        except Exception:
            return
        queues_width = (
            queues_widget.content_region.width or queues_widget.size.width or 52
        )
        queues_render = render_queue_block(
            queues_width,
            self.frame.queue_block,
            outer_width=queues_widget.size.width or queues_width + 4,
            breathe_phase=self._queue_breathe_phase,
        )
        queues_widget.border_title = queues_render.border_title
        queues_widget.border_subtitle = queues_render.border_subtitle
        queues_widget.update(queues_render.body)

    def set_frame(self, frame: CollectorFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            body = self.query_one("#live-overview-body", Static)
            samples_widget = self.query_one("#ov-samples", FieldBlock)
        except Exception:
            return
        frame = self.frame
        stats = frame.sample_stats
        rows = frame.occupancy
        workers = frame.workers
        busy_workers = sum(1 for row in workers if str(row.get("status") or "") == "busy")
        alive_workers = sum(1 for row in workers if row.get("alive") is not False)
        stale = "  [STALE]" if frame.stale else ""
        samples_width = (
            samples_widget.content_region.width
            or samples_widget.size.width
            or 52
        )
        samples_render = render_samples_block(
            samples_width,
            SamplesBlockData(
                completed=int(stats.get("completed", 0) or 0),
                running=int(stats.get("running", 0) or 0),
                failed=int(stats.get("failed", 0) or 0),
                rate=frame.sample_rate,
                sampler=frame.sampler_display,
            ),
        )
        samples_widget.border_title = samples_render.border_title
        samples_widget.border_subtitle = samples_render.border_subtitle
        samples_widget.update(samples_render.body)
        self._paint_queues()
        lines = [
            f"OVERVIEW{stale}  live scan snapshot",
            "",
            "WORKERS",
            f"{alive_workers} / {len(workers)} visible    busy  {busy_workers}",
            "",
            "CALCULATORS  BUSY/FREE/TOTAL",
        ]
        if rows:
            for name, row in rows.items():
                lines.append(
                    f"{name:<20.20} {int(row.get('busy', 0) or 0):>3} / "
                    f"{int(row.get('free', 0) or 0):>3} / {int(row.get('slots', 0) or 0):>3}"
                )
        else:
            lines.append("No calculator pool metadata is available for this scan.")
        lines.extend(
            [
                "",
                "CONTROL PLANE",
                f"core  {str(frame.proc_core.get('status') or 'unknown'):<10}  "
                f"archiver  {str(frame.proc_archiver.get('status') or 'unknown'):<10}  "
                f"redis  {str(frame.proc_redis.get('status') or 'unknown')}",
            ]
        )
        if frame.error:
            lines.extend(["", f"last error  {frame.error}"])
        body.update("\n".join(lines))


class ComingPane(Vertical):
    def __init__(self, title: str, *, id: str) -> None:
        super().__init__(id=id)
        self._title = title

    def compose(self) -> ComposeResult:
        yield Static(
            paint("dim", f"{self._title} — next"),
            id=f"{self.id}-body",
        )


__all__ = [
    "ComingPane",
    "FieldBlock",
    "LiveOverviewPane",
    "OverviewPane",
    "render_adaptive_bridson_samples",
    "render_health_block",
    "render_queue_block",
    "render_spark_block",
]
