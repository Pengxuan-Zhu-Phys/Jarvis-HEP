"""Queue-flow snapshots and the compact Overview QUEUES renderer."""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping
from dataclasses import dataclass, field
import math

from rich.cells import cell_len

from jarvishep2.monitor.styles import paint


_LED_BACKGROUND_RGB = (16, 18, 22)
_LED_COLORS = {
    "green": (53, 201, 138),
    "amber": (246, 211, 63),
    "red": (239, 107, 115),
    "idle": (74, 81, 96),
}
_LED_PHASE_OFFSETS = (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)
_SHARDED_FEEDBACK_SAMPLERS = frozenset({"MCMC", "ToyMCMC", "AMMCMC", "DRAM"})
_SHARED_FEEDBACK_SAMPLERS = frozenset(
    {"EnsembleMCMC", "DEMCMC", "PTMCMC", "PTEnsemble"}
)
_NON_FEEDBACK_SAMPLERS = frozenset(
    {"Random", "Grid", "CSV", "Bridson", "AdaptiveBridson", "Dynesty", "MultiNest"}
)


@dataclass(frozen=True)
class QueueDirection:
    """A bounded Monitor-local reading of one queue depth's movement."""

    state: str = "warming"
    rate_per_min: float | None = None

    def label(self) -> str:
        if self.state == "rising" and self.rate_per_min is not None:
            return f"↑ {max(1, round(abs(self.rate_per_min)))} / min"
        if self.state == "draining" and self.rate_per_min is not None:
            return f"↓ {max(1, round(abs(self.rate_per_min)))} / min"
        if self.state == "flat":
            return "→ flat"
        if self.state == "stale":
            return "— stale"
        if self.state == "unknown":
            return "— unknown"
        return "— warming"


class QueueDirectionWindow:
    """Derive stable direction labels from a short local queue-depth window."""

    def __init__(self, *, window_seconds: float = 8.0, jitter_items: int = 1) -> None:
        self.window_seconds = max(0.1, float(window_seconds))
        self.jitter_items = max(0, int(jitter_items))
        self._observations: dict[str, deque[tuple[float, int]]] = {}

    def reset(self, *lanes: str) -> None:
        """Forget every lane, or only the supplied lanes after a mode change."""
        if not lanes:
            self._observations.clear()
            return
        for lane in lanes:
            self._observations.pop(str(lane), None)

    def observe(
        self,
        lane: str,
        depth: int | None,
        *,
        observed_at: float,
        stale: bool = False,
    ) -> QueueDirection:
        """Record one depth and return its non-flickering local direction."""
        if stale:
            return QueueDirection("stale")
        if depth is None:
            return QueueDirection("unknown")
        name = str(lane)
        now = float(observed_at)
        history = self._observations.setdefault(name, deque())
        history.append((now, max(0, int(depth))))
        cutoff = now - self.window_seconds
        while len(history) > 1 and history[0][0] < cutoff:
            history.popleft()
        if len(history) < 2:
            return QueueDirection("warming")
        started_at, started_depth = history[0]
        elapsed = max(0.0, now - started_at)
        if elapsed <= 0:
            return QueueDirection("warming")
        delta = int(depth) - started_depth
        if abs(delta) <= self.jitter_items:
            return QueueDirection("flat", 0.0)
        rate = delta / elapsed * 60.0
        return QueueDirection("rising" if delta > 0 else "draining", rate)


@dataclass(frozen=True)
class QueueBlockData:
    """Bounded queue snapshot consumed by the renderer, never raw payloads."""

    task_depth: int | None
    archive_depth: int | None
    feedback_depth: int | None = None
    feedback_mode: str = "inactive"
    feedback_shards: Mapping[str, int] = field(default_factory=dict)
    task_direction: QueueDirection = field(default_factory=QueueDirection)
    archive_direction: QueueDirection = field(default_factory=QueueDirection)
    feedback_direction: QueueDirection = field(default_factory=QueueDirection)


@dataclass(frozen=True)
class QueueBlockRender:
    border_title: str
    body: str
    border_subtitle: str


def feedback_queue_plan(
    method: str,
    metrics: Mapping[str, object] | None = None,
) -> tuple[str, tuple[int, ...]]:
    """Return the only feedback lists Monitor may inspect for this sampler.

    Independent-chain shards are listed by already admitted runtime metadata or
    Core telemetry. They are never discovered from the Redis keyspace.
    """
    name = str(method or "Unknown").strip()
    if name in _SHARDED_FEEDBACK_SAMPLERS:
        values = dict(metrics or {})
        raw_count = values.get("chains", values.get("num_chains"))
        try:
            count = max(1, int(raw_count))
        except (TypeError, ValueError, OverflowError):
            return "unknown", ()
        return "sharded", tuple(range(count))
    if name in _SHARED_FEEDBACK_SAMPLERS:
        return "shared", ()
    if name in _NON_FEEDBACK_SAMPLERS:
        return "inactive", ()
    return "unknown", ()


def _clip_cells(text: str, width: int) -> str:
    """Clip plain terminal text by cells, retaining an ellipsis when possible."""
    if width <= 0:
        return ""
    if cell_len(text) <= width:
        return text
    if width == 1:
        return "…"
    chars: list[str] = []
    used = 0
    for char in text:
        char_width = cell_len(char)
        if used + char_width > width - 1:
            break
        chars.append(char)
        used += char_width
    return "".join(chars) + "…"


def _edge_align(left: str, right: str, width: int) -> str:
    """Keep the metric at the right edge while clipping only the left side."""
    if width <= 0:
        return ""
    right = _clip_cells(str(right), width)
    left_width = max(0, width - cell_len(right) - (1 if right else 0))
    left = _clip_cells(str(left), left_width if right else width)
    return left + " " * max(0, width - cell_len(left) - cell_len(right)) + right


def _signal_kind(
    direction: QueueDirection,
    depth: int | None,
    *,
    inactive: bool = False,
) -> str:
    """Map a local flow observation to a traffic-light signal, not HEALTH."""
    if inactive:
        return "idle"
    if direction.state in {"stale", "unknown", "rising"}:
        return "red"
    if direction.state == "draining" or (direction.state == "flat" and not depth):
        return "green"
    return "amber"


def _signal_dot(kind: str, phase: float) -> str:
    """Return a fixed-width, softly breathing traffic-light indicator."""
    rgb = _LED_COLORS.get(kind, _LED_COLORS["idle"])
    if kind == "idle":
        amount = 1.0
    else:
        amount = 0.58 + 0.42 * (1.0 - math.cos(phase)) / 2.0
    mixed = tuple(
        round(_LED_BACKGROUND_RGB[index] + amount * (rgb[index] - _LED_BACKGROUND_RGB[index]))
        for index in range(3)
    )
    return f"[#{mixed[0]:02x}{mixed[1]:02x}{mixed[2]:02x}]●[/]"


def _metric_row(
    label: str,
    depth: int | None,
    metric: str,
    direction: QueueDirection,
    width: int,
    *,
    phase: float,
    inactive: bool = False,
) -> str:
    """Render one fixed-position LED, label, and right-edge primary metric."""
    if width <= 0:
        return ""
    # A half-width Overview column can be 31 cells.  Preserve both endpoints of
    # each path there by collapsing decorative label spacing before clipping.
    if width < 40:
        label = " ".join(label.split()).replace(" → ", "→")
    plain_metric = _clip_cells(metric, width)
    label_width = max(0, width - cell_len(plain_metric) - 1)
    plain_label = _clip_cells(f" {label}", label_width)
    gap = " " * max(0, width - 1 - cell_len(plain_label) - cell_len(plain_metric))
    led = _signal_dot(
        _signal_kind(direction, depth, inactive=inactive),
        phase,
    )
    if inactive or not plain_metric:
        rendered_metric = plain_metric
    elif cell_len(plain_metric) == cell_len(metric):
        depth_text = "—" if depth is None else str(max(0, int(depth)))
        rendered_metric = f"[bold #e6e8eb]{depth_text}[/]  {direction.label()}"
    else:
        rendered_metric = plain_metric
    return led + paint("panel-title", plain_label) + gap + rendered_metric


def _border_title(outer_width: int) -> str:
    base, right = "QUEUES", "FLOW"
    available = max(len(base), int(outer_width) - 6)
    gap = max(1, available - len(base) - len(right) - 2)
    return f"{base} {'─' * gap} {right}"


def _feedback_explanation(data: QueueBlockData) -> str:
    mode = str(data.feedback_mode or "unknown").lower()
    if mode == "sharded":
        shards = {str(key): max(0, int(value)) for key, value in data.feedback_shards.items()}
        if not shards:
            return "  control return · shard metadata unavailable"
        active = sum(1 for value in shards.values() if value > 0)
        return f"  control return · {active}/{len(shards)} shards waiting"
    if mode == "shared":
        return "  control return · shared feedback list"
    if mode == "inactive":
        return "  control return · feedback samplers only"
    return "  control return · unavailable"


def _feedback_metric(data: QueueBlockData) -> str:
    mode = str(data.feedback_mode or "unknown").lower()
    if mode == "inactive":
        return "—  inactive"
    if data.feedback_depth is None:
        return "—  unknown"
    return f"{max(0, int(data.feedback_depth))}  {data.feedback_direction.label()}"


def render_queue_block(
    width: int,
    data: QueueBlockData,
    *,
    outer_width: int | None = None,
    breathe_phase: float = 0.0,
) -> QueueBlockRender:
    """Render the eight-row Core hand-off view without fake capacities or totals."""
    width = max(1, int(width))
    task_depth = "—" if data.task_depth is None else str(max(0, int(data.task_depth)))
    archive_depth = "—" if data.archive_depth is None else str(max(0, int(data.archive_depth)))
    rows = (
        _metric_row(
            "TASK     → WORKERS",
            data.task_depth,
            f"{task_depth}  {data.task_direction.label()}",
            data.task_direction,
            width,
            phase=breathe_phase + _LED_PHASE_OFFSETS[0],
        ),
        paint(
            "dim",
            _edge_align("  dispatch backlog · unclaimed samples", "", width),
        ),
        _metric_row(
            "ARCHIVE  → ARCHIVER",
            data.archive_depth,
            f"{archive_depth}  {data.archive_direction.label()}",
            data.archive_direction,
            width,
            phase=breathe_phase + _LED_PHASE_OFFSETS[1],
        ),
        paint(
            "dim",
            _edge_align("  durability backlog · completed records", "", width),
        ),
        _metric_row(
            "FEEDBACK → SAMPLER",
            data.feedback_depth,
            _feedback_metric(data),
            data.feedback_direction,
            width,
            phase=breathe_phase + _LED_PHASE_OFFSETS[2],
            inactive=str(data.feedback_mode or "").lower() == "inactive",
        ),
        paint("dim", _edge_align(_feedback_explanation(data), "", width)),
    )
    return QueueBlockRender(
        border_title=_border_title(outer_width or width + 4),
        body="\n".join(rows),
        border_subtitle="depths are separate paths",
    )


__all__ = [
    "QueueBlockData",
    "QueueBlockRender",
    "QueueDirection",
    "QueueDirectionWindow",
    "feedback_queue_plan",
    "render_queue_block",
]
