"""Overview pane: rounded fieldset blocks (总 / 分) inside the page frame."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.containers import Horizontal, Vertical
from textual.widgets import Static

from jarvishep2.monitor.bars import bar_rail
from jarvishep2.monitor.simu import OverviewFrame, simu_overview
from jarvishep2.monitor.styles import paint

SPARK_CHARS = "▁▂▃▄▅▆▇█"


def _spark(values: tuple[int, ...], width: int) -> str:
    if width <= 0 or not values:
        return ""
    seq = (list(values) * ((width // len(values)) + 1))[:width]
    lo, hi = min(seq), max(seq)
    span = max(1, hi - lo)
    return "".join(
        SPARK_CHARS[int(round((v - lo) / span * (len(SPARK_CHARS) - 1)))]
        for v in seq
    )


def _layers(total: list[str], detail: list[str]) -> str:
    return "\n".join([*total, "", *detail])


def render_resource_row(
    inner_w: int,
    cpu: float,
    mem_g: float,
    mem_total_g: float,
) -> str:
    """LOCKED RESOURCES inner row. See docs/TUI/STYLES.txt."""
    cpu_sfx = f" {cpu:0.0f}%"
    mem_sfx = f" {mem_g:0.1f}/{mem_total_g:0.0f}G"
    inner_w = max(24, inner_w)
    dot_at = inner_w // 2
    left_w = max(8, dot_at - 2)
    right_w = max(8, inner_w - dot_at - 3)
    cpu_rail = max(4, left_w - len("CPU ") - len(cpu_sfx))
    mem_rail = max(4, right_w - len("MEM ") - len(mem_sfx))
    return (
        "CPU "
        + bar_rail(cpu / 100.0, cpu_rail)
        + cpu_sfx
        + "  ·  "
        + "MEM "
        + bar_rail(mem_g / mem_total_g, mem_rail)
        + mem_sfx
    )


class FieldBlock(Static):
    """Rounded fieldset; ``border_title`` is the label on the top-left edge."""

    def __init__(self, label: str, **kwargs: object) -> None:
        super().__init__(**kwargs)
        self.add_class("ov-block")
        self.border_title = label


class OverviewPane(Vertical):
    def __init__(self, frame: OverviewFrame | None = None, **kwargs: object) -> None:
        kwargs.setdefault("id", "overview")
        super().__init__(**kwargs)
        self.frame = frame or simu_overview()

    def compose(self) -> ComposeResult:
        yield FieldBlock("RESOURCES", id="ov-resources")
        yield FieldBlock("STATUS", id="ov-status")
        with Horizontal(id="ov-row1"):
            yield FieldBlock("SAMPLES", id="ov-samples")
            yield FieldBlock("QUEUES", id="ov-queues")
        with Horizontal(id="ov-row2"):
            yield FieldBlock("WORKERS", id="ov-workers")
            yield FieldBlock("CALCULATORS", id="ov-calcs")
        yield FieldBlock("SPARKS", id="ov-sparks")
        yield FieldBlock("LIVE", id="ov-live")

    def on_mount(self) -> None:
        self.paint()

    def on_resize(self) -> None:
        self.paint()

    def set_frame(self, frame: OverviewFrame) -> None:
        self.frame = frame
        self.paint()

    def paint(self) -> None:
        try:
            self.query_one("#ov-status", FieldBlock)
        except Exception:
            return
        frame = self.frame
        compact = (self.size.width or 80) < 100
        bar_w = 16 if compact else 28
        spark_w = max(20, (self.size.width or 80) - 18)
        remain = max(0, frame.target - frame.done - frame.running)
        done_frac = frame.done / max(1, frame.target)
        pending = frame.task_q + frame.archive_q + frame.feedback_q
        idle = max(0, frame.workers_alive - frame.busy)
        busy_packs = sum(busy for _name, busy, _total in frame.calculators)
        free_packs = sum(total - busy for _name, busy, total in frame.calculators)
        total_packs = busy_packs + free_packs
        live_n = 4

        self.query_one("#ov-status", FieldBlock).update(
            _layers(
                [f"{frame.scan}    ● {frame.mode}    {frame.elapsed}"],
                [f"REF  {frame.ref}     redis  {frame.redis}     hz  {frame.hz}"],
            )
        )
        self.query_one("#ov-samples", FieldBlock).update(
            _layers(
                [
                    f"{frame.done} / {frame.target}  {bar_rail(done_frac, bar_w)}  "
                    f"{done_frac * 100:0.0f}%"
                ],
                [
                    f"done {frame.done}    run {frame.running}    "
                    f"fail {frame.failed}    remain {remain}",
                    f"rate {frame.rate}     eta {frame.eta}     avg {frame.avg}",
                ],
            )
        )
        self.query_one("#ov-queues", FieldBlock).update(
            _layers(
                [
                    f"pending  {pending}     "
                    f"task {frame.task_q} + archive {frame.archive_q} "
                    f"+ feedback {frame.feedback_q}"
                ],
                [
                    f"task        {frame.task_q:>4}  {bar_rail(frame.task_q / 32.0, bar_w)}",
                    f"archive     {frame.archive_q:>4}  {bar_rail(frame.archive_q / 32.0, bar_w)}",
                    f"feedback    {frame.feedback_q:>4}  {bar_rail(0.02, bar_w)}",
                    "chain-0  0      chain-1  0      chain-2  0",
                ],
            )
        )
        self.query_one("#ov-workers", FieldBlock).update(
            _layers(
                [f"{frame.workers_alive} / {frame.workers_total} alive     {frame.stale} stale"],
                [
                    f"busy {frame.busy}     idle {idle}     stale {frame.stale}",
                    f"cpu  {bar_rail(frame.cpu / 100.0, bar_w)}  {frame.cpu:0.0f}%",
                    f"mem  {bar_rail(frame.mem_g / frame.mem_total_g, bar_w)}  "
                    f"{frame.mem_g:0.1f} / {frame.mem_total_g:0.0f} G",
                ],
            )
        )
        calc_detail = [
            f"{name:<12} {busy:>2}/{total:<2}  {bar_rail(busy / max(1, total), 16)}"
            for name, busy, total in frame.calculators
        ]
        self.query_one("#ov-calcs", FieldBlock).update(
            _layers(
                [f"{busy_packs} / {total_packs} busy     {free_packs} free"],
                calc_detail,
            )
        )
        self.query_one("#ov-sparks", FieldBlock).update(
            _layers(
                [f"{frame.rate}     queue {frame.task_q}     cpu {frame.cpu:0.0f}%"],
                [
                    f"samples/min  {_spark(frame.spark_samples, spark_w)}",
                    f"task queue   {_spark(frame.spark_queue, spark_w)}",
                    f"host cpu     {_spark(frame.spark_cpu, spark_w)}",
                ],
            )
        )
        res_widget = self.query_one("#ov-resources", FieldBlock)
        inner_w = max(24, res_widget.size.width or 80)
        res_widget.update(
            render_resource_row(
                inner_w,
                frame.cpu,
                frame.mem_g,
                frame.mem_total_g,
            )
        )
        res_widget.border_subtitle = (
            f"{frame.cpu:0.0f}% · {frame.mem_g:0.1f}G · "
            f"{frame.fds}/{frame.fds_limit}"
        )
        dots = (
            paint("live", "core ●")
            + "   "
            + paint("live", "archiver ●")
            + "   "
            + paint("live", "redis ●")
            + "   "
            + paint("live", "lock ●")
        )
        self.query_one("#ov-live", FieldBlock).update(
            _layers(
                [f"{live_n} / 4 up     {dots}"],
                [
                    f"method  {frame.method}     run  {frame.run_id}",
                    f"host    {frame.host}      started  {frame.started}",
                    f"workers {frame.workers_alive}/{frame.workers_total}               "
                    f"stale {frame.stale}",
                ],
            )
        )


class ComingPane(Vertical):
    def __init__(self, title: str, *, id: str) -> None:
        super().__init__(id=id)
        self._title = title

    def compose(self) -> ComposeResult:
        yield Static(
            paint("dim", f"{self._title} — next"),
            id=f"{self.id}-body",
        )


__all__ = ["ComingPane", "FieldBlock", "OverviewPane"]
