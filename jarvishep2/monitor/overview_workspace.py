"""Overview-only Monitor workspace used by the initial release."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.events import Resize
from textual.screen import Screen
from textual.widgets import ContentSwitcher

from jarvishep2.monitor.chrome import MonitorHint, MonitorTopbar, TabBar
from jarvishep2.monitor.collector import Collector, CollectorFrame
from jarvishep2.monitor.hints import OVERVIEW_ONLY_KEYS
from jarvishep2.monitor.live_overview import LiveOverviewProjector
from jarvishep2.monitor.overview import OverviewPane
from jarvishep2.monitor.scans import ScanChoice, ScanExitWatch
from jarvishep2.monitor.simu import SIMU_INTERVAL_SEC, SimuEngine

_OVERVIEW_PAGES = (("overview", "Overview"),)


class OverviewReleaseScreen(Screen[None]):
    """Monitor v1 workspace exposing only the completed Overview."""

    BINDINGS = [
        Binding("q", "quit_monitor", "Quit", show=True),
        Binding("escape", "back", "Chooser", show=True),
        Binding("r", "refresh", "Refresh", show=False),
        Binding("m", "next_simu_sampler", show=False),
    ]

    def __init__(
        self,
        choice: ScanChoice,
        *,
        collector: Collector | None = None,
    ) -> None:
        super().__init__()
        self.choice = choice
        self._exit_watch = ScanExitWatch(choice)
        self._leaving = False
        self.engine = SimuEngine() if choice.simulated else None
        # Simu must never hold a live Redis collector (FORB-15).
        self.collector = None if choice.simulated else collector
        self.live = (
            None
            if choice.simulated
            else LiveOverviewProjector(
                choice,
                redis_endpoint=(
                    collector.redis_endpoint if collector is not None else "unavailable"
                ),
            )
        )

    def compose(self) -> ComposeResult:
        with Vertical(id="workspace"):
            yield MonitorTopbar()
            yield TabBar(_OVERVIEW_PAGES)
            with Vertical(id="page-frame"):
                with ContentSwitcher(initial="overview", id="pages"):
                    yield OverviewPane(
                        id="overview",
                        on_spark_width_change=self._sync_history,
                    )
            yield MonitorHint(OVERVIEW_ONLY_KEYS)

    def on_mount(self) -> None:
        self.query_one(TabBar).set_active("overview")
        if not self.choice.simulated:
            self.set_interval(0.5, self._check_scan_exit)
        if self.engine is not None:
            pane = self.query_one("#overview", OverviewPane)
            pane.set_frame(self.engine.tick(pane.spark_history_widths()))
            self.set_interval(SIMU_INTERVAL_SEC, self._tick)
        elif self.collector is not None:
            frame = self.collector.tick("overview")
            self._set_live_frame(frame)
            self.set_interval(0.5, self._tick)
        else:
            self._set_live_frame(
                CollectorFrame(
                    page="overview",
                    stale=True,
                    error="Redis connection unavailable",
                )
            )

    def show_page(self, slug: str) -> None:
        if slug != "overview":
            return
        self.query_one("#pages", ContentSwitcher).current = "overview"
        self.query_one(TabBar).set_active("overview")

    def action_refresh(self) -> None:
        self._tick()

    def action_next_simu_sampler(self) -> None:
        if self.engine is None:
            return
        self.engine.cycle_sampler()
        self.query_one("#overview", OverviewPane).set_frame(
            self.engine.current_frame()
        )

    def action_back(self) -> None:
        if self._leaving:
            return
        self._leaving = True
        self._shutdown_collector()
        self.dismiss(None)

    def _check_scan_exit(self) -> None:
        if self._leaving or not self._exit_watch.has_exited():
            return
        self._leaving = True
        self._shutdown_collector()
        # Schedule navigation after the exit-watch callback returns.
        self.app.call_later(self._return_to_chooser)

    def _return_to_chooser(self) -> None:
        """Return without awaiting a screen transition from a timer callback.

        Textual changes the screen stack synchronously, but completes the visual
        removal asynchronously.  Awaiting that removal here used to leave the
        event loop suspended when a scan ended underneath an open help modal.
        Pop one overlay, then re-schedule ourselves for the next UI turn.
        """
        # HEALTH explanations may be open when Core exits. Remove overlays
        # first so dismiss() returns to the existing chooser below this screen.
        if self not in self.app.screen_stack:
            return
        if self.app.screen is not self:
            self.app.pop_screen()
            self.app.call_later(self._return_to_chooser)
            return
        self.dismiss(None)

    def action_quit_monitor(self) -> None:
        self._shutdown_collector()
        self.app.exit()

    def on_unmount(self) -> None:
        self._shutdown_collector()

    def on_resize(self, _event: Resize) -> None:
        try:
            pane = self.query_one("#overview", OverviewPane)
        except Exception:
            return
        self._sync_history(pane.spark_history_widths())

    def _tick(self) -> None:
        if self._leaving:
            return
        if self.engine is not None:
            pane = self.query_one("#overview", OverviewPane)
            pane.set_frame(self.engine.tick(pane.spark_history_widths()))
        elif self.collector is not None:
            frame = self.collector.tick("overview")
            self._set_live_frame(frame)

    def _sync_history(self, widths: tuple[int, int]) -> None:
        try:
            pane = self.query_one("#overview", OverviewPane)
        except Exception:
            return
        if self.engine is not None:
            self.engine.set_history_widths(*widths)
            pane.set_frame(self.engine.current_frame())
        elif self.live is not None:
            self.live.set_history_widths(*widths)

    def _set_live_frame(self, frame: CollectorFrame) -> None:
        if self.live is None:
            return
        pane = self.query_one("#overview", OverviewPane)
        self.live.set_history_widths(*pane.spark_history_widths())
        pane.set_frame(self.live.project(frame))

    def _shutdown_collector(self) -> None:
        collector = self.collector
        self.collector = None
        if collector is not None:
            collector.close()


__all__ = ["OverviewReleaseScreen"]
