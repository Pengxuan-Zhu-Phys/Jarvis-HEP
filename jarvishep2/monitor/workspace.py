"""Future full Monitor workspace with all seven operational pages."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.events import Key, Resize
from textual.screen import Screen
from textual.widgets import ContentSwitcher

from jarvishep2.monitor.chrome import PAGES, MonitorHint, MonitorTopbar, TabBar
from jarvishep2.monitor.collector import Collector
from jarvishep2.monitor.hints import OVERVIEW_KEYS
from jarvishep2.monitor.inflight import CalculatorsPane, SamplesPane
from jarvishep2.monitor.overview import LiveOverviewPane, OverviewPane
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.simu import SIMU_INTERVAL_SEC, SimuEngine
from jarvishep2.monitor.status_pages import FactoryPane, HostPane, SamplerPane, WorkersPane


class WorkspaceScreen(Screen[None]):
    BINDINGS = [
        Binding("q", "quit_monitor", "Quit", show=True),
        Binding("escape", "back", "Chooser", show=True),
        Binding("r", "refresh", "Refresh", show=False),
        Binding("tab", "next_page", show=False, priority=True),
        Binding("shift+tab", "prev_page", show=False, priority=True),
        Binding("left", "prev_page", show=False),
        Binding("right", "next_page", show=False),
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
        self.engine = SimuEngine() if choice.simulated else None
        # Simu must never hold a collector (FORB-15).
        self.collector = None if choice.simulated else collector
        self._page = "overview"

    def compose(self) -> ComposeResult:
        with Vertical(id="workspace"):
            yield MonitorTopbar()
            yield TabBar()
            with Vertical(id="page-frame"):
                with ContentSwitcher(initial="overview", id="pages"):
                    if self.engine is not None:
                        yield OverviewPane(
                            id="overview",
                            on_spark_width_change=self._sync_simu_history,
                        )
                    else:
                        yield LiveOverviewPane()
                    yield WorkersPane()
                    yield FactoryPane()
                    yield SamplerPane()
                    yield CalculatorsPane()
                    yield SamplesPane()
                    yield HostPane()
            yield MonitorHint()

    def on_mount(self) -> None:
        self.query_one(TabBar).set_active("overview")
        if self.engine is not None:
            pane = self.query_one("#overview", OverviewPane)
            pane.set_frame(self.engine.tick(pane.spark_history_widths()))
            self.set_interval(SIMU_INTERVAL_SEC, self._tick)
        elif self.collector is not None:
            self.collector.tick(self._page)
            self.set_interval(0.5, self._tick)

    def show_page(self, slug: str) -> None:
        slugs = {item[0] for item in PAGES}
        if slug not in slugs:
            return
        self._page = slug
        self.query_one("#pages", ContentSwitcher).current = slug
        self.query_one(TabBar).set_active(slug)
        self.query_one(MonitorHint).set_keys(OVERVIEW_KEYS)
        self._tick()

    def action_refresh(self) -> None:
        self._tick()

    def action_next_page(self) -> None:
        self._cycle(1)

    def action_prev_page(self) -> None:
        self._cycle(-1)

    def action_next_simu_sampler(self) -> None:
        """Hidden visual-QA binding: cycle all sampler-specific SAMPLES blocks."""
        if self.engine is None:
            return
        self.engine.cycle_sampler()
        try:
            pane = self.query_one("#overview", OverviewPane)
        except Exception:
            return
        pane.set_frame(self.engine.current_frame())

    def action_back(self) -> None:
        self._shutdown_collector()
        self.dismiss(None)

    def action_quit_monitor(self) -> None:
        self._shutdown_collector()
        self.app.exit()

    def on_unmount(self) -> None:
        self._shutdown_collector()

    def on_key(self, event: Key) -> None:
        if event.key in {"1", "2", "3", "4", "5", "6", "7"}:
            event.stop()
            event.prevent_default()
            self.show_page(PAGES[int(event.key) - 1][0])

    def on_resize(self, _event: Resize) -> None:
        if self.engine is None:
            return
        try:
            pane = self.query_one("#overview", OverviewPane)
        except Exception:
            return
        self._sync_simu_history(pane.spark_history_widths())

    def _sync_simu_history(self, widths: tuple[int, int]) -> None:
        if self.engine is None:
            return
        try:
            pane = self.query_one("#overview", OverviewPane)
        except Exception:
            return
        self.engine.set_history_widths(*widths)
        pane.set_frame(self.engine.current_frame())

    def _cycle(self, delta: int) -> None:
        slugs = [item[0] for item in PAGES]
        index = slugs.index(self._page) if self._page in slugs else 0
        self.show_page(slugs[(index + delta) % len(slugs)])

    def _tick(self) -> None:
        if self.engine is not None:
            pane = self.query_one("#overview", OverviewPane)
            frame = self.engine.tick(pane.spark_history_widths())
            pane.set_frame(frame)
            return
        if self.collector is not None:
            frame = self.collector.tick(self._page)
            if self._page == "overview":
                self.query_one("#overview", LiveOverviewPane).set_frame(frame)
            elif self._page == "workers":
                self.query_one("#workers", WorkersPane).set_frame(frame)
            elif self._page == "factory":
                self.query_one("#factory", FactoryPane).set_frame(frame)
            elif self._page == "sampler":
                self.query_one("#sampler", SamplerPane).set_frame(frame)
            elif self._page == "calculators":
                self.query_one("#calculators", CalculatorsPane).set_frame(frame)
            elif self._page == "samples":
                self.query_one("#samples", SamplesPane).set_frame(frame)
            elif self._page == "host":
                self.query_one("#host", HostPane).set_frame(frame)

    def _shutdown_collector(self) -> None:
        collector = self.collector
        self.collector = None
        if collector is not None:
            collector.close()


__all__ = ["WorkspaceScreen"]
