"""Scan workspace: topbar + tabs + page + hint. Shared by every view."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.events import Key
from textual.screen import Screen
from textual.widgets import ContentSwitcher

from jarvishep2.monitor.chrome import PAGES, MonitorHint, MonitorTopbar, TabBar
from jarvishep2.monitor.hints import OVERVIEW_KEYS
from jarvishep2.monitor.overview import ComingPane, OverviewPane
from jarvishep2.monitor.scans import ScanChoice
from jarvishep2.monitor.simu import SimuEngine


class WorkspaceScreen(Screen[None]):
    BINDINGS = [
        Binding("q", "quit_monitor", "Quit", show=True),
        Binding("escape", "back", "Chooser", show=True),
        Binding("r", "refresh", "Refresh", show=False),
        Binding("tab", "next_page", show=False, priority=True),
        Binding("shift+tab", "prev_page", show=False, priority=True),
        Binding("left", "prev_page", show=False),
        Binding("right", "next_page", show=False),
    ]

    def __init__(self, choice: ScanChoice) -> None:
        super().__init__()
        self.choice = choice
        self.engine = SimuEngine() if choice.simulated else None
        self._page = "overview"

    def compose(self) -> ComposeResult:
        with Vertical(id="workspace"):
            yield MonitorTopbar()
            yield TabBar()
            with Vertical(id="page-frame"):
                with ContentSwitcher(initial="overview", id="pages"):
                    yield OverviewPane(id="overview")
                    yield ComingPane("Workers", id="workers")
                    yield ComingPane("Factory", id="factory")
                    yield ComingPane("Sampler", id="sampler")
                    yield ComingPane("Calculators", id="calculators")
                    yield ComingPane("Samples", id="samples")
                    yield ComingPane("Host", id="host")
            yield MonitorHint()

    def on_mount(self) -> None:
        self.query_one(TabBar).set_active("overview")
        if self.engine is not None:
            pane = self.query_one("#overview", OverviewPane)
            pane.set_frame(self.engine.tick())
            self.set_interval(0.5, self._tick)

    def show_page(self, slug: str) -> None:
        slugs = {item[0] for item in PAGES}
        if slug not in slugs:
            return
        self._page = slug
        self.query_one("#pages", ContentSwitcher).current = slug
        self.query_one(TabBar).set_active(slug)
        self.query_one(MonitorHint).set_keys(OVERVIEW_KEYS)

    def action_refresh(self) -> None:
        self._tick()

    def action_next_page(self) -> None:
        self._cycle(1)

    def action_prev_page(self) -> None:
        self._cycle(-1)

    def action_back(self) -> None:
        self.dismiss(None)

    def action_quit_monitor(self) -> None:
        self.app.exit()

    def on_key(self, event: Key) -> None:
        if event.key in {"1", "2", "3", "4", "5", "6", "7"}:
            event.stop()
            event.prevent_default()
            self.show_page(PAGES[int(event.key) - 1][0])

    def _cycle(self, delta: int) -> None:
        slugs = [item[0] for item in PAGES]
        index = slugs.index(self._page) if self._page in slugs else 0
        self.show_page(slugs[(index + delta) % len(slugs)])

    def _tick(self) -> None:
        if self.engine is None:
            return
        frame = self.engine.tick()
        pane = self.query_one("#overview", OverviewPane)
        pane.set_frame(frame)


__all__ = ["WorkspaceScreen"]
