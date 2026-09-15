"""Textual application: splash chooser, then an attached-scan stub."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from textual.app import App
from textual.screen import Screen

from jarvishep2.monitor.collector import Collector, open_collector
from jarvishep2.monitor.overview_workspace import OverviewReleaseScreen
from jarvishep2.monitor.scans import ScanChoice, resolve_choice
from jarvishep2.monitor.splash import SplashScreen

ScanLister = Callable[[], list]


class MonitorApp(App[None]):
    """Chooser plus the Overview-only Monitor v1 workspace."""

    CSS_PATH = Path(__file__).with_name("theme.tcss")
    TITLE = "Jarvis Monitor"

    def __init__(
        self,
        *,
        scan_ref: str | None = None,
        scan_lister: ScanLister | None = None,
    ) -> None:
        self._scan_ref = str(scan_ref or "").strip() or None
        self._scan_lister = scan_lister
        self._direct: ScanChoice | None = None
        self._splash_notice = ""
        if self._scan_ref:
            try:
                self._direct = resolve_choice(
                    self._scan_ref, [], lister=self._scan_lister
                )
            except ValueError as exc:
                self._splash_notice = str(exc)
        super().__init__()

    def get_default_screen(self) -> Screen:
        return SplashScreen(
            scan_lister=self._scan_lister,
            notice=self._splash_notice,
        )

    def on_mount(self) -> None:
        if self._direct is not None:
            self.attach_scan(self._direct)

    def attach_scan(
        self,
        choice: ScanChoice,
        *,
        collector: Collector | None = None,
    ) -> None:
        if choice.simulated:
            self.push_screen(OverviewReleaseScreen(choice), self._on_session_closed)
            return
        if collector is None:
            collector = open_collector(choice)
        self.push_screen(
            OverviewReleaseScreen(choice, collector=collector),
            self._on_session_closed,
        )

    def _on_session_closed(self, _result: None) -> None:
        # Screen.dismiss invokes its result callback before the Overview is
        # actually popped. Defer the chooser refresh until the next UI turn so
        # it never clears a table while the screen stack is in transition.
        self.call_later(self._refresh_chooser)

    def _refresh_chooser(self) -> None:
        screen = self.screen
        if isinstance(screen, SplashScreen):
            screen.action_refresh()
