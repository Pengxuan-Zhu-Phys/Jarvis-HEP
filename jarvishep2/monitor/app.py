"""Textual application: splash chooser, then an attached-scan stub."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from textual.app import App
from textual.screen import Screen

from jarvishep2.monitor.scans import ScanChoice, resolve_choice
from jarvishep2.monitor.session import SessionScreen
from jarvishep2.monitor.splash import SplashScreen
from jarvishep2.monitor.workspace import WorkspaceScreen

ScanLister = Callable[[], list]


class MonitorApp(App[None]):
    """No composer. Navigation is the chooser table plus later page keys."""

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

    def attach_scan(self, choice: ScanChoice) -> None:
        if choice.simulated:
            self.push_screen(WorkspaceScreen(choice), self._on_session_closed)
            return
        self.push_screen(SessionScreen(choice), self._on_session_closed)

    def _on_session_closed(self, _result: None) -> None:
        screen = self.screen
        if isinstance(screen, SplashScreen):
            screen.action_refresh()
