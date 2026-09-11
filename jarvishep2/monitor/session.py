"""Attached-scan stub. Overview pages land on top of this screen later."""

from __future__ import annotations

from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Vertical
from textual.screen import Screen
from textual.widgets import Static

from jarvishep2.monitor.hints import SESSION_KEYS, render_key_hint
from jarvishep2.monitor.scans import ScanChoice


class SessionScreen(Screen[None]):
    BINDINGS = [
        Binding("q", "quit_monitor", "Quit", show=True),
        Binding("escape", "back", "Chooser", show=True),
    ]

    def __init__(self, choice: ScanChoice) -> None:
        super().__init__()
        self.choice = choice

    def compose(self) -> ComposeResult:
        control = (
            "—" if self.choice.control_pid is None else str(self.choice.control_pid)
        )
        with Vertical(id="session"):
            yield Static(
                f"Jarvis Monitor  {self.choice.name}  ●  {self.choice.reference}",
                id="session-title",
            )
            yield Static(
                "\n".join(
                    [
                        f"control   {control}",
                        f"processes {self.choice.process_count}",
                        "",
                        "Attached. Overview / Workers / Factory / Sampler /",
                        "Calculators / Samples / Host pages come next.",
                        "",
                        "This screen is a stub so the entry path can ship first.",
                    ]
                ),
                id="session-body",
            )
            yield Static(render_key_hint(SESSION_KEYS), id="hint")

    def action_back(self) -> None:
        self.dismiss(None)

    def action_quit_monitor(self) -> None:
        self.app.exit()
