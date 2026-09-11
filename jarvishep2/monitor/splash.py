"""Logo splash + live scan chooser (``Jarvis monitor`` entry)."""

from __future__ import annotations

import os
from typing import Callable

from rich.text import Text
from textual.app import ComposeResult
from textual.binding import Binding
from textual.containers import Horizontal, Vertical
from textual.events import Key
from textual.screen import Screen
from textual.widgets import DataTable, Static

from jarvishep2.monitor.branding import (
    MONITOR_FRAMES,
    REVEAL_FRAMES,
    load_branding,
    monitor_tag_markup,
    render_banner_markup,
    render_logo_monitor_frame,
)
from jarvishep2.monitor.hints import SPLASH_KEYS, render_key_hint
from jarvishep2.monitor.scans import ScanChoice, list_scan_choices, simulated_choice
from jarvishep2.monitor.chrome import MonitorTopbar

CONTROL_COL_WIDTH = 10
PROCS_COL_WIDTH = 5
REF_COL_WIDTH = 4

ScanLister = Callable[[], list]


class SplashScreen(Screen[None]):
    """Agent-style hero plus the ``Jarvis ps`` task list. No composer."""

    BINDINGS = [
        Binding("q", "quit_monitor", "Quit", show=True),
        Binding("enter", "attach", "Attach", show=True),
        Binding("r", "refresh", "Refresh", show=True),
        Binding("j", "cursor_down", "Down", show=False),
        Binding("k", "cursor_up", "Up", show=False),
    ]

    def __init__(
        self,
        *,
        scan_lister: ScanLister | None = None,
        notice: str = "",
    ) -> None:
        super().__init__()
        self._scan_lister = scan_lister
        self._notice = notice
        self._choices: list[ScanChoice] = []
        self._branding = load_branding()
        self._frame = 0
        self._animate = os.environ.get("NO_COLOR", "") == ""
        self._closing = False

    def compose(self) -> ComposeResult:
        with Vertical(id="splash"):
            yield MonitorTopbar()
            with Horizontal(id="hero"):
                yield Static(id="logo-monitor")
                with Vertical(id="home-panel"):
                    yield Static(id="banner")
                    yield Static(id="monitor-tag")
            with Vertical(id="chooser"):
                yield Static("running scans", id="chooser-title")
                yield DataTable(id="scans", cursor_type="row")
            yield Static("", id="notice")
            yield Static(render_key_hint(SPLASH_KEYS), id="hint")

    def on_mount(self) -> None:
        self.query_one("#banner", Static).update(
            render_banner_markup(self._branding.banner_lines)
        )
        self.query_one("#monitor-tag", Static).update(monitor_tag_markup())
        table = self.query_one("#scans", DataTable)
        table.cursor_type = "row"
        table.zebra_stripes = True
        table.add_column("REF", width=REF_COL_WIDTH, key="ref")
        table.add_column("SCAN", width=self._scan_column_width(), key="scan")
        table.add_column(
            Text("CONTROL", justify="right"),
            width=CONTROL_COL_WIDTH,
            key="control",
        )
        table.add_column(
            Text("PROCS", justify="right"),
            width=PROCS_COL_WIDTH,
            key="procs",
        )
        table.cell_padding = 1
        self.action_refresh()
        if self._notice:
            self.set_notice(self._notice)
        self._layout_scan_column()
        self._paint_logo()
        self.call_after_refresh(self._after_first_layout)

    def _after_first_layout(self) -> None:
        self._layout_scan_column()
        if self._animate:
            self.set_interval(0.12, self._tick_logo)

    def on_resize(self) -> None:
        self._layout_scan_column()

    def _scan_column_width(self) -> int:
        try:
            table = self.query_one("#scans", DataTable)
            usable = table.size.width or self.size.width
            pad = int(getattr(table, "cell_padding", 1))
        except Exception:
            usable = self.size.width
            pad = 1
        if usable <= 0:
            usable = 80
        used = REF_COL_WIDTH + CONTROL_COL_WIDTH + PROCS_COL_WIDTH + 8 * pad
        return max(8, usable - used)

    def _layout_scan_column(self) -> None:
        """Pin CONTROL/PROCS to the right; SCAN eats the leftover width."""
        try:
            table = self.query_one("#scans", DataTable)
        except Exception:
            return
        scan_w = self._scan_column_width()
        for column in table.ordered_columns:
            if column.key.value != "scan":
                continue
            column.auto_width = False
            if column.width != scan_w:
                column.width = scan_w
                table.refresh()
            break

    def set_notice(self, message: str) -> None:
        self._notice = message
        try:
            widget = self.query_one("#notice", Static)
            widget.update(message)
            widget.display = bool(message)
        except Exception:
            pass

    def action_refresh(self) -> None:
        self._choices = list_scan_choices(self._scan_lister)
        table = self.query_one("#scans", DataTable)
        cursor = table.cursor_row
        table.clear()
        if not self._notice.startswith("Unknown"):
            self.set_notice("")
        if not self._choices:
            table.add_row(
                "-",
                "-",
                Text("-", justify="right"),
                Text("-", justify="right"),
                key="empty",
            )
            table.focus()
            self.call_after_refresh(self._layout_scan_column)
            return
        for choice in self._choices:
            control = "-" if choice.control_pid is None else str(choice.control_pid)
            table.add_row(
                choice.reference,
                choice.name or "-",
                Text(control, justify="right"),
                Text(str(choice.process_count), justify="right"),
                key=choice.reference,
            )
        if table.row_count:
            table.move_cursor(row=min(max(0, cursor), table.row_count - 1))
        table.focus()
        self.call_after_refresh(self._layout_scan_column)

    def action_cursor_down(self) -> None:
        table = self.query_one("#scans", DataTable)
        if table.row_count:
            table.action_cursor_down()

    def action_cursor_up(self) -> None:
        table = self.query_one("#scans", DataTable)
        if table.row_count:
            table.action_cursor_up()

    def action_attach(self) -> None:
        if self._closing:
            return
        choice = self._selected()
        if choice is None:
            return
        self._closing = True
        attach = getattr(self.app, "attach_scan", None)
        if callable(attach):
            attach(choice)
        self._closing = False

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        event.stop()
        self.action_attach()

    def action_quit_monitor(self) -> None:
        self.app.exit()

    def on_key(self, event: Key) -> None:
        if event.key == "escape":
            event.stop()
            event.prevent_default()
            return
        if event.character not in {"s", "S"} and event.key not in {"s", "S"}:
            return
        event.stop()
        event.prevent_default()
        self.action_simulate()

    def action_simulate(self) -> None:
        if self._closing:
            return
        self._closing = True
        attach = getattr(self.app, "attach_scan", None)
        if callable(attach):
            attach(simulated_choice())
        self._closing = False

    def _selected(self) -> ScanChoice | None:
        if not self._choices:
            return None
        table = self.query_one("#scans", DataTable)
        index = table.cursor_row
        if index < 0 or index >= len(self._choices):
            return None
        return self._choices[index]

    def _tick_logo(self) -> None:
        max_frame = REVEAL_FRAMES + MONITOR_FRAMES
        if self._frame >= max_frame:
            return
        self._frame += 1
        self._paint_logo()

    def _paint_logo(self) -> None:
        self.query_one("#logo-monitor", Static).update(
            render_logo_monitor_frame(
                self._frame,
                self._branding.logo_pattern,
                animate=self._animate,
            )
        )
