#!/usr/bin/env python3
"""Monitor TUI entry: splash chooser, branding, CLI --once."""

from __future__ import annotations

import asyncio
import io
import json
import re
import unittest
from contextlib import redirect_stdout
from unittest import mock

from jarvishep2.client import build_parser, dispatch_monitor, main
from jarvishep2.monitor.branding import (
    load_branding,
    render_banner_markup,
    render_logo_monitor_frame,
)
from jarvishep2.monitor.hints import SPLASH_KEYS, render_key_hint, render_key_hint_plain
from jarvishep2.monitor.bars import bar_rail, eighths, visible_width
from jarvishep2.monitor.overview import render_resource_row
from jarvishep2.monitor.chrome import folder_tab_lines
from jarvishep2.monitor.scans import (
    choice_from_scan,
    list_scan_choices,
    resolve_choice,
    simulated_choice,
)
from jarvishep2.monitor.topbar import (
    clock_label,
    compact_path,
    render_topbar_clock,
    render_topbar_left,
)
from jarvishep2.process_cleanup import JarvisProcess, JarvisScan
from jarvishep2.run_outcome import EXIT_OK, EXIT_USAGE


def _scan(reference: str, name: str, control_pid: int, extra: int = 2) -> JarvisScan:
    processes = [JarvisProcess(control_pid, f"Jarvis:{name}")]
    for index in range(extra):
        processes.append(
            JarvisProcess(
                control_pid + 1 + index,
                f"Jarvis-Worker-{index:02d}:{name}",
            )
        )
    return JarvisScan(
        reference=reference,
        name=name,
        processes=tuple(processes),
    )


class ResourceRowTests(unittest.TestCase):
    def test_cdot_is_on_the_center_column(self) -> None:
        row = render_resource_row(76, 41.0, 18.4, 64.0)
        self.assertEqual(visible_width(row), 76)
        plain = re.sub(r"\[/?[^\]]*\]", "", row)
        self.assertEqual(plain[76 // 2], "·")
        cpu_len = plain.index("┤") - plain.index("├")
        mem_len = plain.rindex("┤") - plain.rindex("├")
        self.assertGreater(cpu_len, mem_len)


class ProgressBarTests(unittest.TestCase):
    def test_rail_gradient_and_visible_width(self) -> None:
        self.assertEqual(len(eighths(0.5, 10, empty="─")), 10)
        rail = bar_rail(0.63, 20)
        shown = visible_width(rail)
        self.assertEqual(shown, 20)
        self.assertTrue(shown >= 2)
        self.assertIn("├", rail)
        self.assertIn("┤", rail)
        self.assertIn("#134a8d", rail)
        self.assertIn("#73b8f4", bar_rail(1.0, 12))


class FolderTabTests(unittest.TestCase):
    def test_active_stays_in_order_and_idles_share_space(self) -> None:
        _top, mid, join, hits = folder_tab_lines(78, 1)
        self.assertIn("│ 1 Overview │", mid)
        self.assertEqual(hits[0][0], "overview")
        self.assertEqual(hits[0][1], 0)
        self.assertEqual(join[0], "│")
        self.assertEqual(join[-1], "╮")

        _top, mid, join, hits = folder_tab_lines(78, 2)
        self.assertIn("│ 2 Workers │", mid)
        self.assertEqual(hits[1][0], "workers")
        self.assertGreater(hits[1][1], 0)
        self.assertLess(hits[0][1], hits[1][1])
        self.assertEqual(join[0], "╭")
        self.assertIn("╯", join)
        self.assertIn("╰", join)
        self.assertNotIn("┌", _top)
        self.assertIn("╭", _top)

    def test_idle_names_shorten_only_when_needed_then_gaps(self) -> None:
        _top, mid, _join, _hits = folder_tab_lines(120, 1)
        self.assertIn("Calculators", mid)
        self.assertNotIn("Cal.", mid)
        without_active = mid.split("│", 2)[-1]
        self.assertIn("     ", without_active)

        _top, mid, _join, _hits = folder_tab_lines(78, 1)
        self.assertIn("Cal.", mid)
        self.assertNotIn("Calculators", mid)
        self.assertNotIn("Cal.6", mid)

        _top, mid, _join, _hits = folder_tab_lines(61, 1)
        self.assertNotIn("Calculators", mid)
        self.assertNotIn("Cal.", mid)
        self.assertIn("Ca.", mid)


class SimulatedScanTests(unittest.TestCase):
    def test_fixture_is_marked_simulated(self) -> None:
        choice = simulated_choice()
        self.assertTrue(choice.simulated)
        self.assertEqual(choice.reference, "SIM")
        self.assertTrue(choice.name.startswith("simu-"))

    def test_engine_ticks_change_the_frame(self) -> None:
        from jarvishep2.monitor.simu import SimuEngine

        engine = SimuEngine()
        first = engine.tick()
        for _ in range(6):
            second = engine.tick()
        self.assertGreaterEqual(second.done, first.done)
        self.assertNotEqual(first.spark_samples, second.spark_samples)


class ScanChoiceTests(unittest.TestCase):
    def test_projects_ps_fields(self) -> None:
        scan = _scan("R1", "iDM_Vector_V1", 44001, extra=3)
        choice = choice_from_scan(scan)
        self.assertEqual(choice.reference, "R1")
        self.assertEqual(choice.name, "iDM_Vector_V1")
        self.assertEqual(choice.control_pid, 44001)
        self.assertEqual(choice.process_count, 4)

    def test_list_and_resolve(self) -> None:
        scans = [
            _scan("R1", "iDM_Vector_V1", 44001),
            _scan("R2", "eggbox-bridson", 51020),
        ]
        choices = list_scan_choices(lambda: scans)
        self.assertEqual([choice.reference for choice in choices], ["R1", "R2"])
        picked = resolve_choice("R2", choices, scans=scans)
        self.assertEqual(picked.name, "eggbox-bridson")
        picked = resolve_choice("44001", choices, scans=scans)
        self.assertEqual(picked.reference, "R1")


class KeyHintTests(unittest.TestCase):
    def test_markup_is_key_colon_meaning_with_bars(self) -> None:
        markup = render_key_hint(SPLASH_KEYS)
        self.assertIn("[bold #f6d33f]Enter[/]", markup)
        self.assertIn("[#8d93a1]: attach[/]", markup)
        self.assertIn("|", markup)
        self.assertIn("J/K", markup)
        self.assertIn("Q", markup)
        self.assertEqual(
            render_key_hint_plain(SPLASH_KEYS),
            "Enter: attach  |  J/K: select  |  R: refresh  |  Q: quit",
        )
        plain = render_key_hint_plain(SPLASH_KEYS)
        self.assertNotIn("S:", plain)
        self.assertNotIn("simu", plain.lower())


class TopbarTests(unittest.TestCase):
    def test_clock_and_left_piece(self) -> None:
        from datetime import datetime

        self.assertEqual(clock_label(datetime(2026, 9, 10, 20, 53)), "8:53 PM")
        left = render_topbar_left(branch="main", path="~/Jarvis-HEP")
        self.assertIn("⎇ main", left)
        self.assertIn("~/Jarvis-HEP", left)
        self.assertIn(" · ", left)
        self.assertFalse(left.endswith("·"))
        bare = render_topbar_left(branch="", path="~/Jarvis-HEP")
        self.assertEqual(bare, "~/Jarvis-HEP")
        self.assertNotIn("⎇", bare)
        self.assertNotIn("·", bare)
        self.assertEqual(
            render_topbar_clock(datetime(2026, 9, 10, 20, 53)),
            "· 8:53 PM",
        )
        self.assertEqual(compact_path("abcdefghij", 7), "...ghij")


class BrandingTests(unittest.TestCase):
    def test_logo_and_banner_come_from_hep_card(self) -> None:
        branding = load_branding()
        self.assertEqual(len(branding.logo_pattern), 8)
        self.assertTrue(any("Version:" in line for line in branding.banner_lines))
        markup = render_banner_markup(branding.banner_lines)
        self.assertIn("2.0.", markup)
        frame = render_logo_monitor_frame(40, branding.logo_pattern, animate=False)
        lines = frame.splitlines()
        self.assertEqual(len(lines), 8)
        self.assertTrue(all("⬤" in line for line in lines))


class CliMonitorOnceTests(unittest.TestCase):
    def test_once_without_ref_prints_scan_table(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        args = build_parser().parse_args(["monitor", "--once"])
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with redirect_stdout(buf):
                code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_OK)
        self.assertIn("iDM_Vector_V1", buf.getvalue())
        self.assertIn("R1", buf.getvalue())

    def test_once_json_emits_rows(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        args = build_parser().parse_args(["monitor", "--once", "--json"])
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with redirect_stdout(buf):
                code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_OK)
        payload = json.loads(buf.getvalue())
        self.assertEqual(payload[0]["reference"], "R1")
        self.assertEqual(payload[0]["control_pid"], 44001)

    def test_unknown_ref_with_once_is_usage(self) -> None:
        args = build_parser().parse_args(["monitor", "R9", "--once"])
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=[],
        ):
            code = dispatch_monitor(args)
        self.assertEqual(code, EXIT_USAGE)

    def test_non_tty_monitor_lists_scans(self) -> None:
        scans = [_scan("R1", "iDM_Vector_V1", 44001)]
        buf = io.StringIO()
        with mock.patch(
            "jarvishep2.process_cleanup.list_active_scans",
            return_value=scans,
        ):
            with mock.patch("sys.stdout.isatty", return_value=False):
                with redirect_stdout(buf):
                    code = main(["monitor"])
        self.assertEqual(code, EXIT_OK)
        self.assertIn("iDM_Vector_V1", buf.getvalue())


class MonitorAppTests(unittest.TestCase):
    def test_splash_lists_ps_rows_and_quit(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [
            _scan("R1", "iDM_Vector_V1", 44001),
            _scan("R2", "eggbox-bridson", 51020),
        ]

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                logo = app.screen.query_one("#logo-monitor", Static)
                self.assertIn("⬤", str(logo.render()))
                table = app.screen.query_one("#scans", DataTable)
                self.assertEqual(table.row_count, 2)
                hint = app.screen.query_one("#hint", Static)
                hint_text = str(hint.render())
                self.assertIn("Enter", hint_text)
                self.assertIn("attach", hint_text)
                self.assertNotIn("simu", hint_text.lower())
                self.assertNotRegex(hint_text, r"\bS:")
                top = str(app.screen.query_one("#topbar-left", Static).render())
                self.assertIn("⎇", top)
                clock = str(app.screen.query_one("#topbar-clock", Static).render())
                self.assertRegex(clock, r"· \d{1,2}:\d{2} (AM|PM)")
                await pilot.press("q")

        asyncio.run(_run())

    def test_escape_does_not_quit_splash(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [_scan("R1", "iDM_Vector_V1", 44001)]

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                await pilot.press("escape")
                await pilot.pause()
                self.assertFalse(app._exit)
                app.screen.query_one("#logo-monitor", Static)
                table = app.screen.query_one("#scans", DataTable)
                scan_col = next(
                    col for col in table.ordered_columns if col.key.value == "scan"
                )
                self.assertFalse(scan_col.auto_width)
                self.assertGreaterEqual(scan_col.width, 8)
                await pilot.press("q")

        asyncio.run(_run())

    def test_empty_scan_list_shows_dash_not_start_hint(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import DataTable, Static

                table = app.screen.query_one("#scans", DataTable)
                self.assertGreaterEqual(table.row_count, 1)
                scan_cell = str(table.get_row_at(0)[1])
                self.assertEqual(scan_cell.strip(), "-")
                notice = str(app.screen.query_one("#notice", Static).render())
                self.assertNotIn("Start one", notice)
                await pilot.press("q")

        asyncio.run(_run())

    def test_direct_ref_opens_session_stub(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        scans = [_scan("R1", "iDM_Vector_V1", 44001)]

        async def _run() -> None:
            app = MonitorApp(scan_ref="R1", scan_lister=lambda: scans)
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import Static

                title = app.screen.query_one("#session-title", Static)
                self.assertIn("iDM_Vector_V1", str(title.render()))
                self.assertIn("R1", str(title.render()))
                await pilot.press("q")

        asyncio.run(_run())

    def test_s_attaches_simulated_scan_without_hint(self) -> None:
        try:
            from jarvishep2.monitor.app import MonitorApp
        except ImportError:
            self.skipTest("textual is not installed")

        async def _run() -> None:
            app = MonitorApp(scan_lister=lambda: [])
            async with app.run_test() as pilot:
                await pilot.pause()
                from textual.widgets import Static

                hint = str(app.screen.query_one("#hint", Static).render())
                self.assertNotIn("simu", hint.lower())
                await pilot.press("s")
                await pilot.pause()
                status = app.screen.query_one("#ov-status")
                self.assertEqual(status.border_title, "STATUS")
                self.assertIn("simu-iDM_Vector_V1", str(status.render()))
                self.assertIn("SIM", str(status.render()))
                samples = app.screen.query_one("#ov-samples")
                self.assertEqual(samples.border_title, "SAMPLES")
                resources = app.screen.query_one("#ov-resources")
                self.assertEqual(resources.border_title, "RESOURCES")
                self.assertIn("CPU", str(resources.render()))
                self.assertIn("MEM", str(resources.render()))
                self.assertIn("/", resources.border_subtitle)
                self.assertIn("├", str(resources.render()))
                app.screen.query_one("#topbar")
                app.screen.query_one("#tab-bar")
                app.screen.query_one("#hint")
                from textual.widgets import ContentSwitcher

                await pilot.press("tab")
                await pilot.pause()
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "workers",
                )
                await pilot.press("left")
                await pilot.pause()
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "overview",
                )
                await pilot.press("right")
                await pilot.pause()
                self.assertEqual(
                    app.screen.query_one("#pages", ContentSwitcher).current,
                    "workers",
                )
                await pilot.press("q")

        asyncio.run(_run())


if __name__ == "__main__":
    unittest.main()
