"""Release workspace regressions for unequal columns and scan completion."""

import asyncio
from dataclasses import replace
import subprocess
import sys
from unittest.mock import Mock, patch

import psutil
import pytest

from jarvishep2.monitor.app import MonitorApp
from jarvishep2.monitor.calculators import CalculatorsBlockData
from jarvishep2.monitor.collector import Collector
from jarvishep2.monitor.overview import HealthDetailModal
from jarvishep2.monitor.scans import ScanChoice, ScanExitWatch
from jarvishep2.monitor.simu import SimuEngine
from jarvishep2.monitor.splash import SplashScreen
from jarvishep2.redis_queue import make_fakeredis_queue


def test_pacman_follows_shorter_column_and_disappears_when_equal():
    async def run():
        app = MonitorApp(scan_lister=lambda: [])
        async with app.run_test(size=(120, 64)) as pilot:
            await pilot.press("s")
            await pilot.pause()
            app.screen.engine = None  # Hold data fixed while testing layout changes.
            pane = app.screen.query_one("#overview")
            frame = pane.frame
            for width in (120, 100):
                await pilot.resize_terminal(width, 64)
                # Sync histories after the width change without recording a sample.
                fixture = SimuEngine()
                fixture.set_history_widths(*pane.spark_history_widths())
                frame = fixture.empty_frame()
                for count in (1, 8, 7, 1):
                    pane.set_frame(replace(
                        frame,
                        calculator_block=CalculatorsBlockData(frame.calculator_block.pools[:count]),
                    ))
                    # Allow the newly visible board's first animation tick.
                    await pilot.pause(0.2)
                    left = pane.query_one("#ov-pacman-left")
                    right = pane.query_one("#ov-pacman")
                    calcs = pane.query_one("#ov-calcs")
                    workers = pane.query_one("#ov-workers")
                    difference = workers.region.bottom - calcs.region.bottom
                    assert left.region.height == max(0, difference)
                    assert right.region.height == max(0, -difference)
                    if difference:
                        game = left if difference > 0 else right
                        block = calcs if difference > 0 else workers
                        assert game.region.y == block.region.bottom
                        assert game.region.bottom == max(calcs.region.bottom, workers.region.bottom)
                        assert game.region.x == block.region.x
                        assert game.region.width == block.region.width
                        assert "·" in str(game.render())
                        assert "👻" in str(game.render())
            await pilot.press("q")

    asyncio.run(run())


@pytest.mark.parametrize("open_help", (False, True))
def test_scan_exit_returns_to_chooser_even_with_health_dialog_open(open_help):
    async def run(process):
        lister = Mock(return_value=[])
        app = MonitorApp(scan_lister=lister)
        collector = Collector(make_fakeredis_queue(), host_snapshotter=lambda: {})
        collector.close = Mock(wraps=collector.close)
        choice = ScanChoice("R1", "test-scan", process.pid, 1, (process.pid,))
        async with app.run_test(size=(120, 64)) as pilot:
            await pilot.pause()
            app.attach_scan(choice, collector=collector)
            await pilot.pause()
            workspace = app.screen
            # A failed Redis poll while Core lives must keep the page open.
            with patch.object(collector._redis, "snapshot_raw", side_effect=TimeoutError):
                workspace._tick()
            workspace._check_scan_exit()
            assert app.screen is workspace
            if open_help:
                await pilot.click("#ov-health-core", offset=(0, 0))
                await pilot.pause()
                assert isinstance(app.screen, HealthDetailModal)
            before_refresh = lister.call_count
            process.stdin.close()
            process.wait(timeout=5)
            # The workspace's timer must continue watching through a modal.
            await pilot.pause(0.8)
            assert isinstance(app.screen, SplashScreen)
            assert lister.call_count > before_refresh
            collector.close.assert_called_once()
            assert workspace not in app.screen_stack
            await pilot.press("q")

    with subprocess.Popen(
        [sys.executable, "-c", "import sys; sys.stdin.read()"],
        stdin=subprocess.PIPE,
    ) as process:
        try:
            asyncio.run(run(process))
        finally:
            if not process.stdin.closed:
                process.stdin.close()
            process.wait(timeout=5)


def test_scan_watch_does_not_treat_access_denied_as_exit():
    process = Mock()
    process.is_running.side_effect = psutil.AccessDenied(123)
    with patch("psutil.Process", return_value=process):
        watch = ScanExitWatch(ScanChoice("R1", "test", 123, 1, (123,)))
        assert not watch.has_exited()
        process.is_running.side_effect = None
        process.is_running.return_value = False  # Includes recycled PID detection.
        assert watch.has_exited()
