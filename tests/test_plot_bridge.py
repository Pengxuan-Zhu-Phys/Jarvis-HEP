#!/usr/bin/env python3
"""Jarvis-PLOT bridge tests — no heavy matplotlib run required."""

from __future__ import annotations

import unittest
from unittest import mock

from jarvishep2.client import main, normalize_argv
from jarvishep2.plot_bridge import PlotBridgeError, run_plot


class _FakePlotter:
    def __init__(self) -> None:
        self.inited = False

    def init(self) -> None:
        self.inited = True


class _FakePlotterSysExit:
    def init(self) -> None:
        raise SystemExit(3)


class PlotBridgeUnitTests(unittest.TestCase):
    def test_run_plot_invokes_plotter_init(self) -> None:
        code = run_plot("scene.yaml", plotter_cls=_FakePlotter)
        self.assertEqual(code, 0)

    def test_run_plot_empty_path_raises(self) -> None:
        with self.assertRaises(PlotBridgeError):
            run_plot("  ", plotter_cls=_FakePlotter)

    def test_run_plot_maps_system_exit(self) -> None:
        code = run_plot("scene.yaml", plotter_cls=_FakePlotterSysExit)
        self.assertEqual(code, 3)

    def test_run_plot_wraps_engine_errors(self) -> None:
        class _Boom:
            def init(self) -> None:
                raise RuntimeError("bad scene")

        with self.assertRaises(PlotBridgeError) as raised:
            run_plot("scene.yaml", plotter_cls=_Boom)
        self.assertIn("bad scene", str(raised.exception))


class PlotCliTests(unittest.TestCase):
    def test_parse_plot_flag(self) -> None:
        self.assertEqual(
            normalize_argv(["scene.yaml", "--plot"]),
            ["plot", "scene.yaml"],
        )

    def test_plot_passthrough_forwards_the_complete_jplot_argv(self) -> None:
        with mock.patch("jarvisplot.client.main", return_value=7) as plot_main:
            self.assertEqual(main(["plot", "scene.yaml", "--dpi", "300"]), 7)
        plot_main.assert_called_once_with(["scene.yaml", "--dpi", "300"])

    def test_legacy_plot_flag_normalizes_to_the_jplot_passthrough(self) -> None:
        with mock.patch("jarvisplot.client.main", return_value=0) as plot_main:
            self.assertEqual(main(["scene.yaml", "--plot", "--dpi", "300"]), 0)
        plot_main.assert_called_once_with(["scene.yaml", "--dpi", "300"])


if __name__ == "__main__":
    unittest.main()
