#!/usr/bin/env python3
"""Jarvis-PLOT bridge tests — no heavy matplotlib run required."""

from __future__ import annotations

import unittest
from unittest import mock

from jarvishep2.client import build_parser, dispatch
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
        args = build_parser().parse_args(["--plot", "scene.yaml"])
        self.assertTrue(args.plot)
        self.assertEqual(args.task_yaml, "scene.yaml")

    def test_dispatch_plot_routes_to_bridge(self) -> None:
        args = build_parser().parse_args(["--plot", "scene.yaml"])
        with mock.patch("jarvishep2.client.run_plot", return_value=0) as run:
            code = dispatch(args)
        self.assertEqual(code, 0)
        run.assert_called_once_with("scene.yaml")

    def test_dispatch_plot_missing_yaml_is_usage(self) -> None:
        args = build_parser().parse_args(["--plot"])
        code = dispatch(args)
        self.assertEqual(code, 2)

    def test_dispatch_plot_missing_package_is_exit_2(self) -> None:
        args = build_parser().parse_args(["--plot", "scene.yaml"])
        with mock.patch(
            "jarvishep2.client.run_plot",
            side_effect=ImportError("install plot"),
        ):
            code = dispatch(args)
        self.assertEqual(code, 2)


if __name__ == "__main__":
    unittest.main()
