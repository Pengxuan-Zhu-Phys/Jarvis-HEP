#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import types
import unittest
from unittest import mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.core import Core  # noqa: E402


class _FakeLogger:
    def __init__(self):
        self.warnings = []

    def warning(self, message):
        self.warnings.append(str(message))


class FlowchartJarvisPlotRenderTests(unittest.TestCase):
    def _core_for_output(self, output_path):
        core = Core.__new__(Core)
        core.info = {"flowchart_path": output_path}
        core.logger = _FakeLogger()
        return core

    def test_render_uses_jarvisplot_public_api_with_semantic_dict(self):
        scene = {"schema": "jarvisplot.scene/v1", "scene_type": "flowchart"}
        calls = []
        fake_jarvisplot = types.ModuleType("jarvisplot")

        def render_flowchart(scene_dict, output_path):
            calls.append((scene_dict, output_path))
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write("rendered")

        fake_jarvisplot.render_flowchart = render_flowchart

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "flowchart.png")
            core = self._core_for_output(output_path)
            with mock.patch.dict(sys.modules, {"jarvisplot": fake_jarvisplot}):
                rendered = core._render_flowchart_with_jarvisplot(scene)

            self.assertTrue(rendered)
            self.assertEqual(calls, [(scene, output_path)])
            self.assertTrue(os.path.exists(output_path))

    def test_missing_jarvisplot_skips_render_without_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "flowchart.png")
            core = self._core_for_output(output_path)
            with mock.patch.dict(sys.modules, {"jarvisplot": None}):
                rendered = core._render_flowchart_with_jarvisplot({})

            self.assertFalse(rendered)
            self.assertFalse(os.path.exists(output_path))
            self.assertTrue(
                any("JarvisPLOT is not installed" in msg for msg in core.logger.warnings)
            )


if __name__ == "__main__":
    unittest.main()
