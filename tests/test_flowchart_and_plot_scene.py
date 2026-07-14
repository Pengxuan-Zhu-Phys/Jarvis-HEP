#!/usr/bin/env python3
"""D11.5 flowchart export + plot scene emitters."""

from __future__ import annotations

import json
import os
import tempfile
import unittest

from jarvishep2.flowchart import (
    build_flowchart_scene_from_config,
    export_flowchart_semantics,
)
from jarvishep2.plot_scene import (
    emit_levelset_overlay_scene,
    emit_scan_scatter_scene_from_hdf5,
)
from jarvishep2.database import SimpleHDF5Writer


class FlowchartExportTests(unittest.TestCase):
    def test_build_scene_from_config_has_layers_and_modules(self) -> None:
        config = {
            "project_name": "demo",
            "Scan": {"name": "demo-scan"},
            "Operas": {
                "Modules": [
                    {
                        "name": "TrivialEggbox",
                        "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                    }
                ]
            },
        }
        scene = build_flowchart_scene_from_config(config)
        self.assertEqual(scene["schema"], "jarvisplot.scene/v1")
        self.assertEqual(scene["scene_type"], "flowchart")
        self.assertTrue(any(n["id"] == "Parameters" for n in scene["nodes"]))
        self.assertTrue(any(n["id"] == "TrivialEggbox" for n in scene["nodes"]))
        self.assertTrue(any(n["id"] == "LogLikelihood" for n in scene["nodes"]))
        self.assertGreaterEqual(len(scene["edges"]), 1)

    def test_export_writes_json(self) -> None:
        config = {
            "Scan": {"name": " cons"},
            "Operas": {"Modules": [{"name": "A", "operator": "x.y"}]},
        }
        scene = build_flowchart_scene_from_config(config, include_likelihood=False)
        with tempfile.TemporaryDirectory() as tmp:
            path = export_flowchart_semantics(scene, os.path.join(tmp, "flowchart.json"))
            self.assertTrue(os.path.isfile(path))
            payload = json.loads(open(path, encoding="utf-8").read())
            self.assertEqual(payload["scene_type"], "flowchart")


class PlotSceneEmitTests(unittest.TestCase):
    def test_levelset_overlay_scene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            levelset = {
                "dim": 2,
                "target_expression": "z",
                "target_value": 1.0,
                "polylines_x": [[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]],
            }
            src = os.path.join(tmp, "levelset.json")
            with open(src, "w", encoding="utf-8") as handle:
                json.dump(levelset, handle)
            out = emit_levelset_overlay_scene(
                src, output_yaml=os.path.join(tmp, "overlay.yaml")
            )
            self.assertIsNotNone(out)
            self.assertTrue(os.path.isfile(str(out)))

    def test_scatter_scene_from_hdf5(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = os.path.join(tmp, "samples.hdf5")
            writer = SimpleHDF5Writer(db)
            writer.add_record({"x": 0.1, "y": 0.2, "LogL": -1.0})
            writer.add_record({"x": 0.3, "y": 0.4, "LogL": -2.0})
            out = emit_scan_scatter_scene_from_hdf5(
                db, output_yaml=os.path.join(tmp, "scatter.yaml")
            )
            self.assertIsNotNone(out)
            self.assertTrue(os.path.isfile(str(out)))


if __name__ == "__main__":
    unittest.main()
