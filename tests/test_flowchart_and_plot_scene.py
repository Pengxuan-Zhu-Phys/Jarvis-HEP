#!/usr/bin/env python3
"""D11.5 flowchart export + plot scene emitters."""

from __future__ import annotations

import json
import os
import tempfile
import unittest

import yaml

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

    def test_eggbox_calculator_io_ports_match_v1_shape(self) -> None:
        """D12.3: Calculator Dump free-symbols + JSON file ports (V1 Eggbox shape)."""
        config = {
            "Scan": {"name": "EggBox_Bridson_V2"},
            "Sampling": {
                "Method": "Bridson",
                "Variables": [
                    {"name": "x", "distribution": {"type": "Flat"}},
                    {"name": "y", "distribution": {"type": "Flat"}},
                ],
                "LogLikelihood": [
                    {"name": "LogL_Z", "expression": "LogGauss(z, 100, 10)"},
                ],
            },
            "Calculators": {
                "Modules": [
                    {
                        "name": "EggBox",
                        "required_modules": [],
                        "execution": {
                            "input": [
                                {
                                    "name": "inpjson",
                                    "path": "input.json",
                                    "type": "JSON",
                                    "actions": [
                                        {
                                            "type": "Dump",
                                            "variables": [
                                                {"name": "xx", "expression": "x * Pi"},
                                                {"name": "yy", "expression": "y * Pi"},
                                            ],
                                        }
                                    ],
                                }
                            ],
                            "output": [
                                {
                                    "name": "oupjson",
                                    "path": "output.json",
                                    "type": "JSON",
                                    "variables": [{"name": "z"}],
                                }
                            ],
                        },
                    }
                ]
            },
        }
        scene = build_flowchart_scene_from_config(config)
        ids = {node["id"]: node for node in scene["nodes"]}
        self.assertIn("Parameters", ids)
        self.assertIn("var::x", ids)
        self.assertIn("var::y", ids)
        self.assertIn("EggBox", ids)
        self.assertIn("file::EggBox::input::inpjson", ids)
        self.assertIn("file::EggBox::output::oupjson", ids)
        self.assertIn("var::z", ids)
        self.assertIn("LogLikelihood", ids)

        egg = ids["EggBox"]
        self.assertEqual(egg["kind"], "module")
        self.assertEqual(egg["role"], "calculator")
        self.assertEqual(egg["layer"], "layer_2")
        self.assertEqual(ids["Parameters"]["layer"], "layer_1")

        edge_pairs = {
            (
                e["source"]["node"],
                e["target"]["node"],
                e["role"],
            )
            for e in scene["edges"]
        }
        self.assertIn(("Parameters", "var::x", "parameterflow"), edge_pairs)
        self.assertIn(("var::x", "file::EggBox::input::inpjson", "fileflow"), edge_pairs)
        self.assertIn(("var::y", "file::EggBox::input::inpjson", "fileflow"), edge_pairs)
        self.assertIn(
            ("file::EggBox::input::inpjson", "EggBox", "fileflow"), edge_pairs
        )
        self.assertIn(
            ("file::EggBox::output::oupjson", "var::z", "fileflow"), edge_pairs
        )

    def test_jarvisplot_renders_exported_scene_without_v2_adapter(self) -> None:
        """Stock jarvisplot.render_flowchart must accept V2-exported scenes."""
        try:
            from jarvisplot import render_flowchart
        except ImportError:
            self.skipTest("JarvisPLOT not installed")

        config = {
            "Scan": {"name": "render-demo"},
            "Sampling": {
                "Variables": [{"name": "x"}, {"name": "y"}],
            },
            "Calculators": {
                "Modules": [
                    {
                        "name": "EggBox",
                        "execution": {
                            "input": [
                                {
                                    "name": "inpjson",
                                    "path": "in.json",
                                    "type": "JSON",
                                    "actions": [
                                        {
                                            "type": "Dump",
                                            "variables": [
                                                {"name": "xx", "expression": "x"},
                                            ],
                                        }
                                    ],
                                }
                            ],
                            "output": [
                                {
                                    "name": "oupjson",
                                    "path": "out.json",
                                    "type": "JSON",
                                    "variables": [{"name": "z"}],
                                }
                            ],
                        },
                    }
                ]
            },
        }
        scene = build_flowchart_scene_from_config(config, include_likelihood=False)
        with tempfile.TemporaryDirectory() as tmp:
            png = os.path.join(tmp, "flowchart.png")
            rendered = render_flowchart(scene, output_path=png)
            self.assertTrue(os.path.isfile(str(rendered)))
            self.assertGreater(os.path.getsize(str(rendered)), 1000)


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

    def test_jplot_levelset_overlay_hook(self) -> None:
        """D10.3: stock jplot YAML; CSV beside HDF5; images under project/images/scan."""
        from jarvishep2.plot_scene import emit_jplot_scan_levelset_yaml, emit_plot_scenes_from_run

        with tempfile.TemporaryDirectory() as project:
            # Layout: project/outputs/scan/DATABASE + project/images/scan
            scan_dir = os.path.join(project, "outputs", "als_demo")
            db_dir = os.path.join(scan_dir, "DATABASE")
            os.makedirs(db_dir)
            db = os.path.join(db_dir, "samples.hdf5")
            writer = SimpleHDF5Writer(db)
            writer.add_record({"x": 0.1, "y": 0.2, "LogL": -1.0})
            writer.add_record({"x": 0.9, "y": 0.8, "LogL": -2.0})

            levelset = {
                "dim": 2,
                "target_expression": "r2",
                "target_value": 0.25,
                "polylines_x": [
                    [[0.0, 0.5], [0.5, 0.0], [1.0, 0.5], [0.5, 1.0], [0.0, 0.5]]
                ],
            }
            with open(os.path.join(scan_dir, "levelset.json"), "w", encoding="utf-8") as handle:
                json.dump(levelset, handle)

            jplot_path = emit_jplot_scan_levelset_yaml(
                scan_dir, scan_name="als_demo", project_root=project
            )
            self.assertIsNotNone(jplot_path)
            self.assertTrue(os.path.isfile(str(jplot_path)))
            # Plot YAML lives under project/images/<scan>/
            self.assertTrue(
                str(jplot_path).startswith(os.path.join(project, "images", "als_demo"))
            )
            # CSV sits next to HDF5
            self.assertTrue(os.path.isfile(os.path.join(db_dir, "samples.csv")))

            with open(str(jplot_path), encoding="utf-8") as handle:
                doc = yaml.safe_load(handle)
            self.assertIn("DataSet", doc)
            self.assertIn("Figures", doc)
            layers = doc["Figures"][0]["layers"]
            methods = {layer.get("method") for layer in layers}
            self.assertIn("scatter", methods)
            self.assertIn("plot", methods)

            written = emit_plot_scenes_from_run(
                scan_dir,
                scan_name="als_demo",
                project_root=project,
                auto_render=False,
            )
            self.assertIn("jplot_levelset", written)
            self.assertIn("levelset_overlay", written)
            self.assertIn("scatter", written)
            self.assertIn("samples_csv", written)
            self.assertEqual(
                written["samples_csv"], os.path.join(db_dir, "samples.csv")
            )


if __name__ == "__main__":
    unittest.main()
