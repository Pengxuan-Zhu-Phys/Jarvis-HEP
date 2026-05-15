#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import json
import os
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.parameters import Parameters  # noqa: E402
from jarvishep.workflow import Workflow  # noqa: E402


class _FakeCalculatorModule:
    def __init__(self):
        self.name = "CalcA"
        self.type = "Calculator"
        self.required_modules = []
        self.input = [
            {
                "name": "input_file",
                "type": "SLHA",
                "path": "/tmp/input.slha",
                "variables": {"x": {"name": "x"}},
            }
        ]
        self.output = [
            {
                "name": "output_file",
                "type": "SLHA",
                "path": "/tmp/output.slha",
                "variables": [{"name": "y"}],
            }
        ]
        self.inputs = {"x": None}
        self.outputs = {"y": None}


class _FakeRelayModule:
    def __init__(self):
        self.name = "CalcB"
        self.type = "Calculator"
        self.required_modules = []
        self.input = ["y"]
        self.output = []
        self.inputs = {"y": None}
        self.outputs = {"z": None}


class _FakeSelectionModule:
    def __init__(self):
        self.name = "SelectedOperas"
        self.type = "Operas"
        self.required_modules = []
        self.input = []
        self.output = []
        self.inputs = {}
        self.outputs = {}
        self.operator = "noop"
        self.call_mode = "call"
        self.selection = "y > 0"
        self._selection_deps = ("y",)


class FlowchartSemanticExportTests(unittest.TestCase):
    def test_export_writes_typed_graph_without_renderer_geometry(self):
        workflow = Workflow()

        parameters = Parameters(
            "Parameters",
            [
                {
                    "name": "x",
                    "description": "x parameter",
                    "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
                }
            ],
        )
        parameters.analyze_ios()

        workflow.add_module(parameters)
        workflow.parameter_module = parameters
        workflow.add_module(_FakeCalculatorModule())
        workflow.add_module(_FakeRelayModule())
        workflow.add_module(_FakeSelectionModule())
        workflow.resolve_dependencies()

        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "flowchart.json")
            graph = asyncio.run(
                workflow.export_flowchart_semantics(
                    save_path=out_path,
                    workflow_name="SemanticExportTest",
                )
            )

            self.assertTrue(os.path.exists(out_path))
            with open(out_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)

        self.assertEqual(payload["schema"], "jarvisplot.scene/v1")
        self.assertEqual(payload["scene_type"], "flowchart")
        self.assertEqual(payload["scene_id"], "workflow_main")
        self.assertEqual(payload["metadata"]["producer"], "Jarvis-HEP")
        self.assertIn("producer_version", payload["metadata"])
        self.assertEqual(payload["metadata"]["workflow_name"], "SemanticExportTest")
        self.assertIn("layers", payload)
        self.assertIn("nodes", payload)
        self.assertIn("edges", payload)
        self.assertEqual(graph["schema"], payload["schema"])

        node_ids = {node["id"] for node in payload["nodes"]}
        self.assertIn("Parameters", node_ids)
        self.assertIn("CalcA", node_ids)
        self.assertIn("CalcB", node_ids)
        self.assertIn("SelectedOperas", node_ids)
        self.assertIn("var::x", node_ids)
        self.assertIn("var::y", node_ids)
        self.assertIn("var::z", node_ids)

        layer_one = next(layer for layer in payload["layers"] if layer["id"] == "layer_1")
        self.assertEqual(layer_one["index"], 1)
        self.assertIn("Parameters", layer_one["nodes"])
        self.assertTrue(all(layer["id"].startswith("layer_") for layer in payload["layers"]))

        edges = payload["edges"]
        self.assertTrue(
            any(
                edge["source"]["node"] == "Parameters"
                and edge["target"]["node"] == "var::x"
                and edge["role"] == "parameterflow"
                for edge in edges
            )
        )
        self.assertTrue(
            any(
                edge["source"]["node"] == "var::x"
                and edge["target"]["node"] == "file::CalcA::input::input_file"
                and edge["role"] == "fileflow"
                for edge in edges
            )
        )
        self.assertTrue(
            any(
                edge["source"]["node"] == "file::CalcA::output::output_file"
                and edge["target"]["node"] == "var::y"
                and edge["role"] == "fileflow"
                for edge in edges
            )
        )
        self.assertTrue(any(edge["role"] == "dataflow" for edge in edges))

        selection_node = next(node for node in payload["nodes"] if node["id"] == "SelectedOperas")
        self.assertEqual(
            selection_node["selection"],
            {"expression": "y > 0", "variables": ["y"]},
        )
        self.assertTrue(
            any(
                port["id"] == "selection::y" and port["role"] == "selection"
                for port in selection_node["in_ports"]
            )
        )
        self.assertTrue(
            any(
                edge["source"]["node"] == "var::y"
                and edge["target"]["node"] == "SelectedOperas"
                and edge["target"]["port"] == "selection::y"
                and edge["role"] == "selectionflow"
                for edge in edges
            )
        )

        for node in payload["nodes"]:
            self.assertNotIn("bp", node)
            self.assertNotIn("width", node)
            self.assertNotIn("ihh", node)
            self.assertNotIn("ohh", node)
            self.assertNotIn("mhh", node)
            self.assertIn("layer", node)


if __name__ == "__main__":
    unittest.main()
