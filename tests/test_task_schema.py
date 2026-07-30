#!/usr/bin/env python3
"""JSON Schema coverage for the stable Jarvis2 task-card surface."""

from __future__ import annotations

import json
import os
import unittest

from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_schema import MANIFEST_PATH, SCHEMA_PATH, task_card_validator
from jarvishep2.task_validation import validate_task_config


def _card() -> dict:
    return {
        "Scan": {"name": "schema-card"},
        "Sampling": {
            "Method": "Random",
            "Point number": 2,
            "Variables": [
                {
                    "name": "x",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                }
            ],
        },
        "Calculators": {
            "Modules": [
                {
                    "name": "Calc",
                    "execution": {
                        "path": ".",
                        "commands": ["true"],
                        "input": [
                            {
                                "name": "input",
                                "path": "input.json",
                                "type": "JSON",
                                "actions": [
                                    {
                                        "type": "Dump",
                                        "variables": [{"name": "x", "expression": "x"}],
                                    }
                                ],
                            }
                        ],
                        "output": [
                            {
                                "name": "output",
                                "path": "output.json",
                                "type": "JSON",
                                "variables": [{"name": "z", "entry": "z"}],
                            }
                        ],
                    },
                }
            ]
        },
        "Operas": {
            "Modules": [
                {
                    "name": "Opera",
                    "operator": "package.module.function",
                    "call_mode": "call",
                    "input": [{"name": "x", "expression": "x"}],
                    "output": [{"name": "z", "entry": "z"}],
                }
            ]
        },
    }


class TaskCardSchemaTests(unittest.TestCase):
    def test_schema_is_bundled_and_compiles(self) -> None:
        self.assertTrue(SCHEMA_PATH.is_file())
        self.assertTrue(MANIFEST_PATH.is_file())
        self.assertEqual(list(task_card_validator().iter_errors(_card())), [])

    def test_manifest_dispatches_sampling_and_io_schemas(self) -> None:
        card = _card()
        del card["Sampling"]["Point number"]
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "UnknownFormat"
        report = validate_task_config(card)
        errors = report.errors()
        self.assertTrue(any(issue.code == "JV2-SCH-001" and issue.path == "$.Sampling" for issue in errors))
        self.assertTrue(any(issue.code == "JV2-SCH-002" and issue.path.endswith("input[0].type") for issue in errors))

    def test_manifest_gives_each_sampling_method_its_own_schema(self) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        method_schemas = manifest["sampling_methods"]
        self.assertEqual(len(method_schemas), len(set(method_schemas.values())))

    def test_json_input_requires_path(self) -> None:
        card = _card()
        del card["Calculators"]["Modules"][0]["execution"]["input"][0]["path"]
        report = validate_task_config(card)
        schema_errors = [issue for issue in report.errors() if issue.code == "JV2-SCH-001"]
        self.assertTrue(schema_errors)
        self.assertTrue(any("input[0]" in issue.path for issue in schema_errors))
        path_error = next(issue for issue in schema_errors if issue.path.endswith("input[0]"))
        self.assertEqual(path_error.suggestion, "Add 'path' at this location using the expected YAML type.")
        self.assertIn("path: input.json", path_error.example or "")

    def test_each_reported_error_has_a_correction_suggestion(self) -> None:
        card = _card()
        del card["Sampling"]["Point number"]
        card["Sampling"]["Variables"] = []
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "UnknownFormat"
        report = validate_task_config(card)
        self.assertTrue(report.errors())
        self.assertTrue(all(issue.suggestion for issue in report.errors()))
        random_error = next(issue for issue in report.errors() if issue.code == "JV2-MTH-020")
        self.assertIn("Point number", random_error.suggestion or "")
        self.assertIn("Point number: 100", random_error.example or "")

    def test_operas_module_rejects_unknown_user_key(self) -> None:
        card = _card()
        card["Operas"]["Modules"][0]["typo"] = True
        report = validate_task_config(card)
        self.assertTrue(any(issue.code == "JV2-SCH-001" for issue in report.errors()))

    def test_existing_parity_card_passes_schema_layer(self) -> None:
        path = os.path.join(os.path.dirname(__file__), "parity_project", "bridson_opera.yaml")
        config = load_task_yaml(path)
        report = validate_task_config(config)
        self.assertFalse(
            [issue for issue in report.errors() if issue.code == "JV2-SCH-001"],
            [issue.format_line() for issue in report.errors()],
        )


if __name__ == "__main__":
    unittest.main()
