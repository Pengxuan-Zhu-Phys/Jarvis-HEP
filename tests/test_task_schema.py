#!/usr/bin/env python3
"""JSON Schema coverage for the stable Jarvis task-card surface."""

from __future__ import annotations

import json
import os
import unittest
from unittest.mock import patch

from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_schema import (
    MANIFEST_PATH, SCHEMA_PATH, schema_catalog_lint_errors, task_card_validator,
)
from jarvishep2.task_validation import validate_task_config
from jarvishep2.io_portal import available_io_formats
from jarvishep2.distributor import Distributor


def _card() -> dict:
    return {
        "Scan": {"name": "schema-card"},
        "Sampling": {
            "Method": "Random",
            "Bounds": {"Point number": 2},
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
        self.assertEqual(schema_catalog_lint_errors(), [])

    def test_manifest_dispatches_sampling_and_io_schemas(self) -> None:
        card = _card()
        del card["Sampling"]["Bounds"]["Point number"]
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "UnknownFormat"
        report = validate_task_config(card)
        errors = report.errors()
        self.assertTrue(
            any(
                issue.code == "JV2-SCH-001"
                and (
                    issue.path == "$.Sampling.Bounds"
                    or issue.path == "$.Sampling"
                    or "Bounds" in issue.path
                )
                for issue in errors
            )
        )
        self.assertTrue(any(issue.code == "JV2-SCH-002" and issue.path.endswith("input[0].type") for issue in errors))

    def test_manifest_gives_each_sampling_method_its_own_schema(self) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        method_schemas = manifest["sampling_methods"]
        self.assertEqual(len(method_schemas), len(set(method_schemas.values())))
        self.assertEqual(set(method_schemas), set(Distributor.available_methods()))

    def test_manifest_matches_builtin_portal_formats(self) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        self.assertTrue(set(manifest["io"]["input"]) <= set(available_io_formats("input")))
        self.assertTrue(set(manifest["io"]["output"]) <= set(available_io_formats("output")))

    def test_portal_format_without_bundled_schema_is_accepted(self) -> None:
        card = _card()
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "FuturePortalFormat"
        real_formats = available_io_formats("input")
        with patch("jarvishep2.io_portal.available_io_formats", return_value=[*real_formats, "FuturePortalFormat"]):
            report = validate_task_config(card)
        self.assertFalse([issue for issue in report.errors() if issue.code == "JV2-SCH-002"])

    def test_every_portal_format_is_accepted_by_its_direction_schema(self) -> None:
        with MANIFEST_PATH.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        for direction, formats in manifest["io"].items():
            for format_name in formats:
                card = _card()
                execution = card["Calculators"]["Modules"][0]["execution"]
                execution[direction] = [{
                    "name": direction,
                    "path": f"{direction}.data",
                    "type": format_name,
                }]
                report = validate_task_config(card)
                self.assertFalse(
                    [issue for issue in report.errors() if issue.code.startswith("JV2-SCH")],
                    f"{direction} {format_name}: {[issue.format_line() for issue in report.errors()]}",
                )

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
        del card["Sampling"]["Bounds"]["Point number"]
        card["Sampling"]["Variables"] = []
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "UnknownFormat"
        report = validate_task_config(card)
        self.assertTrue(report.errors())
        self.assertTrue(all(issue.suggestion for issue in report.errors()))
        random_error = next(issue for issue in report.errors() if issue.code == "JV2-MTH-020")
        self.assertIn("Point number", random_error.suggestion or "")
        self.assertIn("Point number: 100", random_error.example or "")

    def test_strict_user_blocks_reject_typos(self) -> None:
        card = _card()
        card["Scan"]["naem"] = "typo"
        card["Sampling"]["Sead"] = 1
        card["Calculators"]["Modules"][0]["clone_shdow"] = True
        report = validate_task_config(card)
        paths = {issue.path for issue in report.errors() if issue.code == "JV2-SCH-001"}
        self.assertIn("$.Scan", paths)
        self.assertIn("$.Sampling", paths)
        self.assertIn("$.Calculators.Modules[0]", paths)
        suggestion = next(issue.suggestion for issue in report.errors() if issue.path == "$.Scan")
        self.assertIn("Did you mean 'name'?", suggestion or "")

    def test_multiple_closed_block_typos_each_get_a_spelling_hint(self) -> None:
        card = _card()
        card["Sampling"].update({
            "Method": "AdaptiveBridson",
            "Bounds": {
                "target_expression": "x",
                "target_value": 0.5,
                "initial_radiuz": 0.1,
                "max_pointz": 10,
            },
        })
        report = validate_task_config(card)
        issue = next(
            item
            for item in report.errors()
            if item.path in {"$.Sampling.Bounds", "$.Sampling.Bounds"}
            or item.path.startswith("$.Sampling.Bounds")
        )
        self.assertIn("'initial_radiuz' (did you mean 'initial_radius'?)", issue.suggestion or "")
        self.assertIn("'max_pointz' (did you mean 'max_points'?)", issue.suggestion or "")
        self.assertIn("... (+", issue.suggestion or "")

    def test_numeric_union_error_is_actionable(self) -> None:
        card = _card()
        card["Sampling"].update({"Method": "Bridson", "Bounds": {"Radius": "abc", "MaxAttempt": 30}})
        report = validate_task_config(card)
        issue = next(
            item
            for item in report.errors()
            if item.path == "$.Sampling.Bounds.Radius" and item.code == "JV2-SCH-001"
        )
        self.assertEqual(
            issue.message,
            "expected a number (e.g. 0.05 or 1.0e-5), got the string 'abc'",
        )
        self.assertIn("0.05 or 1.0e-5", issue.suggestion or "")

    def test_numeric_string_warns_but_compatibly_passes(self) -> None:
        card = _card()
        card["Sampling"].update({"Method": "Bridson", "Bounds": {"Radius": "1e-1", "MaxAttempt": "30"}})
        report = validate_task_config(card)
        self.assertFalse(
            [issue for issue in report.errors() if issue.path == "$.Sampling.Bounds.Radius"]
        )
        self.assertTrue(any(issue.code == "JV2-SCH-003" for issue in report.warnings()))

    def test_retired_sampling_top_level_knobs_are_rejected(self) -> None:
        card = _card()
        card["Sampling"]["Radius"] = 0.1
        report = validate_task_config(card)
        self.assertTrue(
            any(
                issue.code == "JV2-MTH-001" and "Radius" in issue.path
                for issue in report.errors()
            )
            or any(
                issue.code == "JV2-SCH-001" and issue.path == "$.Sampling"
                for issue in report.errors()
            )
        )

    def test_boolean_where_string_is_required_suggests_quoting(self) -> None:
        card = _card()
        card["Scan"]["name"] = True
        report = validate_task_config(card)
        problem = next(issue for issue in report.errors() if issue.path == "$.Scan.name")
        self.assertIn("Quote YAML boolean-like", problem.suggestion or "")

    def test_open_legacy_sampling_zone_warns_once(self) -> None:
        card = _card()
        card["Sampling"].update({"Control": {"future_key": 1}, "PPO": {"future_key": 2}})
        report = validate_task_config(card)
        self.assertEqual([item.code for item in report.warnings()].count("JV2-SCH-004"), 1)

    def test_adaptive_bridson_lowercase_block_alias_is_valid(self) -> None:
        card = _card()
        card["Sampling"].update({
            "Method": "AdaptiveBridson",
            "Bounds": {"target_expression": "x", "target_value": 0.5},
        })
        report = validate_task_config(card)
        self.assertFalse([item for item in report.errors() if item.code == "JV2-SCH-001"])

    def test_method_and_io_dispatch_match_runtime_whitespace_normalization(self) -> None:
        card = _card()
        card["Sampling"]["Method"] = " Random "
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = " JSON "
        report = validate_task_config(card)
        self.assertFalse([issue for issue in report.errors() if issue.code.startswith("JV2-SCH")])

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

    def test_check_modules_parity_card_passes_closed_block_schema(self) -> None:
        path = os.path.join(os.path.dirname(__file__), "parity_project", "check_modules.yaml")
        report = validate_task_config(load_task_yaml(path), check_modules=True)
        self.assertFalse(
            [issue for issue in report.errors() if issue.code.startswith("JV2-SCH")],
            [issue.format_line() for issue in report.errors()],
        )


if __name__ == "__main__":
    unittest.main()
