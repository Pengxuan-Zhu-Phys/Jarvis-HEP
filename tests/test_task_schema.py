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
            "Bounds": {"point_number": 2},
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
        "EnvReqs": {},
    }


class TaskCardSchemaTests(unittest.TestCase):
    def test_schema_is_bundled_and_compiles(self) -> None:
        self.assertTrue(SCHEMA_PATH.is_file())
        self.assertTrue(MANIFEST_PATH.is_file())
        self.assertEqual(list(task_card_validator().iter_errors(_card())), [])
        self.assertEqual(schema_catalog_lint_errors(), [])

    def test_scan_name_is_required(self) -> None:
        card = _card()
        del card["Scan"]["name"]
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-SCH-001"
                and item.path == "$.Scan"
                and "required" in item.message
                for item in report.errors()
            )
        )

    def test_sampling_method_is_required_for_normal_cards(self) -> None:
        card = _card()
        del card["Sampling"]["Method"]
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-SCH-001"
                and item.path == "$.Sampling"
                and "Method" in item.message
                for item in report.errors()
            )
        )

    def test_envreqs_v2_rejects_unknown_top_level_key(self) -> None:
        card = _card()
        card["EnvReqs"] = {"V2": {"future_setting": True}}
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-SCH-001"
                and item.path == "$.EnvReqs.V2"
                and "future_setting" in item.message
                for item in report.errors()
            )
        )

    def test_envreqs_os_uses_v1_compatible_closed_list_items(self) -> None:
        card = _card()
        card["EnvReqs"] = {
            "OS": [
                {"name": "linux", "version": ">=3.10.0"},
                {"name": "Darwin", "version": ">=10.14"},
            ]
        }
        self.assertFalse(
            [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
        )

        card["EnvReqs"]["OS"][0]["typo"] = True
        self.assertTrue(
            [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
        )

    def test_check_card_may_omit_sampling_method(self) -> None:
        card = _card()
        card["Sampling"] = {}
        card["EnvReqs"] = {"V2": {"check_modules": {"data": "points.csv"}}}
        report = validate_task_config(card, check_modules=True)
        self.assertFalse(
            any(item.code == "JV2-SCH-001" and item.path == "$.Sampling" for item in report.errors()),
            [item.format_line() for item in report.errors()],
        )

    def test_calculator_or_operas_block_is_required(self) -> None:
        card = _card()
        del card["Calculators"]
        del card["Operas"]
        report = validate_task_config(card)
        self.assertTrue(
            any(item.code == "JV2-SCH-001" and item.path == "$" for item in report.errors())
        )

    def test_log_likelihood_has_closed_list_item_schema(self) -> None:
        for invalid in (
            [42],
            [{"name": "LogL"}],
            [{"expression": "x"}],
            [{"name": "LogL", "expression": "x", "typo": True}],
            [],
        ):
            with self.subTest(invalid=invalid):
                card = _card()
                card["Sampling"]["LogLikelihood"] = invalid
                report = validate_task_config(card)
                self.assertTrue(
                    [issue for issue in report.errors() if issue.code == "JV2-SCH-001"]
                )

    def test_log_likelihood_name_and_expression_are_valid(self) -> None:
        card = _card()
        card["Sampling"]["LogLikelihood"] = [
            {"name": "LogL_signal", "expression": "x"},
            {"name": "LogL", "expression": "LogL_signal"},
        ]
        self.assertFalse(
            [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
        )

    def test_json_input_dump_supports_direct_and_expression_variants(self) -> None:
        variants = (
            {"name": "x"},
            {"name": "x", "entry": "parameters.x"},
            {"name": "scaled_x", "expression": "2 * x"},
            {"name": "scaled_x", "expression": "2 * x", "entry": "parameters.x"},
        )
        for variable in variants:
            with self.subTest(variable=variable):
                card = _card()
                card["Calculators"]["Modules"][0]["execution"]["input"][0]["actions"][0]["variables"] = [variable]
                self.assertFalse(
                    [
                        issue
                        for issue in validate_task_config(card).errors()
                        if issue.code == "JV2-SCH-001"
                    ]
                )

    def test_json_input_dump_rejects_unknown_variable_shape(self) -> None:
        card = _card()
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["actions"][0]["variables"] = [
            {"name": "x", "value": 1.0}
        ]
        self.assertTrue(
            [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
        )

    def test_json_output_supports_root_and_entry_variants(self) -> None:
        for variable in ({"name": "z"}, {"name": "log_likelihood", "entry": "fit.log_likelihood"}):
            with self.subTest(variable=variable):
                card = _card()
                card["Calculators"]["Modules"][0]["execution"]["output"][0]["variables"] = [variable]
                self.assertFalse(
                    [
                        issue
                        for issue in validate_task_config(card).errors()
                        if issue.code == "JV2-SCH-001"
                    ]
                )

    def test_json_output_rejects_unknown_variable_shape(self) -> None:
        card = _card()
        card["Calculators"]["Modules"][0]["execution"]["output"][0]["variables"] = [
            {"name": "z", "expression": "z"}
        ]
        self.assertTrue(
            [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
        )

    def test_operas_root_rejects_unknown_key(self) -> None:
        card = _card()
        card["Operas"]["typo"] = True
        report = validate_task_config(card)
        self.assertTrue(
            any(issue.code == "JV2-SCH-001" and issue.path == "$.Operas" for issue in report.errors())
        )

    def test_removed_sampling_check_fields_are_rejected(self) -> None:
        for key, value in (
            ("mode", "check_modules"),
            ("data", "points.csv"),
            ("points_csv", "points.csv"),
        ):
            with self.subTest(key=key):
                card = _card()
                card["Sampling"][key] = value
                self.assertTrue(
                    [issue for issue in validate_task_config(card).errors() if issue.code == "JV2-SCH-001"]
                )

    def test_calculators_archiver_is_not_a_user_yaml_field(self) -> None:
        card = _card()
        card["Calculators"]["Archiver"] = {"mode": "process"}
        report = validate_task_config(card)
        self.assertTrue(
            any(
                issue.code == "JV2-SCH-001"
                and issue.path == "$.Calculators"
                and "Archiver" in issue.message
                for issue in report.errors()
            )
        )

    def test_scan_save_dir_is_not_a_v2_field(self) -> None:
        card = _card()
        card["Scan"]["save_dir"] = "outputs/custom"
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-SCH-001"
                and item.path == "$.Scan"
                and "save_dir" in item.message
                for item in report.errors()
            )
        )

    def test_scan_sample_directory_is_not_a_v2_field(self) -> None:
        card = _card()
        card["Scan"]["sample_directory"] = {"enabled": True, "limit": 10}
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-SCH-001"
                and item.path == "$.Scan"
                and "sample_directory" in item.message
                for item in report.errors()
            )
        )

    def test_manifest_dispatches_sampling_and_io_schemas(self) -> None:
        card = _card()
        del card["Sampling"]["Bounds"]["point_number"]
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

    def test_toymcmc_schema_and_same_scale_contract(self) -> None:
        card = _card()
        card["Sampling"] = {
            "Method": "ToyMCMC",
            "Bounds": {
                "num_chains": 4,
                "num_iters": 20,
                "proposal_scale": 0.2,
                "seed": 7,
            },
            "Variables": card["Sampling"]["Variables"],
        }
        report = validate_task_config(card)
        self.assertTrue(report.ok, report.issues)

        card["Sampling"]["Bounds"]["proposal_scale"] = [0.2, 0.3]
        report = validate_task_config(card)
        self.assertTrue(
            any(
                item.code == "JV2-MTH-054"
                and item.path == "Sampling.Bounds.proposal_scale"
                for item in report.errors()
            )
        )

        card["Sampling"]["Bounds"]["proposal_scale"] = [0.2, 0.2, 0.2]
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-055" for item in report.errors()))

        card["Sampling"]["Bounds"]["proposal_scale"] = 0.2
        card["Sampling"]["Bounds"]["seed"] = -1
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-056" for item in report.errors()))

    def test_ptmcmc_schema_and_ladder_contract(self) -> None:
        card = _card()
        card["Sampling"] = {
            "Method": "PTMCMC",
            "Bounds": {
                "num_chains": 4,
                "num_iters": 20,
                "proposal_scale": 0.2,
                "temperature_ladder": [1.0, 2.0, 4.0, 8.0],
                "exchange_interval": 2,
                "seed": 7,
            },
            "Variables": card["Sampling"]["Variables"],
        }
        report = validate_task_config(card)
        self.assertTrue(report.ok, report.issues)

        card["Sampling"]["Bounds"]["temperature_ladder"] = [1.0, 2.0]
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-067" for item in report.errors()))

        card["Sampling"]["Bounds"]["temperature_ladder"] = [1.5, 2.0, 3.0, 4.0]
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-068" for item in report.errors()))

        card["Sampling"]["Bounds"]["temperature_ladder"] = [1.0, 2.0, 2.0, 4.0]
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-069" for item in report.errors()))

        card["Sampling"]["Bounds"]["temperature_ladder"] = [1.0, 2.0, 4.0, 8.0]
        card["Sampling"]["Bounds"]["proposal_scale"] = [0.1, 0.2]
        report = validate_task_config(card)
        self.assertTrue(any(item.code == "JV2-MTH-065" for item in report.errors()))

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
        del card["Sampling"]["Bounds"]["point_number"]
        card["Sampling"]["Variables"] = []
        card["Calculators"]["Modules"][0]["execution"]["input"][0]["type"] = "UnknownFormat"
        report = validate_task_config(card)
        self.assertTrue(report.errors())
        self.assertTrue(all(issue.suggestion for issue in report.errors()))
        random_error = next(issue for issue in report.errors() if issue.code == "JV2-MTH-020")
        self.assertIn("point_number", random_error.suggestion or "")
        self.assertIn("point_number: 100", random_error.example or "")

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
        card["Sampling"].update({"Method": "Bridson", "Bounds": {"radius": "abc", "max_attempt": 30}})
        report = validate_task_config(card)
        issue = next(
            item
            for item in report.errors()
            if item.path == "$.Sampling.Bounds.radius" and item.code == "JV2-SCH-001"
        )
        self.assertEqual(
            issue.message,
            "expected a number (e.g. 0.05 or 1.0e-5), got the string 'abc'",
        )
        self.assertIn("0.05 or 1.0e-5", issue.suggestion or "")

    def test_numeric_string_warns_but_compatibly_passes(self) -> None:
        card = _card()
        card["Sampling"].update({"Method": "Bridson", "Bounds": {"radius": "1e-1", "max_attempt": "30"}})
        report = validate_task_config(card)
        self.assertFalse(
            [issue for issue in report.errors() if issue.path == "$.Sampling.Bounds.radius"]
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

    def test_removed_yaml_aliases_are_rejected(self) -> None:
        mutations = {
            "Sampling.Bounds.Seed": lambda card: card["Sampling"]["Bounds"].__setitem__("Seed", 2),
            "Sampling.Bounds.Point number": lambda card: card["Sampling"]["Bounds"].__setitem__("Point number", 2),
            "Sampling.Nuisance.TargetMode": lambda card: card["Sampling"].__setitem__("Nuisance", {"TargetMode": "min"}),
            "Operas.Modules.timeout": lambda card: card["Operas"]["Modules"][0].__setitem__("timeout", 1),
            "Calculators.Modules.deps_source": lambda card: card["Calculators"]["Modules"][0].__setitem__("deps_source", "x"),
            "Calculators.Modules.make_paraller": lambda card: card["Calculators"]["Modules"][0].__setitem__("make_paraller", 1),
            "top-level Runtime": lambda card: card.__setitem__("Runtime", {}),
            "top-level Likelihood": lambda card: card.__setitem__("Likelihood", []),
        }
        for path, mutate in mutations.items():
            card = _card()
            mutate(card)
            self.assertTrue(
                validate_task_config(card).errors(),
                msg=f"removed YAML alias unexpectedly accepted: {path}",
            )

    def test_boolean_where_string_is_required_suggests_quoting(self) -> None:
        card = _card()
        card["Scan"]["name"] = True
        report = validate_task_config(card)
        problem = next(issue for issue in report.errors() if issue.path == "$.Scan.name")
        self.assertIn("Quote YAML boolean-like", problem.suggestion or "")

    def test_unmigrated_rltpmcmc_blocks_are_rejected(self) -> None:
        for field in ("Control", "Diagnostics", "PPO", "Reward"):
            with self.subTest(field=field):
                card = _card()
                card["Sampling"][field] = {"future_key": 1}
                report = validate_task_config(card)
                self.assertTrue(
                    any(
                        item.code == "JV2-SCH-001"
                        and item.path == "$.Sampling"
                        and field in item.message
                        for item in report.errors()
                    )
                )

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
