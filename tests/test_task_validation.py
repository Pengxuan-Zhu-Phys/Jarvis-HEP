#!/usr/bin/env python3
"""D14 task YAML validation gate."""

from __future__ import annotations

import contextlib
import io
import json
import os
import tempfile
import unittest
from copy import deepcopy

from jarvishep2.client import dispatch_run, dispatch_validate, main, normalize_argv
from jarvishep2.task_config import TaskCardLoadError, load_task_yaml
from jarvishep2.task_validation import (
    ConfigValidationError,
    format_report,
    validate_or_raise,
    validate_task_config,
)


def _minimal_dynesty_config(**overrides) -> dict:
    base = {
        "Scan": {"name": "valtest"},
        "Sampling": {
            "Method": "Dynesty",
            "Variables": [
                {
                    "name": "x",
                    "description": "dim",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                }
            ],
            "Bounds": {
                "nlive": 50,
                "rseed": 1,
                "dlogz": 0.5,
                "run_nested": {"print_progress": True},
            },
        },
        "EnvReqs": {"V2": {"workers": 2, "batch_size": 16}},
        "Calculators": {},
    }
    cfg = deepcopy(base)
    for key, value in overrides.items():
        cfg[key] = value
    return cfg


class TaskValidationKernelTests(unittest.TestCase):
    def test_many_errors_render_reference_table_before_details(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Scan"]["naem"] = "typo"
        cfg["Sampling"]["Variables"] = []
        report = validate_task_config(cfg)
        rendered = format_report(report)
        self.assertIn("Error summary (reference only):", rendered)
        self.assertIn("| Code", rendered)
        self.assertIn("Detailed diagnostics:", rendered)
        self.assertIn("This table is only a quick reference.", rendered)

    def test_empty_summary_table_is_safe_and_ascii_only(self) -> None:
        from jarvishep2.task_validation import _format_issue_summary_table
        rendered = _format_issue_summary_table([])
        self.assertIn("| YAML path", rendered)
        self.assertTrue(all(ord(char) < 128 for char in rendered))

    def test_non_ascii_key_is_an_encoding_error_before_schema(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Scan"]["\u540d\u79f0"] = "scan"
        report = validate_task_config(cfg)
        self.assertEqual([issue.code for issue in report.errors()], ["JV2-ENC-001"])
        issue = report.errors()[0]
        self.assertIn("$.Scan.\\u540d\\u79f0", issue.path)
        self.assertIn("position(s) 1, 2", issue.message)
        self.assertIn("Chinese text in a # comment", issue.hint or "")

    def test_non_ascii_scan_name_and_description_are_encoding_errors(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Scan"]["name"] = "\u6697\u7269\u8d28"
        cfg["Sampling"]["Variables"][0]["description"] = "\u4e2d\u6587\u8bf4\u660e"
        report = validate_task_config(cfg)
        self.assertEqual([issue.code for issue in report.errors()], ["JV2-ENC-001", "JV2-ENC-001"])
        self.assertEqual({issue.path for issue in report.errors()}, {
            "$.Sampling.Variables[0].description", "$.Scan.name",
        })
        self.assertTrue(all(ord(char) < 128 for char in format_report(report)))

    def test_non_ascii_comment_is_discarded_before_encoding_validation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "comment.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "# \u8fd9\u662f\u5141\u8bb8\u7684\u4e2d\u6587\u6ce8\u91ca\n"
                    "Scan:\n  name: ascii-scan\n"
                    "Sampling:\n  Method: Random\n"
                    "  Bounds:\n    Point number: 1\n"
                    "  Variables:\n    - name: x\n"
                    "      distribution:\n        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                    "EnvReqs: {}\n"
                    "Calculators: {}\n"
                )
            report = validate_task_config(load_task_yaml(path))
        self.assertFalse([issue for issue in report.errors() if issue.code == "JV2-ENC-001"])

    def test_cyrillic_homoglyph_is_an_encoding_error(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Scan"]["n\u0430me"] = "scan"  # Cyrillic a, not ASCII a.
        report = validate_task_config(cfg)
        self.assertEqual([issue.code for issue in report.errors()], ["JV2-ENC-001"])

    def test_non_ascii_key_keeps_unrelated_schema_errors_in_one_report(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Scan"].update({"\u540d\u79f0": "scan", "naem": "scan"})
        report = validate_task_config(cfg)
        encoding = [item for item in report.errors() if item.code == "JV2-ENC-001"]
        schema = [item for item in report.errors() if item.code == "JV2-SCH-001"]
        self.assertEqual(len(encoding), 1)
        self.assertEqual(len(schema), 1)
        self.assertEqual(schema[0].path, "$.Scan")
        self.assertIn("Did you mean 'name'?", schema[0].suggestion or "")
        self.assertNotIn("\\u540d\\u79f0", schema[0].message)

    def test_root_block_typos_are_closed_and_actionable(self) -> None:
        expected = {
            "Calculater": "Calculators",
            "Opera": "Operas",
            "EnvReq": "EnvReqs",
            "Scam": "Scan",
        }
        for typo, canonical in expected.items():
            with self.subTest(typo=typo):
                cfg = _minimal_dynesty_config()
                cfg["Calculators"] = {"Modules": []}
                cfg["Operas"] = {"Modules": []}
                value = cfg.pop(canonical)
                cfg[typo] = value
                schema = [
                    item for item in validate_task_config(cfg).errors()
                    if item.code == "JV2-SCH-001"
                ]
                matching = [
                    item for item in schema
                    if typo in item.message
                    and f"Did you mean {canonical!r}?" in (item.suggestion or "")
                ]
                self.assertEqual(len(matching), 1)
                self.assertEqual(matching[0].path, "$")

    def test_libdeps_is_a_declared_closed_top_level_block(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["LibDeps"] = {
            "path": "deps/library",
            "make_paraller": 2,
            "Modules": [{
                "name": "Tool",
                "installed": False,
                "installation": {
                    "path": "deps/library/Tool",
                    "source": "deps/source/tool.tar.gz",
                    "commands": ["make"],
                },
            }],
            "registered_executables": [{
                "name": "tool",
                "source": "deps/library/Tool/tool",
                "resolution": "direct_path",
            }],
        }
        self.assertFalse([item for item in validate_task_config(cfg).errors() if item.code == "JV2-SCH-001"])

        cfg["LibDeps"]["Moduels"] = []
        schema = [item for item in validate_task_config(cfg).errors() if item.code == "JV2-SCH-001"]
        self.assertEqual(len(schema), 1)
        self.assertIn("Did you mean 'Modules'?", schema[0].suggestion or "")

    def test_removed_v1_root_blocks_explain_the_required_migration(self) -> None:
        expected = {
            "Likelihood": "Move its expressions to Sampling.LogLikelihood",
            "Mapper": "not a V2 task-card interface",
            "project_name": "not a V2 task-card interface",
            "Utils": "Utils.interpolations_1D moved to Jarvis-Operas",
        }
        for key, expected_suggestion in expected.items():
            with self.subTest(key=key):
                cfg = _minimal_dynesty_config()
                cfg[key] = {} if key != "project_name" else "legacy-name"
                schema = [
                    item for item in validate_task_config(cfg).errors()
                    if item.code == "JV2-SCH-001"
                ]
                self.assertEqual(len(schema), 1)
                self.assertEqual(schema[0].path, "$")
                self.assertIn(expected_suggestion, schema[0].suggestion or "")

    def test_valid_dynesty_passes(self) -> None:
        report = validate_task_config(_minimal_dynesty_config())
        self.assertTrue(report.ok)
        self.assertEqual(report.errors(), [])

    def test_unknown_method_errors(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Method"] = "Dynestyy"
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        codes = {i.code for i in report.errors()}
        self.assertIn("JV2-MTH-003", codes)

    def test_empty_method_errors(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Method"] = "  "
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-MTH-002" for i in report.errors()))

    def test_distribution_type_case_sensitive(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Variables"][0]["distribution"]["type"] = "flat"
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-VAR-002" for i in report.errors()))

    def test_nameless_variable_not_dropped(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Variables"][0]["name"] = ""
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-VAR-012" for i in report.errors()))

    def test_missing_flat_params(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Variables"][0]["distribution"]["parameters"] = {"min": 0.0}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-VAR-031" for i in report.errors()))

    def test_min_not_less_than_max(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Variables"][0]["distribution"]["parameters"] = {
            "min": 1.0,
            "max": 0.0,
        }
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-VAR-034" for i in report.errors()))

    def test_unknown_bounds_key(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Bounds"]["live"] = 50
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-BND-001" for i in report.errors()))

    def test_nlive_invalid(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Bounds"]["nlive"] = "auto"
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-BND-010" for i in report.errors()))

    def test_bounds_dynamic_key_rejected_for_any_nested_method(self) -> None:
        for method in ("Dynesty", "MultiNest"):
            cfg = _minimal_dynesty_config()
            cfg["Sampling"]["Method"] = method
            cfg["Sampling"]["Bounds"]["dynamic"] = True
            report = validate_task_config(cfg)
            self.assertFalse(report.ok, method)
            self.assertTrue(
                any(i.code == "JV2-BND-012" for i in report.errors()),
                method,
            )

    def test_archiver_unknown_mode(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Calculators"] = {"Archiver": {"mode": "turbo"}}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-ARC-010" for i in report.errors()))

    def test_archiver_unknown_key(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Calculators"] = {"Archiver": {"mode": "process", "async_io": True}}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-ARC-002" for i in report.errors()))

    def test_archiver_shard_limit_is_validated(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Calculators"] = {"Archiver": {"max_hdf5_bytes": 0}}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-ARC-014" for i in report.errors()))

    def test_workers_illegal_present_value(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["EnvReqs"] = {"V2": {"workers": "many", "batch_size": 16}}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-ENV-010" for i in report.errors()))

    def test_checkpoint_heartbeat_must_be_a_positive_duration(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["EnvReqs"]["V2"]["checkpoint_heartbeat_sec"] = 0
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-ENV-012" for i in report.errors()))

    def test_nested_run_nested_rejects_hep_owned_checkpoint_keys(self) -> None:
        """Checkpoint interval is EnvReqs.V2 only — not Sampling.Bounds.run_nested."""
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Bounds"]["run_nested"] = {
            "dlogz": 0.5,
            "checkpoint_every": 15,
            "checkpoint_file": "foo.pkl",
            "resume": True,
        }
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-BND-021" for i in report.errors()))

    def test_bridson_requires_length_and_radius(self) -> None:
        cfg = {
            "Scan": {"name": "b"},
            "Sampling": {
                "Method": "Bridson",
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0.0, "max": 1.0},
                        },
                    },
                    {
                        "name": "y",
                        "description": "y",
                        "distribution": {
                            "type": "Flat",
                            "parameters": {"min": 0.0, "max": 1.0, "length": 1.0},
                        },
                    },
                ],
            },
        }
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        codes = {i.code for i in report.errors()}
        self.assertIn("JV2-VAR-031", codes)  # missing length on x
        self.assertIn("JV2-MTH-010", codes)  # missing Radius

    def test_check_modules_accepts_csv_only_card(self) -> None:
        cfg = {
            "Scan": {"name": "cm"},
            "Sampling": {
                "mode": "check_modules",
                "Method": "NotARealMethod",
                "data": "points.csv",
            },
            "EnvReqs": {},
            "Calculators": {},
        }
        # Leftover invalid Method is OK when CSV path is configured.
        report = validate_task_config(cfg, check_modules=True)
        self.assertTrue(report.ok)
        report2 = validate_task_config(cfg)
        self.assertTrue(report2.ok)

    def test_check_modules_mode_without_method_ok_via_default_data(self) -> None:
        """Default EnvReqs.V2.check_modules.data is always available as a path spec."""
        cfg = {
            "Scan": {"name": "cm"},
            "Sampling": {"mode": "check_modules"},
            "EnvReqs": {},
            "Calculators": {},
        }
        report = validate_task_config(cfg)
        self.assertTrue(report.ok, [i.format_line() for i in report.errors()])

    def test_check_flag_on_scan_card_allows_method_smoke(self) -> None:
        """CLI `check` on a normal Dynesty/MultiNest card may use sampler smoke points."""
        cfg = _minimal_dynesty_config()
        report = validate_task_config(cfg, check_modules=True)
        self.assertTrue(report.ok, [i.format_line() for i in report.errors()])

    def test_io_save_key_is_valid_in_strict_mode(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Calculators"] = {
            "Modules": [
                {
                    "name": "X",
                    "execution": {
                        "path": ".",
                        "commands": ["true"],
                        "input": [{"name": "a", "path": "input.txt", "type": "Text", "save": True}],
                        "output": [],
                    },
                }
            ]
        }
        report = validate_task_config(cfg)
        self.assertTrue(report.ok)
        self.assertFalse(report.warnings())
        strict = validate_task_config(cfg, strict=True)
        self.assertTrue(strict.ok)

    def test_raise_if_errors(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Sampling"]["Method"] = "Nope"
        with self.assertRaises(ConfigValidationError) as ctx:
            validate_or_raise(cfg)
        self.assertIn("JV2-MTH-003", str(ctx.exception))

    def test_operas_call_mode(self) -> None:
        cfg = _minimal_dynesty_config()
        cfg["Operas"] = {"Modules": [{"name": "T", "operator": "x", "call_mode": "maybe"}]}
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        self.assertTrue(any(i.code == "JV2-OPR-001" for i in report.errors()))

    def test_author_prose_keeps_unicode_in_format_line(self) -> None:
        """D23.14: suggestion/example/hint keep Unicode; path/message still escape users."""
        from jarvishep2.task_validation import issue

        item = issue(
            "error",
            "JV2-MAP-007",
            "Sampling.Mapper[0]",
            "cannot compile Mapper expression",
            hint="See Operas docs → constants",
            suggestion="If this uses pdg.mZ, math.add, … check registration.",
            example='expression: "pdg.mZ * cos(t)"  # ≥ 0',
        )
        line = item.format_line()
        self.assertIn("…", line)
        self.assertIn("→", line)
        self.assertIn("≥", line)
        self.assertNotIn("\\u2026", line)
        self.assertNotIn("\\u2192", line)
        self.assertNotIn("\\u2265", line)
        # Path/message still escape user-provided non-ASCII.
        weird = issue("error", "JV2-X", "Sampling.变量", "value '测' is invalid")
        rendered = weird.format_line()
        self.assertIn("\\u", rendered)
        self.assertNotIn("测", rendered)

    def test_operas_constant_as_module_operator_rejected(self) -> None:
        """D23.4: operator must not name a namespace constant (JV2-OPR-002)."""
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas not installed")
        constants = build_constant_dicts()
        if "pdg.mZ" not in constants:
            self.skipTest("pdg.mZ not registered")

        cfg = _minimal_dynesty_config()
        cfg["Operas"] = {
            "Modules": [
                {
                    "name": "Const",
                    "operator": "pdg.mZ",
                    "call_mode": "call",
                    "input": [],
                    "output": [{"name": "m", "entry": "m"}],
                }
            ]
        }
        report = validate_task_config(cfg)
        self.assertFalse(report.ok)
        errors = [i for i in report.errors() if i.code == "JV2-OPR-002"]
        self.assertEqual(len(errors), 1)
        self.assertIn("pdg.mZ", errors[0].message)
        self.assertIn("expression", errors[0].message)
        # Guidance should point at bare constant form.
        suggestion = (errors[0].suggestion or "").lower()
        example = errors[0].example or ""
        self.assertTrue(
            "constant" in suggestion or "pdg" in example,
            msg=f"suggestion={errors[0].suggestion!r} example={example!r}",
        )


class TaskValidationCliTests(unittest.TestCase):
    def test_normalize_argv_validate(self) -> None:
        self.assertEqual(
            normalize_argv(["task.yaml", "--validate", "--strict"]),
            ["validate", "task.yaml", "--strict"],
        )

    def test_dispatch_validate_json_ok(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: v\n"
                    "Sampling:\n"
                    "  Method: Dynesty\n"
                    "  Variables:\n"
                    "    - name: x\n"
                    "      description: d\n"
                    "      distribution:\n"
                    "        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                    "  Bounds:\n"
                    "    nlive: 20\n"
                    "EnvReqs: {}\n"
                    "Calculators: {}\n"
                )
            code = dispatch_validate(path, as_json=True)
            self.assertEqual(code, 0)

    def test_dispatch_validate_fails_on_typo(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: v\n"
                    "Sampling:\n"
                    "  Method: Dynesty\n"
                    "  Variables:\n"
                    "    - name: x\n"
                    "      description: d\n"
                    "      distribution:\n"
                    "        type: flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                )
            code = dispatch_validate(path, as_json=True)
            self.assertEqual(code, 2)

    def test_dispatch_validate_json_contains_actionable_guidance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("Sampling:\n  Method: Random\n  Variables: []\n")
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = dispatch_validate(path, as_json=True)
            self.assertEqual(code, 2)
            issues = json.loads(stdout.getvalue())["issues"]
            random_issue = next(item for item in issues if item["code"] == "JV2-MTH-020")
            self.assertIn("Point number", random_issue["suggestion"])
            self.assertIn("Point number: 100", random_issue["example"])

    def test_dispatch_validate_reports_yaml_syntax_with_guidance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "broken.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("Sampling:\n  Method: Random\n    Point number: 10\n")
            stdout = io.StringIO()
            with contextlib.redirect_stdout(stdout):
                code = dispatch_validate(path, as_json=True)
            self.assertEqual(code, 2)
            issue = json.loads(stdout.getvalue())["issues"][0]
            self.assertEqual(issue["code"], "JV2-YAML-001")
            self.assertIn("indentation", issue["suggestion"])
            self.assertIn("Method: Random", issue["example"])

    def test_dispatch_run_reports_yaml_syntax_with_guidance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "broken.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write("Sampling:\n  Method: Random\n    Point number: 10\n")
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = dispatch_run(path)
            self.assertEqual(code, 2)
            self.assertIn("JV2-YAML-001", stderr.getvalue())
            self.assertIn("suggestion:", stderr.getvalue())

    def test_dispatch_run_rejects_non_ascii_before_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: \u6697\u7269\u8d28\n"
                    "Sampling:\n  Method: Random\n"
                    "  Bounds:\n    Point number: 1\n"
                    "  Variables:\n    - name: x\n"
                    "      distribution:\n        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                    "EnvReqs: {}\n"
                    "Calculators: {}\n"
                )
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = dispatch_run(path)
            self.assertEqual(code, 2)
            self.assertEqual(stderr.getvalue().count("JV2-ENC-001"), 1)
            self.assertFalse(os.path.exists(os.path.join(tmp, "outputs")))

    def test_dispatch_run_rejects_top_level_calculator_typo_before_creating_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: typo-scan\n"
                    "Sampling:\n  Method: Random\n"
                    "  Bounds:\n    Point number: 1\n"
                    "  Variables:\n    - name: x\n"
                    "      distribution:\n        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                    "Calculater:\n  Modules: []\n"
                    "EnvReqs: {}\n"
                )
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                code = dispatch_run(path)
            self.assertEqual(code, 2)
            self.assertIn("JV2-SCH-001", stderr.getvalue())
            self.assertIn("Did you mean 'Calculators'?", stderr.getvalue())
            self.assertFalse(os.path.exists(os.path.join(tmp, "outputs")))

    def test_referenced_default_yaml_syntax_has_the_same_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            task_path = os.path.join(tmp, "task.yaml")
            default_path = os.path.join(tmp, "environment.yaml")
            with open(default_path, "w", encoding="utf-8") as handle:
                handle.write("EnvReqs:\n  V2:\n    workers: 2\n      batch_size: 1\n")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Check_default_dependencies:\n"
                    "    required: true\n"
                    "    default_yaml_path: environment.yaml\n"
                )
            with self.assertRaises(TaskCardLoadError) as ctx:
                load_task_yaml(task_path)
            self.assertEqual(ctx.exception.code, "JV2-YAML-001")
            self.assertIn("default environment YAML", str(ctx.exception))

    def test_referenced_runtime_yaml_syntax_has_the_same_diagnostic(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            task_path = os.path.join(tmp, "task.yaml")
            default_path = os.path.join(tmp, "runtime.yaml")
            with open(default_path, "w", encoding="utf-8") as handle:
                handle.write("Runtime:\n  workers: 2\n    batch_size: 1\n")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Runtime:\n"
                    "    default_runtime_settings: runtime.yaml\n"
                )
            with self.assertRaises(TaskCardLoadError) as ctx:
                load_task_yaml(task_path)
            self.assertEqual(ctx.exception.code, "JV2-YAML-001")
            self.assertIn("runtime default YAML", str(ctx.exception))

    def test_main_validate_subcommand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: v\n"
                    "Sampling:\n"
                    "  Method: Random\n"
                    "  Bounds:\n"
                    "    Point number: 10\n"
                    "  Variables:\n"
                    "    - name: x\n"
                    "      description: d\n"
                    "      distribution:\n"
                    "        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                    "EnvReqs: {}\n"
                    "Calculators: {}\n"
                )
            code = main(["validate", path, "--json"])
            self.assertEqual(code, 0)

    def test_load_task_yaml_then_gate_on_core(self) -> None:
        from jarvishep2.core import Jarvis2Core

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: v\n"
                    "Sampling:\n"
                    "  Method: Dynestyy\n"
                    "  Variables:\n"
                    "    - name: x\n"
                    "      description: d\n"
                    "      distribution:\n"
                    "        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
                )
            core = Jarvis2Core()
            with self.assertRaises(ConfigValidationError):
                core.load_task_yaml(path)


class SamplingTemplateGateTests(unittest.TestCase):
    """Golden sampling templates must pass when embedded in a minimal task."""

    TEMPLATE_DIR = os.path.join(
        os.path.dirname(__file__),
        "..",
        "jarvishep2",
        "project_template",
        "bin",
        "sampling",
    )

    def test_simple_templates_pass(self) -> None:
        import yaml

        for name in (
            "Sampling_Dynesty_Simple.yaml",
            "Sampling_MultiNest_Simple.yaml",
        ):
            path = os.path.join(self.TEMPLATE_DIR, name)
            if not os.path.isfile(path):
                self.skipTest(f"missing template {name}")
            with open(path, "r", encoding="utf-8") as handle:
                fragment = yaml.safe_load(handle)
            cfg = {
                "Scan": {"name": "tpl"},
                "Sampling": fragment.get("Sampling") or fragment,
                "Calculators": {},
                "Operas": {},
                "EnvReqs": {},
                "LibDeps": {},
            }
            report = validate_task_config(cfg)
            self.assertTrue(
                report.ok,
                msg=f"{name} failed: {[i.format_line() for i in report.errors()]}",
            )


if __name__ == "__main__":
    unittest.main()
