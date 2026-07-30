#!/usr/bin/env python3
"""D14 task YAML validation gate."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from copy import deepcopy

from jarvishep2.client import dispatch_validate, main, normalize_argv
from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_validation import (
    ConfigValidationError,
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
    }
    cfg = deepcopy(base)
    for key, value in overrides.items():
        cfg[key] = value
    return cfg


class TaskValidationKernelTests(unittest.TestCase):
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
                        "input": [{"name": "a", "type": "text", "save": True}],
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

    def test_main_validate_subcommand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Scan:\n  name: v\n"
                    "Sampling:\n"
                    "  Method: Random\n"
                    "  Point number: 10\n"
                    "  Variables:\n"
                    "    - name: x\n"
                    "      description: d\n"
                    "      distribution:\n"
                    "        type: Flat\n"
                    "        parameters: {min: 0.0, max: 1.0}\n"
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
            }
            report = validate_task_config(cfg)
            self.assertTrue(
                report.ok,
                msg=f"{name} failed: {[i.format_line() for i in report.errors()]}",
            )


if __name__ == "__main__":
    unittest.main()
