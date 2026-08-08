#!/usr/bin/env python3
"""Acceptance tests for Jarvis man (DESIGN_YAML_MAN_2.0, D24)."""

from __future__ import annotations

import io
import json
import re
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from jarvishep2.Sampling.dynesty_sampler import NESTED_CONSTRUCTOR_USER_KEYS
from jarvishep2.client import main
from jarvishep2.contracts.nested import BOUNDS_TOP_ALLOWED
from jarvishep2.contracts.variables import PARAMS_REQUIRED
from jarvishep2.core import Jarvis2Core
from jarvishep2.man import normalize_yaml_path, resolve_man_request
from jarvishep2.man_codes import known_diagnostic_codes, man_command_for, man_target_for_code
from jarvishep2.task_schema import schema_by_id, schema_manifest
from jarvishep2.task_validation import issue, validate_task_config


class ManPathTests(unittest.TestCase):
    def test_normalize_path_prefixes(self) -> None:
        self.assertEqual(normalize_yaml_path("Sampling.Bounds"), "$.Sampling.Bounds")
        self.assertEqual(normalize_yaml_path("$.Sampling.Bounds"), "$.Sampling.Bounds")
        self.assertEqual(
            normalize_yaml_path("$.Calculators.Modules[0].execution.input[0]"),
            "$.Calculators.Modules[].execution.input[]",
        )


class ManCliTests(unittest.TestCase):
    def test_center_page_lists_domains(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "--json"])
        self.assertEqual(code, 0)
        payload = json.loads(out.getvalue())
        names = {k["name"] for k in payload["keys"]}
        self.assertTrue(
            {"yaml", "sampler", "calculator", "operas", "tokens", "example"} <= names
        )

    def test_path_prefix_equivalence(self) -> None:
        a = resolve_man_request(["yaml", "Sampling.Bounds"])
        b = resolve_man_request(["yaml", "$.Sampling.Bounds"])
        self.assertEqual(a["path"], b["path"])
        self.assertEqual(a["summary"], b["summary"])

    def test_yaml_prefix_optional_for_calculator(self) -> None:
        a = resolve_man_request(["calculator", "execution"])
        b = resolve_man_request(["yaml", "calculator", "execution"])
        self.assertEqual(a["path"], b["path"])
        self.assertEqual(a["keys"], b["keys"])

    def test_unknown_topic_exit_code(self) -> None:
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            code = main(["man", "yaml", "Samplng"])
        self.assertNotEqual(code, 0)
        self.assertIn("did you mean", err.getvalue().lower())

    def test_sampler_bridson_has_keys_and_descriptions(self) -> None:
        page = resolve_man_request(["sampler", "Bridson"])
        self.assertEqual(page["status"], "stable")
        names = {k["name"] for k in page["keys"]}
        self.assertIn("Radius", names)
        radius = next(k for k in page["keys"] if k["name"] == "Radius")
        self.assertTrue(radius["description"])
        self.assertTrue(radius["description"].isascii())
        self.assertEqual(radius["type"], "number >0")

    def test_mcmc_family_marked_unstable(self) -> None:
        page = resolve_man_request(["sampler", "AMMCMC"])
        self.assertEqual(page["status"], "unstable")
        self.assertIn("not finalised", page["summary"].lower().replace("finalized", "finalised"))
        self.assertEqual(page["keys"], [])

    def test_all_methods_have_page(self) -> None:
        methods = schema_manifest()["sampling_methods"]
        for name in methods:
            page = resolve_man_request(["sampler", name])
            if page["status"] == "unstable":
                self.assertIn("not finalised", page["summary"].lower().replace("finalized", "finalised"))
            else:
                self.assertTrue(page["keys"] or page["examples"] or page["summary"])

    def test_dynesty_bounds_cover_constructor_keys(self) -> None:
        page = resolve_man_request(["sampler", "Dynesty"])
        self.assertEqual(page["status"], "stable")
        names = {k["name"] for k in page["keys"]}
        self.assertTrue(NESTED_CONSTRUCTOR_USER_KEYS <= names)
        # Schema allow-list matches contracts surface (top-level).
        nested = schema_by_id(
            "https://jarvis-hep.org/schema/v2/core/common.json"
        )["$defs"]["nestedBounds"]["properties"]
        self.assertEqual(set(nested), set(BOUNDS_TOP_ALLOWED))

    def test_variables_distribution_params(self) -> None:
        page = resolve_man_request(["yaml", "Sampling.Variables"])
        text = " ".join(k["description"] for k in page["keys"])
        for dtype, req in PARAMS_REQUIRED.items():
            for key in req:
                self.assertIn(key, text, msg=f"{dtype} missing {key}")

    def test_calculator_execution_type_json_self_contained(self) -> None:
        page = resolve_man_request(["calculator", "execution"], io_type="JSON")
        names = {k["name"] for k in page["keys"]}
        self.assertTrue({"name", "path", "type"} <= names)
        self.assertTrue(page["examples"])
        blob = json.dumps(page).lower()
        # A9: no hand-off phrasing that replaces the field table
        self.assertNotIn("run `jarvis portal man`", blob)
        self.assertNotIn("run jarvis portal man", blob)
        self.assertTrue(all(str(k.get("description") or "").isascii() for k in page["keys"]))

    def test_runtime_format_without_schema_is_honest(self) -> None:
        with mock.patch(
            "jarvishep2.io_portal.available_io_formats",
            return_value=["JSON", "FakeFmt"],
        ):
            page = resolve_man_request(["calculator", "execution"], io_type="FakeFmt")
        self.assertIn("FakeFmt", page.get("runtime_formats") or [])
        self.assertIn("no bundled field schema", page["summary"].lower())
        self.assertTrue(page.get("further_reading"))

    def test_tokens_and_operas_pages(self) -> None:
        tokens = resolve_man_request(["tokens"])
        self.assertGreaterEqual(len(tokens["keys"]), 5)
        operas = resolve_man_request(["operas"])
        names = {k["name"] for k in operas["keys"]}
        self.assertIn("operator", names)
        self.assertIn("name", names)

    def test_example_catalog_and_calculator_card_validate(self) -> None:
        page = resolve_man_request(["example"])
        names = {k["name"] for k in page["keys"]}
        self.assertIn("calculator", names)
        self.assertIn("bridson", names)
        calc = resolve_man_request(["example", "calculator"])
        self.assertIn("Calculators:", calc["examples"][0])
        root = Path(__file__).resolve().parents[1] / "jarvishep2" / "project_template"
        for key, rel in {
            "bridson": "bin/quickstart_bridson_operas.yaml",
            "csv": "bin/quickstart_csv_operas.yaml",
            "calculator": "bin/quickstart_calculator.yaml",
            "dynesty": "bin/sampling/Sampling_Dynesty_Simple.yaml",
        }.items():
            path = root / rel
            self.assertTrue(path.is_file(), msg=rel)
            core = Jarvis2Core()
            check = key == "calculator"
            core.load_task_yaml(str(path), validate=True, check_modules=check if check else None)

    def test_code_lookup_and_index_coverage(self) -> None:
        page = resolve_man_request([], code="JV2-MTH-012")
        self.assertEqual(page["context"].get("method"), "Bridson")
        page2 = resolve_man_request([], code="JV2-OPR-002")
        self.assertEqual(page2["path"], "$.Operas.Modules[]")
        for code in known_diagnostic_codes():
            self.assertIsNotNone(man_target_for_code(code), msg=code)
            cmd = man_command_for(code)
            self.assertTrue(cmd.startswith("Run: Jarvis man"), msg=f"{code} -> {cmd}")

    def test_validate_hints_are_executable_man_commands(self) -> None:
        iss = issue("error", "JV2-VAR-001", "Sampling.Variables", "missing variables")
        self.assertIsNotNone(iss.hint)
        self.assertIn("Jarvis man", iss.hint or "")
        self.assertNotIn("docs/", iss.hint or "")
        # Schema path also rewritten
        report = validate_task_config(
            {
                "Scan": {"name": "t"},
                "Sampling": {
                    "Method": "Bridson",
                    "Bounds": {"Radius": 0.1, "MaxAttempt": 30, "totally_bogus": 1},
                    "Variables": [
                        {
                            "name": "x",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0.0, "max": 1.0, "length": 1.0},
                            },
                        }
                    ],
                },
            }
        )
        sch = [i for i in report.issues if i.code == "JV2-SCH-001"]
        self.assertTrue(sch)
        for item in sch:
            self.assertIn("Jarvis man", item.hint or "")
            self.assertNotIn("docs/", item.hint or "")

    def test_json_contract_stable_keys(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            self.assertEqual(main(["man", "sampler", "Bridson", "--json"]), 0)
        payload = json.loads(out.getvalue())
        for key in (
            "path",
            "context",
            "zone",
            "status",
            "summary",
            "keys",
            "examples",
            "diagnostics",
            "see_also",
            "further_reading",
        ):
            self.assertIn(key, payload)

    def test_man_output_is_ascii_prose(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            main(["man", "sampler", "Bridson"])
        text = out.getvalue()
        # Allow box-drawing from Rich; prose letters must be ASCII.
        letters = re.findall(r"[^\x00-\x7F]", text)
        # Rich may emit unicode box chars; filter to letters only
        non_ascii_letters = [c for c in letters if c.isalpha()]
        self.assertEqual(non_ascii_letters, [])


if __name__ == "__main__":
    unittest.main()
