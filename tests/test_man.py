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
from jarvishep2.runtime_config import SUPPORTED_ENVREQS_V2_KEYS


class ManPathTests(unittest.TestCase):
    def test_normalize_path_prefixes(self) -> None:
        self.assertEqual(normalize_yaml_path("Sampling.Bounds"), "$.Sampling.Bounds")
        self.assertEqual(normalize_yaml_path("$.Sampling.Bounds"), "$.Sampling.Bounds")
        self.assertEqual(
            normalize_yaml_path("$.Calculators.Modules.execution.input"),
            "$.Calculators.Modules.execution.input",
        )
        with self.assertRaises(ValueError):
            normalize_yaml_path("$.Calculators.Modules[0].execution.input[0]")

    def test_log_likelihood_page_exposes_list_item_schema(self) -> None:
        sampling_page = resolve_man_request(["yaml.Sampling"])
        log_likelihood_row = next(
            row for row in sampling_page["keys"] if row["name"] == "LogLikelihood"
        )
        self.assertTrue(log_likelihood_row["nav"])
        self.assertEqual(
            log_likelihood_row["man"],
            "Jarvis man yaml.Sampling.LogLikelihood",
        )

        page = resolve_man_request(["yaml.Sampling.LogLikelihood"])
        self.assertEqual(page["path"], "$.Sampling.LogLikelihood")
        self.assertEqual(page["keys"], [])
        self.assertEqual(
            [row["name"] for row in page["list_rows"]],
            ["name", "expression"],
        )
        self.assertEqual(page["item_schema"]["required"], ["name", "expression"])

    def test_envreqs_v2_pages_are_closed(self) -> None:
        self.assertEqual(resolve_man_request(["yaml.EnvReqs.V2"])["zone"], "closed")
        self.assertEqual(
            resolve_man_request(["yaml.EnvReqs.V2.check_modules"])["zone"],
            "closed",
        )
        self.assertEqual(resolve_man_request(["yaml.EnvReqs"])["zone"], "closed")


class ManCliTests(unittest.TestCase):
    def test_man_help_uses_man_cards_and_coding_agent_hint(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "-h"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertIn("╭", rendered)
        self.assertIn("Jarvis man", rendered)
        self.assertIn("Options", rendered)
        self.assertIn("Examples · one dotted path", rendered)
        self.assertIn("Coding agents", rendered)
        self.assertIn("--json", rendered)
        self.assertNotIn("positional arguments:", rendered)

    def test_page_panels_wrap_long_text(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "calculator.execution"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertIn("Keys · $.Calculators.Modules.execution", rendered)
        self.assertIn("Use --type TYPE --direction", rendered)
        self.assertIn("input|output", rendered)
        self.assertIn("file-type", rendered)

    def test_portal_io_actions_are_direction_and_format_specific(self) -> None:
        json_input = resolve_man_request(
            ["calculator.execution.input"], io_type="JSON"
        )
        self.assertTrue(any("Dump:" in line for line in json_input["action_help"]))
        self.assertIn("actions:", json_input["examples"][0])
        variable_schema = json_input["item_schema"]["properties"]["actions"]["items"]["properties"]["variables"]["items"]
        self.assertEqual(len(variable_schema["oneOf"]), 4)
        self.assertIn(
            "entry",
            next(key for key in json_input["keys"] if key["name"] == "actions")["description"],
        )

        json_output = resolve_man_request(
            ["calculator.execution.output"], io_type="JSON"
        )
        self.assertEqual(json_output["action_help"], [])
        output_variable_schema = json_output["item_schema"]["properties"]["variables"]["items"]
        self.assertEqual(len(output_variable_schema["oneOf"]), 2)
        self.assertIn(
            "dotted JSON path",
            next(key for key in json_output["keys"] if key["name"] == "variables")["description"],
        )
        for fmt in ("JSON", "CSV", "DAT", "TSV", "Wolfram", "SLHA", "xSLHA"):
            output_page = resolve_man_request(
                ["calculator.execution.output"], io_type=fmt
            )
            self.assertEqual(len(output_page["examples"]), 1)
            self.assertEqual(output_page["action_rows"], [])
            self.assertTrue(
                all("- {name:" in example for example in output_page["examples"]),
                msg=f"{fmt} output examples should use inline dict variables",
            )

        csv_input = resolve_man_request(
            ["calculator.execution.input"], io_type="CSV"
        )
        csv_key_types = {(key["name"], key["type"]) for key in csv_input["keys"]}
        self.assertIn(("actions", "list"), csv_key_types)
        self.assertIn("column", " ".join(csv_input["action_help"]))

        text_input = resolve_man_request(
            ["calculator.execution.input"], io_type="Text"
        )
        self.assertEqual(len(text_input["action_help"]), 2)
        self.assertIn("Replace", text_input["action_help"][0])
        self.assertIn("File", text_input["action_help"][1])

        slha_input = resolve_man_request(
            ["calculator.execution.input"], io_type="SLHA"
        )
        self.assertIn("SLHA:", " ".join(slha_input["action_help"]))
        self.assertNotIn("operations", {key["name"] for key in slha_input["keys"]})
        self.assertGreaterEqual(len(slha_input["examples"]), 2)
        self.assertIn("{name: M1", slha_input["examples"][0])
        self.assertIn("{name: mix", slha_input["examples"][0])
        self.assertIn("type: File", slha_input["examples"][1])

        text_input = resolve_man_request(
            ["calculator.execution.input"], io_type="Text"
        )
        self.assertNotIn("content", {key["name"] for key in text_input["keys"]})

        with self.assertRaises(KeyError):
            resolve_man_request(["calculator.execution.output"], io_type="Text")
        with self.assertRaises(KeyError):
            resolve_man_request(["calculator.execution.input"], io_type="xSLHA")

    def test_portal_action_table_is_rendered_after_keys(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "calculator.execution.input", "--type", "JSON"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertIn("Keys · $.Calculators.Modules.execution.input", rendered)
        self.assertIn("Actions · $.Calculators.Modules.execution.input", rendered)
        self.assertIn("Dump", rendered)
        self.assertIn("variables list", rendered)
        self.assertLess(
            rendered.index("Keys · $.Calculators.Modules.execution.input"),
            rendered.index("Actions · $.Calculators.Modules.execution.input"),
        )

        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "calculator.execution.output", "--type", "JSON"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertNotIn("Actions · $.Calculators.Modules.execution.output", rendered)
        self.assertNotIn("No input actions on output", rendered)

    def test_execution_direction_lists_file_types_and_type_commands(self) -> None:
        page = resolve_man_request(["calculator.execution.input"])
        file_types = [row["type"] for row in page["type_rows"]]
        self.assertEqual(
            file_types,
            ["CSV", "DAT", "JSON", "SLHA", "TSV", "Text", "Wolfram"],
        )
        self.assertIn(
            "Jarvis man calculator.execution.input --type JSON",
            {row["command"] for row in page["type_rows"]},
        )
        self.assertIn("• Purpose:", page["summary"])
        self.assertIn("• File types:", page["summary"])
        self.assertIn("• Details:", page["summary"])
        self.assertIn("--type TYPE", page["summary"])

        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "calculator.execution.input"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertIn("File types · $.Calculators.Modules.execution.input · input", rendered)
        self.assertIn("Type", rendered)
        self.assertIn("Description", rendered)
        self.assertNotIn("Open with --type", rendered)
        self.assertIn("calculator.execution.input --type TYPE", rendered)

    def test_center_page_lists_domains(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man", "--json"])
        self.assertEqual(code, 0)
        payload = json.loads(out.getvalue())
        names = {k["name"] for k in payload["list_rows"]}
        self.assertTrue(
            {"yaml", "sampler", "calculator", "operas", "tokens", "example"} <= names
        )

    def test_center_page_uses_bulleted_guidance(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            code = main(["man"])
        self.assertEqual(code, 0)
        rendered = out.getvalue()
        self.assertIn("• Purpose:", rendered)
        self.assertIn("• Topic syntax:", rendered)
        self.assertGreaterEqual(rendered.count("•"), 8)

    def test_non_yaml_catalogs_use_list_tables(self) -> None:
        for topic, title in (
            ("tokens", "Tokens · $"),
            ("calculator", "Topics · $.Calculators"),
            ("example", "Examples · $"),
        ):
            out = io.StringIO()
            with redirect_stdout(out):
                code = main(["man", topic])
            self.assertEqual(code, 0)
            rendered = out.getvalue()
            self.assertIn(title, rendered)
            self.assertNotIn("Keys · $", rendered)

    def test_path_prefix_equivalence(self) -> None:
        a = resolve_man_request(["yaml.Sampling.Bounds"])
        b = resolve_man_request(["$.Sampling.Bounds"])
        self.assertEqual(a["path"], b["path"])
        self.assertEqual(a["summary"], b["summary"])

    def test_yaml_root_blocks_are_required(self) -> None:
        page = resolve_man_request(["yaml"])
        required = {key["name"] for key in page["keys"] if key["required"]}
        self.assertEqual(
            required,
            {"Scan", "Sampling", "EnvReqs"},
        )
        conditional = {
            key["name"]: key["requirement"]
            for key in page["keys"]
            if key.get("requirement")
        }
        self.assertEqual(
            conditional,
            {
                "Calculators": "one of Calculators / Operas",
                "Operas": "one of Calculators / Operas",
            },
        )
        self.assertFalse(
            next(key for key in page["keys"] if key["name"] == "LibDeps")["required"]
        )

    def test_sampling_method_is_required_in_man(self) -> None:
        page = resolve_man_request(["yaml.Sampling"])
        method = next(key for key in page["keys"] if key["name"] == "Method")
        self.assertTrue(method["required"])

    def test_envreqs_envelope_is_closed_in_man(self) -> None:
        page = resolve_man_request(["yaml.EnvReqs"])
        self.assertEqual(page["zone"], "closed")
        self.assertIn("OS", {key["name"] for key in page["keys"]})
        self.assertIn("preflight", page["summary"])

    def test_envreqs_os_page_exposes_list_item_schema(self) -> None:
        page = resolve_man_request(["yaml.EnvReqs.OS"])
        self.assertEqual(page["path"], "$.EnvReqs.OS")
        self.assertEqual({row["name"] for row in page["list_rows"]}, {"name", "version"})
        self.assertEqual(page["item_schema"]["required"], ["name", "version"])
        self.assertFalse(page["item_schema"]["additionalProperties"])

    def test_generic_yaml_path_summary_keeps_jsonpath_dot(self) -> None:
        page = resolve_man_request(["yaml.Scan"])
        self.assertEqual(page["path"], "$.Scan")
        self.assertEqual(page["summary"], "YAML path $.Scan")
        self.assertNotIn(".$Scan", page["summary"])
        name = next(key for key in page["keys"] if key["name"] == "name")
        self.assertTrue(name["required"])
        self.assertEqual({key["name"] for key in page["keys"]}, {"name"})

    def test_space_separated_topic_is_rejected(self) -> None:
        with self.assertRaises(KeyError) as ctx:
            resolve_man_request(["calculator", "execution", "input"], io_type="JSON")
        self.assertIn("calculator.execution.input", str(ctx.exception))
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            code = main(["man", "calculator", "execution"])
        self.assertNotEqual(code, 0)
        self.assertIn("calculator.execution", err.getvalue())
        self.assertIn("dotted path", err.getvalue().lower())

    def test_yaml_prefix_optional_for_calculator(self) -> None:
        a = resolve_man_request(["calculator.execution"])
        b = resolve_man_request(["yaml.calculator.execution"])
        self.assertEqual(a["path"], b["path"])
        self.assertEqual(a["keys"], b["keys"])

    def test_unknown_topic_exit_code(self) -> None:
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            code = main(["man", "yaml.Samplng"])
        self.assertNotEqual(code, 0)
        self.assertIn("did you mean", err.getvalue().lower())

    def test_sampler_bridson_has_keys_and_descriptions(self) -> None:
        page = resolve_man_request(["sampler.Bridson"])
        self.assertEqual(page["status"], "stable")
        names = {k["name"] for k in page["keys"]}
        self.assertIn("radius", names)
        radius = next(k for k in page["keys"] if k["name"] == "radius")
        self.assertTrue(radius["description"])
        self.assertTrue(radius["description"].isascii())
        self.assertEqual(radius["type"], "number >0")

    def test_mcmc_family_marked_unstable(self) -> None:
        page = resolve_man_request(["sampler.AMMCMC"])
        self.assertEqual(page["status"], "unstable")
        self.assertIn("not finalised", page["summary"].lower().replace("finalized", "finalised"))
        self.assertEqual(page["keys"], [])

    def test_all_methods_have_page(self) -> None:
        methods = schema_manifest()["sampling_methods"]
        for name in methods:
            page = resolve_man_request([f"sampler.{name}"])
            if page["status"] == "unstable":
                self.assertIn("not finalised", page["summary"].lower().replace("finalized", "finalised"))
            else:
                self.assertTrue(page["keys"] or page["examples"] or page["summary"])

    def test_dynesty_bounds_cover_constructor_keys(self) -> None:
        page = resolve_man_request(["sampler.Dynesty"])
        self.assertEqual(page["status"], "stable")
        names = {k["name"] for k in page["keys"]}
        self.assertTrue(NESTED_CONSTRUCTOR_USER_KEYS <= names)
        # Schema allow-list matches contracts surface (top-level).
        nested = schema_by_id(
            "https://jarvis-hep.org/schema/v2/core/common.json"
        )["$defs"]["nestedBounds"]["properties"]
        self.assertEqual(set(nested), set(BOUNDS_TOP_ALLOWED))

    def test_variables_distribution_params(self) -> None:
        page = resolve_man_request(["yaml.Sampling.Variables"])
        self.assertEqual(page["keys"], [])
        self.assertEqual(page["list_title"], "Distributions")
        text = " ".join(k["description"] for k in page["list_rows"])
        for dtype, req in PARAMS_REQUIRED.items():
            for key in req:
                self.assertIn(key, text, msg=f"{dtype} missing {key}")
        item_schema = page["item_schema"]
        self.assertEqual(item_schema["required"], ["name", "distribution"])
        self.assertEqual(
            set(item_schema["properties"]["distribution"]["properties"]["type"]["enum"]),
            set(PARAMS_REQUIRED),
        )
        variants = {variant["type"]: variant for variant in item_schema["distribution_variants"]}
        self.assertEqual(set(variants), set(PARAMS_REQUIRED))
        for dtype, req in PARAMS_REQUIRED.items():
            params = variants[dtype]["parameters"]
            self.assertEqual(set(params["required"]), set(req), msg=dtype)
            self.assertTrue(set(req) <= set(params["allowed"]), msg=dtype)
            self.assertTrue(params["example"], msg=dtype)
        self.assertEqual(
            item_schema["method_parameter_overrides"]["Bridson"]["required_parameters"],
            ["length"],
        )
        self.assertEqual(
            item_schema["method_parameter_overrides"]["Grid"]["required_parameters"],
            ["num"],
        )

    def test_calculator_execution_type_json_self_contained(self) -> None:
        # Default --type without --direction shows both vocabularies.
        both = resolve_man_request(["calculator.execution"], io_type="JSON")
        both_names = {k["name"] for k in both["keys"]}
        self.assertTrue(
            any(n.startswith("input.") for n in both_names)
            or {"name", "path", "type"} <= both_names
        )
        inp = resolve_man_request(["calculator.execution.input"], io_type="JSON")
        out = resolve_man_request(["calculator.execution.output"], io_type="JSON")
        self.assertEqual((inp.get("context") or {}).get("direction"), "input")
        self.assertEqual((out.get("context") or {}).get("direction"), "output")
        self.assertTrue({"name", "path", "type"} <= {k["name"] for k in inp["keys"]})
        self.assertTrue({"name", "path", "type"} <= {k["name"] for k in out["keys"]})
        # Output vocabulary uses variables/entry, not Load actions.
        out_names = {k["name"] for k in out["keys"]}
        self.assertIn("variables", out_names)
        self.assertTrue(inp["examples"] or out["examples"] or both["examples"])
        blob = json.dumps(inp).lower()
        # A9: field table present; portal CLI is supplementary further_reading only.
        self.assertNotIn("run `jarvis portal man`", blob)
        self.assertTrue(all(str(k.get("description") or "").isascii() for k in inp["keys"]))

    def test_runtime_format_without_schema_is_honest(self) -> None:
        with mock.patch(
            "jarvishep2.io_portal.available_io_formats",
            return_value=["JSON", "FakeFmt"],
        ):
            page = resolve_man_request(
                ["calculator.execution"],
                io_type="FakeFmt",
                direction="input",
            )
        self.assertIn("FakeFmt", page.get("runtime_types") or [])
        self.assertIn("no bundled field schema", page["summary"].lower())
        self.assertTrue(page.get("further_reading"))

    def test_tokens_and_operas_pages(self) -> None:
        tokens = resolve_man_request(["tokens"])
        self.assertEqual(tokens["keys"], [])
        self.assertEqual(tokens["list_title"], "Tokens")
        self.assertGreaterEqual(len(tokens["list_rows"]), 5)
        operas = resolve_man_request(["operas"])
        names = {k["name"] for k in operas["keys"]}
        self.assertIn("operator", names)
        self.assertIn("name", names)
        # Executable mapping-returning example (not math.add scalar).
        example = "\n".join(operas["examples"])
        self.assertIn("helper.eggbox2d", example)
        self.assertIn("entry: z", example)
        self.assertNotIn("math.add", example)
        self.assertTrue(all(str(k.get("description") or "").strip() for k in operas["keys"]))

    def test_calculator_yaml_paths_resolve(self) -> None:
        modules = resolve_man_request(["yaml.Calculators.Modules"])
        self.assertTrue(modules["keys"], msg=modules)
        names = {k["name"] for k in modules["keys"]}
        self.assertIn("name", names)
        self.assertIn("execution", names)
        execution = resolve_man_request(["yaml.Calculators.Modules.execution"])
        self.assertEqual(execution["path"], "$.Calculators.Modules.execution")
        self.assertTrue(execution["keys"])
        # Example must be valid YAML tokens (quoted &J/...) and output variables/entry.
        ex = "\n".join(execution["examples"])
        self.assertIn('path: "&J/', ex)
        self.assertIn("variables:", ex)
        self.assertNotIn("type: Load", ex)
        input_page = resolve_man_request(["yaml.Calculators.Modules.execution.input"])
        self.assertEqual(input_page["path"], "$.Calculators.Modules.execution.input")
        self.assertEqual(input_page["keys"][0]["name"], "type")

    def test_yaml_calculator_paths_preserve_portal_type_details(self) -> None:
        for direction in ("input", "output"):
            with self.subTest(direction=direction):
                direct = resolve_man_request(
                    [f"calculator.execution.{direction}"], io_type="JSON"
                )
                yaml_path = resolve_man_request(
                    [f"yaml.Calculators.Modules.execution.{direction}"],
                    io_type="JSON",
                )
                self.assertEqual(yaml_path["path"], direct["path"])
                self.assertEqual(yaml_path["context"], direct["context"])
                self.assertEqual(yaml_path["keys"], direct["keys"])
                self.assertEqual(
                    yaml_path.get("action_rows"), direct.get("action_rows")
                )

    def test_calculator_module_fields_have_descriptions(self) -> None:
        page = resolve_man_request(["calculator.module"])
        self.assertTrue(page["keys"])
        empty = [k["name"] for k in page["keys"] if not str(k.get("description") or "").strip()]
        self.assertEqual(empty, [])

    def test_keys_nav_marks_expandable_rows(self) -> None:
        """▸ rows carry nav+man; leaf fields like path do not."""
        page = resolve_man_request(["calculator.module"])
        by_name = {k["name"]: k for k in page["keys"]}
        self.assertTrue(by_name["execution"].get("nav"))
        self.assertEqual(
            by_name["execution"].get("man"),
            "Jarvis man calculator.execution",
        )
        # Prefer dotted path (not space-separated).
        self.assertNotIn("calculator execution", by_name["execution"].get("man") or "")
        self.assertTrue(by_name["modes"].get("nav"))
        self.assertFalse(by_name["path"].get("nav"))
        self.assertNotIn("man", by_name["path"])

        center = resolve_man_request([])
        self.assertIn("nav_legend", center)
        self.assertEqual(center["keys"], [])
        self.assertTrue(all(k.get("nav") for k in center["list_rows"]))
        self.assertTrue(all(k.get("man") for k in center["list_rows"]))

    def test_leaf_topics_explain_that_the_parent_page_is_required(self) -> None:
        for topic in ("calculator.module.path", "yaml.Calculators.Modules.path"):
            with self.assertRaises(KeyError) as raised:
                resolve_man_request([topic])
            self.assertIn("leaf field", str(raised.exception))
            self.assertIn("calculator.module", str(raised.exception))

        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            code = main(["man", "calculator.module.path"])
        self.assertEqual(code, 2)
        self.assertIn("Man topic is a leaf field", err.getvalue())
        self.assertNotIn("Unknown man topic: '", err.getvalue())

    def test_example_catalog_and_calculator_card_validate(self) -> None:
        page = resolve_man_request(["example"])
        names = {k["name"] for k in page["list_rows"]}
        self.assertIn("calculator", names)
        self.assertIn("bridson", names)
        calc = resolve_man_request(["example.calculator"])
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
            if key == "dynesty":
                # Sampling files are fragments copied into a full task card.
                continue
            core = Jarvis2Core()
            check = key == "calculator"
            core.load_task_yaml(str(path), validate=True, check_modules=check if check else None)

    def test_code_lookup_and_index_coverage(self) -> None:
        page = resolve_man_request([], code="JV2-MTH-012")
        self.assertEqual(page["context"].get("method"), "Bridson")
        page2 = resolve_man_request([], code="JV2-OPR-002")
        self.assertEqual(page2["path"], "$.Operas.Modules")
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
                    "Bounds": {"radius": 0.1, "max_attempt": 30, "totally_bogus": 1},
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

    def test_envreqs_pages_describe_project_defaults_and_overlays(self) -> None:
        top = resolve_man_request(["yaml.EnvReqs"])
        self.assertEqual(top["path"], "$.EnvReqs")
        self.assertIn("project create", top["summary"])
        self.assertIn("default_yaml", top)
        self.assertIn("EnvReqs:", top["default_yaml"])
        self.assertIn("V2:", top["default_yaml"])
        names = {key["name"] for key in top["keys"]}
        self.assertTrue({"Check_default_dependencies", "V2", "Python", "CERN_ROOT", "OS"} <= names)

        v2 = resolve_man_request(["yaml.EnvReqs.V2"])
        self.assertEqual(v2["path"], "$.EnvReqs.V2")
        self.assertIn("deep-merges", v2["summary"])
        v2_names = {key["name"] for key in v2["keys"]}
        self.assertEqual(v2_names, set(SUPPORTED_ENVREQS_V2_KEYS))
        self.assertEqual(next(key for key in v2["keys"] if key["name"] == "redis")["type"], "mapping")

        check = resolve_man_request(["yaml.EnvReqs.V2.check_modules"])
        self.assertEqual(check["path"], "$.EnvReqs.V2.check_modules")
        self.assertEqual(
            {key["name"] for key in check["keys"]},
            {"data", "n_samples", "timeout_sec"},
        )

        for topic in (
            "sample_directory",
            "archiver",
            "redis",
            "factory",
            "worker",
            "check_modules",
        ):
            nested = resolve_man_request([f"yaml.EnvReqs.V2.{topic}"])
            self.assertTrue(nested["keys"], msg=topic)
            self.assertTrue(all(key.get("description") for key in nested["keys"]))

    def test_array_bracket_topics_are_rejected(self) -> None:
        err = io.StringIO()
        with mock.patch("sys.stderr", err):
            code = main(["man", "yaml.Calculators.Modules[].execution"])
        self.assertEqual(code, 2)
        self.assertIn("brackets are no longer supported", err.getvalue())

    def test_man_uses_list_words_without_array_suffixes(self) -> None:
        pages = [
            resolve_man_request(["calculator.execution"]),
            resolve_man_request(["calculator.module"]),
            resolve_man_request(["operas"]),
            resolve_man_request(["yaml.EnvReqs.V2"]),
        ]
        for page in pages:
            prose = " ".join(
                [
                    str(page.get("path")),
                    str(page.get("summary")),
                    " ".join(str(item) for item in page.get("see_also") or []),
                    " ".join(str(item) for item in page.get("examples") or []),
                    str(page.get("default_yaml") or ""),
                    " ".join(
                        f"{key.get('name')} {key.get('description')}"
                        for key in page.get("keys") or []
                    ),
                ]
            )
            self.assertNotIn("[]", prose)
        execution = resolve_man_request(["calculator.execution"])
        self.assertEqual(
            {key["name"]: key["type"] for key in execution["keys"]}["input"],
            "list",
        )

    def test_json_contract_stable_keys(self) -> None:
        out = io.StringIO()
        with redirect_stdout(out):
            self.assertEqual(main(["man", "sampler.Bridson", "--json"]), 0)
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
            main(["man", "sampler.Bridson"])
        text = out.getvalue()
        # Allow box-drawing from Rich; prose letters must be ASCII.
        letters = re.findall(r"[^\x00-\x7F]", text)
        # Rich may emit unicode box chars; filter to letters only
        non_ascii_letters = [c for c in letters if c.isalpha()]
        self.assertEqual(non_ascii_letters, [])


if __name__ == "__main__":
    unittest.main()
