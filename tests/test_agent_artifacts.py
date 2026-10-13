"""Installed artifact exports and offline YAML-writing workflow for agents."""

from __future__ import annotations

import contextlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from jsonschema import Draft202012Validator
from referencing import Registry, Resource
from referencing.jsonschema import DRAFT202012
import yaml

from jarvishep2.cli.agent_artifacts import (
    SCHEMA_SCOPE,
    dispatch_schema,
    dispatch_template,
    schema_bundle,
    template_catalog,
    template_document,
)
from jarvishep2.cli.man import dispatch_man, resolve_man_request
from jarvishep2.cli.run import dispatch_validate
from jarvishep2.project_scaffold import create_project_scaffold


def _capture(dispatch, argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = dispatch(argv)
    return code, out.getvalue(), err.getvalue()


class AgentSchemaTests(unittest.TestCase):
    def test_default_and_explicit_json_export_complete_local_bundle(self) -> None:
        code, output, errors = _capture(dispatch_schema, [])
        self.assertEqual((code, errors), (0, ""))
        bundle = json.loads(output)
        self.assertEqual(bundle["contract_version"], 1)
        self.assertEqual(bundle["scope"], SCHEMA_SCOPE)
        self.assertEqual(bundle["root_schema_id"], bundle["manifest"]["root"])
        self.assertEqual(len(bundle["schemas"]), len(bundle["manifest"]["schema_files"]))
        self.assertEqual(json.loads(_capture(dispatch_schema, ["--json"])[1]), bundle)

    def test_every_exported_ref_resolves_offline(self) -> None:
        bundle = schema_bundle()

        def forbid_retrieval(uri):
            self.fail(f"Attempted external schema retrieval: {uri}")

        registry = Registry(retrieve=forbid_retrieval).with_resources(
            (schema_id, Resource.from_contents(document, default_specification=DRAFT202012))
            for schema_id, document in bundle["schemas"].items()
        )

        def check_refs(node, resolver):
            if isinstance(node, dict):
                if "$ref" in node:
                    with self.subTest(ref=node["$ref"]):
                        self.assertIsNotNone(resolver.lookup(node["$ref"]).contents)
                for value in node.values():
                    check_refs(value, resolver)
            elif isinstance(node, list):
                for value in node:
                    check_refs(value, resolver)

        for schema_id, document in bundle["schemas"].items():
            check_refs(document, registry.resolver(schema_id))
        card = yaml.safe_load(template_document("random")["body"])
        validator = Draft202012Validator(bundle["schemas"][bundle["root_schema_id"]], registry=registry)
        self.assertEqual(list(validator.iter_errors(card)), [])

    def test_id_selects_exact_bundled_document_without_uri_fetch(self) -> None:
        bundle = schema_bundle()
        schema_id = bundle["manifest"]["sampling_methods"]["Random"]
        code, output, errors = _capture(dispatch_schema, ["--id", schema_id, "--json"])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output), bundle["schemas"][schema_id])
        code, output, errors = _capture(dispatch_schema, ["--id", "https://example.test/schema.json"])
        self.assertEqual((code, errors), (2, ""))
        self.assertFalse(json.loads(output)["ok"])

    def test_small_index_navigates_builtin_methods_and_io_without_schema_bodies(self) -> None:
        code, output, errors = _capture(dispatch_schema, ["--list", "--json"])
        self.assertEqual((code, errors), (0, ""))
        index = json.loads(output)
        bundle = schema_bundle()
        self.assertEqual(set(index["schema_ids"]), set(bundle["schemas"]))
        self.assertEqual(index["sampling_methods"], bundle["manifest"]["sampling_methods"])
        self.assertEqual(index["io"], bundle["manifest"]["io"])
        self.assertNotIn("schemas", index)
        self.assertLess(len(output), len(_capture(dispatch_schema, ["--json"])[1]) // 4)
        code, output, errors = _capture(dispatch_schema, ["--list", "--id", index["root_schema_id"]])
        self.assertEqual((code, errors), (2, ""))
        self.assertFalse(json.loads(output)["ok"])

    def test_local_fragment_pointer_handles_refs_arrays_and_uri_escaping(self) -> None:
        common_id = "https://jarvis-hep.org/schema/v2/core/common.json"
        bundle = schema_bundle()
        for reference in (common_id + "#/$defs/csv", common_id + "#%2F$defs%2Fcsv"):
            with self.subTest(reference=reference):
                code, output, errors = _capture(dispatch_schema, ["--id", reference])
                self.assertEqual((code, errors), (0, ""))
                self.assertEqual(json.loads(output), bundle["schemas"][common_id]["$defs"]["csv"])
        method_id = bundle["manifest"]["sampling_methods"]["Random"]
        code, output, errors = _capture(dispatch_schema, ["--id", method_id + "#/allOf/1/properties/Bounds"])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output), bundle["schemas"][method_id]["allOf"][1]["properties"]["Bounds"])
        escaped = {"schemas": {common_id: {"$defs": {"a/b~c": {"type": "string"}}}}}
        with patch("jarvishep2.cli.agent_artifacts.schema_bundle", return_value=escaped):
            code, output, errors = _capture(dispatch_schema, ["--id", common_id + "#/$defs/a~1b~0c"])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output), {"type": "string"})

    def test_invalid_or_unknown_fragments_fail_without_retrieval(self) -> None:
        bundle = schema_bundle()
        method_id = bundle["manifest"]["sampling_methods"]["Random"]
        for fragment in ("#/missing", "#anchor", "#/allOf/01", "#/allOf/99", "#/allOf/-1", "#/~2", "#/%XX"):
            with self.subTest(fragment=fragment):
                code, output, errors = _capture(dispatch_schema, ["--id", method_id + fragment, "--json"])
                self.assertEqual((code, errors), (2, ""))
                self.assertFalse(json.loads(output)["ok"])
        code, output, errors = _capture(dispatch_schema, ["--id", "https://[bad/schema", "--json"])
        self.assertEqual((code, errors), (2, ""))
        self.assertFalse(json.loads(output)["ok"])

    def test_usage_and_internal_failures_are_single_json_objects(self) -> None:
        for argv in (["--json", "--unknown"], ["--id", "--json"]):
            with self.subTest(argv=argv):
                code, output, errors = _capture(dispatch_schema, argv)
                self.assertEqual((code, errors), (2, ""))
                self.assertFalse(json.loads(output)["ok"])
        with patch("jarvishep2.cli.agent_artifacts.schema_bundle", side_effect=OSError("missing schema")):
            code, output, errors = _capture(dispatch_schema, ["--json"])
        self.assertEqual((code, errors), (1, ""))
        self.assertIn("missing schema", json.loads(output)["error"]["message"])


class BuiltinSamplerSchemaNavigationTests(unittest.TestCase):
    def test_json_renderer_preserves_schema_routes_and_scope(self) -> None:
        code, output, errors = _capture(dispatch_man, ["sampler.CSV", "--json"])
        self.assertEqual((code, errors), (0, ""))
        page = json.loads(output)
        for key in ("schema_id", "schema_argv", "bounds_schema_ref", "bounds_schema_argv", "options_schema", "options_schema_scope"):
            self.assertIn(key, page)
        self.assertIn("Delegated", page["options_schema_scope"])
        code, output, errors = _capture(dispatch_schema, page["bounds_schema_argv"][3:])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output)["required"], ["path"])

    def test_builtin_pages_supply_exact_schema_and_bounds_navigation(self) -> None:
        bundle = schema_bundle()
        for method, schema_id in bundle["manifest"]["sampling_methods"].items():
            with self.subTest(method=method):
                page = resolve_man_request([f"sampler.{method}"])
                self.assertEqual(page["schema_id"], schema_id)
                self.assertEqual(page["schema_argv"], ["Jarvis", "agent", "schema", "--id", schema_id, "--json"])
                self.assertTrue(page["bounds_schema_ref"].startswith("https://jarvis-hep.org/schema/"))
                code, output, errors = _capture(dispatch_schema, page["bounds_schema_argv"][3:])
                self.assertEqual((code, errors), (0, ""))
                self.assertIsInstance(json.loads(output), dict)
                self.assertIsInstance(page["options_schema"], dict)
                Draft202012Validator.check_schema(page["options_schema"])

    def test_copied_bounds_schema_retains_local_ref_base_and_enforces_constraints(self) -> None:
        bundle = schema_bundle()

        def forbid_retrieval(uri):
            self.fail(f"Attempted external schema retrieval: {uri}")

        registry = Registry(retrieve=forbid_retrieval).with_resources(
            (schema_id, Resource.from_contents(document, default_specification=DRAFT202012))
            for schema_id, document in bundle["schemas"].items()
        )
        page = resolve_man_request(["sampler.CSV"])
        self.assertEqual(page["bounds_schema_ref"], "https://jarvis-hep.org/schema/v2/core/common.json#/$defs/csv")
        self.assertEqual(page["options_schema"]["properties"]["path"]["$ref"], "https://jarvis-hep.org/schema/v2/core/common.json#/$defs/nonemptyString")
        validator = Draft202012Validator(page["options_schema"], registry=registry)
        self.assertEqual(list(validator.iter_errors({"path": "points.csv", "uuid_column": "uuid"})), [])
        self.assertTrue(list(validator.iter_errors({"path": ""})))
        self.assertTrue(list(validator.iter_errors({"path": "points.csv", "unknown": True})))
        page = resolve_man_request(["sampler.Random"])
        validator = Draft202012Validator(page["options_schema"], registry=registry)
        self.assertEqual(list(validator.iter_errors({"point_number": 2})), [])
        self.assertTrue(list(validator.iter_errors({"point_number": "many"})))

    def test_dynesty_threshold_description_matches_dynamic_alias_precedence(self) -> None:
        from jarvishep2.sampling.dynesty_sampler import extract_run_nested_kwargs

        dynesty = resolve_man_request(["sampler.Dynesty"])
        descriptions = {row["name"]: row["description"] for row in dynesty["keys"]}
        self.assertIn("dynamic initial-phase", descriptions["dlogz"])
        self.assertIn("Bounds.dlogz_init", descriptions["dlogz"])
        self.assertNotIn("static nested", descriptions["dlogz"])
        self.assertEqual(extract_run_nested_kwargs({"dlogz": 0.3}, dynamic=True)["dlogz_init"], 0.3)
        self.assertEqual(extract_run_nested_kwargs({"dlogz": 0.3, "dlogz_init": 0.2}, dynamic=True)["dlogz_init"], 0.2)
        self.assertEqual(extract_run_nested_kwargs({"dlogz": 0.3, "run_nested": {"dlogz_init": 0.1}}, dynamic=True)["dlogz_init"], 0.1)
        multinest = resolve_man_request(["sampler.MultiNest"])
        descriptions = {row["name"]: row["description"] for row in multinest["keys"]}
        self.assertIn("Static NestedSampler", descriptions["dlogz"])
        self.assertIn("Compatibility alias", descriptions["dlogz_init"])


class AgentTemplateTests(unittest.TestCase):
    def test_catalog_lists_exact_man_catalog_with_fragment_metadata(self) -> None:
        catalog = template_catalog()
        expected = {row["name"] for row in resolve_man_request(["example"])["list_rows"]}
        self.assertEqual({row["name"] for row in catalog["templates"]}, expected)
        for row in catalog["templates"]:
            kind = row["artifact"]["kind"]
            self.assertEqual(kind, "sampling_fragment" if row["name"].startswith(("dynesty", "multinest")) else "complete_task_card")
        code, output, errors = _capture(dispatch_template, ["--json"])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output), catalog)
        code, output, errors = _capture(dispatch_template, [])
        self.assertEqual((code, errors), (0, ""))
        for name in expected:
            self.assertIn(name, output)

    def test_every_template_is_byte_for_byte_shipped_yaml_and_man_body(self) -> None:
        root = Path(__file__).parents[1] / "jarvishep2" / "project_template"
        for entry in template_catalog()["templates"]:
            name = entry["name"]
            with self.subTest(name=name):
                payload = template_document(name)
                shipped = (root / payload["artifact"]["path"]).read_text(encoding="utf-8")
                code, output, errors = _capture(dispatch_template, [name])
                self.assertEqual((code, errors), (0, ""))
                self.assertEqual(output, shipped)
                self.assertEqual(payload["body"], shipped)
                self.assertEqual(payload["body"], resolve_man_request([f"example.{name}"])["body"])
                self.assertTrue(payload["prerequisites"])

    def test_json_prerequisites_cover_actual_paths_csv_calculators_operators(self) -> None:
        csv = template_document("csv")
        reqs = csv["prerequisites"]
        project = next(row for row in reqs if row["kind"] == "project_scaffold")
        self.assertTrue(project["required"])
        self.assertEqual(project["argv"], ["Jarvis", "project", "create", "NAME"])
        paths = {row["value"] for row in reqs if row["kind"] == "project_path"}
        self.assertEqual(paths, {"&J/deps/environment_default.yaml", "&J/data/points.csv"})
        self.assertTrue(all(row["provided_by_scaffold"] for row in reqs if row["kind"] == "project_path"))
        csv_requirement = next(row for row in reqs if row["kind"] == "csv")
        self.assertEqual(csv_requirement["uuid_column"], "uuid")
        self.assertEqual(csv_requirement["input_columns"], ["x", "y"])
        self.assertEqual(csv_requirement["bundled_header"], ["uuid", "x", "y"])
        operator = next(row for row in reqs if row["kind"] == "operator")
        self.assertEqual(operator["operator"], yaml.safe_load(csv["body"])["Operas"]["Modules"][0]["operator"])
        self.assertEqual(operator["reference_kind"], "python_import_path")
        self.assertNotIn("query_argv", operator)
        self.assertFalse(operator["signature_verified"])
        calculator = template_document("random-calculator")
        calc = next(row for row in calculator["prerequisites"] if row["kind"] == "calculator")
        actual = yaml.safe_load(calculator["body"])["Calculators"]["Modules"][0]
        self.assertEqual(calc["path"], actual["path"])
        self.assertEqual(calc["commands"], actual["execution"]["commands"])
        self.assertEqual({row["path"] for row in calc["io_paths"]}, {"@Sdir/input.json", "@Sdir/out.json"})
        code, output, errors = _capture(dispatch_template, ["random-calculator", "--json"])
        self.assertEqual((code, errors), (0, ""))
        self.assertEqual(json.loads(output), calculator)

    def test_fragments_explain_composition_without_rewriting_the_yaml(self) -> None:
        for name in ("dynesty", "dynesty-full", "multinest", "multinest-full"):
            with self.subTest(name=name):
                payload = template_document(name)
                self.assertEqual(set(yaml.safe_load(payload["body"])), {"Sampling"})
                self.assertEqual(payload["artifact"]["kind"], "sampling_fragment")
                composition = next(row for row in payload["prerequisites"] if row["kind"] == "composition")
                self.assertEqual(composition["missing_root_blocks"], ["Scan", "EnvReqs"])
                self.assertEqual(composition["alternative_root_blocks"], ["Calculators", "Operas"])
                self.assertEqual(composition["base_template_argv"], ["Jarvis", "agent", "template", "calculator", "--json"])
                self.assertEqual(composition["operation"], {"action": "replace_mapping", "source_path": "$.Sampling", "target_path": "$.Sampling"})
                self.assertIn("entire root Sampling mapping", " ".join(composition["steps"]))
                validation = next(row for row in payload["prerequisites"] if row["kind"] == "validation")
                self.assertEqual(validation["argv"], ["Jarvis", "validate", "<ComposedTask.yaml>", "--json"])

    def test_operator_source_metadata_never_imports_callables(self) -> None:
        original = template_document("random")
        for operator, kind in (("helper.eggbox2d", "registered_name"), ("private_user.module.function", "registered_or_import_path")):
            with self.subTest(operator=operator):
                document = yaml.safe_load(original["body"])
                document["Operas"]["Modules"][0]["operator"] = operator
                page = {**original, "body": yaml.safe_dump(document)}
                with patch("jarvishep2.cli.agent_artifacts.resolve_man_request", side_effect=[{"list_rows": [{"name": "random"}]}, page]), patch("importlib.import_module") as import_module:
                    payload = template_document("random")
                import_module.assert_not_called()
                requirement = next(row for row in payload["prerequisites"] if row["kind"] == "operator")
                self.assertEqual(requirement["reference_kind"], kind)
                self.assertFalse(requirement["signature_verified"])
                self.assertEqual(requirement["input"], document["Operas"]["Modules"][0]["input"])
                if kind == "registered_name":
                    self.assertEqual(requirement["query_argv"], ["Jarvis", "operas", "info", operator])
                else:
                    self.assertNotIn("query_argv", requirement)

    def test_template_usage_unknown_and_internal_failures_preserve_json(self) -> None:
        for argv in (["missing", "--json"], ["random.extra", "--json"], ["random", "--unknown", "--json"], ["random", "extra", "--json"]):
            with self.subTest(argv=argv):
                code, output, errors = _capture(dispatch_template, argv)
                self.assertEqual((code, errors), (2, ""))
                self.assertFalse(json.loads(output)["ok"])
        with patch("jarvishep2.cli.agent_artifacts.resolve_man_request", side_effect=OSError("missing template")):
            code, output, errors = _capture(dispatch_template, ["random", "--json"])
        self.assertEqual((code, errors), (1, ""))
        self.assertIn("missing template", json.loads(output)["error"]["message"])

    def test_help_supports_human_and_json(self) -> None:
        for dispatch in (dispatch_schema, dispatch_template):
            with self.subTest(dispatch=dispatch.__name__):
                code, output, errors = _capture(dispatch, ["--help"])
                self.assertEqual((code, errors), (0, ""))
                self.assertIn("usage:", output)
                code, output, errors = _capture(dispatch, ["--help", "--json"])
                self.assertEqual((code, errors), (0, ""))
                self.assertIn("usage:", json.loads(output)["help"])

    def test_exports_leave_working_directory_untouched(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            original = os.getcwd()
            try:
                os.chdir(directory)
                with patch("jarvishep2.project_scaffold.create_project_scaffold") as create:
                    for dispatch, argv in ((dispatch_schema, []), (dispatch_template, []), (dispatch_template, ["csv", "--json"])):
                        self.assertEqual(_capture(dispatch, argv)[0], 0)
                    create.assert_not_called()
                    self.assertEqual(list(Path(directory).iterdir()), [])
            finally:
                os.chdir(original)

    def test_exported_cards_and_composed_fragments_validate_in_temporary_scaffold(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(create_project_scaffold("AgentArtifacts", cwd=directory))
            for entry in template_catalog()["templates"]:
                name = entry["name"]
                with self.subTest(name=name):
                    payload = template_document(name)
                    card_path = root / "bin" / f"agent_{name}.yaml"
                    if payload["artifact"]["kind"] == "sampling_fragment":
                        composition = next(row for row in payload["prerequisites"] if row["kind"] == "composition")
                        base_name = composition["base_template_argv"][3]
                        document = yaml.safe_load(template_document(base_name)["body"])
                        document["Scan"]["name"] = f"agent_{name}"
                        document["Sampling"] = yaml.safe_load(payload["body"])["Sampling"]
                        card_path.write_text(yaml.safe_dump(document), encoding="utf-8")
                    else:
                        card_path.write_text(payload["body"], encoding="utf-8")
                    out, err = io.StringIO(), io.StringIO()
                    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                        code = dispatch_validate(str(card_path), as_json=True)
                    result = json.loads(out.getvalue())
                    self.assertEqual(code, 0, msg=result)
                    self.assertTrue(result["ok"], msg=result)


if __name__ == "__main__":
    unittest.main()
