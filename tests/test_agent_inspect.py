"""Loaded-card explanations report real defaults without executing resources."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

import yaml

from jarvishep2.cli.agent_inspect import dispatch_inspect
from jarvishep2.cli.man import dispatch_man
from jarvishep2.task_config import TaskCardLoadError


class AgentInspectTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory(prefix="jarvis-agent-inspect-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        (self.root / "bin").mkdir()
        (self.root / "deps").mkdir()
        (self.root / ".jarvis-project.json").write_text("{}", encoding="utf-8")
        self.card = self.root / "bin" / "task.yaml"
        self.document = {
            "Scan": {"name": "inspection"},
            "Sampling": {
                "Method": "Random", "Bounds": {"point_number": 10},
                "Variables": [{"name": "x", "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}}}],
            },
            "EnvReqs": {"V2": {"workers": 0}},
            "Operas": {"Modules": [{"name": "example", "operator": "user_operator.calculate"}]},
        }
        self.save()
        self.environment = mock.patch.dict("os.environ", {"JARVIS_HEP_TASK_ROOT": "", "JHEP_TASK_ROOT": ""})
        self.environment.start()
        self.addCleanup(self.environment.stop)

    def save(self) -> None:
        self.card.write_text(yaml.safe_dump(self.document, sort_keys=False), encoding="utf-8")

    def inspect(self, argv=None):
        output, errors = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(errors):
            with mock.patch("jarvishep2.core.Jarvis2Core") as core, mock.patch("jarvishep2.redis_queue.RedisQueue") as redis:
                with mock.patch("jarvishep2.command_parser.CommandParser.register", side_effect=AssertionError("no registration")), mock.patch("subprocess.Popen", side_effect=AssertionError("no process execution")):
                    status = dispatch_inspect([str(self.card), "--json"] if argv is None else argv)
                core.assert_not_called()
                redis.assert_not_called()
        payload = json.loads(output.getvalue())
        self.assertNotIn("Traceback", errors.getvalue())
        if "metadata" in payload:
            self.assertFalse(payload["metadata"]["execution_performed"])
            self.assertFalse(payload["metadata"]["resource_preflight_performed"])
            self.assertFalse(payload["metadata"]["user_operator_imports_performed"])
        return status, payload, output.getvalue()

    def test_internal_defaults_and_planned_paths_without_final_resource_claims(self) -> None:
        status, payload, _ = self.inspect()
        self.assertEqual(status, 0)
        self.assertTrue(payload["ok"])
        self.assertEqual(payload["kind"], "agent.inspect")
        self.assertTrue(payload["installed_version"])
        self.assertEqual(payload["raw_task_card"], self.document)
        self.assertEqual(payload["merged_envreqs"], self.document["EnvReqs"])
        self.assertEqual(payload["paths"]["project_root"], str(self.root))
        self.assertIsNone(payload["paths"]["environment_defaults_yaml"])
        self.assertEqual(payload["paths"]["scan_output_root"], str(self.root / "outputs" / "inspection"))
        self.assertEqual(payload["paths"]["logs_dir"], str(self.root / "logs" / "inspection"))
        self.assertEqual(payload["paths"]["declared_checkpoint_file"], str(self.root / "checkpoints" / "inspection" / "Random" / "state.pkl"))
        self.assertEqual(payload["policies"]["checkpoint"], {"enabled": True, "heartbeat": 30.0})
        self.assertEqual(payload["policies"]["redis"]["port"], 6379)
        self.assertEqual(payload["runtime_decisions"]["configured_workers"], 0)
        self.assertIsNone(payload["runtime_decisions"]["selected_workers"])
        self.assertIsNone(payload["runtime_decisions"]["selected_redis_port"])
        self.assertEqual(payload["sources"]["runtime_fields"]["$.Runtime.batch_size"]["layer"], "internal_normalization")
        self.assertFalse(payload["preview"]["truncated"])

    def test_task_overrides_defaults_and_reports_field_sources(self) -> None:
        defaults = self.root / "deps" / "defaults.yaml"
        defaults.write_text(yaml.safe_dump({"EnvReqs": {
            "Python": {"version": ">=3.10"},
            "V2": {"workers": 4, "checkpoint": {"enabled": True, "heartbeat": 90}, "redis": {"port": 6387}},
        }}), encoding="utf-8")
        self.document["EnvReqs"] = {
            "Check_default_dependencies": {"required": True, "default_yaml_path": "../deps/defaults.yaml"},
            "V2": {"workers": 2, "checkpoint": {"heartbeat": 60}, "redis": {"db": 3},
                   "store_samples": False, "check_modules": {"data": "../data/points.csv"}},
        }
        self.save()
        before = {str(path.relative_to(self.root)): (path.read_bytes(), path.stat().st_mtime_ns)
                  for path in self.root.rglob("*") if path.is_file()}
        status, payload, _ = self.inspect()
        self.assertEqual(status, 0)
        self.assertEqual(payload["paths"]["environment_defaults_yaml"], str(defaults))
        self.assertEqual(payload["paths"]["check_modules_data"], str(self.root / "data" / "points.csv"))
        self.assertEqual(payload["runtime"]["workers"], 2)
        self.assertEqual(payload["policies"]["checkpoint"], {"enabled": True, "heartbeat": 60.0})
        self.assertEqual(payload["policies"]["redis"], {"host": "127.0.0.1", "port": 6387, "db": 3})
        self.assertFalse(payload["policies"]["sample_directory"]["enabled"])
        self.assertFalse(payload["path_notes"]["sample_storage_enabled"])
        self.assertFalse(payload["policies"]["archiver"]["pack_buckets"])
        fields = payload["sources"]["runtime_fields"]
        self.assertEqual(fields["$.Runtime.workers"]["layer"], "task")
        self.assertEqual(fields["$.Runtime.checkpoint.enabled"]["layer"], "default_environment")
        self.assertEqual(fields["$.Runtime.checkpoint.heartbeat"]["file"], str(self.card))
        self.assertEqual(payload["sources"]["envreqs_fields"]["$.EnvReqs.Python.version"]["file"], str(defaults))
        after = {str(path.relative_to(self.root)): (path.read_bytes(), path.stat().st_mtime_ns)
                 for path in self.root.rglob("*") if path.is_file()}
        self.assertEqual(before, after)

    def test_j_tokens_and_environment_root_follow_loader_resolution(self) -> None:
        defaults = self.root / "deps" / "defaults.yaml"
        defaults.write_text("EnvReqs:\n  V2:\n    workers: 3\n", encoding="utf-8")
        self.document["EnvReqs"]["Check_default_dependencies"] = {"required": True, "default_yaml_path": "&J/deps/defaults.yaml"}
        self.save()
        with mock.patch.dict("os.environ", {"JARVIS_HEP_TASK_ROOT": str(self.root)}):
            status, payload, _ = self.inspect()
        self.assertEqual(status, 0)
        self.assertEqual(payload["paths"]["environment_defaults_yaml"], str(defaults))
        self.assertEqual(payload["runtime"]["workers"], 0)

    def test_valid_card_does_not_import_user_operator_or_change_project_files(self) -> None:
        sentinel = self.root / "operator-imported"
        module = self.root / "user_operator.py"
        module.write_text(f"from pathlib import Path\nPath({str(sentinel)!r}).touch()\ndef calculate(**kwargs):\n    raise AssertionError('must not execute')\n", encoding="utf-8")
        before = {str(path.relative_to(self.root)): (path.read_bytes(), path.stat().st_mtime_ns)
                  for path in self.root.rglob("*") if path.is_file()}
        with mock.patch.object(sys, "path", [str(self.root), *sys.path]):
            status, payload, _ = self.inspect()
        after = {str(path.relative_to(self.root)): (path.read_bytes(), path.stat().st_mtime_ns)
                 for path in self.root.rglob("*") if path.is_file()}
        self.assertEqual(status, 0)
        self.assertTrue(payload["ok"])
        self.assertNotIn("user_operator", sys.modules)
        self.assertFalse(sentinel.exists())
        self.assertEqual(before, after)

    def test_bad_card_retains_schema_diagnostics_and_callable_navigation(self) -> None:
        self.document["Sampling"] = [1]
        self.save()
        status, payload, _ = self.inspect()
        self.assertEqual(status, 2)
        self.assertFalse(payload["ok"])
        self.assertTrue(payload["metadata"]["validation_completed"])
        self.assertTrue(payload["issues"])
        for diagnostic in payload["issues"]:
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(dispatch_man(diagnostic["manual_argv"][2:]), 0)

    def test_bad_policy_type_does_not_mask_original_validator_errors(self) -> None:
        self.document["EnvReqs"]["V2"]["archiver"] = {"pack_buckets": "maybe"}
        self.save()
        status, payload, _ = self.inspect()
        self.assertEqual(status, 2)
        self.assertFalse(payload["ok"])
        self.assertTrue(payload["issues"])
        self.assertIn("explanation_unavailable", payload)

    def test_strict_warning_and_raw_normalized_values_remain_distinct(self) -> None:
        self.document["Sampling"]["Bounds"]["point_number"] = "10"
        self.save()
        ordinary, payload, _ = self.inspect()
        strict, strict_payload, _ = self.inspect([str(self.card), "--strict", "--json"])
        self.assertEqual(ordinary, 0)
        self.assertTrue(any(item["level"] == "warning" for item in payload["issues"]))
        self.assertEqual(strict, 2)
        self.assertTrue(any(item["level"] == "error" for item in strict_payload["issues"]))

    def test_missing_task_and_syntax_error_use_configuration_exit_code(self) -> None:
        for argv in ([], ["--json"], [str(self.card), "--unknown"]):
            with self.subTest(argv=argv):
                status, payload, _ = self.inspect(argv)
                self.assertEqual(status, 2)
                self.assertFalse(payload["ok"])
        self.card.write_text("Scan: [bad\n", encoding="utf-8")
        status, payload, _ = self.inspect()
        self.assertEqual(status, 2)
        self.assertEqual(payload["issues"][0]["code"], "JV2-YAML-001")

    def test_io_and_validator_internal_failures_return_json_exit_one(self) -> None:
        status, payload, _ = self.inspect([str(self.root / "missing.yaml")])
        self.assertEqual(status, 1)
        self.assertEqual(payload["issues"][0]["code"], "JV2-LOAD-001")
        with mock.patch("jarvishep2.task_config.load_task_yaml", side_effect=PermissionError("unreadable")):
            status, payload, _ = self.inspect()
        self.assertEqual(status, 1)
        self.assertEqual(payload["issues"][0]["code"], "JV2-LOAD-003")
        with mock.patch("jarvishep2.task_validation.validate_task_config", side_effect=RuntimeError("validator unavailable")):
            status, payload, _ = self.inspect()
        self.assertEqual(status, 1)
        self.assertFalse(payload["metadata"]["validation_completed"])
        self.assertEqual(payload["issues"][0]["code"], "JV2-VAL-001")

    def test_load_location_metadata_is_preserved(self) -> None:
        error = TaskCardLoadError("JV2-YAML-003", "duplicate key", "Keep only the intended key.")
        error.path = "$.Scan.name"
        error.location = {"line": 3, "column": 3}
        error.related_locations = [{"line": 2, "column": 3}]
        with mock.patch("jarvishep2.task_config.load_task_yaml", side_effect=error):
            status, payload, _ = self.inspect()
        self.assertEqual(status, 2)
        self.assertEqual(payload["issues"][0]["path"], error.path)
        self.assertEqual(payload["issues"][0]["location"], error.location)
        self.assertEqual(payload["issues"][0]["related_locations"], error.related_locations)

    def test_large_kwargs_have_explicit_preview_omissions(self) -> None:
        self.document["Operas"]["Modules"][0]["kwargs"] = {"large_payload": "a" * 100000}
        self.save()
        status, payload, output = self.inspect()
        self.assertEqual(status, 0)
        self.assertTrue(payload["preview"]["truncated"])
        self.assertTrue(payload["preview"]["omissions"])
        self.assertLess(len(output), 50000)

    def test_installed_template_is_explained_without_duplicate_internal_config(self) -> None:
        template = Path(__file__).parents[1] / "jarvishep2" / "project_template" / "bin" / "quickstart_random.yaml"
        status, payload, _ = self.inspect([str(template)])
        self.assertEqual(status, 0)
        self.assertEqual(payload["runtime"]["workers"], 4)
        self.assertTrue(payload["paths"]["environment_defaults_yaml"].endswith("deps/environment_default.yaml"))
        self.assertNotIn("config", payload)

    def test_sections_keep_only_requested_data_and_preserve_validation(self) -> None:
        status, complete, _ = self.inspect()
        self.assertEqual(status, 0)
        fields = {
            "paths": {"paths", "path_notes"}, "runtime": {"runtime", "runtime_decisions"},
            "sources": {"sources"}, "validation": set(),
            "card": {"raw_task_card", "merged_envreqs"}, "policies": {"policies", "policy_notes"},
        }
        all_data_fields = set().union(*fields.values())
        for section, selected in fields.items():
            with self.subTest(section=section):
                status, payload, _ = self.inspect([str(self.card), "--section", section, "--json"])
                self.assertEqual(status, 0)
                self.assertEqual(set(payload) & all_data_fields, selected)
                self.assertEqual(payload["issues"], complete["issues"])
                self.assertEqual(payload["issue_count"], complete["issue_count"])
                self.assertTrue(payload["metadata"]["validation_completed"])
                self.assertEqual(payload["metadata"]["requested_section"], section)
                self.assertEqual(set(payload["metadata"]["omitted_sections"]), set(fields) - {section})
                self.assertIn("not absent", payload["metadata"]["section_scope"])
                for field in selected:
                    self.assertEqual(payload[field], complete[field])

    def test_section_filters_do_not_hide_configuration_errors_or_change_exit_code(self) -> None:
        self.document["Sampling"]["Bounds"]["point_number"] = 0
        self.save()
        status, complete, _ = self.inspect()
        self.assertEqual(status, 2)
        for section in ("runtime", "paths", "validation", "card", "policies", "sources"):
            with self.subTest(section=section):
                status, payload, _ = self.inspect([str(self.card), "--section", section, "--json"])
                self.assertEqual(status, 2)
                self.assertFalse(payload["ok"])
                self.assertEqual(payload["issues"], complete["issues"])
                self.assertEqual(payload["issue_count"], len(complete["issues"]))

    def test_section_load_failure_and_invalid_section_return_structured_errors(self) -> None:
        status, payload, _ = self.inspect([str(self.root / "missing.yaml"), "--section", "runtime", "--json"])
        self.assertEqual(status, 1)
        self.assertEqual(payload["issues"][0]["code"], "JV2-LOAD-001")
        self.assertEqual(payload["metadata"]["requested_section"], "runtime")
        self.assertFalse(payload["metadata"]["validation_completed"])
        status, payload, _ = self.inspect([str(self.card), "--section", "unknown", "--json"])
        self.assertEqual(status, 2)
        self.assertFalse(payload["ok"])
        self.assertIn("invalid choice", payload["issues"][0]["message"])

    def test_large_card_preview_budget_preserves_effective_runtime_settings(self) -> None:
        self.document["Operas"]["Modules"][0]["kwargs"] = {
            f"payload_{index}": "a" * 4000 for index in range(20)
        }
        self.save()
        status, complete, _ = self.inspect()
        self.assertEqual(status, 0)
        self.assertTrue(complete["preview"]["truncated"])
        self.assertEqual(complete["runtime"]["workers"], 0)
        self.assertEqual(complete["runtime"]["mode"], "redis")
        self.assertEqual(complete["policies"]["redis"]["port"], 6379)
        self.assertEqual(complete["sources"]["runtime_fields"]["$.Runtime.workers"]["layer"], "task")
        status, selected, output = self.inspect([str(self.card), "--section", "runtime"])
        self.assertEqual(status, 0)
        self.assertFalse(selected["preview"]["truncated"])
        self.assertNotIn("raw_task_card", selected)
        self.assertLess(len(output), 5000)


if __name__ == "__main__":
    unittest.main()
