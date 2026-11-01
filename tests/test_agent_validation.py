"""Agent validation responses remain actionable without executing a scan."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from jarvishep2.cli.man import dispatch_man
from jarvishep2.cli.run import dispatch_validate
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE
from jarvishep2.task_validation import ValidationReport, issue


_VALID_CARD = """\
Scan: {name: agent_validation}
Sampling:
  Method: Random
  Bounds: {point_number: 10, seed: 7}
  Variables:
    - name: x
      distribution: {type: Flat, parameters: {min: 0.0, max: 1.0}}
    - name: y
      distribution: {type: Flat, parameters: {min: 0.0, max: 1.0}}
EnvReqs:
  V2: {workers: 1}
Operas:
  Modules:
    - name: eggbox
      operator: jarvishep2.testing.eggbox.eggbox2d_numpy
      input:
        - {name: x, expression: x}
        - {name: y, expression: y}
      output:
        - {name: z, entry: z}
"""


class AgentValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.directory = tempfile.TemporaryDirectory(prefix="jarvis-agent-validation-")
        self.addCleanup(self.directory.cleanup)
        self.card = Path(self.directory.name) / "task.yaml"
        self.card.write_text(_VALID_CARD, encoding="utf-8")

    def validate(self, path: str | None = None, **kwargs):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            with mock.patch("jarvishep2.core.Jarvis2Core") as core:
                with mock.patch("jarvishep2.redis_queue.RedisQueue") as redis:
                    code = dispatch_validate(
                        str(self.card) if path is None else path,
                        as_json=True, **kwargs,
                    )
                    core.assert_not_called()
                    redis.assert_not_called()
        payload = json.loads(stdout.getvalue())
        self.assertNotIn("Traceback", stderr.getvalue())
        self.assertFalse(payload["metadata"]["execution_performed"])
        self.assertFalse(payload["metadata"]["runtime_success_verified"])
        self.assertFalse(payload["metadata"]["physics_correctness_verified"])
        self.assertTrue({"ok", "task_yaml", "scan_name", "issues"} <= payload.keys())
        for diagnostic in payload["issues"]:
            self.assert_navigation(diagnostic)
        return code, payload

    def assert_navigation(self, diagnostic) -> None:
        argv = diagnostic["manual_argv"]
        self.assertEqual(argv[:2], ["Jarvis", "man"])
        self.assertEqual(argv[-1], "--json")
        self.assertTrue(diagnostic["manual_command"].startswith("Jarvis man "))
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            code = dispatch_man(argv[2:])
        self.assertEqual(code, EXIT_OK, msg=(diagnostic, stderr.getvalue()))
        self.assertEqual(json.loads(stdout.getvalue())["diagnostics"], [])

    def test_success_describes_the_limit_of_validation_evidence(self) -> None:
        code, payload = self.validate()
        self.assertEqual(code, EXIT_OK)
        self.assertTrue(payload["ok"])
        self.assertTrue(payload["metadata"]["validation_completed"])
        self.assertEqual(payload["scan_name"], "agent_validation")
        self.assertEqual(payload["issues"], [])

    def test_human_output_is_concise_on_stderr(self) -> None:
        warning = ValidationReport([
            issue("warning", "JV2-FUTURE-001", "$.Sampling.Bounds", "Review the proposed budget")
        ])
        cases = (
            # The base report has no notes API; success still prints.
            ("success", None, EXIT_OK, ["Config validation successful"]),
            (
                "warnings",
                {"return_value": warning},
                EXIT_OK,
                [
                    "Config validation successful",
                    "Config validation warnings",
                    "Review the proposed budget",
                ],
            ),
            (
                "internal error",
                {"side_effect": RuntimeError("validator defect")},
                EXIT_RUN_FAILED,
                ["JV2-VAL-001", "validator defect"],
            ),
        )
        for label, validator, exit_code, fragments in cases:
            with self.subTest(label):
                stdout, stderr = io.StringIO(), io.StringIO()
                patch_validator = (
                    mock.patch("jarvishep2.task_validation.validate_task_config", **validator)
                    if validator
                    else contextlib.nullcontext()
                )
                with (
                    contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr),
                    mock.patch("jarvishep2.core.Jarvis2Core") as core,
                    mock.patch("jarvishep2.redis_queue.RedisQueue") as redis,
                    patch_validator,
                ):
                    code = dispatch_validate(str(self.card))
                self.assertEqual(code, exit_code)
                self.assertEqual(stdout.getvalue(), "")
                for fragment in fragments:
                    self.assertIn(fragment, stderr.getvalue())
                self.assertNotIn("Traceback", stderr.getvalue())
                core.assert_not_called()
                redis.assert_not_called()

    def test_known_schema_and_semantic_issues_have_callable_manual_commands(self) -> None:
        self.card.write_text(
            _VALID_CARD.replace("point_number: 10", "point_number: 0")
            .replace("name: agent_validation", "naem: agent_validation"),
            encoding="utf-8",
        )
        code, payload = self.validate()
        self.assertEqual(code, EXIT_USAGE)
        self.assertFalse(payload["ok"])
        self.assertTrue(payload["metadata"]["validation_completed"])
        self.assertTrue(payload["issues"])
        self.assertTrue(all("--code" in item["manual_argv"] for item in payload["issues"]))

    def test_unloadable_input_keeps_its_exit_code(self) -> None:
        cases = (
            (str(self.card.with_name("missing.yaml")), None, EXIT_RUN_FAILED, "JV2-LOAD-001"),
            (self.directory.name, None, EXIT_RUN_FAILED, "JV2-LOAD-001"),
            (None, "Scan: [broken\n", EXIT_USAGE, "JV2-YAML-001"),
            # A missing argument is structured for direct dispatch too.
            ("", None, EXIT_USAGE, None),
        )
        for path, text, exit_code, issue_code in cases:
            with self.subTest(path=path, text=text):
                if text is not None:
                    self.card.write_text(text, encoding="utf-8")
                code, payload = self.validate(path)
                self.assertEqual(code, exit_code)
                self.assertFalse(payload["ok"])
                self.assertFalse(payload["metadata"]["validation_completed"])
                if issue_code:
                    self.assertEqual(payload["issues"][0]["code"], issue_code)

    def test_loader_and_validator_exceptions_are_run_failures(self) -> None:
        loader = "jarvishep2.task_config.load_task_yaml"
        validator = "jarvishep2.task_validation.validate_task_config"
        cases = (
            (loader, PermissionError("permission denied"), "JV2-LOAD-003", "suggestion", "readable"),
            # A loader defect is distinct from a configuration error.
            (loader, RuntimeError("loader defect"), "JV2-LOAD-004", "message", "loader defect"),
            # Even a valid card fails when the validator itself raises.
            (validator, RuntimeError("registry unavailable"), "JV2-VAL-001", "message", "registry unavailable"),
        )
        for target, error, issue_code, field, fragment in cases:
            with self.subTest(issue_code):
                with mock.patch(target, side_effect=error):
                    code, payload = self.validate()
                self.assertEqual(code, EXIT_RUN_FAILED)
                self.assertFalse(payload["ok"])
                self.assertEqual(payload["issues"][0]["code"], issue_code)
                self.assertIn(fragment, payload["issues"][0][field])
                self.assertFalse(payload["metadata"]["validation_completed"])

    def test_unknown_diagnostic_uses_yaml_overview_without_guessing_a_topic(self) -> None:
        report = ValidationReport([issue("error", "JV2-FUTURE-001", "Future.nested[0].leaf", "future error")])
        with mock.patch("jarvishep2.task_validation.validate_task_config", return_value=report):
            code, payload = self.validate()
        self.assertEqual(code, EXIT_USAGE)
        self.assertEqual(payload["issues"][0]["manual_argv"], ["Jarvis", "man", "yaml", "--json"])


if __name__ == "__main__":
    unittest.main()
