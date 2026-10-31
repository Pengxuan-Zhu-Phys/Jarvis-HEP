#!/usr/bin/env python3
"""Regression coverage for actionable task-card diagnostics."""

from __future__ import annotations

import unittest

import yaml

from jarvishep2.diagnostic_guidance import guidance_for
from jarvishep2.task_validation import ValidationReport, issue


class DiagnosticGuidanceTests(unittest.TestCase):
    def test_every_issue_gets_fallback_guidance(self) -> None:
        diagnostic = issue("error", "JV2-FUTURE-001", "Future.setting", "a future constraint failed")
        self.assertTrue(diagnostic.suggestion)
        self.assertIn("constraint", diagnostic.suggestion or "")

    def test_known_method_guidance_is_copyable(self) -> None:
        suggestion, example = guidance_for(
            "JV2-MTH-020",
            "Sampling.Bounds['point_number']",
            "Random requires a point number",
        )
        self.assertIn("point_number", suggestion)
        self.assertIn("Bounds", example or "")
        self.assertIn("point_number: 100", example or "")

    def test_runtime_guidance_uses_the_relevant_block(self) -> None:
        from jarvishep2.cardload.contracts.operational import (
            _validate_archiver, _validate_sample_directory,
        )
        from jarvishep2.task_validation import validate_task_config

        cases = [
            (_validate_archiver({"unknown": True}, "EnvReqs.V2.archiver")[0], "archiver"),
            (_validate_sample_directory({"unknown": True}, "EnvReqs.V2.sample_directory")[0], "sample_directory"),
        ]
        report = validate_task_config({"EnvReqs": {"V2": {"check_modules": {"unknown": True}}}})
        cases += [(item, "check_modules") for item in report.issues if item.code == "JV2-ENV-061"]
        self.assertEqual(len(cases), 3)
        for diagnostic, block in cases:
            with self.subTest(code=diagnostic.code):
                self.assertNotIn("workers:", diagnostic.example or "")
                self.assertIn(block, diagnostic.example or "")
                self.assertIsInstance(yaml.safe_load(diagnostic.example), dict)

    def test_warning_promotion_keeps_guidance_and_example(self) -> None:
        report = ValidationReport([issue(
            "warning", "JV2-FUTURE-001", "EnvReqs.V2.future", "future setting is not active"
        )])
        before = report.issues[0]
        report.promote_warnings_to_errors()
        after = report.issues[0]
        self.assertEqual(after.level, "error")
        self.assertEqual(after.suggestion, before.suggestion)
        self.assertEqual(after.example, before.example)

    def test_var_031_suggestion_names_missing_distribution_parameters(self) -> None:
        """D21.14: do not fall back to the generic Variables-list suggestion."""
        path = "Sampling.Variables[0].distribution.parameters"
        message = "missing required key(s) for type Poisson: lambda"
        suggestion, _example = guidance_for("JV2-VAR-031", path, message)
        self.assertIn("lambda", suggestion)
        self.assertIn("parameter", suggestion.lower())
        self.assertNotIn("non-empty", suggestion.lower())
        diagnostic = issue("error", "JV2-VAR-031", path, message)
        self.assertIn("lambda", diagnostic.suggestion or "")


if __name__ == "__main__":
    unittest.main()
