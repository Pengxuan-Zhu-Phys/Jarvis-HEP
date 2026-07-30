#!/usr/bin/env python3
"""Regression coverage for actionable task-card diagnostics."""

from __future__ import annotations

import unittest

from jarvishep2.diagnostic_guidance import guidance_for
from jarvishep2.task_validation import ValidationReport, issue


class DiagnosticGuidanceTests(unittest.TestCase):
    def test_every_issue_gets_fallback_guidance(self) -> None:
        diagnostic = issue("error", "JV2-FUTURE-001", "Future.setting", "a future constraint failed")
        self.assertTrue(diagnostic.suggestion)
        self.assertIn("constraint", diagnostic.suggestion or "")

    def test_known_method_guidance_is_copyable(self) -> None:
        suggestion, example = guidance_for(
            "JV2-MTH-020", "Sampling['Point number']", "Random requires a point number"
        )
        self.assertIn("Point number", suggestion)
        self.assertEqual(example, "Point number: 100")

    def test_warning_promotion_keeps_guidance_and_example(self) -> None:
        report = ValidationReport([issue(
            "warning", "JV2-DEAD-001", "Runtime.Subprocess", "key is ignored in V2"
        )])
        before = report.issues[0]
        report.promote_warnings_to_errors()
        after = report.issues[0]
        self.assertEqual(after.level, "error")
        self.assertEqual(after.suggestion, before.suggestion)
        self.assertEqual(after.example, before.example)


if __name__ == "__main__":
    unittest.main()
