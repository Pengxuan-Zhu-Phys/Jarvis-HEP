#!/usr/bin/env python3
"""D11.1 RunOutcome exit-code contract."""

from __future__ import annotations

import unittest

from jarvishep2.run_outcome import (
    EXIT_INTERRUPTED,
    EXIT_OK,
    EXIT_RUN_FAILED,
    RunOutcome,
)


class RunOutcomeTests(unittest.TestCase):
    def test_all_failed_is_non_zero(self) -> None:
        outcome = RunOutcome(submitted=10, completed=0, failed=10, archived=10)
        self.assertEqual(outcome.status, "failed")
        self.assertEqual(outcome.exit_code, EXIT_RUN_FAILED)
        self.assertFalse(outcome.ok)

    def test_partial_failure_is_non_zero(self) -> None:
        outcome = RunOutcome(submitted=10, completed=7, failed=3, archived=10)
        self.assertEqual(outcome.status, "partial_failure")
        self.assertEqual(outcome.exit_code, EXIT_RUN_FAILED)

    def test_all_success_is_zero(self) -> None:
        outcome = RunOutcome(submitted=10, completed=10, failed=0, archived=10)
        self.assertEqual(outcome.status, "success")
        self.assertEqual(outcome.exit_code, EXIT_OK)
        self.assertTrue(outcome.ok)

    def test_archived_rows_do_not_imply_success(self) -> None:
        # Historical bug: CLI used submitted/archived counts and exited 0 on all-fail.
        outcome = RunOutcome(submitted=5, completed=0, failed=5, archived=5)
        self.assertNotEqual(outcome.exit_code, EXIT_OK)

    def test_interrupt_exit_130(self) -> None:
        outcome = RunOutcome.from_counters(
            submitted=3,
            completed=1,
            failed=0,
            interrupted=True,
        )
        self.assertEqual(outcome.status, "interrupted")
        self.assertEqual(outcome.exit_code, EXIT_INTERRUPTED)

    def test_to_dict_includes_exit_code(self) -> None:
        payload = RunOutcome(submitted=1, completed=1, failed=0).to_dict()
        self.assertEqual(payload["exit_code"], EXIT_OK)
        self.assertTrue(payload["ok"])


if __name__ == "__main__":
    unittest.main()
