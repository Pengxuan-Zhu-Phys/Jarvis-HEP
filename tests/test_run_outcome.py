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
    def test_status_and_exit_code(self) -> None:
        cases = (
            # Historical bug: CLI used submitted/archived counts and exited 0 on all-fail.
            (RunOutcome(submitted=10, completed=0, failed=10, archived=10), "failed", EXIT_RUN_FAILED),
            (
                RunOutcome(submitted=10, completed=7, failed=3, archived=10),
                "partial_failure",
                EXIT_RUN_FAILED,
            ),
            (RunOutcome(submitted=10, completed=10, failed=0, archived=10), "success", EXIT_OK),
            (
                RunOutcome.from_counters(submitted=3, completed=1, failed=0, interrupted=True),
                "interrupted",
                EXIT_INTERRUPTED,
            ),
        )
        for outcome, status, exit_code in cases:
            with self.subTest(status=status):
                self.assertEqual(outcome.status, status)
                self.assertEqual(outcome.exit_code, exit_code)
                self.assertEqual(outcome.ok, exit_code == EXIT_OK)
                payload = outcome.to_dict()
                self.assertEqual(payload["exit_code"], exit_code)
                self.assertEqual(payload["ok"], exit_code == EXIT_OK)


if __name__ == "__main__":
    unittest.main()
