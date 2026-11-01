#!/usr/bin/env python3
"""D25.7 DX: naming a slow test file on the CLI still collects those tests."""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

from conftest import drop_default_slow_markexpr

REPO_ROOT = Path(__file__).resolve().parents[1]


class DefaultSlowFilterTests(unittest.TestCase):
    def test_only_a_named_test_file_drops_the_default_filter(self) -> None:
        mcmc = "tests/test_mcmc_sampler.py"
        cases = [
            # Bare pytest keeps the fast gate.
            ("not slow", [], "not slow"),
            ("not slow", ["tests/"], "not slow"),
            ("not slow", ["-q"], "not slow"),
            # A user -m is never overridden.
            ("not slow", ["-m", "not slow", mcmc], "not slow"),
            ("slow", ["-m", "slow", mcmc], "slow"),
            # Option values are not selected test paths.
            ("not slow", ["--deselect=tests/fake.py::test_a", "tests/"], "not slow"),
            ("not slow", ["--deselect", "tests/fake.py::test_a", "tests/real.py"], ""),
            ("not slow", [mcmc], ""),
            ("not slow", ["-q", mcmc + "::MCMCSamplerTests::test_foo"], ""),
        ]
        for option in ("--deselect", "--ignore", "--junitxml", "--override-ini", "-o"):
            cases.append(("not slow", ["tests/", option, "tests/fake.py::test_a"], "not slow"))
        for markexpr, args, expected in cases:
            with self.subTest(args=args):
                self.assertEqual(drop_default_slow_markexpr(markexpr, args), expected)

    def test_named_slow_file_is_collected_under_default_addopts(self) -> None:
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "pytest",
                "--collect-only",
                "-q",
                "tests/test_mcmc_sampler.py",
            ],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
            check=False,
        )
        output = proc.stdout + proc.stderr
        self.assertEqual(proc.returncode, 0, output)
        self.assertNotIn("deselected", output)
        self.assertRegex(output, r"[1-9]\d* test")


if __name__ == "__main__":
    unittest.main()
