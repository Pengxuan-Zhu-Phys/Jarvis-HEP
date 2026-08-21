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
    def test_bare_pytest_keeps_fast_gate(self) -> None:
        self.assertEqual(drop_default_slow_markexpr("not slow", []), "not slow")
        self.assertEqual(drop_default_slow_markexpr("not slow", ["tests/"]), "not slow")
        self.assertEqual(drop_default_slow_markexpr("not slow", ["-q"]), "not slow")

    def test_user_dash_m_is_never_overridden(self) -> None:
        self.assertEqual(
            drop_default_slow_markexpr("not slow", ["-m", "not slow", "tests/test_mcmc_sampler.py"]),
            "not slow",
        )
        self.assertEqual(
            drop_default_slow_markexpr("slow", ["-m", "slow", "tests/test_mcmc_sampler.py"]),
            "slow",
        )

    def test_named_py_file_drops_default_filter(self) -> None:
        self.assertEqual(
            drop_default_slow_markexpr("not slow", ["tests/test_mcmc_sampler.py"]),
            "",
        )
        self.assertEqual(
            drop_default_slow_markexpr(
                "not slow",
                ["-q", "tests/test_mcmc_sampler.py::MCMCSamplerTests::test_foo"],
            ),
            "",
        )

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
