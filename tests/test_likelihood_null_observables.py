#!/usr/bin/env python3
"""Null calculator outputs must yield LogL=-inf, not sample TypeError failures."""

from __future__ import annotations

import math
import unittest

from jarvishep2.inner_func import LogGauss
from jarvishep2.likelihood import LogLikelihoodEvaluator


class LikelihoodNullObservablesTests(unittest.TestCase):
    def test_loggauss_none_is_nan(self) -> None:
        value = LogGauss(None, 0.12, 0.0012)
        self.assertTrue(math.isnan(float(value)))

    def test_calculate_null_omega_returns_neginf_without_raise(self) -> None:
        """micrOMEGAs may emit JSON null for Omega_h2 on unphysical points."""
        evaluator = LogLikelihoodEvaluator(
            [
                {
                    "name": "LogL_Omega",
                    "expression": "LogGauss(Omega_h2, 0.1200, 0.0012)",
                }
            ]
        )
        info = {"observables": {"Omega_h2": None}, "uuid": "null-omega"}
        total = evaluator.calculate(info)
        self.assertEqual(total, float("-inf"))
        self.assertEqual(info["observables"]["LogL"], float("-inf"))
        self.assertEqual(info["observables"]["LogL_Omega"], float("-inf"))
        # Raw portal null is preserved; only LogL terms are written.
        self.assertIsNone(info["observables"]["Omega_h2"])

    def test_evaluate_null_is_neginf(self) -> None:
        evaluator = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "LogGauss(x, 0, 1)"}]
        )
        terms = evaluator.evaluate({"x": None})
        self.assertEqual(terms["LogL"], float("-inf"))


if __name__ == "__main__":
    unittest.main()
