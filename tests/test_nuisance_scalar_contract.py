#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import unittest

import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.nuisance_LogLikelihood import NuisanceExpressionRegistry  # noqa: E402
from jarvishep.Module.nuisance_passCondition import NuisancePassCondition  # noqa: E402


class NuisanceScalarContractTests(unittest.TestCase):
    def test_loglikelihood_scalar_success(self):
        reg = NuisanceExpressionRegistry()
        reg.set_config(name="NLL_base", expression="x + y")

        value = reg.eval({"x": 1.5, "y": 2.0})
        self.assertEqual(value, 3.5)

    def test_loglikelihood_array_input_rejected_with_clear_message(self):
        reg = NuisanceExpressionRegistry()
        reg.set_config(name="NLL_base", expression="x + y")

        with self.assertRaisesRegex(
            ValueError,
            r"Nuisance LogLikelihood 'NLL_base' expects scalar input for 'x'",
        ):
            reg.eval({"x": np.array([1.0, 2.0]), "y": 1.0})

    def test_passcondition_scalar_success(self):
        reg = NuisancePassCondition()
        reg.set_config(name="PC_cut", expression="x > 0")

        self.assertTrue(reg.eval({"x": 0.1}))
        self.assertFalse(reg.eval({"x": -0.1}))

    def test_passcondition_array_input_rejected_with_clear_message(self):
        reg = NuisancePassCondition()
        reg.set_config(name="PC_cut", expression="x > 0")

        with self.assertRaisesRegex(
            ValueError,
            r"Nuisance PassCondition 'PC_cut' expects scalar input for 'x'",
        ):
            reg.eval({"x": np.array([1.0, -1.0])})


if __name__ == "__main__":
    unittest.main()
