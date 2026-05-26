#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import unittest

import numpy as np
import pandas as pd


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.likelihood import LogLikelihood  # noqa: E402


class _NoopLogger:
    def bind(self, **_kwargs):
        return self

    def info(self, *_args, **_kwargs):
        pass

    def warning(self, *_args, **_kwargs):
        pass


class LogLikelihoodTests(unittest.TestCase):
    def _sample_info(self):
        return {
            "logger": _NoopLogger(),
            "logger_name": "Sample@test-likelihood",
        }

    def test_default_total_logl_is_sum_of_named_terms(self):
        values = {"x": 2.0, "y": 3.0}
        result = LogLikelihood(
            [
                {"name": "LogL_x", "expression": "x"},
                {"name": "LogL_y", "expression": "y"},
            ]
        ).calculate(values, self._sample_info())

        self.assertEqual(result["LogL_x"], 2.0)
        self.assertEqual(result["LogL_y"], 3.0)
        self.assertEqual(result["LogL"], 5.0)
        self.assertEqual(values["LogL_x"], 2.0)
        self.assertEqual(values["LogL_y"], 3.0)
        self.assertEqual(values["LogL"], 5.0)

    def test_explicit_logl_expression_defines_total_and_keeps_terms(self):
        values = {"x": 2.0, "y": 3.0}
        result = LogLikelihood(
            [
                {"name": "LogL", "expression": "LogL_x - 2 * LogL_y"},
                {"name": "LogL_x", "expression": "x"},
                {"name": "LogL_y", "expression": "y"},
            ]
        ).calculate(values, self._sample_info())

        self.assertEqual(result["LogL_x"], 2.0)
        self.assertEqual(result["LogL_y"], 3.0)
        self.assertEqual(result["LogL"], -4.0)
        self.assertEqual(values["LogL_x"], 2.0)
        self.assertEqual(values["LogL_y"], 3.0)
        self.assertEqual(values["LogL"], -4.0)

    def test_calculate4dnn_uses_explicit_logl_total(self):
        loglike = LogLikelihood(
            [
                {"name": "LogL", "expression": "LogL_x - LogL_y"},
                {"name": "LogL_x", "expression": "x"},
                {"name": "LogL_y", "expression": "y"},
            ]
        )

        actual = loglike.calculate4dnn(pd.DataFrame([{"x": 2.0, "y": 3.0}]))

        np.testing.assert_allclose(actual, np.array([-1.0]))


if __name__ == "__main__":
    unittest.main()
