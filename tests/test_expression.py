#!/usr/bin/env python3
"""Tests for the shared compile-once expression runtime."""

from __future__ import annotations

import math
import unittest

import numpy as np

from jarvishep2 import ExpressionContext as PublicExpressionContext
from jarvishep2.expression import (
    CompiledExpression,
    ExpressionContext,
    MissingExpressionVariablesError,
)
from jarvishep2.inner_func import (
    EXPRESSION_CONSTANTS,
    EXPRESSION_FUNCTION_NAMES,
)
from jarvishep2.likelihood import LogLikelihoodEvaluator


class ExpressionContextTests(unittest.TestCase):
    def test_context_is_part_of_the_package_public_api(self) -> None:
        self.assertIs(PublicExpressionContext, ExpressionContext)

    def test_compile_caches_by_expression_and_symbol_contract(self) -> None:
        context = ExpressionContext()
        first = context.compile("x + y", symbols=("x", "y"))
        second = context.compile("x + y", symbols=("y", "x"))

        self.assertIs(first, second)
        self.assertIsInstance(first, CompiledExpression)
        self.assertEqual(first.variable_names, ("x", "y"))
        self.assertEqual(first.evaluate({"x": 1.5, "y": 2.5}), 4.0)
        self.assertEqual(context.cache_info().misses, 1)
        self.assertEqual(context.cache_info().hits, 1)

    def test_v1_constants_and_numeric_functions_share_one_language(self) -> None:
        context = ExpressionContext()
        compiled = context.compile("LogGauss(sin(x * Pi), 1, 2)", symbols=("x",))

        self.assertAlmostEqual(compiled.evaluate({"x": 0.5}), 0.0, places=12)

    def test_missing_variables_are_structured(self) -> None:
        compiled = ExpressionContext().compile("x + y", symbols=("x", "y"))

        with self.assertRaises(MissingExpressionVariablesError) as raised:
            compiled.evaluate({"x": 1.0})

        self.assertEqual(raised.exception.expression, "x + y")
        self.assertEqual(raised.exception.missing, ("y",))

    def test_function_update_invalidates_cache(self) -> None:
        context = ExpressionContext(functions={"twice": lambda x: 2 * x})
        first = context.compile("twice(x)", symbols=("x",))
        self.assertEqual(first.evaluate({"x": 3}), 6)

        context.update_functions({"twice": lambda x: 3 * x})
        second = context.compile("twice(x)", symbols=("x",))
        self.assertIsNot(first, second)
        self.assertEqual(second.evaluate({"x": 3}), 9)

    def test_pi_variants_are_equivalent(self) -> None:
        context = ExpressionContext()
        value = context.evaluate("Pi + pi + PI", {})
        self.assertAlmostEqual(float(value), 3.0 * math.pi, places=12)


class V1LightweightFunctionCompatibilityTests(unittest.TestCase):
    V1_FUNCTION_NAMES = {
        "log",
        "exp",
        "ln",
        "sin",
        "cos",
        "tan",
        "sec",
        "csc",
        "cot",
        "sinc",
        "asin",
        "acos",
        "atan",
        "asec",
        "acsc",
        "acot",
        "atan2",
        "sinh",
        "cosh",
        "tanh",
        "sech",
        "csch",
        "coth",
        "asinh",
        "acosh",
        "atanh",
        "acoth",
        "asech",
        "acsch",
        "sqrt",
        "Min",
        "Max",
        "root",
        "Abs",
        "Gauss",
        "Normal",
        "LogGauss",
        "Heaviside",
    }

    def test_complete_v1_function_and_constant_surface_is_present(self) -> None:
        self.assertEqual(EXPRESSION_FUNCTION_NAMES, self.V1_FUNCTION_NAMES)
        self.assertTrue({"Pi", "E", "Inf"}.issubset(EXPRESSION_CONSTANTS))

    def test_every_v1_lightweight_function_evaluates(self) -> None:
        cases = [
            ("log(x)", math.e, 1.0),
            ("exp(x)", 1.0, math.e),
            ("ln(x)", math.e, 1.0),
            ("sin(x)", 0.3, math.sin(0.3)),
            ("cos(x)", 0.3, math.cos(0.3)),
            ("tan(x)", 0.3, math.tan(0.3)),
            ("sec(x)", 0.3, 1.0 / math.cos(0.3)),
            ("csc(x)", 0.7, 1.0 / math.sin(0.7)),
            ("cot(x)", 0.7, 1.0 / math.tan(0.7)),
            ("sinc(x)", 0.2, math.sin(0.2) / 0.2),
            ("asin(x)", 0.3, math.asin(0.3)),
            ("acos(x)", 0.3, math.acos(0.3)),
            ("atan(x)", 0.3, math.atan(0.3)),
            ("asec(x)", 2.0, math.acos(0.5)),
            ("acsc(x)", 2.0, math.asin(0.5)),
            ("acot(x)", 2.0, math.atan(0.5)),
            ("atan2(x, 1)", 0.3, math.atan2(0.3, 1.0)),
            ("sinh(x)", 0.3, math.sinh(0.3)),
            ("cosh(x)", 0.3, math.cosh(0.3)),
            ("tanh(x)", 0.3, math.tanh(0.3)),
            ("sech(x)", 0.3, 1.0 / math.cosh(0.3)),
            ("csch(x)", 0.7, 1.0 / math.sinh(0.7)),
            ("coth(x)", 0.7, 1.0 / math.tanh(0.7)),
            ("asinh(x)", 0.3, math.asinh(0.3)),
            ("acosh(x)", 2.0, math.acosh(2.0)),
            ("atanh(x)", 0.3, math.atanh(0.3)),
            ("acoth(x)", 2.0, math.atanh(0.5)),
            ("asech(x)", 0.5, math.acosh(2.0)),
            ("acsch(x)", 2.0, math.asinh(0.5)),
            ("sqrt(x)", 4.0, 2.0),
            ("Min(x, 3)", 2.0, 2.0),
            ("Max(x, 3)", 2.0, 3.0),
            ("root(x, 3)", 8.0, 2.0),
            ("Abs(x)", -3.0, 3.0),
            ("Gauss(x, 0, 1)", 0.5, math.exp(-0.125)),
            (
                "Normal(x, 0, 1)",
                0.5,
                math.exp(-0.125) / math.sqrt(2.0 * math.pi),
            ),
            ("LogGauss(x, 0, 1)", 0.5, -0.125),
            ("Heaviside(x)", 0.0, 0.5),
        ]
        context = ExpressionContext()
        for expression, x, expected in cases:
            with self.subTest(expression=expression):
                value = context.evaluate(expression, {"x": x})
                self.assertAlmostEqual(float(value), expected, places=12)

    def test_v1_constants_and_vector_gaussians(self) -> None:
        context = ExpressionContext()
        self.assertEqual(context.evaluate("log(E) + Min(Inf, 2)", {}), 3)
        self.assertAlmostEqual(
            float(context.evaluate("E", {"E": 123.0})),
            math.e,
            places=12,
        )

        x = np.asarray([-1.0, 0.0, 1.0])
        log_gauss = context.evaluate("LogGauss(x, 0, 1)", {"x": x})
        np.testing.assert_allclose(log_gauss, [-0.5, 0.0, -0.5])

    def test_released_gmfit_acceptance_expression(self) -> None:
        expression = (
            "flag * Heaviside(vac) * "
            "Max(Heaviside(3 - Abs(mh - 125)), Heaviside(3 - Abs(mH - 125)))"
        )
        value = ExpressionContext().evaluate(
            expression,
            {"flag": 1, "vac": 1, "mh": 126.0, "mH": 140.0},
        )
        self.assertEqual(float(value), 1.0)

        likelihood = LogLikelihoodEvaluator(
            [{"name": "LogL_accept", "expression": expression}]
        ).evaluate({"flag": 1, "vac": 1, "mh": 126.0, "mH": 140.0})
        self.assertEqual(likelihood, {"LogL_accept": 1.0, "LogL": 1.0})


class SharedConsumerTests(unittest.TestCase):
    def test_likelihood_uses_compiled_expressions_with_ordered_dependencies(self) -> None:
        evaluator = LogLikelihoodEvaluator(
            [
                {"name": "term", "expression": "x + 1"},
                {"name": "LogL", "expression": "term * 2"},
            ]
        )

        values = evaluator.evaluate({"x": 2.0})

        self.assertEqual(values, {"term": 3.0, "LogL": 6.0})
        self.assertTrue(
            all(
                isinstance(compiled, CompiledExpression)
                for _, compiled in evaluator._compiled
            )
        )


if __name__ == "__main__":
    unittest.main()
