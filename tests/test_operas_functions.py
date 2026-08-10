#!/usr/bin/env python3
"""Dynamic Jarvis-Operas discovery with the released V1 Modules YAML surface."""

from __future__ import annotations

import asyncio
import unittest
from unittest.mock import patch

from jarvishep2.io_portal import evaluate_io_expression
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.operas import resolve_operator
from jarvishep2.operas_functions import (
    build_operas_expression_context,
    discover_operas_expression_functions,
    expression_uses_operas_function,
    expression_uses_operas_reference,
    operas_expression_functions_required,
    operas_expression_references_required,
    reset_operas_discovery_for_tests,
)
from jarvishep2.Sampling.sampling_utils import evaluate_selection
from jarvishep2.worker_config import build_worker_config

try:
    import jarvis_operas
    from jarvis_operas import OperasRegistry, oper
except ImportError:  # pragma: no cover - optional-dependency test environment
    jarvis_operas = None
    OperasRegistry = None
    oper = None


@unittest.skipIf(jarvis_operas is None, "Jarvis-Operas is not installed")
class DynamicOperasFunctionTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_operas_discovery_for_tests()
        self.registry = OperasRegistry()

        @oper("double", namespace="hep2test", registry=self.registry)
        def double(x):
            return 2 * x

    def test_all_registered_functions_are_discovered_by_qualified_name(self) -> None:
        context = build_operas_expression_context(
            registry=self.registry,
            discover_plugins=False,
            required=True,
        )

        self.assertEqual(context.evaluate("hep2test.double(x)", {"x": 4}), 8)

    def test_entrypoint_discovery_runs_once_per_process_registry(self) -> None:
        def discover(target_registry):
            @oper("plugin_square", namespace="hep2plugin", registry=target_registry)
            def plugin_square(x):
                return x * x

            return ["hep2plugin.plugin_square"]

        with patch.object(jarvis_operas, "discover_entrypoints", side_effect=discover) as mocked:
            first = discover_operas_expression_functions(registry=self.registry, required=True)
            second = discover_operas_expression_functions(registry=self.registry, required=True)

        self.assertEqual(mocked.call_count, 1)
        self.assertEqual(
            first.create_context().evaluate("hep2plugin.plugin_square(x)", {"x": 5}),
            25,
        )
        self.assertIn("hep2plugin.plugin_square", second.registered_names)

    def test_registered_module_operator_triggers_entrypoint_discovery(self) -> None:
        def discover(target_registry):
            @oper("plugin_square", namespace="hep2plugin", registry=target_registry)
            def plugin_square(x):
                return {"z": x * x}

            return ["hep2plugin.plugin_square"]

        with patch.object(
            jarvis_operas,
            "get_global_operas_registry",
            return_value=self.registry,
        ), patch.object(jarvis_operas, "discover_entrypoints", side_effect=discover) as mocked:
            operator = resolve_operator("hep2plugin.plugin_square")

        self.assertEqual(mocked.call_count, 1)
        with patch.object(
            self.registry,
            "resolve_name",
            side_effect=AssertionError("registry queried in module hot path"),
        ), patch.object(
            self.registry,
            "get",
            side_effect=AssertionError("registry queried in module hot path"),
        ):
            self.assertEqual(operator(x=4), {"z": 16})

    def test_acall_wraps_a_snapshotted_sync_operator_without_registry_dispatch(self) -> None:
        with patch.object(
            jarvis_operas,
            "get_global_operas_registry",
            return_value=self.registry,
        ):
            operator = resolve_operator("hep2test.double", call_mode="acall")

        with patch.object(
            self.registry,
            "resolve_name",
            side_effect=AssertionError("registry queried in module hot path"),
        ), patch.object(
            self.registry,
            "get",
            side_effect=AssertionError("registry queried in module hot path"),
        ):
            self.assertEqual(asyncio.run(operator(x=4)), 8)

    def test_evaluation_hot_path_does_not_query_registry(self) -> None:
        context = build_operas_expression_context(
            registry=self.registry,
            discover_plugins=False,
            required=True,
        )
        qualified = context.compile("hep2test.double(x)", symbols=("x",))

        with patch.object(
            self.registry,
            "resolve_name",
            side_effect=AssertionError("registry queried in hot path"),
        ), patch.object(
            self.registry,
            "get",
            side_effect=AssertionError("registry queried in hot path"),
        ):
            for value in range(5):
                self.assertEqual(qualified.evaluate({"x": value}), 2 * value)

    def test_same_context_serves_likelihood_selection_and_portal_dump(self) -> None:
        context = build_operas_expression_context(
            registry=self.registry,
            discover_plugins=False,
            required=True,
        )
        likelihood = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "-hep2test.double(x)"}],
            expression_context=context,
        )

        self.assertEqual(likelihood.evaluate({"x": 3})["LogL"], -6.0)
        self.assertTrue(
            evaluate_selection(
                "hep2test.double(x) > 5",
                {"x": 3},
                context=context,
            )
        )
        self.assertEqual(
            evaluate_io_expression(
                "hep2test.double(x)",
                {"x": 3},
                context=context,
            ),
            6.0,
        )

    def test_likelihood_auto_discovers_qualified_global_function(self) -> None:
        likelihood = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "math.add(x, 2)"}]
        )

        self.assertEqual(likelihood.evaluate({"x": 3})["LogL"], 5.0)


class OperasModulesYamlCompatibilityTests(unittest.TestCase):
    def test_startup_gate_skips_internal_only_expressions(self) -> None:
        self.assertFalse(
            operas_expression_functions_required({"expression": "LogGauss(x, 0, 1)"})
        )
        self.assertTrue(
            operas_expression_functions_required({"expression": "user.external(x)"})
        )
        # D23 aliases stay wired.
        self.assertIs(expression_uses_operas_function, expression_uses_operas_reference)
        self.assertIs(
            operas_expression_functions_required,
            operas_expression_references_required,
        )

    def test_startup_gate_accepts_bare_qualified_constants(self) -> None:
        """D23 §8.1 — bare constants trip the discovery gate; path-like dotted strings do not."""
        cases = [
            ("(mz - pdg.mZ)/pdg.mZ", True),
            ("pdg.hbarc", True),
            ("math.add(x, 2)", True),
            ("sqrt(pdg.mZ**2)", True),
            ("/usr/local/bin/foo.sh", False),
            ("./run.sh", False),
            ("results/out.dat", False),
            (r"C:\tools\a.exe", False),
            ("0.5", False),
            ("1.5e3", False),
            ("x*1.0", False),
            ("sin(t)*2.5", False),
            ("Pi", False),
            ("LogGauss(x, 0, 1)", False),
        ]
        for text, expected in cases:
            with self.subTest(text=text):
                self.assertEqual(
                    expression_uses_operas_reference(text),
                    expected,
                    msg=text,
                )
                self.assertEqual(
                    operas_expression_references_required({"expression": text}),
                    expected,
                    msg=f"nested {text}",
                )

    def test_startup_gate_ignores_calculator_cmd_and_path_strings(self) -> None:
        """Qualified calls in shell/cmd must not force Operas expression discovery."""
        calculator_payload = {
            "calculator_modules": [
                {
                    "name": "EggBox",
                    "source": "&J/deps/assets/inertial/EggBox",
                    "path": "&J/calculators/runtime/program/EggBox/@PackID",
                    "installation": ["cp -r ${source}/* ${path}"],
                    "execution": {
                        "commands": [
                            {"cmd": './eggbox.py --out numpy.savetxt("z.npy")'},
                            "python3 -c \"import numpy; numpy.savetxt('a', [])\"",
                        ],
                        "input": [
                            {
                                "type": "JSON",
                                "actions": [
                                    {
                                        "type": "Dump",
                                        "variables": [
                                            {"name": "xx", "expression": "x * Pi"},
                                        ],
                                    }
                                ],
                            }
                        ],
                    },
                }
            ],
            "likelihood_expressions": [{"name": "LogL", "expression": "LogGauss(z, 1, 0.1)"}],
        }
        self.assertFalse(operas_expression_functions_required(calculator_payload))

        with_qualified_dump = {
            "calculator_modules": [
                {
                    "execution": {
                        "commands": ['python3 -c "numpy.savetxt(x)"'],
                        "input": [
                            {
                                "actions": [
                                    {
                                        "variables": [
                                            {
                                                "name": "xx",
                                                "expression": "user.external(x)",
                                            }
                                        ]
                                    }
                                ]
                            }
                        ],
                    }
                }
            ]
        }
        self.assertTrue(operas_expression_functions_required(with_qualified_dump))

    def test_startup_gate_reads_selection_and_target_expression_keys(self) -> None:
        self.assertTrue(
            operas_expression_functions_required(
                {"Sampling": {"selection": "helper.eggbox2d(x) > 0"}}
            )
        )
        self.assertTrue(
            operas_expression_functions_required(
                {
                    "Sampling": {
                        "Bounds": {
                            "target_expression": "ns.level(f)",
                        }
                    }
                }
            )
        )
        self.assertFalse(
            operas_expression_functions_required(
                {"Sampling": {"selection": "x + y < 1.5", "Bounds": {"radius": 0.3}}}
            )
        )

    def test_v1_operas_modules_shape_is_forwarded_without_new_functions_block(self) -> None:
        module = {
            "name": "EggBox",
            "operator": "helper.eggbox2d",
            "call_mode": "call",
            "required_modules": [],
            "input": [
                {"name": "x", "expression": "xx * Pi"},
                {"name": "y", "expression": "yy * Pi"},
            ],
            "output": [{"name": "z", "entry": "z"}],
        }
        worker = build_worker_config(
            {"Operas": {"Modules": [module]}},
            task_result_dir="/tmp/result",
        )

        self.assertEqual(worker["opera_modules"], [module])
        self.assertNotIn("operas_functions", worker)


@unittest.skipIf(jarvis_operas is None, "Jarvis-Operas is not installed")
class NamespaceConstantExpressionTests(unittest.TestCase):
    """D23: bare pdg.* constants fold into expressions (JO + HEP2 wiring)."""

    def setUp(self) -> None:
        reset_operas_discovery_for_tests()

    def test_snapshot_includes_constant_table(self) -> None:
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas build_constant_dicts not available")
        snap = discover_operas_expression_functions(required=True)
        expected = build_constant_dicts()
        self.assertIn("pdg.mZ", expected)
        self.assertEqual(snap.constants.get("pdg.mZ"), expected["pdg.mZ"])
        self.assertGreaterEqual(len(snap.constants), 1)

    def test_bare_pdg_constant_folds_and_coexists_with_observable(self) -> None:
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas build_constant_dicts not available")
        if "pdg.mZ" not in build_constant_dicts():
            self.skipTest("pdg.mZ not registered")

        context = build_operas_expression_context(required=True)
        compiled = context.compile("(mz - pdg.mZ)/pdg.mZ", symbols=("mz",))
        self.assertEqual(compiled.variable_names, ("mz",))
        # mz = 92 → (92 - 91.1876)/91.1876
        val = compiled.evaluate({"mz": 92.0})
        self.assertAlmostEqual(val, (92.0 - 91.1876) / 91.1876, places=12)

        # Bare observable mZ and qualified pdg.mZ are distinct free-symbol namespaces.
        co = context.compile("mZ + pdg.mZ", symbols=("mZ",))
        self.assertEqual(co.variable_names, ("mZ",))
        self.assertAlmostEqual(co.evaluate({"mZ": 1.0}), 1.0 + 91.1876, places=12)

    def test_constant_subexpression_algebraically_folds(self) -> None:
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas build_constant_dicts not available")
        constants = build_constant_dicts()
        if "pdg.mZ" not in constants or "pdg.mW" not in constants:
            self.skipTest("pdg masses not registered")

        context = build_operas_expression_context(required=True)
        compiled = context.compile(
            "sqrt(pdg.mZ**2 + pdg.mW**2) * x",
            symbols=("x",),
        )
        self.assertEqual(compiled.variable_names, ("x",))
        expected = (constants["pdg.mZ"] ** 2 + constants["pdg.mW"] ** 2) ** 0.5
        self.assertAlmostEqual(compiled.evaluate({"x": 1.0}), expected, places=9)
        # Algebraic fold: only free symbol is x; lambdify consts carry the float.
        consts = compiled._callable.__code__.co_consts
        self.assertTrue(
            any(isinstance(c, float) and abs(c - expected) < 1e-6 for c in consts),
            msg=f"expected folded float ~{expected} in co_consts={consts}",
        )

    def test_likelihood_and_selection_see_bare_constants(self) -> None:
        try:
            from jarvis_operas import build_constant_dicts
        except ImportError:
            self.skipTest("Jarvis-Operas build_constant_dicts not available")
        if "pdg.mZ" not in build_constant_dicts():
            self.skipTest("pdg.mZ not registered")

        likelihood = LogLikelihoodEvaluator(
            [{"name": "LogL", "expression": "-(mz - pdg.mZ)**2"}]
        )
        # Near the pole → near zero logL contribution.
        near = likelihood.evaluate({"mz": 91.1876})["LogL"]
        far = likelihood.evaluate({"mz": 100.0})["LogL"]
        self.assertAlmostEqual(near, 0.0, places=10)
        self.assertLess(far, near)

        self.assertTrue(
            evaluate_selection("mz > pdg.mZ", {"mz": 92.0})
        )
        self.assertFalse(
            evaluate_selection("mz > pdg.mZ", {"mz": 90.0})
        )


if __name__ == "__main__":
    unittest.main()
