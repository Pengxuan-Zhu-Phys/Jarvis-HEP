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
    operas_expression_functions_required,
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
                {"Sampling": {"selection": "x + y < 1.5", "Bounds": {"Radius": 0.3}}}
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
            {"Operas": {"make_paraller": 16, "Modules": [module]}},
            task_result_dir="/tmp/result",
        )

        self.assertEqual(worker["opera_modules"], [module])
        self.assertNotIn("operas_functions", worker)


if __name__ == "__main__":
    unittest.main()
