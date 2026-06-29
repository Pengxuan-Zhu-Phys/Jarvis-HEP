#!/usr/bin/env python3
"""Workflow execution-plan unit tests (WP-D2.2 extras)."""

from __future__ import annotations

import pickle
import unittest

from jarvishep2.mp_context import get_spawn_context
from jarvishep2.runtime_config import workflow_has_calculator, workflow_references_sdir
from jarvishep2.workflow import (
    build_execution_plan,
    concurrency_groups,
    execution_plan_template,
    resolve_module_layers,
)

from test_layer_concurrency import SLOW_A_MODULE, SLOW_B_MODULE
from test_worker_calculator import EGGBOX_CALC_MODULE


class WorkflowExecutionPlanTests(unittest.TestCase):
    def test_plan_determinism(self) -> None:
        modules = [SLOW_A_MODULE, SLOW_B_MODULE]
        first = build_execution_plan(calculator_modules=modules, include_likelihood=False)
        second = build_execution_plan(calculator_modules=modules, include_likelihood=False)
        self.assertEqual(
            [(step.type, step.name, step.layer) for step in first],
            [(step.type, step.name, step.layer) for step in second],
        )

    def test_layering_correctness_fixture(self) -> None:
        dependent_b = dict(SLOW_B_MODULE)
        dependent_b["required_modules"] = ["SlowA"]
        layers = resolve_module_layers([SLOW_A_MODULE, dependent_b])
        plan = build_execution_plan(
            calculator_modules=[SLOW_A_MODULE, dependent_b],
            include_likelihood=False,
        )
        self.assertEqual(layers["SlowA"], 0)
        self.assertEqual(layers["SlowB"], 1)
        self.assertEqual(concurrency_groups(plan), [["SlowA"], ["SlowB"]])

    def test_concurrency_groups_match_layer_membership(self) -> None:
        modules = [SLOW_A_MODULE, SLOW_B_MODULE]
        plan = build_execution_plan(calculator_modules=modules, include_likelihood=False)
        layer_map = resolve_module_layers(modules)
        grouped = concurrency_groups(plan)
        layer_buckets: dict[int, list[str]] = {}
        for name, layer in layer_map.items():
            layer_buckets.setdefault(layer, []).append(name)
        expected = [sorted(layer_buckets[layer]) for layer in sorted(layer_buckets)]
        self.assertEqual([sorted(group) for group in grouped], expected)

    def test_workflow_flags_on_calculator_vs_opera_configs(self) -> None:
        calc_config = {
            "Calculators": {
                "Modules": [
                    {
                        "name": "EggBox",
                        "execution": {"commands": [{"cmd": "run", "cwd": "@Sdir"}]},
                    }
                ]
            }
        }
        opera_config = {
            "Operas": {
                "Modules": [{"name": "Trivial", "function": "lambda x: x"}],
            }
        }
        self.assertTrue(workflow_has_calculator(calc_config))
        self.assertFalse(workflow_has_calculator(opera_config))
        self.assertTrue(workflow_references_sdir(calc_config))
        self.assertFalse(workflow_references_sdir(opera_config))

    def test_execution_plan_template_pickles_under_spawn(self) -> None:
        template = execution_plan_template(
            calculator_modules=[EGGBOX_CALC_MODULE],
            include_likelihood=True,
        )
        blob = pickle.dumps(template, protocol=pickle.HIGHEST_PROTOCOL)
        restored = pickle.loads(blob)
        self.assertEqual(restored, template)
        self.assertEqual(get_spawn_context().get_start_method(), "spawn")

    @unittest.skip("flowchart export not implemented in jarvishep2 yet (V1 parity deferred)")
    def test_flowchart_parity_golden(self) -> None:
        raise NotImplementedError


if __name__ == "__main__":
    unittest.main()