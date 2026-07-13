#!/usr/bin/env python3
"""Operas operator resolution — Jarvis-Operas registry + importlib fallback."""

from __future__ import annotations

import math
import os
import tempfile
import unittest
from unittest import mock

import jarvishep2.expression as expression_backend
from jarvishep2.operas import OperasModule, resolve_operator
from jarvishep2.sample import Sample
from jarvishep2.testing.eggbox import eggbox2d_numpy


def _operas_available() -> bool:
    try:
        from jarvis_operas import get_global_operas_registry

        registry = get_global_operas_registry()
        registry.resolve_name("helper.eggbox2d")
        return True
    except Exception:
        return False


class ResolveOperatorTests(unittest.TestCase):
    def test_importlib_dotted_path(self) -> None:
        func = resolve_operator("jarvishep2.testing.eggbox.eggbox2d_numpy")
        result = func(x=0.25, y=0.25)
        self.assertIn("z", result)
        expected = eggbox2d_numpy(0.25, 0.25)["z"]
        self.assertAlmostEqual(float(result["z"]), float(expected), places=12)

    def test_unknown_operator_raises(self) -> None:
        with self.assertRaises(ValueError) as raised:
            resolve_operator("not.a.real.operator.anywhere")
        self.assertIn("Cannot resolve", str(raised.exception))

    @unittest.skipUnless(_operas_available(), "Jarvis-Operas not installed")
    def test_jarvis_operas_registry_name(self) -> None:
        func = resolve_operator("helper.eggbox2d")
        result = func(x=0.25, y=0.25)
        self.assertIn("z", result)
        # Match local reference implementation on radians-style inputs.
        expected = float((math.sin(0.25) * math.cos(0.25) + 2.0) ** 5)
        self.assertAlmostEqual(float(result["z"]), expected, places=6)


class OperasModuleBridgeTests(unittest.TestCase):
    def test_module_importlib_operator_still_works(self) -> None:
        module = OperasModule(
            "TrivialEggbox",
            {
                "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                "input": [
                    {"name": "x", "expression": "x"},
                    {"name": "y", "expression": "y"},
                ],
                "output": [{"name": "z", "entry": "z"}],
            },
        )
        module.preload()
        out = module.execute({"x": 0.3, "y": 0.4}, {"uuid": "u1"})
        self.assertIn("z", out)
        expected = eggbox2d_numpy(0.3, 0.4)["z"]
        self.assertAlmostEqual(float(out["z"]), float(expected), places=12)

    @unittest.skipUnless(_operas_available(), "Jarvis-Operas not installed")
    def test_module_registry_operator_end_to_end(self) -> None:
        module = OperasModule(
            "RegistryEggbox",
            {
                "operator": "helper.eggbox2d",
                "call_mode": "call",
                "input": [
                    {"name": "x", "expression": "x"},
                    {"name": "y", "expression": "y"},
                ],
                "output": [{"name": "z", "entry": "z"}],
            },
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 0.2, "y": 0.3, "uuid": "reg-egg"})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                }
            )
            sample.materialize()
            module.preload()
            out = module.execute(sample.observables, sample.info)
            expected = float((math.sin(0.2) * math.cos(0.3) + 2.0) ** 5)
            self.assertAlmostEqual(float(out["z"]), expected, places=6)

    @unittest.skipUnless(_operas_available(), "Jarvis-Operas not installed")
    def test_v1_eggbox_pi_input_expressions(self) -> None:
        module = OperasModule(
            "EggBox",
            {
                "operator": "helper.eggbox2d",
                "call_mode": "call",
                "input": [
                    {"name": "x", "expression": "xx * Pi"},
                    {"name": "y", "expression": "yy * Pi"},
                ],
                "output": [{"name": "z", "entry": "z"}],
            },
        )
        module.preload()
        out = module.execute({"xx": 0.5, "yy": 0.0}, {"uuid": "v1-eggbox"})
        expected = float((math.sin(0.5 * math.pi) * math.cos(0.0) + 2.0) ** 5)
        self.assertAlmostEqual(float(out["z"]), expected, places=6)

    def test_input_expressions_compile_once_and_reuse(self) -> None:
        module = OperasModule(
            "EggBox",
            {
                "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                "input": [
                    {"name": "x", "expression": "xx * Pi"},
                    {"name": "y", "expression": "yy * Pi"},
                ],
                "output": [{"name": "z", "entry": "z"}],
            },
        )
        original_lambdify = expression_backend.lambdify
        with mock.patch(
            "jarvishep2.expression.lambdify",
            wraps=original_lambdify,
        ) as compile_spy:
            module.preload()
            self.assertEqual(compile_spy.call_count, 2)
            first = module.execute({"xx": 0.25, "yy": 0.5}, {"uuid": "first"})
            second = module.execute({"xx": 0.5, "yy": 0.25}, {"uuid": "second"})

        self.assertEqual(compile_spy.call_count, 2)
        self.assertNotEqual(float(first["z"]), float(second["z"]))

    def test_sample_logger_records_full_opera_call(self) -> None:
        module = OperasModule(
            "EggBox",
            {
                "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                "input": [
                    {"name": "x", "expression": "xx * Pi"},
                    {"name": "y", "expression": "yy * Pi"},
                ],
                "output": [{"name": "z", "entry": "z"}],
            },
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"xx": 0.25, "yy": 0.5})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                }
            )
            sample.materialize()
            module.preload()
            output = module.execute(
                sample.observables,
                sample.info,
                sample_logger=sample.child_logger(module="Sample@logger-test (Operas:EggBox)"),
            )
            sample.merge_observables(output)
            sample.close()
            with open(sample.info["run_log"], encoding="utf-8") as handle:
                log = handle.read()

        self.assertIn("Evaluating   x:", log)
        self.assertIn("Operas input dispatch:", log)
        self.assertIn("operator \t-> jarvishep2.testing.eggbox.eggbox2d_numpy", log)
        self.assertIn("EggBox input loaded: x ->", log)
        self.assertIn("Operas output:", log)
        self.assertIn("mapped \t-> [z :", log)


if __name__ == "__main__":
    unittest.main()
