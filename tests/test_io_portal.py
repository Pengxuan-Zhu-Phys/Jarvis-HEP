#!/usr/bin/env python3
"""Portal IO bridge tests — registry, unknown types, JSON/CSV via HEP facade."""

from __future__ import annotations

import asyncio
import json
import math
import os
import tempfile
import unittest

from jarvishep2.io_portal import (
    UnsupportedIOTypeError,
    available_io_formats,
    build_io_context,
    evaluate_io_expression,
    get_io_registry,
    read_io_output,
    read_io_output_sync,
    reset_io_registry_for_tests,
    write_io_input,
    write_io_input_sync,
)
from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig
from jarvishep2.command_parser import CommandParser
from jarvishep2.sample import Sample


def _attach_calc_runtime(module: CalculatorModule, project_root: str) -> AsyncSubprocessScheduler:
    scheduler = AsyncSubprocessScheduler(
        SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
    )
    scheduler.start()
    module.attach_scheduler(scheduler)
    module.attach_command_parser(CommandParser(project_root=project_root, scan_name="t"))
    return scheduler


class PortalRegistryTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_registry_exposes_json_and_csv(self) -> None:
        formats = {name.casefold() for name in available_io_formats()}
        self.assertIn("json", formats)
        self.assertIn("csv", formats)

    def test_unknown_type_raises_with_registered_list(self) -> None:
        context = build_io_context(runtime_values={"x": 1.0})
        with self.assertRaises(UnsupportedIOTypeError) as raised:
            write_io_input_sync(
                {"name": "bad", "path": "x.slha", "type": "NotARealFormat", "actions": []},
                {"x": 1.0},
                context=context,
                path="/tmp/x.slha",
            )
        message = str(raised.exception)
        self.assertIn("NotARealFormat", message)
        self.assertIn("registered:", message)
        self.assertIn("JSON", message)

    def test_empty_type_raises(self) -> None:
        context = build_io_context()
        with self.assertRaises(ValueError):
            write_io_input_sync(
                {"name": "bad", "path": "x.json", "type": "  ", "actions": []},
                {},
                context=context,
                path="/tmp/x.json",
            )


class ExpressionAndJsonBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_evaluate_io_expression_pi(self) -> None:
        value = evaluate_io_expression("x * Pi", {"x": 0.5})
        self.assertAlmostEqual(float(value), 0.5 * math.pi, places=12)

    def test_json_write_read_roundtrip_sync(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.json")
            out_path = os.path.join(tmpdir, "output.json")
            context = build_io_context(
                sample_info={"uuid": "s1", "save_dir": tmpdir},
                module="Unit",
                runtime_values={"x": 0.25, "y": 0.5},
            )
            write_obs = write_io_input_sync(
                {
                    "name": "inp",
                    "path": in_path,
                    "type": "JSON",
                    "actions": [
                        {
                            "type": "Dump",
                            "variables": [
                                {"name": "xx", "expression": "x * Pi"},
                                {"name": "yy", "expression": "y * Pi"},
                            ],
                        }
                    ],
                },
                {"x": 0.25, "y": 0.5},
                context=context,
                path=in_path,
            )
            self.assertAlmostEqual(float(write_obs["xx"]), 0.25 * math.pi, places=12)
            payload = json.loads(open(in_path, encoding="utf-8").read())
            self.assertAlmostEqual(float(payload["xx"]), 0.25 * math.pi, places=12)

            with open(out_path, "w", encoding="utf-8") as handle:
                json.dump({"z": 1.5, "fit": {"loglike": -2.0}}, handle)

            read_obs = read_io_output_sync(
                {
                    "name": "out",
                    "path": out_path,
                    "type": "JSON",
                    "variables": [
                        {"name": "z", "entry": "z"},
                        {"name": "loglike", "entry": "fit.loglike"},
                    ],
                },
                context=context,
                path=out_path,
            )
            self.assertEqual(read_obs["z"], 1.5)
            self.assertEqual(read_obs["loglike"], -2.0)

    def test_json_async_path(self) -> None:
        async def _run() -> None:
            with tempfile.TemporaryDirectory() as tmpdir:
                path = os.path.join(tmpdir, "in.json")
                context = build_io_context(runtime_values={"x": 2.0})
                obs = await write_io_input(
                    {
                        "name": "inp",
                        "path": path,
                        "type": "json",
                        "actions": [{"type": "Dump", "variables": [{"name": "x"}]}],
                    },
                    {"x": 2.0},
                    context=context,
                    path=path,
                )
                self.assertEqual(json.loads(open(path, encoding="utf-8").read())["x"], 2.0)
                self.assertIsInstance(obs, dict)

        asyncio.run(_run())


class CsvBridgeTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_csv_write_and_read_via_bridge(self) -> None:
        formats = {name.casefold() for name in available_io_formats()}
        if "csv" not in formats:
            self.skipTest("CSV adapter not exposed by installed Portal")

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.csv")
            out_path = os.path.join(tmpdir, "output.csv")
            context = build_io_context(runtime_values={"x": 1.5, "y": 2.5})

            write_io_input_sync(
                {
                    "name": "params",
                    "path": in_path,
                    "type": "CSV",
                    "header": True,
                    "actions": [
                        {
                            "type": "Dump",
                            "variables": [
                                {"name": "var_mass", "expression": "x", "column": "mass"},
                                {"name": "var_coupling", "expression": "y", "column": "coupling"},
                            ],
                        }
                    ],
                },
                {"x": 1.5, "y": 2.5},
                context=context,
                path=in_path,
            )
            text = open(in_path, encoding="utf-8").read()
            self.assertIn("mass", text)
            self.assertIn("1.5", text)

            with open(out_path, "w", encoding="utf-8") as handle:
                handle.write("mass,fit_loglike\n")
                handle.write("10.0,-3.25\n")

            obs = read_io_output_sync(
                {
                    "name": "observables",
                    "path": out_path,
                    "type": "CSV",
                    "header": True,
                    "variables": [
                        {"name": "best_mass", "column": "mass", "row": 0},
                        {"name": "loglike", "column": "fit_loglike", "row": 0},
                    ],
                },
                context=context,
                path=out_path,
            )
            self.assertEqual(float(obs["best_mass"]), 10.0)
            self.assertEqual(float(obs["loglike"]), -3.25)


class CalculatorPortalIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_load_input_does_not_clobber_params_with_dump_observables(self) -> None:
        """Dump expression observables must not overwrite physical params of the same name."""
        with tempfile.TemporaryDirectory() as tmpdir:
            module = CalculatorModule(
                "ClobberGuard",
                {
                    "name": "ClobberGuard",
                    "execution": {
                        "commands": [],
                        "input": [
                            {
                                "name": "inp",
                                "path": "@Sdir/input.json",
                                "type": "JSON",
                                "actions": [
                                    {
                                        "type": "Dump",
                                        "variables": [
                                            {"name": "x", "expression": "x * Pi"},
                                            {"name": "xx", "expression": "x * Pi"},
                                        ],
                                    }
                                ],
                            }
                        ],
                        "output": [],
                    },
                },
            )
            sample = Sample.from_params({"x": 0.5, "y": 0.1, "uuid": "clobber"})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": True,
                }
            )
            sample.materialize()
            scheduler = _attach_calc_runtime(module, tmpdir)
            try:
                module.preload_templates()
                module.acquire_pack_id("pack-clobber")
                result = module.execute(sample.info)
            finally:
                scheduler.shutdown(wait=True)

            self.assertAlmostEqual(float(result["x"]), 0.5, places=12)
            self.assertAlmostEqual(float(result["xx"]), 0.5 * math.pi, places=12)
            payload = json.loads(open(os.path.join(sample.save_dir, "input.json"), encoding="utf-8").read())
            self.assertAlmostEqual(float(payload["x"]), 0.5 * math.pi, places=12)

    def test_calculator_rejects_unknown_io_type(self) -> None:
        module = CalculatorModule(
            "BadIO",
            {
                "name": "BadIO",
                "execution": {
                    "commands": [],
                    "input": [
                        {
                            "name": "inp",
                            "path": "@Sdir/x.notreal",
                            "type": "NotARealFormat",
                            "actions": [],
                        }
                    ],
                    "output": [],
                },
            },
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0, "uuid": "bad-io"})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": True,
                }
            )
            sample.materialize()
            scheduler = _attach_calc_runtime(module, tmpdir)
            try:
                module.preload_templates()
                module.acquire_pack_id("pack-1")
                with self.assertRaises(UnsupportedIOTypeError):
                    module.execute(sample.info)
            finally:
                scheduler.shutdown(wait=True)

    def test_calculator_json_still_works_through_portal(self) -> None:
        eggbox_dir = os.path.join(
            os.path.dirname(__file__),
            "parity_project",
            "calculators",
            "runtime",
            "program",
            "eggbox",
            "001",
        )
        script = os.path.join(eggbox_dir, "eggbox.py")
        module = CalculatorModule(
            "EggBox",
            {
                "name": "EggBox",
                "path": eggbox_dir,
                "execution": {
                    "path": eggbox_dir,
                    "commands": [{"cmd": f"python3 {script}", "cwd": "@Sdir"}],
                    "input": [
                        {
                            "name": "inpjson",
                            "path": "@Sdir/input.json",
                            "type": "JSON",
                            "actions": [
                                {
                                    "type": "Dump",
                                    "variables": [
                                        {"name": "xx", "expression": "x * Pi"},
                                        {"name": "yy", "expression": "y * Pi"},
                                    ],
                                }
                            ],
                        }
                    ],
                    "output": [
                        {
                            "name": "oupjson",
                            "path": "@Sdir/output.json",
                            "type": "JSON",
                            "variables": [{"name": "z", "entry": "z"}],
                        }
                    ],
                },
            },
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 0.25, "y": 0.25, "uuid": "portal-eggbox"})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": True,
                }
            )
            sample.materialize()
            scheduler = _attach_calc_runtime(module, tmpdir)
            try:
                module.preload_templates()
                module.acquire_pack_id("pack-egg")
                result = module.execute(sample.info)
            finally:
                scheduler.shutdown(wait=True)
            from numpy import cos, sin

            expected_z = float((sin(0.25 * math.pi) * cos(0.25 * math.pi) + 2) ** 5)
            self.assertAlmostEqual(float(result["z"]), expected_z, places=6)
            self.assertIn("x", result)
            self.assertTrue(os.path.exists(os.path.join(sample.save_dir, "input.json")))


class RegistryContractTests(unittest.TestCase):
    def test_get_io_registry_is_cached(self) -> None:
        reset_io_registry_for_tests()
        first = get_io_registry()
        second = get_io_registry()
        self.assertIs(first, second)


if __name__ == "__main__":
    unittest.main()
