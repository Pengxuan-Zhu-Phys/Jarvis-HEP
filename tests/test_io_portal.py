#!/usr/bin/env python3
"""Portal IO contracts: registry, format roundtrips and HEP-owned saved files."""

from __future__ import annotations

import asyncio
import json
import math
import os
import tempfile
import unittest
from pathlib import Path

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
from jarvishep2.file_operation_service import FileOperationService


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

    def test_evaluate_io_expression_keeps_numeric_string_coercion(self) -> None:
        value = evaluate_io_expression("x * Pi", {"x": "0.5"})
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
                module.acquire_pack_id("001")
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
                module.acquire_pack_id("001")
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
                module.acquire_pack_id("001")
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


MINIMAL_SLHA = """\
Block MASS
   25  125.0  # higgs
Block MINPAR
   1   0.1
"""


class PortalV2SlhaTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def tearDown(self) -> None:
        reset_io_registry_for_tests()

    def test_v2_registry_exposes_slha(self) -> None:
        formats = available_io_formats("input")
        self.assertIn("SLHA", formats)
        self.assertIn("JSON", formats)
        out_formats = available_io_formats("output")
        self.assertIn("SLHA", out_formats)
        self.assertIn("xSLHA", out_formats)

    def test_slha_write_and_read_end_to_end(self) -> None:
        registry = get_io_registry()
        adapter = registry.get("SLHA", "input")
        self.assertEqual(adapter.format_name, "SLHA")

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "input.slha")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(MINIMAL_SLHA)

            context = build_io_context(sample_info={"uuid": "slha-test"}, pack_id="001")

            async def _run():
                await write_io_input(
                    {
                        "name": "card",
                        "type": "SLHA",
                        "path": path,
                        "actions": [
                            {
                                "type": "SLHA",
                                "variables": [
                                    {
                                        "name": "mH",
                                        "block": "MASS",
                                        "entry": 25,
                                        "value": 130.0,
                                    }
                                ],
                            }
                        ],
                    },
                    {"mH": 130.0},
                    context=context,
                    path=path,
                )
                return await read_io_output(
                    {
                        "name": "card_out",
                        "type": "SLHA",
                        "path": path,
                        "variables": [
                            {"name": "mH", "block": "MASS", "entry": 25},
                        ],
                    },
                    context=context,
                    path=path,
                )

            observables = asyncio.run(_run())
            self.assertAlmostEqual(float(observables["mH"]), 130.0, places=6)


class IOSaveTrueTests(unittest.TestCase):
    def setUp(self) -> None:
        reset_io_registry_for_tests()

    def _assert_eggbox_save(
        self,
        *,
        file_ops: FileOperationService | None,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            proj = Path(tmp)
            runtime = proj / "calc" / "001"
            runtime.mkdir(parents=True)
            (runtime / "input.json").write_text('{"xx": 0, "yy": 0}\n', encoding="utf-8")
            (runtime / "eggbox.py").write_text(
                "import json\n"
                "d = json.load(open('input.json'))\n"
                "json.dump({'z': float(d.get('xx', 0)) + float(d.get('yy', 0))}, "
                "open('output.json', 'w'))\n",
                encoding="utf-8",
            )
            sample_root = proj / "SAMPLE"
            sample_root.mkdir()
            cfg = {
                "name": "EggBox",
                "path": str(proj / "calc" / "@PackID"),
                "execution": {
                    "path": str(proj / "calc" / "@PackID"),
                    "commands": [
                        {"cmd": "python3 eggbox.py", "cwd": str(runtime)}
                    ],
                    "input": [
                        {
                            "name": "inpjson",
                            "type": "JSON",
                            "path": str(runtime / "input.json"),
                            "save": True,
                            "actions": [
                                {
                                    "type": "Dump",
                                    "variables": [
                                        {"name": "xx", "expression": "x"},
                                        {"name": "yy", "expression": "y"},
                                    ],
                                }
                            ],
                        }
                    ],
                    "output": [
                        {
                            "name": "oupjson",
                            "type": "JSON",
                            "path": str(runtime / "output.json"),
                            "save": True,
                            "variables": [{"name": "z"}],
                        }
                    ],
                },
            }
            mod = CalculatorModule("EggBox", cfg)
            mod.attach_scheduler(AsyncSubprocessScheduler())
            if file_ops is not None:
                mod.attach_file_ops(file_ops)
            mod.acquire_pack_id("001")
            sample = Sample(uuid="save-uuid")
            sample.info = {
                "uuid": "save-uuid",
                "params": {"x": 1.5, "y": 2.5},
                "observables": {"x": 1.5, "y": 2.5, "uuid": "save-uuid"},
                "_bucket_parent": str(sample_root),
                "sample_artifacts": "always",
                "workflow_has_calculator": True,
            }
            sample.materialize()
            result = mod.execute(sample.info)
            self.assertAlmostEqual(float(result["z"]), 4.0)
            self.assertIn("SAMPLE", str(result["inpjson"]))
            self.assertIn("SAMPLE", str(result["oupjson"]))
            self.assertTrue(os.path.isfile(result["inpjson"]))
            self.assertTrue(os.path.isfile(result["oupjson"]))
            self.assertTrue(result["inpjson"].endswith("input.json@EggBox"))
            self.assertTrue(result["oupjson"].endswith("output.json@EggBox"))

    def test_json_save_true_copies_into_sample_dir(self) -> None:
        """Fallback path: Calculator apply_hep_io_save without a service client."""
        self._assert_eggbox_save(file_ops=None)

    def test_json_save_true_via_file_operation_service(self) -> None:
        service = FileOperationService.start(mode="inline")
        try:
            self._assert_eggbox_save(file_ops=service)
        finally:
            service.shutdown()


if __name__ == "__main__":
    unittest.main()
