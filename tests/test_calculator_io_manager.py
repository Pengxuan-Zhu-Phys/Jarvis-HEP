#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import json
import math
import os
import sys
import tempfile
import threading
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.IOs import registry as io_registry  # noqa: E402
from jarvishep.IOs.IOs import IOfile  # noqa: E402
from jarvishep.IOs.Output import FileOutput  # noqa: E402
from jarvishep.IOs.portal import PortalInputFile, PortalOutputFile  # noqa: E402
from jarvishep.Module.calculator import CalculatorModule  # noqa: E402
import jarvis_portal  # noqa: E402


class _NoopLogger:
    def __init__(self):
        self._lock = threading.Lock()
        self.messages = []

    def _record(self, level, message):
        with self._lock:
            self.messages.append((level, str(message)))

    def debug(self, message, *args, **kwargs):
        self._record("DEBUG", message)

    def info(self, message, *args, **kwargs):
        self._record("INFO", message)

    def warning(self, message, *args, **kwargs):
        self._record("WARNING", message)

    def error(self, message, *args, **kwargs):
        self._record("ERROR", message)


class _CountingIOManager:
    def __init__(self):
        self.run_blocking_calls = 0

    async def run_blocking(self, fn, *args, **kwargs):
        self.run_blocking_calls += 1
        return fn(*args, **kwargs)

    async def read_text(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def write_text(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def make_dirs(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def exists(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")


class _FakeTableAdapter:
    format_name = "FAKE_TABLE"
    direction = "both"

    def __init__(self):
        self.write_context = None
        self.read_context = None

    async def write_input(self, context, spec, data):
        self.write_context = context
        rows = []
        observables = {}
        for action in spec.get("actions", []):
            if action.get("type") != "Dump":
                continue
            for variable in action.get("variables", []):
                name = variable["name"]
                if "expression" in variable:
                    value = context.evaluate_expression(variable["expression"], data)
                    observables[name] = value
                else:
                    value = data.get(name)
                rows.append(f"{name}={value}")
        path = context.path(spec["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(rows), encoding="utf-8")
        return observables

    async def read_output(self, context, spec):
        self.read_context = context
        path = context.path(spec["path"])
        payload = {}
        if path.exists():
            for line in path.read_text(encoding="utf-8").splitlines():
                if "=" not in line:
                    continue
                key, value = line.split("=", 1)
                payload[key] = float(value)
        return {variable["name"]: payload.get(variable["name"]) for variable in spec.get("variables", [])}


class CalculatorIOManagerTests(unittest.TestCase):
    def test_registry_exposes_builtin_calculator_io_formats(self):
        self.assertIn("JSON", io_registry.available_formats("input"))
        self.assertIn("CSV", io_registry.available_formats("input"))
        self.assertIn("TSV", io_registry.available_formats("input"))
        self.assertIn("DAT", io_registry.available_formats("input"))
        self.assertIn("SLHA", io_registry.available_formats("input"))
        self.assertIn("JSON", io_registry.available_formats("output"))
        self.assertIn("CSV", io_registry.available_formats("output"))
        self.assertIn("TSV", io_registry.available_formats("output"))
        self.assertIn("DAT", io_registry.available_formats("output"))
        self.assertIn("SLHA", io_registry.available_formats("output"))
        self.assertIn("xSLHA", io_registry.available_formats("output"))
        self.assertIn("File", io_registry.available_formats("output"))

    def test_iofile_dispatch_uses_case_insensitive_registry_lookup(self):
        logger = _NoopLogger()

        input_handler = IOfile.create(
            "input_cfg",
            "/tmp/input.json",
            "json",
            [{"type": "Dump", "variables": [{"name": "alpha"}]}],
            False,
            logger,
            None,
            "S-1",
            "/tmp/sample",
            "DemoModule",
            {},
        )
        output_handler = IOfile.load(
            "output_cfg",
            "/tmp/output.json",
            "json",
            [{"name": "alpha"}],
            False,
            logger,
            None,
            "S-1",
            "/tmp/sample",
            "DemoModule",
            {},
        )

        self.assertIsInstance(input_handler, PortalInputFile)
        self.assertIsInstance(output_handler, PortalOutputFile)

    def test_iofile_dispatch_exposes_generic_portal_handlers_for_table_formats(self):
        logger = _NoopLogger()

        input_handler = IOfile.create(
            "input_cfg",
            "/tmp/input.csv",
            "csv",
            [{"type": "Dump", "variables": [{"name": "alpha", "column": "alpha"}]}],
            False,
            logger,
            None,
            "S-1",
            "/tmp/sample",
            "DemoModule",
            {},
            header=True,
        )
        output_handler = IOfile.load(
            "output_cfg",
            "/tmp/output.csv",
            "csv",
            [{"name": "alpha"}],
            False,
            logger,
            None,
            "S-1",
            "/tmp/sample",
            "DemoModule",
            {},
            header=True,
        )

        self.assertIsInstance(input_handler, PortalInputFile)
        self.assertIsInstance(output_handler, PortalOutputFile)

    def test_table_formats_dispatch_through_portal_registry_only(self):
        logger = _NoopLogger()

        for format_name in ("CSV", "TSV", "DAT"):
            with self.subTest(format_name=format_name):
                input_handler = IOfile.create(
                    f"{format_name.lower()}_input",
                    f"/tmp/input.{format_name.lower()}",
                    format_name,
                    [{"type": "Dump", "variables": [{"name": "alpha", "column": "alpha"}]}],
                    False,
                    logger,
                    None,
                    "S-1",
                    "/tmp/sample",
                    "DemoModule",
                    {},
                    header=True,
                )
                output_handler = IOfile.load(
                    f"{format_name.lower()}_output",
                    f"/tmp/output.{format_name.lower()}",
                    format_name,
                    [{"name": "alpha"}],
                    False,
                    logger,
                    None,
                    "S-1",
                    "/tmp/sample",
                    "DemoModule",
                    {},
                    header=True,
                )

                self.assertIsInstance(input_handler, PortalInputFile)
                self.assertIsInstance(output_handler, PortalOutputFile)
                self.assertEqual(input_handler.portal_adapter.format_name, format_name)
                self.assertEqual(output_handler.portal_adapter.format_name, format_name)

    def test_missing_calculator_io_format_reports_available_formats(self):
        logger = _NoopLogger()

        with self.assertRaisesRegex(
            ValueError,
            r"Unsupported IO format 'ROOT' for output.*Available output formats:.*JSON",
        ):
            IOfile.load(
                "root_output",
                "/tmp/output.root",
                "ROOT",
                [],
                False,
                logger,
                None,
                "S-1",
                "/tmp/sample",
                "DemoModule",
                {},
            )

    def test_json_input_file_runs_whole_write_in_io_manager(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "input.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"alpha": 1}, handle)

            io_manager = _CountingIOManager()
            logger = _NoopLogger()
            sample_dir = os.path.join(tmpdir, "sample")
            handler = IOfile.create(
                "input_cfg",
                path,
                "JSON",
                [{"type": "Dump", "variables": [{"name": "alpha"}]}],
                False,
                logger,
                None,
                "S-1",
                sample_dir,
                "DemoModule",
                {},
                io_manager,
            )

            observables = asyncio.run(handler.write({"alpha": 7}))

            self.assertEqual(io_manager.run_blocking_calls, 1)
            self.assertEqual(observables, {})
            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["alpha"], 7)

    def test_json_input_file_save_copies_written_payload_to_sample_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime_dir = os.path.join(tmpdir, "runtime")
            sample_dir = os.path.join(tmpdir, "outputs", "SAMPLE", "S-1")
            os.makedirs(runtime_dir, exist_ok=True)
            path = os.path.join(runtime_dir, "input.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"alpha": 1}, handle)

            logger = _NoopLogger()
            handler = IOfile.create(
                "input_cfg",
                path,
                "JSON",
                [{"type": "Dump", "variables": [{"name": "alpha"}]}],
                True,
                logger,
                None,
                "S-1",
                sample_dir,
                "DemoModule",
                {},
            )

            observables = asyncio.run(handler.write({"alpha": 7}))

            saved_path = os.path.join(sample_dir, "input.json@DemoModule")
            self.assertEqual(observables["input_cfg"], os.path.realpath(path))
            self.assertTrue(os.path.exists(saved_path))
            with open(saved_path, "r", encoding="utf-8") as handle:
                saved_payload = json.load(handle)
            self.assertEqual(saved_payload["alpha"], 7)

    def test_json_input_file_expression_writes_nested_entry_and_returns_observable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "input.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"test": {"config": {"x": 0}}}, handle)

            logger = _NoopLogger()
            handler = IOfile.create(
                "input_cfg",
                path,
                "JSON",
                [
                    {
                        "type": "Dump",
                        "variables": [
                            {"name": "x"},
                            {
                                "name": "cx",
                                "expression": "(x + y) * Pi",
                                "entry": "test.config.x",
                            },
                        ],
                    }
                ],
                False,
                logger,
                None,
                "S-1",
                os.path.join(tmpdir, "sample"),
                "DemoModule",
                {},
            )

            observables = asyncio.run(handler.write({"x": 1.0, "y": 2.0}))

            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["x"], 1.0)
            self.assertAlmostEqual(payload["test"]["config"]["x"], 3.0 * math.pi)
            self.assertNotIn("x", observables)
            self.assertAlmostEqual(observables["cx"], 3.0 * math.pi)
            self.assertTrue(any("Evaluating: expression" in msg for _, msg in logger.messages))

    def test_json_output_file_runs_whole_read_in_io_manager(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "output.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"alpha": 9}, handle)

            io_manager = _CountingIOManager()
            logger = _NoopLogger()
            sample_dir = os.path.join(tmpdir, "sample")
            handler = IOfile.load(
                "output_cfg",
                path,
                "JSON",
                [{"name": "alpha"}],
                True,
                logger,
                None,
                "S-1",
                sample_dir,
                "DemoModule",
                {},
                io_manager,
            )

            observables = asyncio.run(handler.read())

            self.assertEqual(io_manager.run_blocking_calls, 1)
            self.assertEqual(observables["alpha"], 9)
            self.assertEqual(observables["output_cfg"], os.path.realpath(path))
            self.assertTrue(os.path.exists(os.path.join(sample_dir, "output.json@DemoModule")))

    def test_json_output_file_reads_nested_entry(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "output.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"z": 9, "fit": {"loglike": -3.5}}, handle)

            logger = _NoopLogger()
            handler = IOfile.load(
                "output_cfg",
                path,
                "JSON",
                [{"name": "z"}, {"name": "likelihood", "entry": "fit.loglike"}],
                False,
                logger,
                None,
                "S-1",
                os.path.join(tmpdir, "sample"),
                "DemoModule",
                {},
            )

            observables = asyncio.run(handler.read())

            self.assertEqual(observables["z"], 9)
            self.assertEqual(observables["likelihood"], -3.5)

    def test_csv_input_and_output_delegate_to_portal_adapter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.csv")
            output_path = os.path.join(tmpdir, "output.csv")
            with open(input_path, "w", encoding="utf-8") as handle:
                handle.write("mass,coupling,config\n0,0,0\n")
            with open(output_path, "w", encoding="utf-8") as handle:
                handle.write("chi2,mass,fit_loglike\n3.42,125.1,-1.25\n1.08,125.3,-0.87\n")

            logger = _NoopLogger()
            sample_dir = os.path.join(tmpdir, "sample")
            input_handler = IOfile.create(
                "params",
                input_path,
                "CSV",
                [
                    {
                        "type": "Dump",
                        "variables": [
                            {"name": "mass", "expression": "x * Pi", "column": "mass"},
                            {"name": "plain", "column": "coupling"},
                        ],
                    }
                ],
                False,
                logger,
                None,
                "S-1",
                sample_dir,
                "DemoModule",
                {},
                header=True,
            )
            observables = asyncio.run(input_handler.write({"x": 2.0, "plain": 7}))
            self.assertAlmostEqual(observables["mass"], 2.0 * math.pi)

            with open(input_path, "r", encoding="utf-8") as handle:
                self.assertIn(str(2.0 * math.pi), handle.read())

            output_handler = IOfile.load(
                "observables",
                output_path,
                "CSV",
                [
                    {"name": "chi2", "row": 0},
                    {"name": "best_mass", "column": "mass", "row": 0},
                    {"name": "loglike", "column": "fit_loglike"},
                ],
                False,
                logger,
                None,
                "S-1",
                sample_dir,
                "DemoModule",
                {},
                header=True,
            )
            result = asyncio.run(output_handler.read())

            self.assertEqual(result["chi2"], 3.42)
            self.assertEqual(result["best_mass"], 125.1)
            self.assertEqual(result["loglike"], [-1.25, -0.87])

    def test_fake_portal_adapter_runs_through_calculator_without_hep_format_hardcode(self):
        fake_adapter = _FakeTableAdapter()
        jarvis_portal.register("FAKE_TABLE", fake_adapter, "both", override=True)

        with tempfile.TemporaryDirectory() as tmpdir:
            module = CalculatorModule(
                "FakeCalc",
                {
                    "modes": False,
                    "required_modules": [],
                    "clone_shadow": False,
                    "installation": [],
                    "initialization": [],
                    "path": tmpdir,
                    "execution": {
                        "commands": [],
                        "input": [
                            {
                                "name": "params",
                                "path": os.path.join(tmpdir, "input.fake"),
                                "type": "FAKE_TABLE",
                                "save": False,
                                "actions": [
                                    {
                                        "type": "Dump",
                                        "variables": [
                                            {"name": "mass", "expression": "x * Pi"},
                                            {"name": "plain"},
                                        ],
                                    }
                                ],
                            }
                        ],
                        "output": [
                            {
                                "name": "observables",
                                "path": os.path.join(tmpdir, "output.fake"),
                                "type": "FAKE_TABLE",
                                "save": False,
                                "variables": [{"name": "chi2"}],
                            }
                        ],
                    },
                },
            )
            module.logger = _NoopLogger()
            module.PackID = "P-1"
            module.sample_info = {"uuid": "S-1", "save_dir": os.path.join(tmpdir, "sample")}
            os.makedirs(module.sample_info["save_dir"], exist_ok=True)
            with open(os.path.join(tmpdir, "output.fake"), "w", encoding="utf-8") as handle:
                handle.write("chi2=3.5\n")

            input_observables = asyncio.run(module.load_input({"x": 2.0, "plain": 7.0}))
            output_observables = asyncio.run(module.read_output())

            self.assertIn("FAKE_TABLE", io_registry.available_formats("input"))
            self.assertAlmostEqual(input_observables["mass"], 2.0 * math.pi)
            self.assertEqual(output_observables["chi2"], 3.5)
            self.assertEqual(fake_adapter.write_context.sample_uuid, "S-1")
            self.assertEqual(fake_adapter.write_context.pack_id, "P-1")
            self.assertEqual(fake_adapter.write_context.sample_save_dir, module.sample_info["save_dir"])
            self.assertEqual(fake_adapter.write_context.module, "FakeCalc")
            self.assertEqual(fake_adapter.write_context.runtime_values["plain"], 7.0)
            self.assertIsNotNone(fake_adapter.write_context.resolve_path)
            self.assertIsNotNone(fake_adapter.write_context.evaluate_expression)
            self.assertIsNotNone(fake_adapter.read_context.resolve_path)

    def test_calculator_analyze_config_accepts_table_dump_actions(self):
        module = CalculatorModule(
            "CsvCalc",
            {
                "modes": False,
                "required_modules": [],
                "clone_shadow": False,
                "installation": [],
                "initialization": [],
                "path": "/tmp/csv-calc",
                "execution": {
                    "commands": [],
                    "input": [
                        {
                            "name": "params",
                            "path": "input.csv",
                            "type": "CSV",
                            "header": True,
                            "save": False,
                            "actions": [
                                {
                                    "type": "Dump",
                                    "variables": [
                                        {"name": "mass", "expression": "x * Pi", "column": "mass"},
                                        {"name": "plain", "column": "coupling"},
                                    ],
                                }
                            ],
                        }
                    ],
                    "output": [
                        {
                            "name": "observables",
                            "path": "output.csv",
                            "type": "CSV",
                            "header": True,
                            "save": False,
                            "variables": [{"name": "chi2"}],
                        }
                    ],
                },
            },
        )

        self.assertIn("x", module.inputs)
        self.assertIn("plain", module.inputs)
        self.assertIn("chi2", module.outputs)

    def test_file_output_still_uses_internal_handler(self):
        logger = _NoopLogger()
        handler = IOfile.load(
            "raw_file",
            "/tmp/output.txt",
            "File",
            [],
            True,
            logger,
            None,
            "S-1",
            "/tmp/sample",
            "DemoModule",
            {},
        )

        self.assertIsInstance(handler, FileOutput)


if __name__ == "__main__":
    unittest.main()
