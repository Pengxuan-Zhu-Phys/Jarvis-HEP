#!/usr/bin/env python3
"""WP-D1.2 calculator-in-Worker tests — fakeredis + spawn parity path."""

from __future__ import annotations

import csv
import glob
import json
import os
import tempfile
import threading
import time
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.sample import ExecutionStep, Sample
from jarvishep2.workflow import build_execution_plan


TESTS_ROOT = os.path.dirname(__file__)
PARITY_PROJECT = os.path.join(TESTS_ROOT, "parity_project")
FIXTURES = os.path.join(TESTS_ROOT, "fixtures", "parity_m1")
ECHO_CALC_DIR = os.path.join(
    PARITY_PROJECT,
    "calculators",
    "runtime",
    "program",
    "echo_calc",
    "001",
)

ECHO_CALC_MODULE = {
    "name": "EchoCalc",
    "required_modules": ["Parameters"],
    "clone_shadow": False,
    "path": ECHO_CALC_DIR,
    "installation": [],
    "initialization": [],
    "execution": {
        "path": ECHO_CALC_DIR,
        "commands": [
            {
                "cmd": (
                    "python3 -c 'import json; p=json.load(open(\"@Sdir/in.json\")); "
                    "json.dump({\"calc_z\": float(p[\"x\"])*2.0}, open(\"@Sdir/out.json\",\"w\"))'"
                ),
                "cwd": ECHO_CALC_DIR,
            },
        ],
        "input": [
            {
                "name": "params",
                "path": "@Sdir/in.json",
                "type": "JSON",
                "save": False,
                "actions": [
                    {
                        "type": "Dump",
                        "variables": [{"name": "x"}],
                    }
                ],
            }
        ],
        "output": [
            {
                "name": "observables",
                "path": "@Sdir/out.json",
                "type": "JSON",
                "save": False,
                "variables": [{"name": "calc_z", "entry": "calc_z"}],
            }
        ],
    },
}

LIKELIHOOD_EXPRESSIONS = [{"name": "LogL", "expression": "LogGauss(calc_z, calc_z, 1.0)"}]


def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict[str, Any]]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0}


def _worker_config(tmpdir: str) -> dict[str, Any]:
    return {
        "sample_config": {
            "task_result_dir": tmpdir,
            "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
            "sample_artifacts": "always",
            "workflow_has_calculator": True,
            "workflow_references_sdir": True,
        },
        "mapper": {"type": "identity", "keys": ["x"]},
        "calculator_modules": [ECHO_CALC_MODULE],
        "likelihood_expressions": LIKELIHOOD_EXPRESSIONS,
        "pull_timeout": 1,
    }


def _load_csv_points() -> list[dict[str, str]]:
    csv_path = os.path.join(PARITY_PROJECT, "data", "check_modules_points.csv")
    with open(csv_path, "r", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _normalize_database_records(records: list[dict[str, Any]]) -> list[dict[str, float]]:
    normalized = []
    for row in records:
        normalized.append(
            {
                "x": float(row["x"]),
                "calc_z": float(row["calc_z"]),
                "LogL": float(row["LogL"]),
            }
        )
    return sorted(normalized, key=lambda item: item["x"])


def _sample_tree_file_sets(sample_root: str) -> list[list[str]]:
    manifests: list[list[str]] = []
    if not os.path.isdir(sample_root):
        return manifests
    for child in sorted(os.listdir(sample_root)):
        child_path = os.path.join(sample_root, child)
        if not os.path.isdir(child_path):
            continue
        files = sorted(
            os.path.relpath(path, child_path)
            for path in glob.glob(os.path.join(child_path, "**", "*"), recursive=True)
            if os.path.isfile(path)
        )
        manifests.append(files)
    return sorted(manifests)


class WorkflowCalculatorTests(unittest.TestCase):
    def test_build_execution_plan_orders_calculator_before_likelihood(self) -> None:
        plan = build_execution_plan(
            calculator_modules=[ECHO_CALC_MODULE],
            include_likelihood=True,
        )
        self.assertEqual([step.type for step in plan], ["calculator", "likelihood"])
        self.assertEqual(plan[0].name, "EchoCalc")
        self.assertEqual(plan[0].layer, 0)
        # calculator-only plans reserve layer 1 for a future opera layer
        self.assertEqual(plan[1].layer, 2)


class CalculatorModuleUnitTests(unittest.TestCase):
    def test_preload_templates_is_idempotent(self) -> None:
        module = CalculatorModule("EchoCalc", ECHO_CALC_MODULE)
        module.preload_templates()
        module.preload_templates()
        self.assertTrue(module._templates_loaded)

    def test_execute_produces_expected_observables(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 0.25, "uuid": "unit-sample"})
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
            module = CalculatorModule("EchoCalc", ECHO_CALC_MODULE)
            module.preload_templates()
            module.acquire_pack_id("pack-test")
            result = module.execute(sample.info)
            self.assertAlmostEqual(float(result["calc_z"]), 0.5)
            self.assertTrue(os.path.exists(os.path.join(sample.save_dir, "in.json")))
            self.assertTrue(os.path.exists(os.path.join(sample.save_dir, "out.json")))


class WorkerCalculatorIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_single_worker_calculator_database_and_sample_parity(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(os.path.join(FIXTURES, "expected_calculator_records.json"), encoding="utf-8") as handle:
                    expected_records = json.load(handle)
                with open(os.path.join(FIXTURES, "expected_sample_files.json"), encoding="utf-8") as handle:
                    expected_files = json.load(handle)

                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 1,
                            "redis": redis_config,
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                core.init_factory(_worker_config(tmpdir))
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[ECHO_CALC_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                samples = []
                for row in _load_csv_points():
                    sample = sampler._build_sample(np.array([float(row["x"])], dtype=np.float64))
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)

                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=60.0)
                core.shutdown()

                from jarvishep2.database import SimpleHDF5Writer

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))

                sample_root = os.path.join(tmpdir, "SAMPLE")
                tree = _sample_tree_file_sets(sample_root)
                self.assertEqual(len(tree), 10)
                for files in tree:
                    self.assertEqual(files, expected_files)
        finally:
            server.shutdown()
            server.server_close()

    def test_worker_calculator_pack_id_traceability_in_result(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_worker_config(tmpdir))
                assert factory.redis is not None

                plan = build_execution_plan(
                    calculator_modules=[ECHO_CALC_MODULE],
                    include_likelihood=True,
                )
                sample = Sample(
                    uuid="pack-trace-sample",
                    u_coords=np.array([0.25], dtype=np.float64),
                    execution_plan=plan,
                )
                factory.redis.push_task(sample.to_task_dict())

                deadline = time.monotonic() + 20.0
                result = None
                while time.monotonic() < deadline:
                    result = factory.redis.pull_result(timeout=1)
                    if result is not None:
                        break
                factory.shutdown()

                self.assertIsNotNone(result)
                assert result is not None
                self.assertEqual(result["status"], "Completed")
                self.assertIn("pack_id", result)
                self.assertTrue(str(result["pack_id"]).strip())
                self.assertAlmostEqual(float(result["observables"]["calc_z"]), 0.5)
                sample_dir = os.path.join(tmpdir, "SAMPLE", "pack-trace-sample")
                self.assertTrue(os.path.isdir(sample_dir))
                self.assertTrue(os.path.exists(os.path.join(sample_dir, "in.json")))
                self.assertTrue(os.path.exists(os.path.join(sample_dir, "out.json")))
        finally:
            server.shutdown()
            server.server_close()

    def test_calculator_failure_submits_failed_result(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_worker_config(tmpdir))

                sample = Sample(
                    uuid="calc-failed",
                    u_coords=np.array([0.5]),
                    execution_plan=[
                        ExecutionStep(type="calculator", name="MissingCalc", layer=0),
                    ],
                )
                sample.observables = {"x": 0.5, "uuid": "calc-failed"}
                assert factory.redis is not None
                factory.redis.push_task(sample.to_task_dict())

                deadline = time.monotonic() + 15.0
                result = None
                while time.monotonic() < deadline:
                    result = factory.redis.pull_result(timeout=1)
                    if result is not None:
                        break
                factory.shutdown()

                self.assertIsNotNone(result)
                assert result is not None
                self.assertEqual(result["status"], "Failed")
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()