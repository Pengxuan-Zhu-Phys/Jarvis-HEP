#!/usr/bin/env python3
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

import h5py  # noqa: E402
import yaml  # noqa: E402

from jarvishep.Module.calculator import CalculatorModule  # noqa: E402
from jarvishep.Module.likelihood import LogLikelihood  # noqa: E402
from jarvishep.runtime_config import (  # noqa: E402
    normalize_runtime_block,
    should_eager_materialize,
    workflow_references_sdir,
)
from jarvishep.sample import Sample, ensure_sample_materialized, materialize_failure_artifacts  # noqa: E402


BENCHMARK_PROJECT = os.path.join(PROJECT_ROOT, "tests", "benchmark_project")
BENCHMARK_TASK = os.path.join(BENCHMARK_PROJECT, "bin", "benchmark_random_operas.yaml")
BENCHMARK_OUTPUT = os.path.join(
    BENCHMARK_PROJECT,
    "outputs",
    "benchmark_random_operas",
)


def _run_benchmark_task(
    *,
    runtime_overrides: dict | None = None,
    seconds: str = "0.5",
) -> subprocess.CompletedProcess[str]:
    with open(BENCHMARK_TASK, "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)

    if runtime_overrides is not None:
        runtime = dict(payload.get("Runtime", {}) or {})
        runtime.update(runtime_overrides)
        payload["Runtime"] = runtime

    task_path = os.path.join(BENCHMARK_PROJECT, "bin", "_lazy_materialization_task.yaml")
    with open(task_path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)

    env = dict(os.environ)
    env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "jarvishep",
            task_path,
            "--benchmark",
            seconds,
            "--skip-draw-flowchart",
        ],
        cwd=BENCHMARK_PROJECT,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=60,
    )


def _list_sample_uuid_dirs(sample_root: str) -> list[str]:
    found = []
    if not os.path.isdir(sample_root):
        return found
    for bucket_name in os.listdir(sample_root):
        bucket_path = os.path.join(sample_root, bucket_name)
        if not os.path.isdir(bucket_path):
            continue
        for child in os.listdir(bucket_path):
            child_path = os.path.join(bucket_path, child)
            if os.path.isdir(child_path):
                found.append(child_path)
    for child in os.listdir(sample_root):
        child_path = os.path.join(sample_root, child)
        if os.path.isdir(child_path) and child not in {"archive"}:
            # Direct uuid dirs created by failure materialization.
            if os.path.exists(os.path.join(child_path, "Sample_running.log")) or len(child) >= 32:
                found.append(child_path)
    return sorted(set(found))


def _read_database_records(output_dir: str) -> list[dict]:
    records = []
    db_dir = os.path.join(output_dir, "DATABASE")
    h5_paths = sorted(glob.glob(os.path.join(db_dir, "samples*.hdf5")))
    for h5_path in h5_paths:
        with h5py.File(h5_path, "r") as handle:
            if "records" not in handle:
                continue
            dataset = handle["records"]
            for item in dataset[()]:
                if isinstance(item, bytes):
                    item = json.loads(item.decode("utf-8"))
                if isinstance(item, dict):
                    records.append(dict(item))
    return records


class RuntimeConfigTests(unittest.TestCase):
    def test_runtime_defaults_normalize(self):
        runtime = normalize_runtime_block(None)
        self.assertEqual(runtime["mode"], "auto")
        self.assertEqual(runtime["workers"], 0)
        self.assertEqual(runtime["batch_size"], 256)
        self.assertEqual(runtime["sample_artifacts"], "auto")

    def test_workflow_references_sdir_detects_operas_config(self):
        config = {
            "Operas": {
                "Modules": [
                    {
                        "name": "Dump",
                        "operator": "helper.echo",
                        "kwargs": {"path": "@Sdir/result.json"},
                    }
                ]
            }
        }
        self.assertTrue(workflow_references_sdir(config))


class SampleMaterializationUnitTests(unittest.TestCase):
    def test_lazy_sample_exposes_logger_name_for_likelihood(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"xx": 1.0, "yy": 2.0, "z": 100.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            self.assertIn("logger_name", sample.info)
            result = LogLikelihood(
                [{"name": "LogL_Z", "expression": "LogGauss(z, 100, 10)"}]
            ).calculate(sample.observables, sample.info)
            self.assertNotEqual(result["LogL"], float("-inf"))

    def test_auto_operas_only_defers_directory_creation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            self.assertFalse(sample.info["_materialized"])
            self.assertIsNone(sample.info["save_dir"])
            self.assertFalse(os.path.exists(os.path.join(tmpdir, sample.uuid)))

    def test_always_mode_materializes_immediately(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                }
            )
            self.assertTrue(sample.info["_materialized"])
            self.assertTrue(os.path.isdir(sample.info["save_dir"]))
            self.assertTrue(os.path.exists(sample.info["run_log"]))

    def test_never_mode_skips_failure_materialization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "dead-beef",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "never",
                "_materialized": False,
            }
            result = materialize_failure_artifacts(sample_info, error="boom")
            self.assertIsNone(result)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, "dead-beef")))

    def test_failure_materialization_writes_minimal_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "dead-beef",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "auto",
                "workflow_has_calculator": False,
                "workflow_references_sdir": False,
                "_materialized": False,
            }
            save_dir = materialize_failure_artifacts(sample_info, error="boom")
            self.assertIsNotNone(save_dir)
            log_path = os.path.join(save_dir, "Sample_running.log")
            self.assertTrue(os.path.exists(log_path))
            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertIn("Sample failed", handle.read())

    def test_sdir_resolution_triggers_materialization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "sdir-case",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "auto",
                "workflow_has_calculator": False,
                "workflow_references_sdir": False,
                "_materialized": False,
            }
            save_dir = ensure_sample_materialized(sample_info)
            self.assertIsNotNone(save_dir)
            self.assertTrue(os.path.isdir(save_dir))

    def test_calculator_workflow_is_eager_in_auto_mode(self):
        self.assertTrue(
            should_eager_materialize(
                {
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": False,
                }
            )
        )

    def test_operas_sdir_reference_is_eager_in_auto_mode(self):
        self.assertTrue(
            should_eager_materialize(
                {
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": True,
                }
            )
        )

    def test_calculator_workflow_materializes_per_sample_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "save_dir": os.path.join(tmpdir, "000001"),
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": False,
                }
            )
            self.assertTrue(sample.info["_materialized"])
            self.assertTrue(os.path.isdir(sample.info["save_dir"]))
            self.assertTrue(os.path.exists(sample.info["run_log"]))

    def test_calculator_sdir_token_materializes_lazy_sample(self):
        module = CalculatorModule(
            "DemoCalc",
            {
                "modes": False,
                "required_modules": [],
                "clone_shadow": False,
                "installation": [],
                "initialization": [],
                "execution": {"commands": [], "input": [], "output": []},
                "path": "&J/calculators/runtime/program/demo",
            },
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_info = {
                "uuid": "calc-sdir",
                "params": {"x": 1.0},
                "sample_dirs": tmpdir,
                "task_result_dir": tmpdir,
                "sample_artifacts": "auto",
                "workflow_has_calculator": False,
                "workflow_references_sdir": False,
                "_materialized": False,
            }
            module.sample_info = sample_info
            resolved = module._resolve_sample_runtime_tokens(
                "echo @Sdir/out.txt",
                stage="execution",
                field="cmd",
            )
            self.assertTrue(resolved.endswith("/out.txt"))
            self.assertTrue(os.path.isdir(sample_info["save_dir"]))


class LazyMaterializationIntegrationTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()
        task_path = os.path.join(BENCHMARK_PROJECT, "bin", "_lazy_materialization_task.yaml")
        if os.path.exists(task_path):
            os.remove(task_path)

    def test_operas_auto_scan_leaves_sample_tree_empty_on_success(self):
        proc = _run_benchmark_task(runtime_overrides=None)
        self.assertEqual(proc.returncode, 0, msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

        sample_root = os.path.join(BENCHMARK_OUTPUT, "SAMPLE")
        self.assertEqual(_list_sample_uuid_dirs(sample_root), [])

    def test_operas_always_scan_materializes_sample_dirs(self):
        proc = _run_benchmark_task(runtime_overrides={"sample_artifacts": "always"})
        self.assertEqual(proc.returncode, 0, msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

        sample_root = os.path.join(BENCHMARK_OUTPUT, "SAMPLE")
        self.assertGreater(len(_list_sample_uuid_dirs(sample_root)), 0)

    def _assert_database_contract(self, *, sample_artifacts: str) -> None:
        proc = _run_benchmark_task(
            runtime_overrides={"sample_artifacts": sample_artifacts},
            seconds="0.35",
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        records = _read_database_records(BENCHMARK_OUTPUT)
        benchmark = json.load(
            open(os.path.join(BENCHMARK_OUTPUT, "benchmark.json"), "r", encoding="utf-8")
        )
        self.assertGreater(len(records), 0)
        self.assertEqual(len(records), benchmark["total_completed"])
        expected_keys = {"x", "y", "shift", "uuid", "z", "LogL"}
        self.assertTrue(expected_keys.issubset(set(records[0].keys())))
        self.assertNotEqual(records[0]["LogL"], float("-inf"))

    def test_database_contract_auto_mode(self):
        self._assert_database_contract(sample_artifacts="auto")

    def test_database_contract_always_mode(self):
        self._assert_database_contract(sample_artifacts="always")


if __name__ == "__main__":
    unittest.main()