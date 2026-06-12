#!/usr/bin/env python3
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import unittest

import h5py
import yaml


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

PARITY_PROJECT = os.path.join(PROJECT_ROOT, "tests", "parity_project")
PARITY_TASK = os.path.join(PARITY_PROJECT, "bin", "parity_calculator_check.yaml")
PARITY_OUTPUT = os.path.join(PARITY_PROJECT, "outputs", "parity_calculator_check")
BENCHMARK_PROJECT = os.path.join(PROJECT_ROOT, "tests", "benchmark_project")
BENCHMARK_TASK = os.path.join(BENCHMARK_PROJECT, "bin", "benchmark_random_operas.yaml")
BENCHMARK_OUTPUT = os.path.join(BENCHMARK_PROJECT, "outputs", "benchmark_random_operas")
FIXTURES = os.path.join(PROJECT_ROOT, "tests", "fixtures", "parity_m1")

# WP-0.1 baseline on this machine (ledger); M1 scaled gate is >=2.5×.
WP01_BASELINE_SAMPLES_PER_SEC = 311.88
M1_SCALED_MIN_SAMPLES_PER_SEC = WP01_BASELINE_SAMPLES_PER_SEC * 2.5


def _run_task(
    *,
    project_root: str,
    task_path: str,
    extra_args: list[str] | None = None,
    runtime_overrides: dict | None = None,
) -> subprocess.CompletedProcess[str]:
    with open(task_path, "r", encoding="utf-8") as handle:
        payload = yaml.safe_load(handle)

    if runtime_overrides is not None:
        runtime = dict(payload.get("Runtime", {}) or {})
        runtime.update(runtime_overrides)
        payload["Runtime"] = runtime

    generated = os.path.join(project_root, "bin", "_output_parity_task.yaml")
    with open(generated, "w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False)

    env = dict(os.environ)
    env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
    return subprocess.run(
        [sys.executable, "-m", "jarvishep", generated, *(extra_args or [])],
        cwd=project_root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=120,
    )


def _read_database_records(output_dir: str) -> list[dict]:
    records = []
    for h5_path in sorted(glob.glob(os.path.join(output_dir, "DATABASE", "samples*.hdf5"))):
        with h5py.File(h5_path, "r") as handle:
            if "records" not in handle:
                continue
            for item in handle["records"][()]:
                if isinstance(item, bytes):
                    item = json.loads(item.decode("utf-8"))
                if isinstance(item, dict):
                    records.append(dict(item))
    return sorted(records, key=lambda row: float(row.get("x", 0)))


def _normalize_calculator_records(records: list[dict]) -> list[dict]:
    normalized = []
    for row in records:
        normalized.append(
            {
                "x": float(row["x"]),
                "calc_z": float(row["calc_z"]),
                "LogL": float(row["LogL"]),
            }
        )
    return normalized


def _sample_tree_file_sets(sample_root: str) -> list[list[str]]:
    if not os.path.isdir(sample_root):
        return []
    manifests = []
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


class CalculatorOutputParityTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(PARITY_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(PARITY_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(PARITY_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(PARITY_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()
        generated = os.path.join(PARITY_PROJECT, "bin", "_output_parity_task.yaml")
        if os.path.exists(generated):
            os.remove(generated)

    def _run_check_modules(self, *, sample_artifacts: str) -> None:
        proc = _run_task(
            project_root=PARITY_PROJECT,
            task_path=PARITY_TASK,
            extra_args=["--check-modules", "--skip-draw-flowchart"],
            runtime_overrides={"sample_artifacts": sample_artifacts},
        )
        self.assertEqual(
            proc.returncode,
            0,
            msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
        )

    def test_calculator_always_matches_pinned_expectation(self):
        with open(
            os.path.join(FIXTURES, "expected_calculator_records.json"),
            "r",
            encoding="utf-8",
        ) as handle:
            expected_records = json.load(handle)
        with open(
            os.path.join(FIXTURES, "expected_sample_files.json"),
            "r",
            encoding="utf-8",
        ) as handle:
            expected_files = json.load(handle)

        self._run_check_modules(sample_artifacts="always")
        records = _normalize_calculator_records(_read_database_records(PARITY_OUTPUT))
        self.assertEqual(records, expected_records)

        sample_root = os.path.join(PARITY_OUTPUT, "SAMPLE", "tests")
        tree = _sample_tree_file_sets(sample_root)
        self.assertEqual(len(tree), 10)
        for files in tree:
            self.assertEqual(files, expected_files)

    def test_calculator_auto_matches_always_outputs(self):
        self._run_check_modules(sample_artifacts="always")
        always_records = _normalize_calculator_records(_read_database_records(PARITY_OUTPUT))
        always_tree = _sample_tree_file_sets(os.path.join(PARITY_OUTPUT, "SAMPLE", "tests"))

        shutil.rmtree(os.path.join(PARITY_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(PARITY_PROJECT, "logs"), ignore_errors=True)

        self._run_check_modules(sample_artifacts="auto")
        auto_records = _normalize_calculator_records(_read_database_records(PARITY_OUTPUT))
        auto_tree = _sample_tree_file_sets(os.path.join(PARITY_OUTPUT, "SAMPLE", "tests"))

        self.assertEqual(auto_records, always_records)
        self.assertEqual(auto_tree, always_tree)


class OperaOutputParityTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()
        generated = os.path.join(BENCHMARK_PROJECT, "bin", "_output_parity_task.yaml")
        if os.path.exists(generated):
            os.remove(generated)

    def _assert_operas_database_contract(self, *, sample_artifacts: str | None) -> int:
        proc = _run_task(
            project_root=BENCHMARK_PROJECT,
            task_path=BENCHMARK_TASK,
            extra_args=["--benchmark", "0.35", "--skip-draw-flowchart"],
            runtime_overrides=(
                {"sample_artifacts": sample_artifacts}
                if sample_artifacts is not None
                else None
            ),
        )
        self.assertEqual(proc.returncode, 0, msg=proc.stderr)
        records = _read_database_records(BENCHMARK_OUTPUT)
        benchmark = json.load(
            open(os.path.join(BENCHMARK_OUTPUT, "benchmark.json"), "r", encoding="utf-8")
        )
        expected_keys = {"x", "y", "shift", "uuid", "z", "LogL"}
        self.assertGreater(len(records), 0)
        self.assertEqual(len(records), benchmark["total_completed"])
        self.assertTrue(expected_keys.issubset(records[0].keys()))
        self.assertNotEqual(records[0]["LogL"], float("-inf"))
        return len(records)

    def test_operas_auto_and_always_share_database_contract(self):
        auto_count = self._assert_operas_database_contract(sample_artifacts=None)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        always_count = self._assert_operas_database_contract(sample_artifacts="always")
        self.assertGreater(auto_count, 0)
        self.assertGreater(always_count, 0)


class M1BenchmarkGateTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()

    def test_m1_throughput_meets_scaled_gate(self):
        proc = _run_task(
            project_root=BENCHMARK_PROJECT,
            task_path=BENCHMARK_TASK,
            extra_args=["--benchmark", "10", "--skip-draw-flowchart"],
            runtime_overrides={"sample_artifacts": "auto"},
        )
        self.assertEqual(proc.returncode, 0, msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

        benchmark_json = os.path.join(BENCHMARK_OUTPUT, "benchmark.json")
        with open(benchmark_json, "r", encoding="utf-8") as handle:
            payload = json.load(handle)

        self.assertGreaterEqual(
            payload["samples_per_sec"],
            M1_SCALED_MIN_SAMPLES_PER_SEC,
            msg=(
                f"M1 scaled gate requires >= {M1_SCALED_MIN_SAMPLES_PER_SEC:.2f} samples/s "
                f"(2.5× WP-0.1 baseline {WP01_BASELINE_SAMPLES_PER_SEC:.2f}); "
                f"measured {payload['samples_per_sec']:.2f}"
            ),
        )
        self.assertEqual(payload["total_failed"], 0)


if __name__ == "__main__":
    unittest.main()