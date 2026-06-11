#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BENCHMARK_PROJECT = os.path.join(PROJECT_ROOT, "tests", "benchmark_project")
BENCHMARK_TASK = os.path.join(BENCHMARK_PROJECT, "bin", "benchmark_random_operas.yaml")
BENCHMARK_OUTPUT = os.path.join(
    BENCHMARK_PROJECT,
    "outputs",
    "benchmark_random_operas",
)


class BenchmarkModeTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()

    def test_benchmark_mode_writes_positive_throughput_json(self):
        env = dict(os.environ)
        env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "jarvishep",
                BENCHMARK_TASK,
                "--benchmark",
                "0.25",
                "--skip-draw-flowchart",
            ],
            cwd=BENCHMARK_PROJECT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
        )
        self.assertEqual(
            proc.returncode,
            0,
            msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}",
        )

        benchmark_json = os.path.join(BENCHMARK_OUTPUT, "benchmark.json")
        self.assertTrue(os.path.exists(benchmark_json), benchmark_json)
        with open(benchmark_json, "r", encoding="utf-8") as handle:
            payload = json.load(handle)

        self.assertGreater(payload["samples_per_sec"], 0)
        self.assertGreater(payload["total_submitted"], 0)
        self.assertGreater(payload["total_completed"], 0)
        self.assertEqual(payload["total_failed"], 0)
        self.assertIn("run_summary", payload)
        self.assertEqual(
            payload["run_summary"]["total_points_finished"],
            payload["total_completed"],
        )
        self.assertIn("stage_timers", payload)
        stages = payload["stage_timers"]["stages"]
        for stage in {
            "sample_setup",
            "submit_future_locks",
            "module_dispatch",
            "module_execution",
            "hdf5_enqueue",
            "sampler_loop",
        }:
            self.assertIn(stage, stages)
            self.assertGreater(stages[stage]["total_sec"], 0.0)


if __name__ == "__main__":
    unittest.main()
