#!/usr/bin/env python3
"""check-modules must pin calculator PackID pool size to 1 (no per-point install)."""

from __future__ import annotations

import os
import tempfile
import unittest

from jarvishep2.calculator_pools import resolve_calculator_pools
from jarvishep2.core import Jarvis2Core
from jarvishep2.worker_config import build_worker_config


class CheckModulesPoolPolicyTests(unittest.TestCase):
    def test_check_policy_pins_make_parallel_and_pools_to_one(self) -> None:
        core = Jarvis2Core.__new__(Jarvis2Core)
        core.config = {
            "Runtime": {"workers": 8},
            "EnvReqs": {
                "V2": {
                    "workers": 8,
                    "sample_directory": {"pack": True, "enabled": True},
                }
            },
            "Scan": {"name": "smoke"},
            "Calculators": {
                "make_parallel": 8,
                "Pools": {"MicroOMEGAs_Vector": 8},
                "Modules": [
                    {
                        "name": "MicroOMEGAs_Vector",
                        "make_parallel": 8,
                        "path": "&J/calculators/microOMEGAs/@PackID",
                    }
                ],
                "Archiver": {"pack_buckets": True},
            },
            "Operas": {"make_parallel": 4},
        }
        core.runtime = {}
        core.info = {}
        with tempfile.TemporaryDirectory() as tmpdir:
            core.info["task_result_dir"] = tmpdir
            core.config["task_result_dir"] = tmpdir
            core._apply_check_modules_runtime_policy()

            calc = core.config["Calculators"]
            self.assertEqual(calc["make_parallel"], 1)
            self.assertEqual(calc["Modules"][0]["make_parallel"], 1)
            self.assertEqual(calc["Pools"]["MicroOMEGAs_Vector"], 1)
            self.assertFalse(calc["Archiver"]["pack_buckets"])
            self.assertEqual(core.config["Runtime"]["workers"], 1)
            self.assertEqual(core.config["EnvReqs"]["V2"]["workers"], 1)
            self.assertEqual(core.config["Operas"]["make_parallel"], 1)

            worker = build_worker_config(core.config, task_result_dir=tmpdir)
            pools = resolve_calculator_pools(worker)
            self.assertEqual(pools.get("MicroOMEGAs_Vector"), 1)
            self.assertTrue(
                os.path.isdir(os.path.join(tmpdir, "SAMPLE", "test"))
                or core.info.get("sample_root", "").endswith(
                    os.path.join("SAMPLE", "test")
                )
            )


if __name__ == "__main__":
    unittest.main()
