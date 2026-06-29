#!/usr/bin/env python3
"""WP-D2.3 clone_shadow isolation + LibDeps symlink tests."""

from __future__ import annotations

import json
import os
import tempfile
import unittest

import numpy as np

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.sample import ExecutionStep, Sample
from jarvishep2.workflow import build_execution_plan

from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    FIXTURES,
    _load_csv_points,
    _normalize_database_records,
    _start_tcp_fakeredis,
    _worker_config,
)

TESTS_ROOT = os.path.dirname(__file__)
SHADOW_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "shadow_calc", "source")
SAFE_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "safe_calc", "source")
STATEFUL_SCRIPT = os.path.join(SHADOW_SOURCE, "stateful.py")
SAFE_SCRIPT = os.path.join(SAFE_SOURCE, "safe.py")


def _shadow_calc_module(runtime_root: str) -> dict:
    slot = os.path.join(runtime_root, "ShadowCalc", "@PackID")
    return {
        "name": "ShadowCalc",
        "required_modules": [],
        "clone_shadow": True,
        "path": slot,
        "source": SHADOW_SOURCE,
        "installation": [],
        "initialization": [],
        "execution": {
            "path": slot,
            "commands": [{"cmd": "python3 stateful.py output.json", "cwd": slot}],
            "input": [],
            "output": [
                {
                    "name": "oupjson",
                    "path": os.path.join(slot, "output.json"),
                    "type": "JSON",
                    "save": False,
                    "variables": [{"name": "isolation_count", "entry": "isolation_count"}],
                }
            ],
        },
    }


def _safe_calc_module() -> dict:
    return {
        "name": "SafeCalc",
        "required_modules": [],
        "clone_shadow": False,
        "path": SAFE_SOURCE,
        "source": SAFE_SOURCE,
        "symlink_name": "SafeCalc",
        "installation": [],
        "initialization": [],
        "execution": {
            "path": SAFE_SOURCE,
            "commands": [
                {
                    "cmd": "python3 @Sdir/SafeCalc/safe.py @Sdir/output_safe.json",
                    "cwd": "@Sdir",
                }
            ],
            "input": [],
            "output": [
                {
                    "name": "oupjson",
                    "path": "@Sdir/output_safe.json",
                    "type": "JSON",
                    "save": False,
                    "variables": [{"name": "safe", "entry": "safe"}],
                }
            ],
        },
    }


class CalculatorShadowUnitTests(unittest.TestCase):
    def test_decode_shadow_path_replaces_pack_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            module = CalculatorModule(
                "ShadowCalc",
                _shadow_calc_module(os.path.join(tmpdir, "runtime")),
            )
            module.acquire_pack_id("pack-abc")
            decoded = module.decode_shadow_path(module.basepath)
            self.assertIn("pack-abc", decoded)
            self.assertNotIn("@PackID", decoded)

    def test_decode_shadow_commands_replaces_pack_id(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            module = CalculatorModule(
                "ShadowCalc",
                _shadow_calc_module(os.path.join(tmpdir, "runtime")),
            )
            module.acquire_pack_id("pack-xyz")
            decoded = module.decode_shadow_commands({"cmd": "cd @PackID", "cwd": "@PackID"})
            self.assertEqual(decoded["cmd"], "cd pack-xyz")
            self.assertEqual(decoded["cwd"], "pack-xyz")


class CloneShadowIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def _run_samples(
        self,
        modules: list[dict],
        *,
        workers: int = 1,
        sample_count: int = 1,
    ) -> list[dict]:
        server, redis_config = _start_tcp_fakeredis()
        records: list[dict] = []
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                runtime_root = os.path.join(tmpdir, "runtime")
                resolved_modules = []
                for item in modules:
                    module = dict(item)
                    if module.get("name") == "ShadowCalc":
                        module = _shadow_calc_module(runtime_root)
                    resolved_modules.append(module)

                plan = build_execution_plan(
                    calculator_modules=resolved_modules,
                    include_likelihood=False,
                )
                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": workers,
                            "redis": redis_config,
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _worker_config(tmpdir)
                worker_config["calculator_modules"] = resolved_modules
                worker_config["likelihood_expressions"] = []
                pools = {str(m["name"]): max(2, workers) for m in resolved_modules}
                worker_config["calculator_pools"] = pools
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler._execution_plan_template = [step.to_dict() for step in plan]
                core.set_sampler(sampler)

                samples = []
                for _ in range(sample_count):
                    samples.append(sampler._build_sample([0.0, 0.0]))
                core.submit_samples(samples)
                core.wait_for_results(sample_count, timeout=60.0)

                from jarvishep2.database import SimpleHDF5Writer

                records = SimpleHDF5Writer(
                    os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                ).read_records()
                core.shutdown()
        finally:
            server.shutdown()
            server.server_close()
        return records

    def test_shadow_calculator_isolated_under_concurrency(self) -> None:
        records = self._run_samples(
            [_shadow_calc_module("/tmp/unused")],
            workers=2,
            sample_count=2,
        )
        self.assertEqual(len(records), 2)
        for record in records:
            self.assertEqual(int(record.get("isolation_count", 0)), 1)

    def test_safe_calculator_uses_symlink_without_copy(self) -> None:
        marker = os.path.join(SAFE_SOURCE, ".safe_marker")
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
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
                safe_module = _safe_calc_module()
                worker_config = _worker_config(tmpdir)
                worker_config["calculator_modules"] = [safe_module]
                worker_config["calculator_pools"] = {"SafeCalc": 1}
                worker_config["likelihood_expressions"] = []
                core.init_factory(worker_config)

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[safe_module],
                    include_likelihood=False,
                )
                core.set_sampler(sampler)
                sample = sampler._build_sample([0.0, 0.0])
                sample.set_config(worker_config["sample_config"])
                sample.start()
                sample.materialize(worker_id="0")

                module = CalculatorModule("SafeCalc", safe_module)
                module.acquire_pack_id("pack-safe")
                link_path = module.ensure_symlink_runtime(sample.info)
                assert link_path is not None
                self.assertTrue(os.path.islink(link_path))
                with open(marker, "w", encoding="utf-8") as handle:
                    handle.write("shared")
                self.assertTrue(os.path.exists(os.path.join(os.path.realpath(link_path), ".safe_marker")))
                core.shutdown()
        finally:
            server.shutdown()
            server.server_close()
            if os.path.exists(marker):
                os.remove(marker)

    def test_eggbox_parity_preserved_with_clone_shadow_false(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(
                    os.path.join(FIXTURES, "expected_calculator_records.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_records = json.load(handle)

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
                worker_config = _worker_config(tmpdir)
                worker_config["calculator_modules"] = [EGGBOX_CALC_MODULE]
                worker_config["calculator_pools"] = {"EggBox": 1}
                core.init_factory(worker_config)
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[EGGBOX_CALC_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                samples = []
                for row in _load_csv_points():
                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)

                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=90.0)
                core.shutdown()

                from jarvishep2.database import SimpleHDF5Writer

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()