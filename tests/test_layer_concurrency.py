#!/usr/bin/env python3
"""WP-D2.2 layer-internal calculator concurrency tests."""

from __future__ import annotations

import os
import tempfile
import threading
import time
import unittest

from fakeredis import TcpFakeServer

from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.sample import ExecutionStep
from jarvishep2.workflow import (
    build_execution_plan,
    concurrency_groups,
    resolve_module_layers,
)

from test_worker_calculator import _start_tcp_fakeredis, _worker_config

TESTS_ROOT = os.path.dirname(__file__)
SLOW_DIR = os.path.join(TESTS_ROOT, "fixtures", "slow_calculators")
SLOW_A_SCRIPT = os.path.join(SLOW_DIR, "slow_a.py")
SLOW_B_SCRIPT = os.path.join(SLOW_DIR, "slow_b.py")
SLEEP_SEC = 0.35


def _slow_calc_module(name: str, script_path: str, output_key: str) -> dict:
    output_path = f"@Sdir/output_{output_key}.json"
    return {
        "name": name,
        "required_modules": [],
        "clone_shadow": False,
        "path": SLOW_DIR,
        "installation": [],
        "initialization": [],
        "execution": {
            "path": SLOW_DIR,
            "commands": [{"cmd": f"python3 {script_path} {output_path}", "cwd": "@Sdir"}],
            "input": [
                {
                    "name": "inpjson",
                    "path": "@Sdir/input.json",
                    "type": "JSON",
                    "save": False,
                    "actions": [
                        {
                            "type": "Dump",
                            "variables": [{"name": "seed", "expression": "1"}],
                        }
                    ],
                }
            ],
            "output": [
                {
                    "name": "oupjson",
                    "path": output_path,
                    "type": "JSON",
                    "save": False,
                    "variables": [{"name": output_key, "entry": output_key}],
                }
            ],
        },
    }


SLOW_A_MODULE = _slow_calc_module("SlowA", SLOW_A_SCRIPT, "a")
SLOW_B_MODULE = _slow_calc_module("SlowB", SLOW_B_SCRIPT, "b")


class WorkflowLayerDerivationTests(unittest.TestCase):
    def test_independent_calculators_share_layer(self) -> None:
        layers = resolve_module_layers([SLOW_A_MODULE, SLOW_B_MODULE])
        self.assertEqual(layers["SlowA"], layers["SlowB"])

    def test_dependent_calculator_lands_in_later_layer(self) -> None:
        dependent = dict(SLOW_B_MODULE)
        dependent["required_modules"] = ["SlowA"]
        layers = resolve_module_layers([SLOW_A_MODULE, dependent])
        self.assertLess(layers["SlowA"], layers["SlowB"])

    def test_concurrency_groups_match_layer_membership(self) -> None:
        plan = build_execution_plan(
            calculator_modules=[SLOW_A_MODULE, SLOW_B_MODULE],
            include_likelihood=False,
        )
        groups = concurrency_groups(plan)
        self.assertEqual(groups[0], ["SlowA", "SlowB"])


class LayerConcurrencyIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def _run_dual_calculator_sample(
        self,
        modules: list[dict],
        *,
        serial_layers: bool,
        force_serial_layers: bool = False,
    ) -> tuple[dict[str, object], float]:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                if serial_layers:
                    plan = [
                        ExecutionStep(type="calculator", name="SlowA", layer=0),
                        ExecutionStep(type="calculator", name="SlowB", layer=1),
                    ]
                else:
                    plan = build_execution_plan(
                        calculator_modules=modules,
                        include_likelihood=False,
                    )

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
                worker_config["calculator_modules"] = modules
                worker_config["calculator_pools"] = {"SlowA": 2, "SlowB": 2}
                worker_config["likelihood_expressions"] = []
                worker_config["force_serial_layers"] = force_serial_layers
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler._calculator_modules = list(modules)
                sampler._execution_plan_template = [step.to_dict() for step in plan]
                core.set_sampler(sampler)

                sample = sampler._build_sample([0.0, 0.0])
                started = time.monotonic()
                core.submit_samples([sample])
                core.wait_for_results(1, timeout=30.0)
                elapsed = time.monotonic() - started
                from jarvishep2.database import SimpleHDF5Writer

                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                records = SimpleHDF5Writer(db_path).read_records()
                core.shutdown()
                self.assertEqual(len(records), 1)
                return {"status": "Completed", "observables": records[0]}, elapsed
        finally:
            server.shutdown()
            server.server_close()

    def test_same_layer_runs_faster_than_serial(self) -> None:
        modules = [SLOW_A_MODULE, SLOW_B_MODULE]
        parallel_payload, parallel_elapsed = self._run_dual_calculator_sample(
            modules,
            serial_layers=False,
        )
        serial_payload, serial_elapsed = self._run_dual_calculator_sample(
            modules,
            serial_layers=True,
        )

        self.assertEqual(parallel_payload["status"], "Completed")
        self.assertEqual(serial_payload["status"], "Completed")

        def _calc_observables(payload: dict[str, object]) -> dict[str, object]:
            obs = dict(payload.get("observables") or {})
            obs.pop("uuid", None)
            return obs

        self.assertEqual(
            _calc_observables(parallel_payload),
            _calc_observables(serial_payload),
        )
        # Wall-clock includes spawn/archiver overhead; compare relative durations.
        self.assertGreater(serial_elapsed, parallel_elapsed)
        self.assertGreater(serial_elapsed - parallel_elapsed, SLEEP_SEC * 0.35)
        self.assertLess(parallel_elapsed, serial_elapsed * 0.92)

    def test_force_serial_layers_matches_explicit_serial_layers(self) -> None:
        modules = [SLOW_A_MODULE, SLOW_B_MODULE]
        forced_payload, forced_elapsed = self._run_dual_calculator_sample(
            modules,
            serial_layers=False,
            force_serial_layers=True,
        )
        explicit_payload, explicit_elapsed = self._run_dual_calculator_sample(
            modules,
            serial_layers=True,
        )

        def _calc_observables(payload: dict[str, object]) -> dict[str, object]:
            obs = dict(payload.get("observables") or {})
            obs.pop("uuid", None)
            return obs

        self.assertEqual(
            _calc_observables(forced_payload),
            _calc_observables(explicit_payload),
        )
        # Both paths run the calculators serially; allow spawn/archiver jitter.
        min_serial = SLEEP_SEC * 2 * 0.75
        self.assertGreater(forced_elapsed, min_serial)
        self.assertGreater(explicit_elapsed, min_serial)
        self.assertLess(abs(forced_elapsed - explicit_elapsed), SLEEP_SEC * 1.6)


if __name__ == "__main__":
    unittest.main()