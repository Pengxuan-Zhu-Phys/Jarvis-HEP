#!/usr/bin/env python3
"""WP-D2.1 multi-Worker pool + Redis calculator free-pool tests."""

from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.calculator_pools import resolve_calculator_pools
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import calc_status_busy_field, calc_status_free_field

from test_layer_concurrency import SLEEP_SEC, SLOW_A_MODULE
from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    FIXTURES,
    LIKELIHOOD_EXPRESSIONS,
    PARITY_PROJECT,
    _load_csv_points,
    _normalize_database_records,
    _start_tcp_fakeredis,
    _worker_config,
)

FAIL_CALC_MODULE = {
    "name": "FailCalc",
    "required_modules": [],
    "clone_shadow": False,
    "path": ".",
    "installation": [],
    "initialization": [],
    "execution": {
        "path": ".",
        "commands": [{"cmd": "python3 -c 'import sys; sys.exit(42)'", "cwd": "@Sdir"}],
        "input": [],
        "output": [],
    },
}


class CalculatorPoolConfigTests(unittest.TestCase):
    def test_resolve_calculator_pools_from_modules(self) -> None:
        config = {
            "calculator_modules": [
                {"name": "EggBox", "make_paraller": 3},
            ]
        }
        self.assertEqual(resolve_calculator_pools(config), {"EggBox": 3})

    def test_explicit_calculator_pools_override(self) -> None:
        config = {
            "calculator_pools": {"EggBox": 2},
            "calculator_modules": [{"name": "EggBox", "make_paraller": 9}],
        }
        self.assertEqual(resolve_calculator_pools(config), {"EggBox": 2})


class WorkerPoolIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def _run_slow_samples(self, *, workers: int, sample_count: int) -> float:
        server, redis_config = _start_tcp_fakeredis()
        started = time.monotonic()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
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
                worker_config["calculator_modules"] = [SLOW_A_MODULE]
                worker_config["calculator_pools"] = {"SlowA": workers}
                worker_config["likelihood_expressions"] = []
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[SLOW_A_MODULE],
                    include_likelihood=False,
                )
                core.set_sampler(sampler)

                samples = [sampler._build_sample([0.0, 0.0]) for _ in range(sample_count)]
                core.submit_samples(samples)
                core.wait_for_results(sample_count, timeout=120.0)
                core.shutdown()
        finally:
            server.shutdown()
            server.server_close()
        return time.monotonic() - started

    def test_throughput_scales_from_one_to_two_workers(self) -> None:
        sample_count = 12
        min_gain = SLEEP_SEC * (sample_count / 2) * 0.25

        def _measure() -> tuple[float, float]:
            one_elapsed = self._run_slow_samples(workers=1, sample_count=sample_count)
            two_elapsed = self._run_slow_samples(workers=2, sample_count=sample_count)
            return one_elapsed, two_elapsed

        one_worker_elapsed, two_worker_elapsed = _measure()
        gain = one_worker_elapsed - two_worker_elapsed
        if gain < min_gain:
            # Retry once when the host is busy (full-suite load).
            one_worker_elapsed, two_worker_elapsed = _measure()
            gain = one_worker_elapsed - two_worker_elapsed

        min_per_sample = SLEEP_SEC * 0.75
        self.assertGreater(one_worker_elapsed, min_per_sample * sample_count * 0.5)
        self.assertGreater(two_worker_elapsed, min_per_sample * (sample_count / 2) * 0.5)
        self.assertLess(two_worker_elapsed, one_worker_elapsed)
        self.assertGreater(gain, min_gain)

    def test_slots_released_after_forced_calculator_failure(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                pool_size = 1
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
                worker_config["calculator_modules"] = [FAIL_CALC_MODULE]
                worker_config["calculator_pools"] = {"FailCalc": pool_size}
                worker_config["likelihood_expressions"] = []
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[FAIL_CALC_MODULE],
                    include_likelihood=False,
                )
                core.set_sampler(sampler)

                sample = sampler._build_sample([0.0, 0.0])
                core.submit_samples([sample])
                core.wait_for_results(1, timeout=30.0)

                assert core.redis is not None
                final_status = core.redis.fetch_calculator_status()
                core.shutdown()

                free_field = calc_status_free_field("FailCalc")
                busy_field = calc_status_busy_field("FailCalc")
                self.assertEqual(int(final_status.get(free_field, 0) or 0), pool_size)
                self.assertEqual(int(final_status.get(busy_field, 0) or 0), 0)
        finally:
            server.shutdown()
            server.server_close()

    def test_two_workers_eggbox_database_parity(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(os.path.join(FIXTURES, "expected_calculator_records.json"), encoding="utf-8") as handle:
                    expected_records = json.load(handle)

                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 2,
                            "redis": redis_config,
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _worker_config(tmpdir)
                worker_config["calculator_pools"] = {"EggBox": 2}
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

    def test_calc_pool_busy_never_exceeds_configured_slots(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                pool_size = 1
                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 2,
                            "redis": redis_config,
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _worker_config(tmpdir)
                worker_config["calculator_pools"] = {"EggBox": pool_size}
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[EGGBOX_CALC_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                assert core.redis is not None
                polling = True
                peak_busy = 0
                lock = threading.Lock()

                def _poll_busy() -> None:
                    nonlocal peak_busy
                    busy_field = calc_status_busy_field("EggBox")
                    while polling:
                        status = core.redis.fetch_calculator_status()
                        busy = int(status.get(busy_field, 0) or 0)
                        with lock:
                            peak_busy = max(peak_busy, busy)
                        time.sleep(0.01)

                poll_thread = threading.Thread(target=_poll_busy, daemon=True)
                poll_thread.start()

                samples = []
                for row in _load_csv_points():
                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)
                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=90.0)
                polling = False
                poll_thread.join(timeout=2.0)
                assert core.redis is not None
                final_status = core.redis.fetch_calculator_status()
                core.shutdown()

                self.assertLessEqual(peak_busy, pool_size)
                free_field = calc_status_free_field("EggBox")
                self.assertEqual(int(final_status.get(free_field, 0) or 0), pool_size)
                self.assertEqual(int(final_status.get(calc_status_busy_field("EggBox"), 0) or 0), 0)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()