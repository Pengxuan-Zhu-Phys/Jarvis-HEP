#!/usr/bin/env python3
"""WP-D7.1 slow-regime distributed acceptance gates."""

from __future__ import annotations

import json
import os
import platform
import signal
import tempfile
import time
import unittest
from datetime import datetime, timezone
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.archive_handoff import resolve_staging_dir
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.factory import TaskFactory

from test_layer_concurrency import _slow_calc_module
from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    FIXTURES,
    LIKELIHOOD_EXPRESSIONS,
    _load_csv_points,
    _normalize_database_records,
    _sample_tree_file_sets,
    _start_tcp_fakeredis,
    _worker_config,
)
from test_worker_failure import _wait_until, _worker_is_busy
from test_worker_mvp import BENCHMARK_OPERA_MODULE, SAMPLING_VARIABLES


TESTS_ROOT = os.path.dirname(__file__)
BENCHMARK_BIN = os.path.join(TESTS_ROOT, "benchmark_project", "bin")
ACCEPTANCE_SLOW_SCRIPT = os.path.join(BENCHMARK_BIN, "acceptance_slow.py")
SLOW_SLEEP_SEC = 0.35

# Machine-relative gates (scale to the host running pytest).
MIN_MONITOR_HZ = 60
MIN_WORKER_SPEEDUP_RATIO = 1.25
MAX_STAGING_BACKLOG = 12
HIGH_WORKER_CAP = min(4, os.cpu_count() or 2)


def _acceptance_slow_module() -> dict[str, Any]:
    return _slow_calc_module("AcceptanceSlow", ACCEPTANCE_SLOW_SCRIPT, "z")


def _slow_worker_config(tmpdir: str, *, workers: int) -> dict[str, Any]:
    module = _acceptance_slow_module()
    return {
        "sample_config": {
            "task_result_dir": tmpdir,
            "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
            "sample_artifacts": "always",
            "workflow_has_calculator": True,
            "workflow_references_sdir": True,
        },
        "mapper": {"type": "identity", "keys": ["x", "y"]},
        "calculator_modules": [module],
        "calculator_pools": {"AcceptanceSlow": max(1, workers)},
        "likelihood_expressions": [],
        "pull_timeout": 1,
        "handoff_to_staging": False,
    }


def _count_staging_entries(task_result_dir: str) -> int:
    staging_root = resolve_staging_dir(task_result_dir)
    if not os.path.isdir(staging_root):
        return 0
    return sum(
        1
        for name in os.listdir(staging_root)
        if os.path.isdir(os.path.join(staging_root, name))
    )


def _records_to_csv_rows(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in records:
        rows.append(
            {
                "x": f"{float(record['x']):.12g}",
                "y": f"{float(record['y']):.12g}",
                "z": f"{float(record['z']):.12g}",
                "LogL": f"{float(record['LogL']):.12g}",
            }
        )
    return sorted(rows, key=lambda item: (item["x"], item["y"]))


def _golden_csv_rows(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    return _records_to_csv_rows(records)


def _stop_factory_workers() -> None:
    factory = TaskFactory._instance
    if factory is not None:
        try:
            factory.shutdown(wait=False)
        except Exception:
            pass
    TaskFactory.reset_instance()


class DistributedAcceptanceTests(unittest.TestCase):
    """Slow-regime acceptance gates for the distributed V2 runtime (WP-D7.1)."""

    metrics: dict[str, Any] = {}

    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def _run_slow_calculator_batch(
        self,
        *,
        workers: int,
        sample_count: int,
    ) -> dict[str, Any]:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": workers,
                            "redis": redis_config,
                            "Watchdog": {"enabled": False},
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _slow_worker_config(tmpdir, workers=workers)
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[_acceptance_slow_module()],
                    include_likelihood=False,
                )
                core.set_sampler(sampler)

                samples = []
                for index in range(sample_count):
                    sample = sampler._build_sample(
                        np.array([0.1 * (index + 1), 0.2 * (index + 1)], dtype=np.float64)
                    )
                    sample.uuid = f"scale-{workers}-{index}"
                    samples.append(sample)

                started = time.monotonic()
                core.submit_samples(samples)
                core.wait_for_results(sample_count, timeout=120.0)
                elapsed = max(1e-6, time.monotonic() - started)
                core.shutdown(wait=True)
                return {
                    "workers": workers,
                    "sample_count": sample_count,
                    "elapsed_sec": elapsed,
                    "samples_per_sec": sample_count / elapsed,
                }
        finally:
            server.shutdown()
            server.server_close()

    def test_gate_worker_scaling_near_linear(self) -> None:
        sample_count = 4
        one = self._run_slow_calculator_batch(workers=1, sample_count=sample_count)
        two = self._run_slow_calculator_batch(workers=2, sample_count=sample_count)
        speedup = two["samples_per_sec"] / one["samples_per_sec"]
        self.assertGreaterEqual(
            speedup,
            MIN_WORKER_SPEEDUP_RATIO,
            f"expected near-linear scaling 1→2 workers; got speedup={speedup:.2f}",
        )
        DistributedAcceptanceTests.metrics["worker_scaling"] = {
            "workers_1": one,
            "workers_2": two,
            "speedup_ratio": speedup,
            "min_required_ratio": MIN_WORKER_SPEEDUP_RATIO,
        }

    def test_gate_high_worker_count_remains_stable(self) -> None:
        if HIGH_WORKER_CAP < 2:
            self.skipTest("single-core host; high-worker gate not applicable")
        sample_count = HIGH_WORKER_CAP * 2
        result = self._run_slow_calculator_batch(
            workers=HIGH_WORKER_CAP,
            sample_count=sample_count,
        )
        self.assertGreater(result["samples_per_sec"], 0.0)
        DistributedAcceptanceTests.metrics["high_worker_stability"] = {
            "workers": HIGH_WORKER_CAP,
            "cpu_count": os.cpu_count(),
            **result,
        }

    def test_gate_archive_backlog_stays_bounded(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 2,
                            "redis": redis_config,
                            "Watchdog": {"enabled": False},
                        },
                        "Calculators": {
                            "Cleanup": {"strategy": "mv_to_staging"},
                            "Archiver": {
                                "mode": "thread",
                                "batch_size": 6,
                                "flush_interval_sec": 30.0,
                            },
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _worker_config(tmpdir)
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
                for row in _load_csv_points()[:6]:
                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)

                peak_backlog = 0
                core.submit_samples(samples)

                deadline = time.monotonic() + 90.0
                while time.monotonic() < deadline:
                    peak_backlog = max(peak_backlog, _count_staging_entries(tmpdir))
                    if core._archiver_records_written() >= len(samples):
                        break
                    time.sleep(0.05)

                core.wait_for_results(len(samples), timeout=30.0)
                core.shutdown()

                self.assertLessEqual(
                    peak_backlog,
                    MAX_STAGING_BACKLOG,
                    f"staging backlog grew without bound (peak={peak_backlog})",
                )
                DistributedAcceptanceTests.metrics["archive_backlog"] = {
                    "peak_staging_entries": peak_backlog,
                    "max_allowed": MAX_STAGING_BACKLOG,
                    "samples_submitted": len(samples),
                }
        finally:
            server.shutdown()
            server.server_close()

    def test_gate_chaos_sigkill_preserves_parity(self) -> None:
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
                            "Watchdog": {
                                "enabled": True,
                                "stale_sec": 30.0,
                                "poll_interval_sec": 0.1,
                                "max_sample_retries": 3,
                            },
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = _worker_config(tmpdir)
                worker_config["test_process_delay_sec"] = 2.0
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

                points = _load_csv_points()
                samples = []
                for row in points:
                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)
                core.submit_samples(samples)

                assert core.factory is not None
                _wait_until(lambda: _worker_is_busy(core.factory), timeout=20.0)
                worker = core.factory.workers[0]
                self.assertIsNotNone(worker.pid)
                os.kill(worker.pid, signal.SIGKILL)
                worker.join(timeout=5.0)

                core.wait_for_results(len(samples), timeout=120.0)
                respawn_count = int(getattr(core.factory, "_respawn_count", 0) or 0)
                core.shutdown()

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))
                DistributedAcceptanceTests.metrics["chaos"] = {
                    "samples": len(samples),
                    "respawn_count": respawn_count,
                    "parity": "ok",
                }
        finally:
            server.shutdown()
            server.server_close()

    def test_gate_calculator_parity_matches_v1_golden(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(
                    os.path.join(FIXTURES, "expected_calculator_records.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_records = json.load(handle)
                with open(
                    os.path.join(FIXTURES, "expected_sample_files.json"),
                    encoding="utf-8",
                ) as handle:
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

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))

                tree = _sample_tree_file_sets(os.path.join(tmpdir, "SAMPLE"))
                self.assertEqual(len(tree), 10)
                for files in tree:
                    self.assertEqual(files, expected_files)

                csv_rows = _records_to_csv_rows(records)
                golden_csv_rows = _golden_csv_rows(expected_records)
                self.assertEqual(csv_rows, golden_csv_rows)
                DistributedAcceptanceTests.metrics["parity"] = {
                    "database_rows": len(records),
                    "sample_dirs": len(tree),
                    "csv_rows": len(csv_rows),
                }
        finally:
            server.shutdown()
            server.server_close()

    def test_gate_monitor_snapshot_sustains_60hz(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 1,
                            "redis": redis_config,
                            "Watchdog": {"enabled": False},
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                worker_config = {
                    **_worker_config(tmpdir),
                    "opera_modules": [BENCHMARK_OPERA_MODULE],
                    "mapper": {"type": "flat", "variables": SAMPLING_VARIABLES},
                    "likelihood_expressions": LIKELIHOOD_EXPRESSIONS,
                }
                core.init_factory(worker_config)
                assert core.factory is not None

                ready_deadline = time.monotonic() + 2.0
                while time.monotonic() < ready_deadline:
                    if core.get_monitor_snapshot().get("workers_total") is not None:
                        break
                    time.sleep(0.05)
                self.assertIsNotNone(core.get_monitor_snapshot().get("workers_total"))

                started = time.monotonic()
                snapshot_calls = 0
                while time.monotonic() - started < 1.0:
                    snapshot = core.get_monitor_snapshot()
                    self.assertIn("workers_total", snapshot)
                    snapshot_calls += 1

                self.assertGreaterEqual(snapshot_calls, MIN_MONITOR_HZ)

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    opera_modules=[BENCHMARK_OPERA_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                samples = []
                for index in range(3):
                    sample = sampler._build_sample(
                        np.array([0.1 * index, 0.2 * index, 0.3 * index], dtype=np.float64)
                    )
                    sample.uuid = f"monitor-gate-{index}"
                    samples.append(sample)

                workload_started = time.monotonic()
                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=30.0)
                workload_elapsed = time.monotonic() - workload_started
                core.shutdown()

                self.assertLess(workload_elapsed, 25.0)
                DistributedAcceptanceTests.metrics["monitor"] = {
                    "snapshot_hz": snapshot_calls,
                    "min_required_hz": MIN_MONITOR_HZ,
                    "workload_elapsed_sec": workload_elapsed,
                }
        finally:
            server.shutdown()
            server.server_close()


def _write_benchmark_report() -> None:
    """Persist gate metrics for docs/benchmarks (invoked after the suite)."""
    if not DistributedAcceptanceTests.metrics:
        return
    report_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "docs", "benchmarks")
    os.makedirs(report_dir, exist_ok=True)
    payload = {
        "milestone": "D7.1",
        "recorded_at_utc": datetime.now(timezone.utc).isoformat(),
        "machine": {
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "python": platform.python_version(),
        },
        "gates": DistributedAcceptanceTests.metrics,
    }
    path = os.path.join(report_dir, "d7_1_acceptance.json")
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


    @classmethod
    def tearDownClass(cls) -> None:
        _write_benchmark_report()


if __name__ == "__main__":
    unittest.main()