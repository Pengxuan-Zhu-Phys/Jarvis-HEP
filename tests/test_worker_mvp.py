#!/usr/bin/env python3
"""WP-D1.1 Worker MVP tests — fakeredis + spawn end-to-end opera path."""

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
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.redis_queue import RedisQueue, make_fakeredis_queue
from jarvishep2.sample import ExecutionStep, Sample
from jarvishep2.worker import Worker


FIXTURES = os.path.join(os.path.dirname(__file__), "fixtures")
GOLDEN_PATH = os.path.join(FIXTURES, "golden_operas_records.json")

BENCHMARK_OPERA_MODULE = {
    "name": "TrivialEggbox",
    "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
    "call_mode": "call",
    "input": [
        {"name": "x", "expression": "x + shift * 0"},
        {"name": "y", "expression": "y"},
    ],
    "output": [{"name": "z", "entry": "z"}],
}

SAMPLING_VARIABLES = [
    {"name": "x", "distribution": {"parameters": {"min": 0, "max": 1}}},
    {"name": "y", "distribution": {"parameters": {"min": 0, "max": 1}}},
    {"name": "shift", "distribution": {"parameters": {"min": 0, "max": 1}}},
]

LIKELIHOOD_EXPRESSIONS = [{"name": "LogL_Z", "expression": "z"}]


def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict[str, Any]]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0}


def _normalize_database_records(records: list[dict[str, Any]]) -> list[dict[str, float]]:
    keys = ("x", "y", "shift", "z", "LogL")
    normalized = []
    for row in records:
        normalized.append({key: float(row[key]) for key in keys})
    return sorted(normalized, key=lambda item: (item["x"], item["y"], item["shift"]))


def _spawn_build_worker_label(args: tuple[int, dict[str, Any], dict[str, Any]]) -> str:
    worker_id, redis_cfg, worker_cfg = args
    worker = Worker(worker_id, redis_cfg, worker_cfg)
    return f"{worker.worker_id}:{worker.worker_config['pull_timeout']}"


def _spawn_extract_redis_config(args: tuple[int, dict[str, Any], dict[str, Any]]) -> str:
    worker_id, redis_cfg, worker_cfg = args
    queue = RedisQueue(redis_cfg)
    worker = Worker(worker_id, queue, worker_cfg)
    return f"{worker.worker_id}:{worker.redis_config['port']}"


def _worker_config(tmpdir: str) -> dict[str, Any]:
    return {
        "sample_config": {
            "task_result_dir": tmpdir,
            "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
            "sample_artifacts": "auto",
            "workflow_has_calculator": False,
            "workflow_references_sdir": False,
        },
        "mapper": {
            "type": "flat",
            "variables": SAMPLING_VARIABLES,
        },
        "opera_modules": [BENCHMARK_OPERA_MODULE],
        "likelihood_expressions": LIKELIHOOD_EXPRESSIONS,
        "pull_timeout": 1,
    }


class WorkerMVPTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_start_workers_refuses_duplicate_start(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_worker_config(tmpdir))
                with self.assertRaises(RuntimeError):
                    factory.start_workers(1, **_worker_config(tmpdir))
                factory.shutdown()
        finally:
            server.shutdown()
            server.server_close()

    def test_factory_shutdown_closes_redis_connection(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            factory = TaskFactory.get_instance(redis_config)
            queue = factory.init_redis()
            self.assertIsNotNone(queue.r)
            factory.shutdown()
            self.assertIsNone(factory.redis)
            self.assertIsNone(queue.r)
        finally:
            server.shutdown()
            server.server_close()

    def test_worker_submits_failed_result_on_execution_error(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_worker_config(tmpdir))
                assert factory.redis is not None

                sample = Sample(
                    uuid="failed-sample",
                    u_coords=np.array([0.1, 0.2, 0.3], dtype=np.float64),
                    execution_plan=[
                        ExecutionStep(type="opera", name="UnknownOpera", layer=0),
                    ],
                )
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
                self.assertEqual(result["uuid"], "failed-sample")
                self.assertEqual(result["status"], "Failed")
        finally:
            server.shutdown()
            server.server_close()

    def test_task_factory_starts_worker_and_monitor_snapshot(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                factory = TaskFactory.get_instance(redis_config)
                factory.init_redis()
                factory.start_workers(1, **_worker_config(tmpdir))
                factory.start_monitor(update_hz=5.0)
                time.sleep(0.3)

                snapshot = factory.get_monitor_snapshot()
                self.assertIn("workers", snapshot)
                self.assertEqual(snapshot["workers_total"], 1)
                self.assertEqual(snapshot["workers_alive"], 1)
                self.assertEqual(len(snapshot["workers"]), 1)
                self.assertTrue(snapshot["workers"][0]["alive"])
                self.assertIn("task_queue_length", snapshot)
                self.assertIn("op_counts", snapshot)

                factory.shutdown()
        finally:
            server.shutdown()
            server.server_close()

    def test_core_skips_factory_when_mode_not_redis(self) -> None:
        core = Jarvis2Core({"Runtime": {"mode": "auto"}})
        self.assertFalse(core.is_redis_runtime())
        self.assertIsNone(core.init_factory())

    def test_worker_config_pickles_under_spawn_context(self) -> None:
        config = _worker_config(tempfile.mkdtemp())
        redis_config = {"host": "127.0.0.1", "port": 6379}
        payload = (0, redis_config, config)
        ctx = get_spawn_context()

        with ctx.Pool(1) as pool:
            label = pool.apply(_spawn_build_worker_label, (payload,))
        self.assertEqual(label, "0:1")

    def test_worker_accepts_live_redis_queue_for_config_extraction(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            queue = RedisQueue(redis_config)
            queue.connect()
            payload = (0, redis_config, _worker_config(tempfile.mkdtemp()))
            ctx = get_spawn_context()
            with ctx.Pool(1) as pool:
                label = pool.apply(_spawn_extract_redis_config, (payload,))
            self.assertEqual(label, f"0:{redis_config['port']}")
            queue.close()
        finally:
            server.shutdown()
            server.server_close()

    def test_end_to_end_opera_worker_archiver_database_parity(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(GOLDEN_PATH, "r", encoding="utf-8") as handle:
                    expected = json.load(handle)

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
                core.init_logger()
                core.init_redis()
                worker_config = _worker_config(tmpdir)
                core.init_factory(worker_config)
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template([BENCHMARK_OPERA_MODULE])
                core.set_sampler(sampler)

                samples = []
                for index, point in enumerate(expected):
                    sample = sampler._build_sample(
                        np.array([point["x"], point["y"], point["shift"]], dtype=np.float64)
                    )
                    sample.uuid = f"golden-sample-{index}"
                    samples.append(sample)
                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=30.0)
                core.shutdown()

                from jarvishep2.database import SimpleHDF5Writer

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected))
        finally:
            server.shutdown()
            server.server_close()

    def test_get_monitor_snapshot_does_not_touch_redis(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        factory._snapshot = {
            "timestamp": time.time(),
            "workers": [],
            "workers_alive": 0,
            "workers_total": 0,
            "task_queue_length": 3,
            "archive_queue_length": 1,
            "sample_stats": {"completed": 2},
            "calculator_status": {},
            "op_counts": {"task": 1, "worker": 0, "calculator": 0, "sample": 2},
        }

        def _forbidden_redis_call(*_args: Any, **_kwargs: Any) -> None:
            raise AssertionError("get_monitor_snapshot must not call Redis")

        assert queue.r is not None
        queue.r.get = _forbidden_redis_call  # type: ignore[method-assign]
        queue.r.hgetall = _forbidden_redis_call  # type: ignore[method-assign]
        queue.r.llen = _forbidden_redis_call  # type: ignore[method-assign]

        snapshot = factory.get_monitor_snapshot()
        self.assertEqual(snapshot["task_queue_length"], 3)
        self.assertEqual(snapshot["sample_stats"]["completed"], 2)

    def test_collect_latest_status_gates_hgetall_on_idle_ticks(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue
        assert queue.r is not None

        factory._snapshot = factory._collect_latest_status()
        hgetall = queue.r.hgetall
        calls = {"count": 0}

        def _counting_hgetall(*args: Any, **kwargs: Any) -> Any:
            calls["count"] += 1
            return hgetall(*args, **kwargs)

        queue.r.hgetall = _counting_hgetall  # type: ignore[method-assign]
        factory._collect_latest_status()
        self.assertEqual(calls["count"], 0)

    def test_collect_latest_status_refreshes_sample_stats_on_op_count_bump(self) -> None:
        factory = TaskFactory.get_instance()
        queue = make_fakeredis_queue()
        factory.redis = queue

        factory._snapshot = factory._collect_latest_status()
        queue.submit_result({"uuid": "sample-1", "status": "Completed"})
        updated = factory._collect_latest_status()

        self.assertEqual(updated["sample_stats"].get("completed"), 1)
        self.assertGreaterEqual(updated["op_counts"]["sample"], 1)

    def test_shutdown_stops_monitor_thread(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            factory = TaskFactory.get_instance(redis_config)
            factory.init_redis()
            factory.start_monitor(update_hz=20.0)
            time.sleep(0.1)
            self.assertIsNotNone(factory._updater_thread)
            assert factory._updater_thread is not None
            self.assertTrue(factory._updater_thread.is_alive())

            factory.shutdown()
            self.assertFalse(factory._running)
            self.assertIsNone(factory._updater_thread)
            self.assertIsNone(factory.redis)
        finally:
            server.shutdown()
            server.server_close()

    def test_graceful_sigterm_finishes_inflight_sample(self) -> None:
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
                worker_config = _worker_config(tmpdir)
                core.init_factory(worker_config)
                core.init_archiver(os.path.join(tmpdir, "DATABASE", "samples.hdf5"))

                assert core.redis is not None
                sample = Sample(
                    uuid="sigterm-sample",
                    u_coords=np.array([0.2, 0.3, 0.4], dtype=np.float64),
                    execution_plan=[
                        ExecutionStep(type="opera", name="TrivialEggbox", layer=0),
                        ExecutionStep(type="likelihood", name="LogLikelihood", layer=1),
                    ],
                )
                core.redis.push_task(sample.to_task_dict())
                time.sleep(0.5)
                assert core.factory is not None
                core.factory.request_worker_shutdown()
                core.factory.stop_all_workers(graceful=True, join_timeout=15.0)
                core.archiver.stop(drain=True)
                self.assertEqual(core.archiver.records_written, 1)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()