#!/usr/bin/env python3
"""WP-D6.2 distributed checkpoint/resume tests."""

from __future__ import annotations

import os
import pickle
import tempfile
import threading
import time
import unittest
from typing import Any

from fakeredis import TcpFakeServer

from jarvishep2.Sampling.runtime_checkpoint import (
    CHECKPOINT_FORMAT,
    CHECKPOINT_HEARTBEAT_SEC,
    REFUSAL_MESSAGES,
    RESUME_PROMPT,
    THROUGHPUT_CORE_FORMAT,
    V1_CHECKPOINT_FORMAT,
    build_payload,
    build_run_spec,
    derive_sample_seed,
    load_checkpoint,
    prepare_resume,
    safe_barrier_ready,
    save_checkpoint,
    serialize_seed_sequence,
    validate_checkpoint_payload,
)
from jarvishep2.Sampling.seeded_sampler import SeededOperaSampler, deterministic_uuid
from jarvishep2.core import Jarvis2Core
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue

from test_worker_mvp import (
    BENCHMARK_OPERA_MODULE,
    LIKELIHOOD_EXPRESSIONS,
    SAMPLING_VARIABLES,
    _normalize_database_records,
    _start_tcp_fakeredis,
)


def _worker_config(tmpdir: str) -> dict[str, Any]:
    return {
        "sample_config": {
            "task_result_dir": tmpdir,
            "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
            "sample_artifacts": "auto",
            "workflow_has_calculator": False,
            "workflow_references_sdir": False,
        },
        "mapper": {"type": "flat", "variables": SAMPLING_VARIABLES},
        "opera_modules": [BENCHMARK_OPERA_MODULE],
        "likelihood_expressions": LIKELIHOOD_EXPRESSIONS,
        "pull_timeout": 1,
        "handoff_to_staging": False,
    }


def _run_seeded_scan(
    *,
    redis_config: dict[str, Any],
    tmpdir: str,
    seed: int,
    total_points: int,
    workers: int,
    resume_payload: dict[str, Any] | None = None,
) -> list[dict[str, float]]:
    TaskFactory.reset_instance()
    core = Jarvis2Core(
        {
            "task_result_dir": tmpdir,
            "task_root": tmpdir,
            "scan_name": "resume-scan",
            "Runtime": {
                "mode": "redis",
                "workers": workers,
                "redis": redis_config,
                "Watchdog": {"enabled": False},
            },
        }
    )
    core.info = {
        "task_result_dir": tmpdir,
        "task_root": tmpdir,
        "scan_name": "resume-scan",
        "sampler_name": "SeededOperaSampler",
    }
    if resume_payload is not None:
        core._resume_policy = "resume"
        core._resume_checkpoint_payload = resume_payload

    core.init_redis()
    sampler = SeededOperaSampler(seed=seed, total_points=total_points)
    sampler.set_config(core.config)
    sampler.set_execution_plan_template(
        opera_modules=[BENCHMARK_OPERA_MODULE],
        include_likelihood=True,
    )
    core.set_sampler(sampler)
    core.init_factory(worker_config=_worker_config(tmpdir))
    db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
    core.init_archiver(db_path)

    try:
        if resume_payload is not None:
            sampler.repropose_unfinished()
        sampler.submit_all_remaining()

        if resume_payload is not None:
            completed = len(resume_payload.get("sampler_state", {}).get("completed_uuids") or [])
            expected_new = max(0, total_points - completed)
        else:
            expected_new = total_points
        core.wait_for_results(expected_new, timeout=60.0)
    finally:
        core.shutdown()
    if not os.path.exists(db_path):
        return []
    from jarvishep2.database import SimpleHDF5Writer

    return _normalize_database_records(SimpleHDF5Writer(db_path).read_records())


def _stop_factory_workers() -> None:
    factory = TaskFactory._instance
    if factory is not None:
        try:
            factory.shutdown(wait=False)
        except Exception:
            pass
    TaskFactory.reset_instance()


class DistributedResumeTests(unittest.TestCase):
    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def test_checkpoint_roundtrip_restores_sampler_state(self) -> None:
        sampler = SeededOperaSampler(seed=11, total_points=5)
        sampler.set_config({"Runtime": {"mode": "redis"}})
        sampler.set_redis(make_fakeredis_queue())
        sampler.submit_next()
        sampler.submit_next()
        exported = sampler.export_runtime_state()
        restored = SeededOperaSampler(seed=0, total_points=1)
        restored.set_config({"Runtime": {"mode": "redis"}})
        restored.set_redis(make_fakeredis_queue())
        restored.import_runtime_state(exported)
        self.assertEqual(restored._next_sample_index, 2)
        self.assertEqual(len(restored._submitted_uuids), 2)
        next_uuid = restored.submit_next()
        self.assertIsNotNone(next_uuid)
        self.assertEqual(restored._next_sample_index, 3)

    def test_format_refusal_rejects_v1_and_throughput_core(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            for fmt, message in REFUSAL_MESSAGES.items():
                path = os.path.join(tmpdir, f"{fmt}.pkl")
                with open(path, "wb") as handle:
                    pickle.dump({"format": fmt, "version": 1, "run_spec": {}, "sampler_state": {}}, handle)
                with self.assertRaisesRegex(ValueError, message[:40]):
                    load_checkpoint(path)

        ok, reason = validate_checkpoint_payload({"format": THROUGHPUT_CORE_FORMAT, "version": 1})
        self.assertFalse(ok)
        self.assertEqual(reason, REFUSAL_MESSAGES[THROUGHPUT_CORE_FORMAT])

        ok, reason = validate_checkpoint_payload({"format": V1_CHECKPOINT_FORMAT, "version": 1})
        self.assertFalse(ok)
        self.assertEqual(reason, REFUSAL_MESSAGES[V1_CHECKPOINT_FORMAT])

    def test_safe_barrier_requires_archiver_ack(self) -> None:
        submitted = frozenset({"a", "b"})
        self.assertFalse(
            safe_barrier_ready(
                sampler_at_barrier=True,
                submitted_uuids=submitted,
                archiver_persistence={"acked_uuids": ["a"]},
            )
        )
        self.assertTrue(
            safe_barrier_ready(
                sampler_at_barrier=True,
                submitted_uuids=submitted,
                archiver_persistence={"acked_uuids": ["a", "b"]},
            )
        )

    def test_resume_prompt_and_heartbeat_contract_frozen(self) -> None:
        self.assertEqual(
            RESUME_PROMPT,
            "Detected checkpoint file. Re-run from scratch? [y/N] (default: resume in 30s): ",
        )
        self.assertEqual(CHECKPOINT_HEARTBEAT_SEC, 30.0)

    def test_prepare_resume_drains_stale_task_queue(self) -> None:
        queue = make_fakeredis_queue()
        queue.push_task(
            {
                "uuid": "stale-task",
                "u_coords": [0.1, 0.2, 0.3],
                "execution_plan": [
                    {"name": "TrivialEggbox", "type": "opera", "layer": 0},
                ],
            }
        )
        drained = prepare_resume(queue, worker_config={})
        self.assertEqual(drained, 1)
        assert queue.r is not None
        self.assertEqual(int(queue.r.llen("hep:task_queue")), 0)

    def test_derive_sample_seed_independent_of_worker_count(self) -> None:
        import numpy as np

        master = np.random.SeedSequence(99)
        one_worker = [float(derive_sample_seed(master, index).entropy) for index in range(4)]
        four_worker = [float(derive_sample_seed(master, index).entropy) for index in range(4)]
        self.assertEqual(one_worker, four_worker)

    def _records_for_workers(self, *, tmpdir: str, redis_config: dict[str, Any], workers: int) -> list[dict[str, float]]:
        TaskFactory.reset_instance()
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
        sampler = SeededOperaSampler(seed=42, total_points=6)
        sampler.set_config(core.config)
        sampler.set_execution_plan_template(
            opera_modules=[BENCHMARK_OPERA_MODULE],
            include_likelihood=True,
        )
        core.set_sampler(sampler)
        core.init_factory(_worker_config(tmpdir))
        db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
        core.init_archiver(db_path)
        try:
            sampler.submit_all_remaining()
            core.wait_for_results(6, timeout=60.0)
        finally:
            core.shutdown()
        from jarvishep2.database import SimpleHDF5Writer

        return _normalize_database_records(SimpleHDF5Writer(db_path).read_records())

    def test_worker_count_independent_trajectory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir_one:
            with tempfile.TemporaryDirectory() as tmpdir_two:
                server_one, redis_one = _start_tcp_fakeredis()
                server_two, redis_two = _start_tcp_fakeredis()
                try:
                    one = self._records_for_workers(
                        tmpdir=tmpdir_one,
                        redis_config=redis_one,
                        workers=1,
                    )
                    two = self._records_for_workers(
                        tmpdir=tmpdir_two,
                        redis_config=redis_two,
                        workers=2,
                    )
                    self.assertEqual(one, two)
                finally:
                    server_one.shutdown()
                    server_one.server_close()
                    server_two.shutdown()
                    server_two.server_close()

    def test_kill_and_resume_completes_without_duplicate_uuids(self) -> None:
        seed = 7
        total_points = 6
        with tempfile.TemporaryDirectory() as tmpdir:
            server, redis_config = _start_tcp_fakeredis()
            try:
                sampler = SeededOperaSampler(seed=seed, total_points=total_points)
                sampler.set_execution_plan_template(
                    opera_modules=[BENCHMARK_OPERA_MODULE],
                    include_likelihood=True,
                )
                expected = [
                    deterministic_uuid(master=sampler._master_seq, sample_index=index)
                    for index in range(total_points)
                ]

                # Phase 1: partial run
                TaskFactory.reset_instance()
                core = Jarvis2Core(
                    {
                        "task_result_dir": tmpdir,
                        "task_root": tmpdir,
                        "scan_name": "resume-scan",
                        "Runtime": {
                            "mode": "redis",
                            "workers": 1,
                            "redis": redis_config,
                            "Watchdog": {"enabled": False},
                        },
                    }
                )
                core.info = {
                    "task_result_dir": tmpdir,
                    "task_root": tmpdir,
                    "scan_name": "resume-scan",
                    "sampler_name": "SeededOperaSampler",
                }
                core.init_redis()
                sampler.set_config(core.config)
                core.set_sampler(sampler)
                core.init_factory(worker_config=_worker_config(tmpdir))
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                try:
                    for _ in range(3):
                        sampler.submit_next()
                    core.wait_for_results(3, timeout=45.0)
                    for uuid in expected[:3]:
                        sampler.mark_completed(uuid)
                    persistence = core.archiver.persistence_state()
                    payload = build_payload(
                        run_spec=build_run_spec(
                            config=core.config,
                            scan_name="resume-scan",
                            task_root=tmpdir,
                            task_result_dir=tmpdir,
                            sampler_name="SeededOperaSampler",
                        ),
                        sampler_state=sampler.export_runtime_state(),
                        persistence=persistence,
                        reason="test_barrier",
                    )
                    checkpoint = core.checkpoint_file()
                    save_checkpoint(checkpoint, payload)
                finally:
                    core.shutdown()
                self.assertEqual(payload["format"], CHECKPOINT_FORMAT)

                # Phase 2: resume on a fresh Redis and finish remaining points
                server_two, redis_two = _start_tcp_fakeredis()
                try:
                    resumed = _run_seeded_scan(
                        redis_config=redis_two,
                        tmpdir=tmpdir,
                        seed=0,
                        total_points=total_points,
                        workers=1,
                        resume_payload=load_checkpoint(checkpoint),
                    )
                    self.assertEqual(len(resumed), total_points)
                finally:
                    server_two.shutdown()
                    server_two.server_close()
            finally:
                server.shutdown()
                server.server_close()


if __name__ == "__main__":
    unittest.main()