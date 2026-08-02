#!/usr/bin/env python3
"""WP-D6.2 distributed checkpoint/resume tests."""

from __future__ import annotations

import os
import pickle
import signal
import tempfile
import time
import unittest
from typing import Any

from jarvishep2.Sampling.runtime_checkpoint import (
    CHECKPOINT_FORMAT,
    CHECKPOINT_HEARTBEAT_SEC,
    REFUSAL_MESSAGES,
    RESUME_PROMPT,
    THROUGHPUT_CORE_FORMAT,
    V1_CHECKPOINT_FORMAT,
    derive_sample_seed,
    load_checkpoint,
    prepare_resume,
    safe_barrier_ready,
    validate_checkpoint_payload,
)
from jarvishep2.Sampling.seeded_sampler import SeededOperaSampler, deterministic_uuid
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer, read_persisted_uuids
from jarvishep2.factory import TaskFactory
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.process_cleanup import (
    kill_jarvis_processes,
    list_running_scans,
)
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
    persisted = read_persisted_uuids(os.path.dirname(db_path))
    sampler.set_persisted_uuids(persisted)

    try:
        if resume_payload is not None:
            sampler.repropose_unfinished()
        sampler.submit_all_remaining()

        expected_new = max(0, total_points - len(persisted))
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


def _run_crashable_control(
    redis_config: dict[str, Any],
    tmpdir: str,
    scan_name: str,
    resume: bool,
    outcome_queue: Any,
    total_points: int = 200,
) -> None:
    config = {
        "task_root": tmpdir,
        "task_result_dir": tmpdir,
        "scan_name": scan_name,
        "Scan": {"name": scan_name},
        "Runtime": {
            "mode": "redis",
            "workers": 2,
            "batch_size": 4,
            "max_inflight": 8,
            "redis": dict(redis_config, managed=False),
            "Watchdog": {"enabled": False},
            "Archiver": {"mode": "process", "batch_size": 4},
        },
        "Sampling": {
            "Method": "Random",
            "Point number": int(total_points),
            "Seed": 29,
            "Variables": [
                {
                    "name": "x",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                },
                {
                    "name": "y",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                },
            ],
            "LogLikelihood": [{"name": "LogL_Z", "expression": "z"}],
        },
        "Operas": {
            "Modules": [
                {
                    "name": "SlowResume",
                    "operator": "jarvishep2.testing.resume.slow_sum",
                    "call_mode": "call",
                    "input": [
                        {"name": "x", "expression": "x"},
                        {"name": "y", "expression": "y"},
                    ],
                    "output": [{"name": "z", "entry": "z"}],
                }
            ]
        },
    }
    core = Jarvis2Core(config)
    outcome = core.run(resume=resume, write_run_summary=False)
    outcome_queue.put(
        {
            "ok": outcome.ok,
            "submitted": outcome.submitted,
            "archived": outcome.archived,
            "exit_code": outcome.exit_code,
        }
    )
    TaskFactory.reset_instance()


class DistributedResumeTests(unittest.TestCase):
    @staticmethod
    def _database_records(database_dir: str) -> list[dict[str, Any]]:
        records: list[dict[str, Any]] = []
        for path in sorted(
            os.path.join(database_dir, name)
            for name in os.listdir(database_dir)
            if name.endswith(".hdf5")
        ):
            records.extend(SimpleHDF5Writer(path).read_records())
        return records

    def test_ctrl_c_resume_matches_uninterrupted_4000_point_run(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        ctx = get_spawn_context()
        baseline = interrupted = resumed = None
        baseline_queue = ctx.Queue()
        interrupted_queue = ctx.Queue()
        resumed_queue = ctx.Queue()
        point_count = 4000
        scan_prefix = f"resume-ctrl-c-{os.getpid()}"
        try:
            with (
                tempfile.TemporaryDirectory() as baseline_dir,
                tempfile.TemporaryDirectory() as resume_dir,
            ):
                baseline = ctx.Process(
                    target=_run_crashable_control,
                    args=(
                        redis_config,
                        baseline_dir,
                        f"{scan_prefix}-baseline",
                        False,
                        baseline_queue,
                        point_count,
                    ),
                )
                baseline.start()
                baseline.join(timeout=120.0)
                self.assertFalse(baseline.is_alive(), "uninterrupted baseline did not finish")
                self.assertEqual(baseline.exitcode, 0)
                self.assertTrue(baseline_queue.get(timeout=5.0)["ok"])
                baseline_records = self._database_records(
                    os.path.join(baseline_dir, "DATABASE")
                )

                interrupted = ctx.Process(
                    target=_run_crashable_control,
                    args=(
                        redis_config,
                        resume_dir,
                        f"{scan_prefix}-resume",
                        False,
                        interrupted_queue,
                        point_count,
                    ),
                )
                interrupted.start()
                database_dir = os.path.join(resume_dir, "DATABASE")
                deadline = time.monotonic() + 60.0
                durable_at_interrupt: set[str] = set()
                while time.monotonic() < deadline:
                    durable_at_interrupt = read_persisted_uuids(database_dir)
                    if 16 <= len(durable_at_interrupt) < point_count:
                        break
                    if not interrupted.is_alive():
                        self.fail("control process exited before Ctrl+C checkpoint")
                    time.sleep(0.05)
                self.assertGreaterEqual(len(durable_at_interrupt), 16)
                os.kill(interrupted.pid, signal.SIGINT)
                interrupted.join(timeout=60.0)
                self.assertFalse(interrupted.is_alive(), "Ctrl+C shutdown did not finish")

                checkpoint = os.path.join(
                    resume_dir,
                    "checkpoints",
                    f"{scan_prefix}-resume",
                    "Random",
                    "state.pkl",
                )
                self.assertTrue(os.path.isfile(checkpoint))
                durable_before_resume = read_persisted_uuids(database_dir)
                self.assertLess(len(durable_before_resume), point_count)

                resumed = ctx.Process(
                    target=_run_crashable_control,
                    args=(
                        redis_config,
                        resume_dir,
                        f"{scan_prefix}-resume",
                        True,
                        resumed_queue,
                        point_count,
                    ),
                )
                resumed.start()
                resumed.join(timeout=120.0)
                self.assertFalse(resumed.is_alive(), "Ctrl+C resume did not finish")
                self.assertEqual(resumed.exitcode, 0)
                outcome = resumed_queue.get(timeout=5.0)
                self.assertTrue(outcome["ok"], outcome)
                self.assertLessEqual(
                    int(outcome["submitted"]),
                    point_count - len(durable_at_interrupt) + 12,
                )

                resumed_records = self._database_records(database_dir)
                physics_keys = ("uuid", "x", "y", "z", "LogL", "LogL_Z", "status")
                baseline_by_uuid = {
                    str(row["uuid"]): {key: row.get(key) for key in physics_keys}
                    for row in baseline_records
                }
                resumed_by_uuid = {
                    str(row["uuid"]): {key: row.get(key) for key in physics_keys}
                    for row in resumed_records
                }
                self.assertEqual(len(baseline_records), point_count)
                self.assertEqual(len(resumed_records), point_count)
                self.assertEqual(len(resumed_by_uuid), point_count)
                self.assertEqual(resumed_by_uuid, baseline_by_uuid)
        finally:
            for proc in (baseline, interrupted, resumed):
                if proc is not None and proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=5.0)
            for scan in list_running_scans():
                if scan.name.startswith(scan_prefix):
                    kill_jarvis_processes(list(scan.processes), force=True)
            for queue in (baseline_queue, interrupted_queue, resumed_queue):
                queue.close()
                queue.join_thread()
            server.shutdown()
            server.server_close()

    def test_sigkill_control_then_database_only_resume_is_unique(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        scan_name = f"resume-e2e-{os.getpid()}"
        ctx = get_spawn_context()
        outcome_queue = ctx.Queue()
        first = None
        second = None
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                first = ctx.Process(
                    target=_run_crashable_control,
                    args=(redis_config, tmpdir, scan_name, False, outcome_queue),
                )
                first.start()
                db_dir = os.path.join(tmpdir, "DATABASE")
                deadline = time.monotonic() + 45.0
                durable_before_kill: set[str] = set()
                while time.monotonic() < deadline:
                    try:
                        durable_before_kill = read_persisted_uuids(db_dir)
                    except OSError:
                        durable_before_kill = set()
                    if 1 <= len(durable_before_kill) < 200:
                        break
                    if not first.is_alive():
                        self.fail("control process exited before SIGKILL checkpoint")
                    time.sleep(0.1)
                self.assertTrue(durable_before_kill)
                os.kill(first.pid, signal.SIGKILL)
                first.join(timeout=10.0)
                self.assertFalse(first.is_alive())

                # No heartbeat checkpoint is expected this early.  Resume must
                # reconcile from DATABASE and clean the exact-name orphan fleet.
                checkpoint = os.path.join(
                    tmpdir, "checkpoints", scan_name, "Random", "state.pkl"
                )
                self.assertFalse(os.path.exists(checkpoint))
                second = ctx.Process(
                    target=_run_crashable_control,
                    args=(redis_config, tmpdir, scan_name, True, outcome_queue),
                )
                second.start()
                second.join(timeout=120.0)
                self.assertFalse(second.is_alive(), "resume control did not finish")
                self.assertEqual(second.exitcode, 0)
                outcome = outcome_queue.get(timeout=5.0)
                self.assertTrue(outcome["ok"], outcome)

                records = self._database_records(db_dir)
                uuids = [str(row.get("uuid") or "") for row in records]
                self.assertEqual(len(records), 200)
                self.assertEqual(len(set(uuids)), len(records))
                self.assertLessEqual(
                    int(outcome["submitted"]),
                    200 - len(durable_before_kill) + 12,
                )
        finally:
            for proc in (first, second):
                if proc is not None and proc.is_alive():
                    proc.terminate()
                    proc.join(timeout=5.0)
            scans = [scan for scan in list_running_scans() if scan.name == scan_name]
            for scan in scans:
                kill_jarvis_processes(list(scan.processes), force=True)
            outcome_queue.close()
            outcome_queue.join_thread()
            server.shutdown()
            server.server_close()

    def test_database_only_resume_without_checkpoint_skips_durable_uuids(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sampler = SeededOperaSampler(seed=11, total_points=6)
            sampler.set_config({"Runtime": {"mode": "redis"}})
            sampler.set_execution_plan_template(include_likelihood=False)
            queue = make_fakeredis_queue()
            queue.connect()
            sampler.set_redis(queue)
            durable = {
                deterministic_uuid(master=sampler._master_seq, sample_index=index)
                for index in range(3)
            }
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            SimpleHDF5Writer(db_path).add_records(
                [{"uuid": uuid, "status": "Completed"} for uuid in sorted(durable)]
            )

            sampler.set_persisted_uuids(read_persisted_uuids(os.path.dirname(db_path)))
            submitted = sampler.submit_all_remaining()

            self.assertEqual(len(submitted), 3)
            self.assertTrue(set(submitted).isdisjoint(durable))
            self.assertEqual(int(queue.r.llen("hep:task_queue")), 3)

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

    def test_prepare_resume_seeds_failed_counts_from_database(self) -> None:
        """D21.13: resume must not re-label durable Failed rows as completed."""
        queue = make_fakeredis_queue()
        queue.connect()
        prepare_resume(
            queue,
            worker_config={},
            persisted_count=40,
            persisted_failed=20,
        )
        stats = queue.fetch_sample_stats()
        self.assertEqual(int(stats.get("completed") or 0), 20)
        self.assertEqual(int(stats.get("failed") or 0), 20)
        self.assertEqual(int(stats.get("running") or 0), 0)

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
                    checkpoint = core.save_runtime_checkpoint(reason="test_interrupt")
                    assert checkpoint is not None
                    payload = load_checkpoint(checkpoint)
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
