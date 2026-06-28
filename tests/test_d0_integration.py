#!/usr/bin/env python3
"""D0 integration tests: Sample ↔ RedisQueue end-to-end (fakeredis)."""

from __future__ import annotations

import multiprocessing as mp
import os
import pickle
import tempfile
import unittest
from uuid import uuid4

import numpy as np

from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import ExecutionStep, Sample, materialize_failure_artifacts


def _pickle_round_trip(obj: object) -> object:
    return pickle.loads(pickle.dumps(obj))


def _rebuild_sample_from_task_dict(task: dict) -> str:
    sample = Sample.from_task_dict(task)

    class _Mapper:
        def map(self, u_coords: np.ndarray) -> dict[str, float]:
            return {"x": float(u_coords[0]), "y": float(u_coords[1])}

    sample.bind_params(_Mapper())
    return f"{sample.uuid}:{sample.params['x']}"


class D0IntegrationTests(unittest.TestCase):
    def test_sample_task_dict_redis_round_trip_and_bind_params(self):
        queue = make_fakeredis_queue()
        original = Sample(
            uuid=str(uuid4()),
            u_coords=np.array([0.2, 0.8], dtype=np.float64),
            execution_plan=[
                ExecutionStep(type="calculator", name="DemoCalc", layer=0),
                ExecutionStep(type="opera", name="LogL", layer=1),
            ],
            opera_params={"shift": 0.1},
        )

        queue.push_task(original.to_task_dict())
        wire = queue.pull_task(timeout=1)
        self.assertIsNotNone(wire)

        rebuilt = Sample.from_task_dict(wire)
        self.assertFalse(rebuilt._materialized)
        self.assertIsNone(rebuilt._logger)

        class _Mapper:
            def map(self, u_coords: np.ndarray) -> dict[str, float]:
                return {"x": float(u_coords[0]), "y": float(u_coords[1])}

        rebuilt.bind_params(_Mapper())
        self.assertEqual(rebuilt.params["x"], 0.2)
        self.assertEqual(rebuilt.observables["uuid"], rebuilt.uuid)

    def test_worker_path_materialize_and_submit_result(self):
        queue = make_fakeredis_queue()
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample(
                uuid=str(uuid4()),
                u_coords=np.array([0.5]),
                execution_plan=[ExecutionStep(type="opera", name="LogL", layer=0)],
            )
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            queue.push_task(sample.to_task_dict())
            wire = queue.pull_task(timeout=1)
            worker = Sample.from_task_dict(wire)
            worker.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            worker.start()
            worker.observables.update({"x": 1.0, "LogL": -2.0})
            worker.status = "Completed"
            worker.info["elapsed_s"] = 0.42

            result = worker.to_info_dict()
            queue.submit_result(result)
            archived = queue.pull_result(timeout=1)
            self.assertEqual(archived["uuid"], worker.uuid)
            self.assertEqual(archived["status"], "Completed")
            self.assertEqual(archived["timings"]["elapsed_s"], 0.42)
            self.assertFalse(os.path.exists(os.path.join(tmpdir, worker.uuid)))

    def test_failure_replay_after_redis_round_trip(self):
        queue = make_fakeredis_queue()
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample.from_params({"x": 1.0})
            sample.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            queue.push_task(sample.to_task_dict())
            wire = queue.pull_task(timeout=1)
            worker = Sample.from_task_dict(wire)
            worker.set_config(
                {
                    "sample_dirs": tmpdir,
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "auto",
                    "workflow_has_calculator": False,
                    "workflow_references_sdir": False,
                }
            )
            worker.start()
            worker.info["logger"].info("workflow failed inside worker")
            worker.status = "Failed"
            save_dir = materialize_failure_artifacts(worker.info, error="boom")
            self.assertIsNotNone(save_dir)
            log_path = os.path.join(save_dir, "Sample_running.log")
            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("workflow failed inside worker", text)
            self.assertIn("Sample failed -> boom", text)

            queue.submit_result(worker.to_info_dict())
            archived = queue.pull_result(timeout=1)
            self.assertEqual(archived["status"], "Failed")


class SpawnPicklingTests(unittest.TestCase):
    def test_execution_step_pickles_under_spawn_context(self):
        step = ExecutionStep(type="opera", name="LogL", layer=1, params={"k": 1})
        ctx = mp.get_context("spawn")
        with ctx.Pool(1) as pool:
            restored = pool.apply(_pickle_round_trip, (step,))
        self.assertIsInstance(restored, ExecutionStep)
        self.assertEqual(restored.name, "LogL")

    def test_sample_pickles_under_spawn_context(self):
        sample = Sample(
            uuid="spawn-sample",
            u_coords=np.array([0.1, 0.9], dtype=np.float64),
            execution_plan=[ExecutionStep(type="calculator", name="DemoCalc", layer=0)],
        )
        ctx = mp.get_context("spawn")
        with ctx.Pool(1) as pool:
            restored = pool.apply(_pickle_round_trip, (sample,))
        self.assertIsInstance(restored, Sample)
        self.assertTrue(np.allclose(restored.u_coords, sample.u_coords))

    def test_task_dict_rebuild_pickles_under_spawn_context(self):
        task = Sample(
            uuid="spawn-task",
            u_coords=np.array([0.3, 0.7]),
            execution_plan=[ExecutionStep(type="opera", name="LogL", layer=0)],
        ).to_task_dict()
        ctx = mp.get_context("spawn")
        with ctx.Pool(1) as pool:
            label = pool.apply(_rebuild_sample_from_task_dict, (task,))
        self.assertEqual(label, "spawn-task:0.3")


if __name__ == "__main__":
    unittest.main()