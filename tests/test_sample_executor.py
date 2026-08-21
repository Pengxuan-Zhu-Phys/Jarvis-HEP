#!/usr/bin/env python3
"""D25.5: Worker is a Redis/FileOperation process; physics lives on SampleExecutor."""

from __future__ import annotations

import pickle
import unittest

from jarvishep2.sample_executor import SampleExecutor, stamp_task_sample_bucket
from jarvishep2.worker import Worker


class SampleExecutorSplitTests(unittest.TestCase):
    def test_worker_lazily_holds_executor(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        self.assertFalse(hasattr(worker, "_executor"))
        executor = worker._get_executor()
        self.assertIsInstance(executor, SampleExecutor)
        self.assertIs(executor._worker, worker)
        self.assertIs(worker._get_executor(), executor)

    def test_spawn_boundary_does_not_pickle_live_runtime(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        self.assertIsNone(worker._scheduler)
        self.assertIsNone(worker._observables_lock)
        self.assertIsNone(worker._redis)
        self.assertFalse(hasattr(worker, "_executor"))
        restored = pickle.loads(pickle.dumps(worker.redis_config))
        self.assertEqual(restored["host"], "127.0.0.1")

    def test_stamp_helper_skips_when_disabled_or_already_allocated(self) -> None:
        payload = {"uuid": "u1", "bucket_dir": "/tmp/SAMPLE/000001"}
        self.assertIs(stamp_task_sample_bucket(object(), payload, enabled=True), payload)
        self.assertEqual(
            stamp_task_sample_bucket(object(), {"uuid": "u1"}, enabled=False),
            {"uuid": "u1"},
        )


if __name__ == "__main__":
    unittest.main()
