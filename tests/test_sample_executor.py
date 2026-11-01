#!/usr/bin/env python3
"""SampleExecutor ownership, bucket allocation and calculator release."""

from __future__ import annotations

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

    def test_stamp_helper_skips_when_disabled_or_already_allocated(self) -> None:
        payload = {"uuid": "u1", "bucket_dir": "/tmp/SAMPLE/000001"}
        self.assertIs(stamp_task_sample_bucket(object(), payload, enabled=True), payload)
        self.assertEqual(
            stamp_task_sample_bucket(object(), {"uuid": "u1"}, enabled=False),
            {"uuid": "u1"},
        )


class _FailingRedis:
    def force_release_calc(self, _name: str, _pack: str) -> bool:
        raise ConnectionError("redis unavailable")


class _SuccessfulRedis:
    def force_release_calc(self, _name: str, _pack: str) -> bool:
        return True


class WorkerReleaseTests(unittest.TestCase):
    def test_release_error_keeps_pack_visible_to_watchdog(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        worker._redis = _FailingRedis()
        worker._held_calc_packs["DemoCalc"] = "001"
        worker._force_release_pack("DemoCalc", "001")
        self.assertEqual(worker._held_calc_packs, {"DemoCalc": "001"})

    def test_release_success_removes_pack_after_confirmation(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        worker._redis = _SuccessfulRedis()
        worker._held_calc_packs["DemoCalc"] = "001"
        worker._force_release_pack("DemoCalc", "001")
        self.assertEqual(worker._held_calc_packs, {})


if __name__ == "__main__":
    unittest.main()
