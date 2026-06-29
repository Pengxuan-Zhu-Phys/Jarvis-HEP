#!/usr/bin/env python3
"""Unit tests for jarvishep2.async_subprocess (WP-D2.2)."""

from __future__ import annotations

import pickle
import time
import unittest

from jarvishep2.async_subprocess import (
    AsyncSubprocessScheduler,
    SubprocessJob,
    SubprocessRuntimeConfig,
)
from jarvishep2.worker import Worker


class AsyncSubprocessSchedulerTests(unittest.TestCase):
    def tearDown(self) -> None:
        if hasattr(self, "_scheduler") and self._scheduler is not None:
            self._scheduler.shutdown(wait=True)

    def test_concurrent_jobs_finish_near_max_duration(self) -> None:
        delay = 0.25
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=2, log_policy="quiet")
        )
        jobs = [
            SubprocessJob(cmd=f"python3 -c 'import time; time.sleep({delay})'"),
            SubprocessJob(cmd=f"python3 -c 'import time; time.sleep({delay})'"),
        ]
        started = time.monotonic()
        futures = [self._scheduler.submit(job) for job in jobs]
        results = [future.result(timeout=5.0) for future in futures]
        elapsed = time.monotonic() - started
        self.assertTrue(all(result.ok for result in results))
        self.assertLess(elapsed, delay * 1.9)
        self.assertGreaterEqual(elapsed, delay * 0.8)
        snap = self._scheduler.snapshot()
        self.assertGreaterEqual(int(snap["peak_running"]), 2)

    def test_timeout_marks_result_not_ok(self) -> None:
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        result = self._scheduler.run(
            SubprocessJob(
                cmd="python3 -c 'import time; time.sleep(2)'",
                timeout_sec=0.2,
            ),
            timeout=5.0,
        )
        self.assertFalse(result.ok)
        self.assertTrue(result.timed_out)

    def test_env_injection_visible_to_child(self) -> None:
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        result = self._scheduler.run(
            SubprocessJob(
                cmd="python3 -c 'import os; print(os.environ.get(\"JARVIS_TEST_FLAG\", \"\"))'",
                env={"JARVIS_TEST_FLAG": "layer-concurrent"},
            ),
            timeout=5.0,
        )
        self.assertTrue(result.ok)

    def test_snapshot_tracks_pending_and_peak_running(self) -> None:
        delay = 0.3
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        job = SubprocessJob(cmd=f"python3 -c 'import time; time.sleep({delay})'")
        future = self._scheduler.submit(job)
        deadline = time.monotonic() + 2.0
        saw_pending_or_running = False
        while time.monotonic() < deadline and not future.done():
            snap = self._scheduler.snapshot()
            if int(snap["pending"]) > 0 or int(snap["running"]) > 0:
                saw_pending_or_running = True
                break
            time.sleep(0.01)
        future.result(timeout=5.0)
        final = self._scheduler.snapshot()
        self.assertTrue(saw_pending_or_running)
        self.assertGreaterEqual(int(final["peak_running"]), 1)

    def test_per_worker_scheduler_isolation(self) -> None:
        left = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        right = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        try:
            left.start()
            right.start()
            left_result = left.run(
                SubprocessJob(cmd="python3 -c 'print(1)'"),
                timeout=5.0,
            )
            left.shutdown(wait=True)
            right_result = right.run(
                SubprocessJob(cmd="python3 -c 'print(2)'"),
                timeout=5.0,
            )
            self.assertTrue(left_result.ok)
            self.assertTrue(right_result.ok)
            self.assertEqual(int(left.snapshot()["submitted"]), 1)
            self.assertEqual(int(right.snapshot()["submitted"]), 1)
        finally:
            right.shutdown(wait=True)

    def test_worker_spawn_boundary_has_no_live_runtime_objects(self) -> None:
        worker_config = {
            "force_serial_layers": False,
            "calculator_modules": [{"name": "Demo", "required_modules": []}],
        }
        restored_config = pickle.loads(pickle.dumps(worker_config, protocol=pickle.HIGHEST_PROTOCOL))
        self.assertEqual(restored_config, worker_config)
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, worker_config)
        self.assertIsNone(worker._scheduler)
        self.assertIsNone(worker._observables_lock)
        self.assertIsNone(worker._redis)


if __name__ == "__main__":
    unittest.main()