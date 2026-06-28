#!/usr/bin/env python3
"""Unit tests for jarvishep2.async_subprocess (WP-D2.2)."""

from __future__ import annotations

import time
import unittest

from jarvishep2.async_subprocess import (
    AsyncSubprocessScheduler,
    SubprocessJob,
    SubprocessRuntimeConfig,
)


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


if __name__ == "__main__":
    unittest.main()