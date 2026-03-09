#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.async_subprocess import (  # noqa: E402
    AsyncSubprocessScheduler,
    SubprocessJob,
    SubprocessRuntimeConfig,
)


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _CaptureLogger:
    def __init__(self):
        self.info_lines = []
        self.warning_lines = []

    def info(self, message, *_args, **_kwargs):
        self.info_lines.append(str(message))

    def warning(self, message, *_args, **_kwargs):
        self.warning_lines.append(str(message))

    def error(self, *_args, **_kwargs):
        return None


class TestAsyncSubprocessScheduler(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="jarvis-async-subprocess-")
        self.tmpdir = Path(self._tmp.name)
        self.scheduler = None

    def tearDown(self):
        if self.scheduler is not None:
            self.scheduler.shutdown(wait=True, timeout=20.0)
            self.scheduler = None
        self._tmp.cleanup()

    def _new_scheduler(self, **overrides):
        cfg = SubprocessRuntimeConfig(
            max_concurrency=overrides.get("max_concurrency", 2),
            max_pending=overrides.get("max_pending", 32),
            queue_put_timeout_sec=overrides.get("queue_put_timeout_sec", 5.0),
            per_task_timeout_sec=overrides.get("per_task_timeout_sec", None),
            terminate_grace_sec=overrides.get("terminate_grace_sec", 0.2),
            log_policy=overrides.get("log_policy", "file"),
            progress_interval_sec=overrides.get("progress_interval_sec", 999.0),
            diagnostics_enabled=overrides.get("diagnostics_enabled", False),
            diagnostics_interval_sec=overrides.get("diagnostics_interval_sec", 0.1),
        )
        self.scheduler = AsyncSubprocessScheduler(
            config=cfg,
            logger=_NoopLogger(),
            status_path=str(self.tmpdir / "status.jsonl"),
        )
        return self.scheduler

    def test_concurrency_cap_enforced(self):
        scheduler = self._new_scheduler(max_concurrency=3, max_pending=8)

        futures = []
        for idx in range(24):
            task_dir = self.tmpdir / "tasks" / f"cap-{idx:03}"
            job = SubprocessJob(
                cmd=[sys.executable, "-c", "import time; time.sleep(0.06)"],
                shell=False,
                log_dir=str(task_dir),
                task_id=f"cap-{idx:03}",
            )
            futures.append(scheduler.submit(job))

        results = [f.result(timeout=15.0) for f in futures]
        self.assertTrue(all(r.ok for r in results))

        snap = scheduler.snapshot()
        self.assertLessEqual(int(snap.get("peak_running", 0)), 3)
        self.assertEqual(int(snap.get("completed", 0)), len(results))

    def test_timeout_terminate_with_fallback(self):
        scheduler = self._new_scheduler(
            max_concurrency=1,
            per_task_timeout_sec=0.2,
            terminate_grace_sec=0.1,
            max_pending=2,
        )

        result = scheduler.run(
            SubprocessJob(
                cmd=[sys.executable, "-c", "import time; time.sleep(2.0)"],
                shell=False,
                task_id="timeout-001",
                log_dir=str(self.tmpdir / "tasks" / "timeout-001"),
            ),
            timeout=10.0,
        )
        self.assertTrue(result.timed_out)
        self.assertNotEqual(int(result.returncode), 0)

    def test_large_stdout_stderr_streaming(self):
        scheduler = self._new_scheduler(max_concurrency=2, max_pending=4)

        payload = (
            "import sys;"
            "chunk='x'*65536;"
            "[sys.stdout.write(chunk) for _ in range(24)];"
            "sys.stdout.flush();"
            "[sys.stderr.write(chunk) for _ in range(24)];"
            "sys.stderr.flush()"
        )

        jobs = []
        for idx in range(8):
            jobs.append(
                SubprocessJob(
                    cmd=[sys.executable, "-c", payload],
                    shell=False,
                    task_id=f"large-{idx:03}",
                    log_dir=str(self.tmpdir / "tasks" / f"large-{idx:03}"),
                )
            )

        futures = [scheduler.submit(job) for job in jobs]
        results = [f.result(timeout=20.0) for f in futures]
        self.assertTrue(all(r.ok for r in results))

        for result in results:
            self.assertGreater(result.stdout_bytes, 0)
            self.assertGreater(result.stderr_bytes, 0)
            self.assertTrue(Path(result.stdout_path).exists())
            self.assertTrue(Path(result.stderr_path).exists())

    def test_diagnostics_writes_shutdown_final_snapshot(self):
        scheduler = self._new_scheduler(
            max_concurrency=2,
            max_pending=8,
            diagnostics_enabled=True,
            diagnostics_interval_sec=999.0,
        )

        futures = []
        for idx in range(12):
            futures.append(
                scheduler.submit(
                    SubprocessJob(
                        cmd=[sys.executable, "-c", "import sys; sys.stdout.write('ok\\n')"],
                        shell=False,
                        task_id=f"diag-{idx:03}",
                        log_dir=str(self.tmpdir / "tasks" / f"diag-{idx:03}"),
                    )
                )
            )
        results = [f.result(timeout=20.0) for f in futures]
        self.assertTrue(all(r.ok for r in results))

        scheduler.shutdown(wait=True, timeout=20.0)
        self.scheduler = None

        status_path = self.tmpdir / "status.jsonl"
        self.assertTrue(status_path.exists())
        lines = [json.loads(line) for line in status_path.read_text().splitlines() if line.strip()]
        self.assertGreaterEqual(len(lines), 1)

        final = lines[-1]
        self.assertEqual(final.get("reason"), "shutdown-final")
        self.assertEqual(int(final.get("submitted", -1)), len(results))
        self.assertEqual(int(final.get("completed", -1)), len(results))
        self.assertEqual(int(final.get("failed", -1)), 0)
        self.assertEqual(int(final.get("running", -1)), 0)

    def test_logger_policy_streams_output_without_file_artifacts(self):
        cfg = SubprocessRuntimeConfig(
            max_concurrency=1,
            max_pending=4,
            queue_put_timeout_sec=5.0,
            per_task_timeout_sec=None,
            terminate_grace_sec=0.2,
            log_policy="logger",
            progress_interval_sec=999.0,
            diagnostics_enabled=False,
            diagnostics_interval_sec=1.0,
        )
        caplog = _CaptureLogger()
        self.scheduler = AsyncSubprocessScheduler(
            config=cfg,
            logger=caplog,
            status_path=str(self.tmpdir / "status.jsonl"),
        )

        task_dir = self.tmpdir / "tasks" / "logger-001"
        result = self.scheduler.run(
            SubprocessJob(
                cmd=[sys.executable, "-c", "import sys;print('hello');print('err', file=sys.stderr)"],
                shell=False,
                log_dir=str(task_dir),
                task_id="logger-001",
            ),
            timeout=10.0,
        )
        self.assertTrue(result.ok)
        self.assertIsNone(result.stdout_path)
        self.assertIsNone(result.stderr_path)
        self.assertFalse((task_dir / "stdout.log").exists())
        self.assertFalse((task_dir / "stderr.log").exists())
        self.assertFalse((task_dir / "meta.json").exists())
        joined = "\n".join(caplog.info_lines)
        self.assertIn("[Subprocess][logger-001][stdout] hello", joined)
        self.assertIn("[Subprocess][logger-001][stderr] err", joined)


if __name__ == "__main__":
    unittest.main()
