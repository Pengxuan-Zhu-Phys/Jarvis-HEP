#!/usr/bin/env python3
"""Control-plane completion progress + duration helpers."""

from __future__ import annotations

import threading
import time
import unittest
from unittest import mock

from jarvishep2.archiver import ArchiverProcess, SimpleArchiver
from jarvishep2.core import Jarvis2Core
from jarvishep2.log_kv import PermilleProgress, format_duration
from jarvishep2.mp_context import get_spawn_context


class FormatDurationTests(unittest.TestCase):
    def test_shape(self) -> None:
        self.assertRegex(format_duration(0.192), r"^\d{2}:\d{2}:\d{2}\.\d{3}$")
        self.assertTrue(format_duration(0.0).startswith("00:00:00"))


class PermilleProgressTests(unittest.TestCase):
    def test_logs_each_new_permille_and_warns_on_percent(self) -> None:
        info: list[str] = []
        warn: list[str] = []

        class _Log:
            def info(self, msg, *a, **k):  # noqa: ANN001
                info.append(str(msg) % a if a else str(msg))

            def warning(self, msg, *a, **k):  # noqa: ANN001
                warn.append(str(msg) % a if a else str(msg))

        progress = PermilleProgress(_Log(), total=1000, label="samples finished", t0=time.time())
        progress.update(0, force=True)
        progress.update(1)  # 1‰ → info
        progress.update(1)  # no-op same ‰
        progress.update(10)  # 10‰ → warning (exact 1%)

        self.assertTrue(any("0‰ of 0/1000 samples finished" in line for line in info))
        self.assertTrue(any("1‰ of 1/1000" in line for line in info))
        self.assertTrue(any(line.startswith("10‰") for line in warn))


class WaitForResultsProgressTests(unittest.TestCase):
    def test_wait_for_results_emits_completion_heartbeats(self) -> None:
        core = Jarvis2Core()
        debug_messages: list[str] = []
        info_messages: list[str] = []

        class _Log:
            def debug(self, msg, *a, **k):  # noqa: ANN001
                debug_messages.append(str(msg) % a if a else str(msg))

            def info(self, msg, *a, **k):  # noqa: ANN001
                info_messages.append(str(msg) % a if a else str(msg))

            def warning(self, msg, *a, **k):  # noqa: ANN001
                info_messages.append(str(msg) % a if a else str(msg))

        core._logger = _Log()  # type: ignore[assignment]
        # Fake archiver that reaches target after a few polls.
        state = {"n": 0}

        class _Archiver:
            @property
            def records_written(self):
                state["n"] += 1
                return min(state["n"], 5)

        core.archiver = _Archiver()  # type: ignore[assignment]
        core.redis = None
        core.wait_for_results(
            5,
            timeout=2.0,
            poll_interval=0.01,
            progress_total=5,
            progress_base=0,
        )
        # ‰ archive heartbeats are DEBUG (avoid duplicating DataRecorder noise).
        self.assertTrue(any("samples archived" in line for line in debug_messages))
        # Final drain line remains INFO on Jarvis-HEP.
        self.assertTrue(any("sample drain complete" in line for line in info_messages))

    def test_wait_for_results_uses_archive_baseline_and_worker_completion(self) -> None:
        core = Jarvis2Core()
        messages: list[str] = []

        class _Log:
            def debug(self, msg, *a, **k):  # noqa: ANN001
                messages.append(str(msg) % a if a else str(msg))

            def info(self, msg, *a, **k):  # noqa: ANN001
                messages.append(str(msg) % a if a else str(msg))

            def warning(self, msg, *a, **k):  # noqa: ANN001
                messages.append(str(msg) % a if a else str(msg))

        core._logger = _Log()  # type: ignore[assignment]
        state = {"polls": 0}

        class _Archiver:
            @property
            def records_written(self):
                state["polls"] += 1
                # Twenty rows already existed before this ten-sample batch.
                return 20 if state["polls"] < 4 else 30

        class _Redis:
            def fetch_sample_stats(self):
                if state["polls"] < 4:
                    return {"completed": 0, "failed": 0, "running": 1}
                return {"completed": 10, "failed": 0, "running": 0}

            def get_queue_lengths(self):
                if state["polls"] < 4:
                    return {"task_queue_length": 9, "archive_queue_length": 0}
                return {"task_queue_length": 0, "archive_queue_length": 0}

        core.archiver = _Archiver()  # type: ignore[assignment]
        core.redis = _Redis()  # type: ignore[assignment]
        core.wait_for_results(
            30,
            timeout=2.0,
            poll_interval=0.01,
            progress_total=10,
            progress_base=20,
            require_worker_completion=True,
        )

        self.assertGreaterEqual(state["polls"], 4)
        self.assertTrue(any("sample drain complete" in line for line in messages))

    def test_wait_for_results_waits_out_a_slow_pending_batch(self) -> None:
        core = Jarvis2Core()

        class _Log:
            def debug(self, msg, *a, **k):  # noqa: ANN001
                return None

            def info(self, msg, *a, **k):  # noqa: ANN001
                return None

            def warning(self, msg, *a, **k):  # noqa: ANN001
                return None

        core._logger = _Log()  # type: ignore[assignment]
        clock = _VirtualClock()

        class _Archiver:
            def persistence_busy(self) -> bool:
                return clock.now < 8.0

            @property
            def records_written(self) -> int:
                return 5 if clock.now >= 8.0 else 4

        class _Redis:
            def fetch_sample_stats(self):
                return {"completed": 5, "failed": 0, "running": 0}

            def get_queue_lengths(self):
                return {"task_queue_length": 0, "archive_queue_length": 0}

        core.archiver = _Archiver()  # type: ignore[assignment]
        core.redis = _Redis()  # type: ignore[assignment]
        with (
            mock.patch(
                "jarvishep2.runtime._scan_driver.time.monotonic",
                clock.monotonic,
            ),
            mock.patch("jarvishep2.runtime._scan_driver.time.sleep", clock.sleep),
        ):
            core.wait_for_results(5, timeout=30.0, poll_interval=0.1)
        self.assertGreaterEqual(clock.now, 8.0)

    def test_wait_for_results_stalls_when_archiver_is_idle(self) -> None:
        core = Jarvis2Core()

        class _Log:
            def debug(self, msg, *a, **k):  # noqa: ANN001
                return None

            def info(self, msg, *a, **k):  # noqa: ANN001
                return None

            def warning(self, msg, *a, **k):  # noqa: ANN001
                return None

        core._logger = _Log()  # type: ignore[assignment]
        clock = _VirtualClock()

        class _Archiver:
            def persistence_busy(self) -> bool:
                return False

            @property
            def records_written(self) -> int:
                return 4

        class _Redis:
            def fetch_sample_stats(self):
                return {"completed": 5, "failed": 0, "running": 0}

            def get_queue_lengths(self):
                return {"task_queue_length": 0, "archive_queue_length": 0}

        core.archiver = _Archiver()  # type: ignore[assignment]
        core.redis = _Redis()  # type: ignore[assignment]
        with (
            mock.patch(
                "jarvishep2.runtime._scan_driver.time.monotonic",
                clock.monotonic,
            ),
            mock.patch("jarvishep2.runtime._scan_driver.time.sleep", clock.sleep),
        ):
            with self.assertRaisesRegex(TimeoutError, "archive drain stalled"):
                core.wait_for_results(5, timeout=30.0, poll_interval=0.1)
        self.assertLess(clock.now, 30.0)
        self.assertGreaterEqual(clock.now, 5.0)


class PersistenceBusyTests(unittest.TestCase):
    def test_thread_archiver_stays_busy_while_flush_holds_the_lock(self) -> None:
        archiver = SimpleArchiver.__new__(SimpleArchiver)
        archiver.processor = mock.Mock()
        archiver.processor._lock = threading.Lock()
        archiver.processor._batch = [{"uuid": "tail"}]
        archiver._thread = mock.Mock()
        archiver._thread.is_alive.return_value = True

        self.assertTrue(archiver.persistence_busy())
        self.assertTrue(archiver.processor._lock.acquire(blocking=False))
        try:
            self.assertTrue(archiver.persistence_busy())
        finally:
            archiver.processor._lock.release()
        archiver.processor._batch.clear()
        self.assertFalse(archiver.persistence_busy())
        archiver.processor._batch.append({"uuid": "tail"})
        archiver._thread.is_alive.return_value = False
        self.assertFalse(archiver.persistence_busy())

    def test_process_archiver_busy_requires_live_loop_and_pending_batch(self) -> None:
        proc = ArchiverProcess.__new__(ArchiverProcess)
        ctx = get_spawn_context()
        proc._loop_alive = ctx.Value("i", 1)
        proc.pending_batch = ctx.Value("i", 49)
        proc.is_alive = mock.Mock(return_value=True)  # type: ignore[method-assign]

        self.assertTrue(proc.persistence_busy())
        proc.pending_batch.value = 0
        self.assertFalse(proc.persistence_busy())
        proc.pending_batch.value = 49
        with proc._loop_alive.get_lock():
            proc._loop_alive.value = 0
        self.assertFalse(proc.persistence_busy())
        with proc._loop_alive.get_lock():
            proc._loop_alive.value = 1
        proc.is_alive.return_value = False
        self.assertFalse(proc.persistence_busy())


class _VirtualClock:
    def __init__(self) -> None:
        self.now = 0.0

    def monotonic(self) -> float:
        return self.now

    def sleep(self, seconds: float) -> None:
        self.now += float(seconds)


if __name__ == "__main__":
    unittest.main()
