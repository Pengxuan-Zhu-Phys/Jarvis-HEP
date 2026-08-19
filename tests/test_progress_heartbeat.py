#!/usr/bin/env python3
"""Control-plane completion progress + duration helpers."""

from __future__ import annotations

import time
import unittest
from unittest import mock

from jarvishep2.core import Jarvis2Core
from jarvishep2.log_kv import PermilleProgress, format_duration


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


if __name__ == "__main__":
    unittest.main()
