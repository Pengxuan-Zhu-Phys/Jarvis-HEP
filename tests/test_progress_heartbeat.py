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
        messages: list[str] = []

        class _Log:
            def info(self, msg, *a, **k):  # noqa: ANN001
                messages.append(str(msg) % a if a else str(msg))

            def warning(self, msg, *a, **k):  # noqa: ANN001
                messages.append(str(msg) % a if a else str(msg))

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
        joined = "\n".join(messages)
        self.assertIn("samples archived", joined)
        self.assertIn("sample drain complete", joined)


if __name__ == "__main__":
    unittest.main()
