#!/usr/bin/env python3
"""D8.3 human path: interrupt → force runtime checkpoint (agent pieces parked)."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.core import Jarvis2Core
from jarvishep2.Sampling.runtime_checkpoint import load_checkpoint


class InterruptCheckpointTests(unittest.TestCase):
    def test_save_interrupt_checkpoint_writes_payload(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            core = Jarvis2Core(
                {
                    "project_name": "interrupt-ckpt",
                    "Scan": {"name": "scan"},
                    "task_result_dir": tmpdir,
                    "Runtime": {"mode": "redis", "workers": 0},
                }
            )
            # Minimal sampler with export_runtime_state.
            sampler = mock.Mock()
            sampler.export_runtime_state.return_value = {
                "submitted_uuids": ["a", "b"],
                "completed_uuids": ["a"],
                "cursor": 1,
            }
            core.sampler = sampler
            core.info["task_result_dir"] = tmpdir
            core.info["scan_name"] = "scan"
            core.info["sampler_name"] = "MockSampler"

            path = core._save_interrupt_checkpoint()
            self.assertIsNotNone(path)
            assert path is not None
            self.assertTrue(os.path.isfile(path))
            payload = load_checkpoint(path)
            integrity = dict(payload.get("integrity") or {})
            self.assertEqual(str(integrity.get("checkpoint_reason") or ""), "interrupt")
            state = dict(payload.get("sampler_state") or {})
            self.assertEqual(state.get("submitted_uuids"), ["a", "b"])
            self.assertEqual(state.get("completed_uuids"), ["a"])

    def test_keyboard_interrupt_triggers_checkpoint_before_shutdown(self) -> None:
        core = Jarvis2Core(
            {
                "task_result_dir": "/tmp",
                "Runtime": {"mode": "redis", "workers": 0},
            }
        )
        calls: list[str] = []

        def _bootstrap() -> None:
            calls.append("bootstrap")

        def _scan() -> int:
            calls.append("scan")
            raise KeyboardInterrupt("test interrupt")

        def _ckpt() -> str | None:
            calls.append("checkpoint")
            return "/tmp/fake.pkl"

        def _shutdown(*, wait: bool = True, write_run_summary: bool = False) -> None:
            calls.append(f"shutdown:summary={write_run_summary}")

        core.bootstrap_distributed_runtime = _bootstrap  # type: ignore[method-assign]
        core.run_distributed_scan = _scan  # type: ignore[method-assign]
        core._save_interrupt_checkpoint = _ckpt  # type: ignore[method-assign]
        core.shutdown = _shutdown  # type: ignore[method-assign]
        core.prepare_resume = lambda **_k: None  # type: ignore[method-assign]
        core._install_control_signal_handlers = lambda: None  # type: ignore[method-assign]
        core._restore_control_signal_handlers = lambda: None  # type: ignore[method-assign]

        with mock.patch(
            "jarvishep2.core.sampling_method",
            return_value="Bridson",
        ), mock.patch(
            "jarvishep2.core.STATELESS_METHODS",
            frozenset({"Bridson"}),
        ), mock.patch(
            "jarvishep2.core.is_check_modules_task",
            return_value=False,
        ):
            outcome = core.run(write_run_summary=True)

        self.assertEqual(outcome.status, "interrupted")
        self.assertEqual(
            calls,
            ["bootstrap", "scan", "checkpoint", "shutdown:summary=False"],
        )


if __name__ == "__main__":
    unittest.main()
