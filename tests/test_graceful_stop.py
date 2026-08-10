#!/usr/bin/env python3
"""D8.3 human path: interrupt → force runtime checkpoint (agent pieces parked)."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.core import Jarvis2Core
from jarvishep2.run_outcome import RunOutcome
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
                "cursor": 1,
                "ready_queue": [],
            }
            sampler.submitted_uuids = frozenset({"a", "b"})
            sampler.at_safe_barrier.return_value = False
            sampler.checkpoint_runtime_state.return_value = (
                sampler.export_runtime_state.return_value
            )
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
            self.assertEqual(state.get("cursor"), 1)
            self.assertNotIn("completed_uuids", state)

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
        ):
            outcome = core.run(write_run_summary=True)

        self.assertEqual(outcome.status, "interrupted")
        self.assertEqual(
            calls,
            ["bootstrap", "scan", "checkpoint", "shutdown:summary=False"],
        )


class FinalArchiveVerificationTests(unittest.TestCase):
    def _make_core(self) -> Jarvis2Core:
        core = Jarvis2Core({"Runtime": {"mode": "redis", "workers": 0}})
        core.prepare_resume = lambda **_k: None  # type: ignore[method-assign]
        core._install_control_signal_handlers = lambda: None  # type: ignore[method-assign]
        core._restore_control_signal_handlers = lambda: None  # type: ignore[method-assign]
        core.bootstrap_distributed_runtime = lambda: None  # type: ignore[method-assign]
        core.shutdown = lambda **_k: None  # type: ignore[method-assign]
        return core

    def test_partial_failure_still_emits_plot_scenes(self) -> None:
        """D19.2: mixed completed/failed must emit plots, not silent-skip."""
        core = self._make_core()
        core.run_distributed_scan = mock.Mock(return_value=40)  # type: ignore[method-assign]
        core._wait_for_archive_caught_up = mock.Mock()  # type: ignore[method-assign]
        core._finalize_nested_result_csv = mock.Mock()  # type: ignore[method-assign]
        core._emit_plot_scenes = mock.Mock()  # type: ignore[method-assign]
        partial = RunOutcome(submitted=40, completed=20, failed=20, archived=40)
        core._capture_run_outcome = mock.Mock(return_value=partial)  # type: ignore[method-assign]

        with (
            mock.patch("jarvishep2.core.sampling_method", return_value="Bridson"),
            mock.patch("jarvishep2.core.STATELESS_METHODS", frozenset({"Bridson"})),
        ):
            outcome = core.run(write_run_summary=False)

        self.assertEqual(outcome.status, "partial_failure")
        core._finalize_nested_result_csv.assert_called_once()
        core._emit_plot_scenes.assert_called_once()

    def test_full_failure_skips_plot_scenes_with_reason(self) -> None:
        core = self._make_core()
        core.run_distributed_scan = mock.Mock(return_value=10)  # type: ignore[method-assign]
        core._wait_for_archive_caught_up = mock.Mock()  # type: ignore[method-assign]
        core._finalize_nested_result_csv = mock.Mock()  # type: ignore[method-assign]
        core._emit_plot_scenes = mock.Mock()  # type: ignore[method-assign]
        failed = RunOutcome(submitted=10, completed=0, failed=10, archived=10)
        core._capture_run_outcome = mock.Mock(return_value=failed)  # type: ignore[method-assign]

        with (
            mock.patch("jarvishep2.core.sampling_method", return_value="Bridson"),
            mock.patch("jarvishep2.core.STATELESS_METHODS", frozenset({"Bridson"})),
        ):
            outcome = core.run(write_run_summary=False)

        self.assertEqual(outcome.status, "failed")
        core._finalize_nested_result_csv.assert_not_called()
        core._emit_plot_scenes.assert_not_called()

    def test_normal_run_verifies_archive_and_recaptures_outcome(self) -> None:
        core = self._make_core()
        core.run_distributed_scan = mock.Mock(return_value=3)  # type: ignore[method-assign]
        core._wait_for_archive_caught_up = mock.Mock()  # type: ignore[method-assign]
        core._finalize_nested_result_csv = mock.Mock()  # type: ignore[method-assign]
        core._emit_plot_scenes = mock.Mock()  # type: ignore[method-assign]
        core._capture_run_outcome = mock.Mock(  # type: ignore[method-assign]
            side_effect=[
                RunOutcome(submitted=3, completed=3, archived=0),
                RunOutcome(submitted=3, completed=3, archived=3),
            ]
        )

        with (
            mock.patch("jarvishep2.core.sampling_method", return_value="Bridson"),
            mock.patch("jarvishep2.core.STATELESS_METHODS", frozenset({"Bridson"})),
        ):
            outcome = core.run()

        core._wait_for_archive_caught_up.assert_called_once_with(timeout=120.0)
        self.assertEqual(outcome.archived, 3)
        core._finalize_nested_result_csv.assert_called_once_with()
        core._emit_plot_scenes.assert_called_once_with()

    def test_check_run_skips_normal_exit_archive_verification(self) -> None:
        core = self._make_core()
        core._apply_check_modules_runtime_policy = mock.Mock()  # type: ignore[method-assign]
        core.run_check_modules = mock.Mock(return_value=1)  # type: ignore[method-assign]
        core._wait_for_archive_caught_up = mock.Mock()  # type: ignore[method-assign]
        core._capture_run_outcome = mock.Mock(  # type: ignore[method-assign]
            return_value=RunOutcome(submitted=1, completed=1, archived=1)
        )

        outcome = core.run(check_modules=True)

        self.assertTrue(outcome.ok)
        core._wait_for_archive_caught_up.assert_not_called()


if __name__ == "__main__":
    unittest.main()
