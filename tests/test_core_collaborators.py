#!/usr/bin/env python3
"""D25.3: Jarvis2Core is a façade over private runtime/scan/resume collaborators."""

from __future__ import annotations

import inspect
import logging
import threading
import unittest
from unittest import mock

from jarvishep2._resume_service import _ResumeService
from jarvishep2._runtime_supervisor import _RuntimeSupervisor
from jarvishep2._scan_driver import _ScanDriver
from jarvishep2.core import Jarvis2Core


class CoreCollaboratorTests(unittest.TestCase):
    def test_core_holds_private_collaborators(self) -> None:
        core = Jarvis2Core()
        self.assertIsInstance(core._runtime, _RuntimeSupervisor)
        self.assertIsInstance(core._scan, _ScanDriver)
        self.assertIsInstance(core._resume, _ResumeService)
        self.assertIs(core._runtime._core, core)
        self.assertIs(core._scan._core, core)
        self.assertIs(core._resume._core, core)
        self.assertIsInstance(core._proc_board_lock, type(threading.Lock()))

    def test_new_without_init_still_lazily_builds_collaborators(self) -> None:
        core = Jarvis2Core.__new__(Jarvis2Core)
        self.assertIsInstance(core._get_scan(), _ScanDriver)
        self.assertIsInstance(core._get_runtime(), _RuntimeSupervisor)
        self.assertIsInstance(core._get_resume(), _ResumeService)

    def test_run_and_shutdown_signatures_unchanged(self) -> None:
        run_params = inspect.signature(Jarvis2Core.run).parameters
        self.assertEqual(
            list(run_params),
            [
                "self",
                "resume",
                "check_modules",
                "verify_golden",
                "write_run_summary",
                "check_timeout",
            ],
        )
        shutdown_params = inspect.signature(Jarvis2Core.shutdown).parameters
        self.assertEqual(list(shutdown_params), ["self", "wait", "write_run_summary"])

    def test_dead_archiver_process_is_restarted(self) -> None:
        class _DeadArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 4242

            def is_alive(self) -> bool:
                return False

            def join(self, timeout=None) -> None:
                return None

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.archiver_restart")
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _DeadArchiver
        ):
            core.archiver = _DeadArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_called_once_with("/tmp/samples.hdf5")
        self.assertFalse(core._interrupt_requested)

    def test_archiver_restart_failure_requests_shutdown(self) -> None:
        class _DeadArchiver:
            db_path = "/tmp/samples.hdf5"
            pid = 7

            def is_alive(self) -> bool:
                return False

            def join(self, timeout=None) -> None:
                return None

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core._logger = logging.getLogger("test.archiver_restart")
        core.init_archiver = mock.Mock(side_effect=RuntimeError("hdf5 locked"))
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _DeadArchiver
        ):
            core.archiver = _DeadArchiver()
            core._runtime._ensure_archiver_alive()
        self.assertTrue(core._interrupt_requested)

    def test_archiver_liveness_skips_live_process(self) -> None:
        class _LiveArchiver:
            db_path = "/tmp/samples.hdf5"

            def is_alive(self) -> bool:
                return True

        core = Jarvis2Core()
        core._shutdown_done = False
        core._interrupt_requested = False
        core.init_archiver = mock.Mock()
        with mock.patch(
            "jarvishep2.runtime._runtime_supervisor.ArchiverProcess", _LiveArchiver
        ):
            core.archiver = _LiveArchiver()
            core._runtime._ensure_archiver_alive()
        core.init_archiver.assert_not_called()


if __name__ == "__main__":
    unittest.main()
