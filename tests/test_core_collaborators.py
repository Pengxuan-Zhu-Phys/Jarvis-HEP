#!/usr/bin/env python3
"""D25.3: Jarvis2Core is a façade over private runtime/scan/resume collaborators."""

from __future__ import annotations

import inspect
import unittest

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


if __name__ == "__main__":
    unittest.main()
