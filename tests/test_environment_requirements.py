#!/usr/bin/env python3
"""V1 EnvReqs Python/CERN_ROOT compatibility tests."""

from __future__ import annotations

import subprocess
import unittest
from unittest import mock

from jarvishep2.core import Jarvis2Core
from jarvishep2.environment_requirements import check_environment_requirements


class EnvironmentRequirementsTests(unittest.TestCase):
    def test_python_version_and_package_requirement_pass(self) -> None:
        config = {
            "EnvReqs": {
                "Python": {
                    "version": ">=3.10",
                    "Dependencies": [{"name": "example", "required": True, "version": ">=2.0"}],
                }
            }
        }
        with mock.patch("jarvishep2.environment_requirements.importlib_metadata.version", return_value="2.1"):
            report = check_environment_requirements(config)
        self.assertTrue(report.ok)
        self.assertIn("Python version", report.summary)

    def test_required_python_package_failure_is_reported(self) -> None:
        config = {"EnvReqs": {"Python": {"Dependencies": [{"name": "missing", "required": True, "version": ">=1"}]}}}
        with mock.patch(
            "jarvishep2.environment_requirements.importlib_metadata.version",
            side_effect=__import__("importlib").metadata.PackageNotFoundError,
        ):
            report = check_environment_requirements(config)
        self.assertFalse(report.ok)
        self.assertIn("Python package missing>=1 is not installed", report.errors)

    def test_root_version_and_feature_are_checked(self) -> None:
        config = {
            "EnvReqs": {
                "CERN_ROOT": {
                    "required": True,
                    "version": ">=6.20",
                    "get_path_command": "root-config --prefix",
                    "Dependencies": [{"name": "minuit2", "required": True, "check_command": "root-config --has-minuit2", "expected_output": "yes"}],
                }
            }
        }
        completed = [
            subprocess.CompletedProcess([], 0, "6.30/04\n", ""),
            subprocess.CompletedProcess([], 0, "/opt/root\n", ""),
            subprocess.CompletedProcess([], 0, "yes\n", ""),
        ]
        with mock.patch("jarvishep2.environment_requirements._run_command", side_effect=completed):
            report = check_environment_requirements(config)
        self.assertTrue(report.ok)
        self.assertEqual(report.summary["ROOT version"], "6.30/04")
        self.assertTrue(report.summary["ROOT-minuit2"])

    def test_missing_required_root_is_reported(self) -> None:
        config = {"EnvReqs": {"CERN_ROOT": {"required": True}}}
        with mock.patch(
            "jarvishep2.environment_requirements._run_command",
            side_effect=FileNotFoundError("root-config"),
        ):
            report = check_environment_requirements(config)
        self.assertFalse(report.ok)
        self.assertIn("CERN ROOT is not installed", report.errors[0])

    def test_optional_root_path_is_available_to_libdeps_tokens(self) -> None:
        report = check_environment_requirements(
            {"EnvReqs": {"CERN_ROOT": {"required": False, "path": "/opt/root"}}}
        )
        self.assertEqual(report.summary["ROOT path"], "/opt/root")

    def test_core_logs_environment_preflight_success(self) -> None:
        core = Jarvis2Core(
            {
                "EnvReqs": {
                    "Python": {
                        "version": ">=3.10",
                        "Dependencies": [{"name": "example", "required": True, "version": ">=2.0"}],
                    }
                }
            }
        )
        core._logger = mock.Mock()

        with mock.patch("jarvishep2.environment_requirements.importlib_metadata.version", return_value="2.1"):
            core.check_environment_requirements()

        core._logger.warning.assert_called_once()
        message = core._logger.warning.call_args.args[1]
        self.assertIn("Python                     | >=3.10", message)
        self.assertIn("package example", message)
        self.assertIn("Result                     | -", message)
        rows = [line for line in message.splitlines() if "|" in line and "---" not in line]
        self.assertEqual({row.index("|") for row in rows}, {29})


if __name__ == "__main__":
    unittest.main()
