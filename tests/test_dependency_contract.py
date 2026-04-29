#!/usr/bin/env python3
from __future__ import annotations

import os
import re
import sys
import unittest
from unittest.mock import Mock, patch

import yaml


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.config import ConfigLoader  # noqa: E402


def _normalize_requirement_name(requirement: str) -> str:
    name = re.split(r"[<>=!~\[\]\s]", requirement, maxsplit=1)[0]
    return name.strip().lower().replace("_", "-")


def _load_pyproject_dependency_names(pyproject_path: str) -> set[str]:
    with open(pyproject_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    in_dependencies = False
    dep_entries = []
    for line in lines:
        stripped = line.strip()
        if not in_dependencies:
            if stripped.startswith("dependencies") and stripped.endswith("["):
                in_dependencies = True
            continue
        if stripped.startswith("]"):
            break
        match = re.search(r'"([^"]+)"', line)
        if match:
            dep_entries.append(match.group(1))

    return {_normalize_requirement_name(item) for item in dep_entries}


def _load_env_required_dependency_names(env_yaml_path: str) -> set[str]:
    with open(env_yaml_path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    deps = data["EnvReqs"]["Python"]["Dependencies"]
    required = [item["name"] for item in deps if item.get("required", False)]
    return {_normalize_requirement_name(item) for item in required}


class DependencyContractTests(unittest.TestCase):
    def test_pyproject_and_env_default_required_dependencies_are_aligned(self):
        pyproject_path = os.path.join(PROJECT_ROOT, "pyproject.toml")
        env_yaml_path = os.path.join(PROJECT_ROOT, "jarvishep", "card", "environment_default.yaml")

        pyproject_names = _load_pyproject_dependency_names(pyproject_path)
        env_required_names = _load_env_required_dependency_names(env_yaml_path)

        self.assertEqual(
            pyproject_names,
            env_required_names,
            msg=(
                f"Dependency mismatch.\n"
                f"Only in pyproject: {sorted(pyproject_names - env_required_names)}\n"
                f"Only in environment_default: {sorted(env_required_names - pyproject_names)}"
            ),
        )

    def test_legacy_optional_dependencies_not_in_required_env_defaults(self):
        env_yaml_path = os.path.join(PROJECT_ROOT, "jarvishep", "card", "environment_default.yaml")
        env_required_names = _load_env_required_dependency_names(env_yaml_path)
        self.assertFalse({"dynesty", "pyhf", "sqlalchemy"} & env_required_names)

    def test_missing_required_dependency_reports_error_without_auto_install(self):
        loader = ConfigLoader()
        loader.logger = Mock()
        loader.config = {
            "EnvReqs": {
                "Python": {
                    "version": ">=3.10",
                    "Dependencies": [
                        {
                            "name": "__jarvishep_missing_dependency__",
                            "required": True,
                            "version": ">=1.0.0",
                        }
                    ],
                }
            }
        }

        with patch("jarvishep.config.subprocess.run") as subprocess_run:
            with self.assertRaises(SystemExit) as cm:
                loader.check_PYTHON_env()

            self.assertEqual(cm.exception.code, 2)
            subprocess_run.assert_not_called()

        error_text = " ".join(
            str(call.args[0]) for call in loader.logger.error.call_args_list if call.args
        )
        self.assertIn("Install/upgrade required dependencies via:", error_text)
        self.assertIn("python3 -m pip install -U", error_text)

    def test_directory_setting_archive_samples_default_true(self):
        loader = ConfigLoader()
        loader.logger = Mock()
        loader.config = {}
        loader._normalize_optional_sections()
        self.assertIn("Directory_Setting", loader.config)
        self.assertTrue(loader.config["Directory_Setting"].get("archive_samples"))

    def test_directory_setting_archive_samples_respects_user_false(self):
        loader = ConfigLoader()
        loader.logger = Mock()
        loader.config = {"Directory_Setting": {"archive_samples": False}}
        loader._normalize_optional_sections()
        self.assertFalse(loader.config["Directory_Setting"].get("archive_samples"))


if __name__ == "__main__":
    unittest.main()
