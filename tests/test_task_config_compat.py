#!/usr/bin/env python3
"""V1 task-YAML compatibility at the V2 loader boundary."""

from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.task_config import load_task_yaml


class TaskConfigCompatibilityTests(unittest.TestCase):
    def test_v1_yaml_without_runtime_uses_redis(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            task_path = os.path.join(tmpdir, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write("Scan:\n  name: legacy\nSampling:\n  Method: Bridson\n")

            config = load_task_yaml(task_path)

        self.assertEqual(config["Runtime"]["mode"], "redis")

    def test_top_level_runtime_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            task_path = os.path.join(tmpdir, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "Runtime:\n  mode: auto\nScan:\n  name: explicit\n"
                    "Sampling:\n  Method: Bridson\n"
                )

            with self.assertRaisesRegex(ValueError, "top-level Runtime"):
                load_task_yaml(task_path)


    def test_v1_envreqs_path_supplies_v2_defaults_and_task_overrides_them(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "bin"))
            os.makedirs(os.path.join(root, "deps"))
            with open(os.path.join(root, "jarvis.project.yaml"), "w", encoding="utf-8") as handle:
                handle.write("project: compatibility-test\n")
            with open(
                os.path.join(root, "deps", "environment_default.yaml"), "w", encoding="utf-8"
            ) as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    workers: 2\n"
                    "    batch_size: 32\n"
                )
            task_path = os.path.join(root, "bin", "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Check_default_dependencies:\n"
                    "    required: true\n"
                    "    default_yaml_path: '&J/deps/environment_default.yaml'\n"
                    "  V2:\n"
                    "    worker: 4\n"
                )

            config = load_task_yaml(task_path)

        self.assertEqual(config["Runtime"]["mode"], "redis")
        self.assertEqual(config["Runtime"]["workers"], 4)
        self.assertEqual(config["Runtime"]["batch_size"], 32)
        self.assertEqual(config["EnvReqs"]["V2"], {"workers": 4, "batch_size": 32})
        self.assertTrue(config["EnvReqs"]["Check_default_dependencies"]["required"])

    def test_required_v1_default_yaml_must_exist(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            task_path = os.path.join(root, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Check_default_dependencies:\n"
                    "    required: true\n"
                    "    default_yaml_path: '&J/deps/missing.yaml'\n"
                )

            with self.assertRaisesRegex(FileNotFoundError, "default YAML not found"):
                load_task_yaml(task_path)

    def test_task_root_environment_override_is_honored(self) -> None:
        with tempfile.TemporaryDirectory() as yaml_dir, tempfile.TemporaryDirectory() as root:
            task_path = os.path.join(yaml_dir, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write("Scan:\n  name: isolated\nSampling:\n  Method: Bridson\n")

            with mock.patch.dict(os.environ, {"JARVIS_HEP_TASK_ROOT": root}, clear=False):
                config = load_task_yaml(task_path)

        self.assertEqual(config["project_root"], root)
        self.assertEqual(config["task_root"], root)
        self.assertEqual(config["task_result_dir"], os.path.join(root, "outputs", "isolated"))


if __name__ == "__main__":
    unittest.main()
