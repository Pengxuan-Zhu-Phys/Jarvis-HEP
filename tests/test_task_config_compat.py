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

    def test_envreqs_v2_redis_factory_worker_groups(self) -> None:
        """D12.4: optional redis/factory/worker groups merge into Runtime."""
        from jarvishep2.runtime_config import (
            get_factory_config,
            get_redis_config,
            get_watchdog_config,
        )
        from jarvishep2.worker_config import build_worker_config

        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "bin"))
            os.makedirs(os.path.join(root, "deps"))
            with open(os.path.join(root, "jarvis.project.yaml"), "w", encoding="utf-8") as handle:
                handle.write("project: envreqs-groups\n")
            with open(
                os.path.join(root, "deps", "environment_default.yaml"), "w", encoding="utf-8"
            ) as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    workers: 2\n"
                    "    redis:\n"
                    "      host: 10.0.0.5\n"
                    "      port: 6381\n"
                    "    factory:\n"
                    "      monitor_hz: 60\n"
                    "      watchdog:\n"
                    "        enabled: false\n"
                    "        stale_sec: 12\n"
                    "    worker:\n"
                    "      force_serial_layers: true\n"
                    "      sample_artifacts: always\n"
                )
            task_path = os.path.join(root, "bin", "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Check_default_dependencies:\n"
                    "    required: true\n"
                    "    default_yaml_path: '&J/deps/environment_default.yaml'\n"
                    "  V2:\n"
                    "    redis:\n"
                    "      db: 3\n"
                    "Scan:\n"
                    "  name: groups\n"
                    "Sampling:\n"
                    "  Method: Bridson\n"
                )

            config = load_task_yaml(task_path)

        redis = get_redis_config(config)
        self.assertEqual(redis["host"], "10.0.0.5")
        self.assertEqual(redis["port"], 6381)
        self.assertEqual(redis["db"], 3)
        self.assertEqual(config["Runtime"]["redis"]["host"], "10.0.0.5")
        self.assertEqual(config["Runtime"]["sample_artifacts"], "always")
        self.assertTrue(config["Runtime"]["force_serial_layers"])

        factory = get_factory_config(config)
        self.assertEqual(factory["monitor_hz"], 60.0)
        watchdog = get_watchdog_config(config)
        self.assertFalse(watchdog["enabled"])
        self.assertEqual(watchdog["stale_sec"], 12.0)

        worker_cfg = build_worker_config(config, task_result_dir="/tmp/out")
        self.assertTrue(worker_cfg["force_serial_layers"])
        self.assertEqual(worker_cfg["sample_config"]["sample_artifacts"], "always")

    def test_redis_partial_override_keeps_internal_defaults(self) -> None:
        from jarvishep2.runtime_config import get_redis_config, normalize_redis_config
        from jarvishep2.redis_queue import INTERNAL_REDIS_CONFIG

        only_port = normalize_redis_config({"port": 6390})
        self.assertEqual(only_port["host"], INTERNAL_REDIS_CONFIG["host"])
        self.assertEqual(only_port["port"], 6390)
        self.assertEqual(only_port["db"], INTERNAL_REDIS_CONFIG["db"])

        with tempfile.TemporaryDirectory() as tmpdir:
            task_path = os.path.join(tmpdir, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    redis:\n"
                    "      port: 6399\n"
                    "Scan:\n"
                    "  name: port-only\n"
                    "Sampling:\n"
                    "  Method: Bridson\n"
                )
            config = load_task_yaml(task_path)

        redis = get_redis_config(config)
        self.assertEqual(redis["port"], 6399)
        self.assertEqual(redis["host"], INTERNAL_REDIS_CONFIG["host"])
        self.assertEqual(redis["db"], INTERNAL_REDIS_CONFIG["db"])

    def test_unknown_envreqs_v2_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            task_path = os.path.join(tmpdir, "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    not_a_real_knob: 1\n"
                    "Scan:\n"
                    "  name: bad\n"
                    "Sampling:\n"
                    "  Method: Bridson\n"
                )
            with self.assertRaisesRegex(ValueError, "unsupported.*not_a_real_knob"):
                load_task_yaml(task_path)

    def test_legacy_runtime_defaults_file_strips_unknown_keys(self) -> None:
        """EnvReqs.Runtime defaults may carry V1 Runtime keys; only V2 knobs apply."""
        with tempfile.TemporaryDirectory() as root:
            os.makedirs(os.path.join(root, "bin"))
            os.makedirs(os.path.join(root, "deps"))
            with open(os.path.join(root, "jarvis.project.yaml"), "w", encoding="utf-8") as handle:
                handle.write("project: legacy-runtime\n")
            defaults = os.path.join(root, "deps", "runtime_defaults.yaml")
            with open(defaults, "w", encoding="utf-8") as handle:
                handle.write(
                    "Runtime:\n"
                    "  mode: redis\n"
                    "  workers: 6\n"
                    "  batch_size: 64\n"
                    "  Subprocess:\n"
                    "    timeout: 1\n"
                )
            task_path = os.path.join(root, "bin", "task.yaml")
            with open(task_path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  Runtime:\n"
                    "    default_runtime_settings: '&J/deps/runtime_defaults.yaml'\n"
                    "Scan:\n"
                    "  name: legacy\n"
                    "Sampling:\n"
                    "  Method: Bridson\n"
                )
            config = load_task_yaml(task_path)

        self.assertEqual(config["Runtime"]["workers"], 6)
        self.assertEqual(config["Runtime"]["batch_size"], 64)
        self.assertEqual(config["Runtime"]["mode"], "redis")


if __name__ == "__main__":
    unittest.main()
