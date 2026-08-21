#!/usr/bin/env python3
"""D25.8: EnvReqs.V2 is the user surface; get_runtime_block is the only reader."""

from __future__ import annotations

import os
import tempfile
import unittest

from jarvishep2.core import Jarvis2Core
from jarvishep2.runtime_config import (
    RUNTIME_DEFAULTS,
    get_runtime_block,
    get_watchdog_config,
    normalize_runtime_block,
)
from jarvishep2.task_config import load_task_yaml


class RuntimeReadModelTests(unittest.TestCase):
    def test_defaults_are_redis_not_auto(self) -> None:
        self.assertEqual(RUNTIME_DEFAULTS["mode"], "redis")
        self.assertEqual(normalize_runtime_block(None)["mode"], "redis")
        self.assertEqual(normalize_runtime_block({"mode": "auto"})["mode"], "redis")
        self.assertEqual(get_runtime_block(None)["mode"], "redis")
        self.assertEqual(get_runtime_block({})["mode"], "redis")

    def test_skip_loader_envreqs_v2_matches_loader(self) -> None:
        skip_loader = {
            "EnvReqs": {
                "V2": {
                    "workers": 4,
                    "batch_size": 16,
                    "checkpoint_heartbeat_sec": 12,
                    "worker": {
                        "force_serial_layers": True,
                        "sample_artifacts": "always",
                    },
                    "factory": {"watchdog": {"enabled": False, "stale_sec": 9}},
                    "redis": {"host": "10.1.2.3", "port": 6388, "db": 2},
                }
            }
        }
        runtime = get_runtime_block(skip_loader)
        self.assertEqual(runtime["mode"], "redis")
        self.assertEqual(runtime["workers"], 4)
        self.assertEqual(runtime["batch_size"], 16)
        self.assertEqual(runtime["checkpoint_heartbeat_sec"], 12.0)
        self.assertTrue(runtime["force_serial_layers"])
        self.assertEqual(runtime["sample_artifacts"], "always")
        self.assertEqual(runtime["redis"]["host"], "10.1.2.3")
        self.assertFalse(get_watchdog_config(skip_loader)["enabled"])

        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "task.yaml")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(
                    "EnvReqs:\n"
                    "  V2:\n"
                    "    workers: 4\n"
                    "    batch_size: 16\n"
                    "    checkpoint_heartbeat_sec: 12\n"
                    "    worker:\n"
                    "      force_serial_layers: true\n"
                    "      sample_artifacts: always\n"
                    "    factory:\n"
                    "      watchdog:\n"
                    "        enabled: false\n"
                    "        stale_sec: 9\n"
                    "    redis:\n"
                    "      host: 10.1.2.3\n"
                    "      port: 6388\n"
                    "      db: 2\n"
                    "Scan:\n"
                    "  name: read-model\n"
                    "Sampling:\n"
                    "  Method: Bridson\n"
                )
            loaded = load_task_yaml(path)

        loaded_runtime = get_runtime_block(loaded)
        for key in (
            "mode",
            "workers",
            "batch_size",
            "checkpoint_heartbeat_sec",
            "force_serial_layers",
            "sample_artifacts",
        ):
            self.assertEqual(loaded_runtime[key], runtime[key], key)
        self.assertEqual(loaded_runtime["redis"]["host"], runtime["redis"]["host"])
        self.assertEqual(loaded["Runtime"], loaded_runtime)

    def test_core_raw_dict_is_redis(self) -> None:
        core = Jarvis2Core({"EnvReqs": {"V2": {"workers": 3}}})
        self.assertTrue(core.is_redis_runtime())
        self.assertEqual(core.runtime["workers"], 3)
        self.assertEqual(core.runtime["mode"], "redis")


if __name__ == "__main__":
    unittest.main()
