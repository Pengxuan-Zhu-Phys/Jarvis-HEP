#!/usr/bin/env python3
"""WP-D3.2 env_setup capture-from-source + cache tests."""

from __future__ import annotations

import os
import subprocess
import tempfile
import unittest

from unittest.mock import patch

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig
from jarvishep2.env_setup import EnvCapture, resolve_env_setup_sources
from jarvishep2.worker import Worker

TESTS_ROOT = os.path.dirname(__file__)
ENV_FIXTURES = os.path.join(TESTS_ROOT, "fixtures", "env_setup")
EXPORT_FOO = os.path.join(ENV_FIXTURES, "export_foo.sh")
FIRST_SCRIPT = os.path.join(ENV_FIXTURES, "first.sh")
SECOND_SCRIPT = os.path.join(ENV_FIXTURES, "second.sh")


class EnvCaptureTests(unittest.TestCase):
    def setUp(self) -> None:
        EnvCapture.clear_cache()
        EnvCapture.set_runner(None)

    def tearDown(self) -> None:
        EnvCapture.clear_cache()
        EnvCapture.set_runner(None)

    def test_capture_from_source_exports_variable(self) -> None:
        captured = EnvCapture.capture_from_source(EXPORT_FOO)
        self.assertEqual(captured.get("JARVIS_ENV_FOO"), "layer-concurrent")
        self.assertIn("PATH", captured)

    def test_merged_env_sources_script_once(self) -> None:
        calls: list[list[str]] = []

        def _spy_runner(command, *, env, check=False, **kwargs):
            calls.append(list(command))
            return subprocess.run(
                command,
                capture_output=True,
                text=True,
                env=env,
                check=check,
            )

        EnvCapture.set_runner(_spy_runner)
        first = EnvCapture.merged_env([EXPORT_FOO])
        second = EnvCapture.merged_env([EXPORT_FOO])
        self.assertEqual(first.get("JARVIS_ENV_FOO"), "layer-concurrent")
        self.assertEqual(second.get("JARVIS_ENV_FOO"), "layer-concurrent")
        self.assertEqual(len(calls), 1)

    def test_merge_order_later_script_overrides(self) -> None:
        merged = EnvCapture.merged_env([FIRST_SCRIPT, SECOND_SCRIPT])
        self.assertEqual(merged.get("JARVIS_ENV_FIRST"), "one")
        self.assertEqual(merged.get("JARVIS_ENV_SECOND"), "two")
        self.assertEqual(merged.get("JARVIS_ENV_ORDER"), "second")

    def test_missing_script_raises_clear_error(self) -> None:
        with self.assertRaises(FileNotFoundError):
            EnvCapture.capture_from_source("/tmp/does-not-exist-env-setup.sh")

    def test_failing_script_raises_clear_error(self) -> None:
        with tempfile.NamedTemporaryFile("w", suffix=".sh", delete=False) as handle:
            handle.write("exit 9\n")
            bad_script = handle.name
        try:
            with self.assertRaises(RuntimeError) as ctx:
                EnvCapture.capture_from_source(bad_script)
            self.assertIn("env_setup failed", str(ctx.exception))
        finally:
            os.remove(bad_script)

    def test_resolve_env_setup_sources_extracts_paths(self) -> None:
        sources = resolve_env_setup_sources(
            [{"source": EXPORT_FOO}, {"source": ""}, {"other": "x"}]
        )
        self.assertEqual(sources, [EXPORT_FOO])


class CalculatorEnvSetupTests(unittest.TestCase):
    def setUp(self) -> None:
        EnvCapture.clear_cache()
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
        )
        self._scheduler.start()

    def tearDown(self) -> None:
        EnvCapture.clear_cache()
        self._scheduler.shutdown(wait=True)

    def test_calculator_subprocess_sees_bound_env(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            marker = os.path.join(tmpdir, "env-marker.txt")
            module = CalculatorModule(
                "EnvCalc",
                {
                    "name": "EnvCalc",
                    "env_setup": [{"source": EXPORT_FOO}],
                    "execution": {
                        "path": tmpdir,
                        "commands": [
                            {
                                "cmd": (
                                    "python3 -c "
                                    f"\"import os; open('{marker}', 'w').write("
                                    "os.environ.get('JARVIS_ENV_FOO', ''))\""
                                ),
                                "cwd": tmpdir,
                            }
                        ],
                        "input": [],
                        "output": [],
                    },
                },
            )
            module.attach_scheduler(self._scheduler)
            module.bind_env(EnvCapture.merged_env([EXPORT_FOO]))
            module.acquire_pack_id("env-pack")
            module.execute({"observables": {}, "params": {}})
            with open(marker, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "layer-concurrent")

    def test_calculator_without_env_setup_unchanged(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "marker.txt")
            module = CalculatorModule(
                "PlainCalc",
                {
                    "name": "PlainCalc",
                    "execution": {
                        "path": tmpdir,
                        "commands": [
                            {
                                "cmd": f"python3 -c \"open('{output_path}', 'w').write('ok')\"",
                                "cwd": tmpdir,
                            }
                        ],
                        "input": [],
                        "output": [],
                    },
                },
            )
            module.attach_scheduler(self._scheduler)
            module.acquire_pack_id("plain-pack")
            module.execute({"observables": {}, "params": {}})
            self.assertTrue(os.path.exists(output_path))


class WorkerEnvSetupTests(unittest.TestCase):
    def setUp(self) -> None:
        EnvCapture.clear_cache()

    def tearDown(self) -> None:
        EnvCapture.clear_cache()

    def test_worker_init_runtime_binds_env_setup_once(self) -> None:
        merge_calls: list[list[str]] = []
        original = EnvCapture.merged_env

        def _counting_merge(scripts, base_env=None):
            merge_calls.append(list(scripts))
            return original(scripts, base_env=base_env)

        worker_config = {
            "calculator_modules": [
                {
                    "name": "EnvCalc",
                    "env_setup": [{"source": EXPORT_FOO}],
                    "execution": {
                        "path": ".",
                        "commands": [{"cmd": "true", "cwd": "."}],
                        "input": [],
                        "output": [],
                    },
                }
            ],
            "likelihood_expressions": [],
        }
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, worker_config)
        with patch.object(EnvCapture, "merged_env", side_effect=_counting_merge):
            worker._init_runtime()
            try:
                module = worker._calculators["EnvCalc"]
                self.assertEqual(len(merge_calls), 1)
                self.assertEqual(merge_calls[0], [EXPORT_FOO])
                self.assertEqual(
                    module._subprocess_env.get("JARVIS_ENV_FOO"),
                    "layer-concurrent",
                )
            finally:
                if worker._scheduler is not None:
                    worker._scheduler.shutdown(wait=True)


if __name__ == "__main__":
    unittest.main()