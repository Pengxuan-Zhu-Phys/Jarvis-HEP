#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.calculator import CalculatorModule  # noqa: E402


def _run_async(coro):
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(coro)
    finally:
        loop.close()


def _minimal_calc_config() -> dict:
    return {
        "modes": False,
        "required_modules": [],
        "clone_shadow": False,
        "installation": [],
        "initialization": [],
        "execution": {"commands": [], "input": [], "output": []},
        "path": "&J/calculators/runtime/program/demo",
    }


class CalculatorRuntimeSampleIDTests(unittest.TestCase):
    def test_module_selection_expression_evaluates_against_observables(self):
        config = _minimal_calc_config()
        config["selection"] = "x > 0.5"
        module = CalculatorModule("DemoCalc", config)
        module.set_funcs({})

        self.assertTrue(module.evaluate_selection({"x": 0.6}))
        self.assertFalse(module.evaluate_selection({"x": 0.2}))

    def test_execution_stage_replaces_runtime_tokens_in_cmd_and_cwd(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "outputs", "SAMPLE", "S-1001")
            module.sample_info = {"uuid": "S-1001", "save_dir": sample_dir}
            captured = {}

            async def _fake_local(command, stage, command_index, timeout_sec=None):
                captured["command"] = dict(command)
                captured["stage"] = stage
                captured["command_index"] = command_index

            module._run_command_local = _fake_local  # type: ignore[method-assign]
            _run_async(
                module.run_command(
                    {
                        "cmd": "echo @SampleID > @Sdir/out_@SampleID.txt",
                        "cwd": os.path.join(tmpdir, "@SampleID"),
                    },
                    stage="execution",
                    command_index=9,
                )
            )

            self.assertEqual(captured["stage"], "execution")
            self.assertEqual(captured["command_index"], 9)
            self.assertEqual(captured["command"]["cmd"], f"echo S-1001 > {sample_dir}/out_S-1001.txt")
            self.assertEqual(captured["command"]["cwd"], os.path.join(tmpdir, "S-1001"))

    def test_execution_stage_without_uuid_raises(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        module.sample_info = {"save_dir": "/tmp/none"}
        with self.assertRaises(RuntimeError):
            _run_async(
                module.run_command(
                    {"cmd": "echo @SampleID", "cwd": "."},
                    stage="execution",
                    command_index=1,
                )
            )

    def test_install_stage_skips_runtime_sample_id_replacement(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        module.sample_info = {}
        captured = {}

        async def _fake_local(command, stage, command_index, timeout_sec=None):
            captured["command"] = dict(command)
            captured["stage"] = stage
            captured["command_index"] = command_index

        module._run_command_local = _fake_local  # type: ignore[method-assign]
        _run_async(
            module.run_command(
                {"cmd": "echo @SampleID @Sdir", "cwd": "."},
                stage="install",
                command_index=2,
            )
        )
        self.assertEqual(captured["stage"], "install")
        self.assertEqual(captured["command"]["cmd"], "echo @SampleID @Sdir")

    def test_execution_stage_without_save_dir_raises_for_sdir(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        module.sample_info = {"uuid": "S-1001"}
        with self.assertRaises(RuntimeError):
            _run_async(
                module.run_command(
                    {"cmd": "echo @Sdir", "cwd": "."},
                    stage="execution",
                    command_index=3,
                )
            )

    def test_scheduler_command_submission_includes_stream_logger(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        module.sample_info = {"uuid": "S-42", "save_dir": "/tmp/s42"}

        class _StreamLogger:
            def __init__(self):
                self.messages = []

            def bind(self, **_kwargs):
                return self

            def info(self, message, *_args, **_kwargs):
                self.messages.append(str(message))

        class _FakeScheduler:
            def __init__(self):
                self.jobs = []

            async def arun(self, job):
                self.jobs.append(job)
                return SimpleNamespace(
                    returncode=0,
                    timed_out=False,
                    duration_sec=0.01,
                    stdout_bytes=12,
                    stderr_bytes=0,
                    ok=True,
                )

        fake_scheduler = _FakeScheduler()
        module.subprocess_scheduler = fake_scheduler
        module.logger = _StreamLogger()

        _run_async(
            module.run_command(
                {"cmd": "echo @SampleID", "cwd": "."},
                stage="execution",
                command_index=11,
            )
        )
        self.assertEqual(len(fake_scheduler.jobs), 1)
        submitted = fake_scheduler.jobs[0]
        self.assertEqual(submitted.log_policy, "logger")
        self.assertIsNone(submitted.log_dir)
        self.assertIs(submitted.stream_logger, module.logger)

    def test_install_stage_uses_install_log_even_when_scheduler_exists(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        with tempfile.TemporaryDirectory() as tmpdir:
            module.basepath = tmpdir
            module.assign_ID("007")
            module.create_basic_logger()

            local_calls = []

            async def _fake_local(command, stage, command_index, timeout_sec=None):
                local_calls.append((dict(command), stage, command_index))
                module.logger.bind(raw=True).info("install output line")

            class _FakeScheduler:
                def __init__(self):
                    self.jobs = []

                async def arun(self, job):
                    self.jobs.append(job)
                    return SimpleNamespace(
                        returncode=0,
                        timed_out=False,
                        duration_sec=0.01,
                        stdout_bytes=0,
                        stderr_bytes=0,
                        ok=True,
                    )

            fake_scheduler = _FakeScheduler()
            module._run_command_local = _fake_local  # type: ignore[method-assign]
            module.subprocess_scheduler = fake_scheduler

            try:
                _run_async(
                    module.run_command(
                        {"cmd": "echo install-output", "cwd": tmpdir},
                        stage="install",
                        command_index=1,
                    )
                )
            finally:
                module.logger.remove(module.handlers["install"])
                del module.handlers["install"]

            install_log = os.path.join(
                tmpdir,
                f"Installation_{module.name}-{module.PackID}.log",
            )
            self.assertEqual(len(local_calls), 1)
            self.assertEqual(fake_scheduler.jobs, [])
            with open(install_log, "r", encoding="utf-8") as handle:
                contents = handle.read()
            self.assertIn("install output line", contents)
            self.assertNotIn("[install#", contents)

    def test_install_local_output_is_written_raw_without_formatter_prefix(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        with tempfile.TemporaryDirectory() as tmpdir:
            module.basepath = tmpdir
            module.assign_ID("009")
            module.create_basic_logger()

            try:
                _run_async(
                    module._run_command_local(
                        {
                            "cmd": "printf 'install output line\\n'",
                            "cwd": tmpdir,
                        },
                        stage="install",
                        command_index=3,
                    )
                )
            finally:
                module.logger.remove(module.handlers["install"])
                del module.handlers["install"]

            install_log = os.path.join(
                tmpdir,
                f"Installation_{module.name}-{module.PackID}.log",
            )
            with open(install_log, "r", encoding="utf-8") as handle:
                lines = handle.read().splitlines()

            output_lines = [line for line in lines if "install output line" in line]
            self.assertTrue(output_lines)
            self.assertTrue(any(line == "install output line" for line in output_lines))
            self.assertTrue(all("·•·" not in line for line in output_lines))
            self.assertTrue(all("[install#" not in line for line in output_lines))

            done_lines = [line for line in lines if line.startswith("Command done [install#")]
            self.assertTrue(done_lines)
            self.assertTrue(all("·•·" not in line for line in done_lines))

    def test_calculator_timeout_is_normalized_from_module_config(self):
        config = _minimal_calc_config()
        config["timeout"] = "2.5"

        module = CalculatorModule("DemoCalc", config)

        self.assertEqual(module.timeout, 2.5)

    def test_execution_timeout_budget_is_passed_to_commands(self):
        config = _minimal_calc_config()
        config["timeout"] = 0.5
        config["execution"]["commands"] = [
            {"cmd": "echo one", "cwd": "."},
            {"cmd": "echo two", "cwd": "."},
        ]
        module = CalculatorModule("DemoCalc", config)
        captured_timeouts = []

        async def _fake_run_command(command, stage="execution", command_index=0, timeout_sec=None):
            captured_timeouts.append(timeout_sec)
            await asyncio.sleep(0.01)

        module.run_command = _fake_run_command  # type: ignore[method-assign]

        _run_async(module.execute_commands())

        self.assertEqual(len(captured_timeouts), 2)
        self.assertTrue(all(timeout is not None for timeout in captured_timeouts))
        self.assertTrue(all(0 < timeout <= 0.5 for timeout in captured_timeouts))
        self.assertLess(captured_timeouts[1], captured_timeouts[0])

    def test_execution_timeout_stops_before_next_command_after_budget_expires(self):
        config = _minimal_calc_config()
        config["timeout"] = 0.01
        config["execution"]["commands"] = [
            {"cmd": "echo one", "cwd": "."},
            {"cmd": "echo two", "cwd": "."},
        ]
        module = CalculatorModule("DemoCalc", config)
        calls = []

        async def _fake_run_command(command, stage="execution", command_index=0, timeout_sec=None):
            calls.append(command["cmd"])
            await asyncio.sleep(0.02)

        module.run_command = _fake_run_command  # type: ignore[method-assign]

        with self.assertRaises(RuntimeError) as ctx:
            _run_async(module.execute_commands())

        self.assertEqual(calls, ["echo one"])
        self.assertIn("Calculator execution timed out after", str(ctx.exception))

    def test_local_execution_command_timeout_terminates_process(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(RuntimeError) as ctx:
                _run_async(
                    module.run_command(
                        {
                            "cmd": f'"{sys.executable}" -c "import time; time.sleep(2)"',
                            "cwd": tmpdir,
                        },
                        stage="execution",
                        command_index=4,
                        timeout_sec=0.1,
                    )
                )

        self.assertIn("timeout=True", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
