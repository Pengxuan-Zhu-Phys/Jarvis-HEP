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
    def test_execution_stage_replaces_sample_id_in_cmd_and_cwd(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "samples")
            module.sample_info = {"uuid": "S-1001", "save_dir": sample_dir}
            captured = {}

            async def _fake_local(command, stage, command_index):
                captured["command"] = dict(command)
                captured["stage"] = stage
                captured["command_index"] = command_index

            module._run_command_local = _fake_local  # type: ignore[method-assign]
            _run_async(
                module.run_command(
                    {
                        "cmd": "echo @SampleID > out_@SampleID.txt",
                        "cwd": os.path.join(tmpdir, "@SampleID"),
                    },
                    stage="execution",
                    command_index=9,
                )
            )

            self.assertEqual(captured["stage"], "execution")
            self.assertEqual(captured["command_index"], 9)
            self.assertEqual(captured["command"]["cmd"], "echo S-1001 > out_S-1001.txt")
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

        async def _fake_local(command, stage, command_index):
            captured["command"] = dict(command)
            captured["stage"] = stage
            captured["command_index"] = command_index

        module._run_command_local = _fake_local  # type: ignore[method-assign]
        _run_async(
            module.run_command(
                {"cmd": "echo @SampleID", "cwd": "."},
                stage="install",
                command_index=2,
            )
        )
        self.assertEqual(captured["stage"], "install")
        self.assertEqual(captured["command"]["cmd"], "echo @SampleID")

    def test_scheduler_command_submission_is_loguru_only(self):
        module = CalculatorModule("DemoCalc", _minimal_calc_config())
        module.sample_info = {"uuid": "S-42", "save_dir": "/tmp/s42"}

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


if __name__ == "__main__":
    unittest.main()
