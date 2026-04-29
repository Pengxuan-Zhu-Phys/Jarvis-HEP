#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import os
import sys
import tempfile
import unittest
from unittest.mock import Mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.library import LibraryModule  # noqa: E402
from jarvishep.config import ConfigLoader  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None

    def bind(self, **_kwargs):
        return self


def _build_loader(cwd: str) -> ConfigLoader:
    loader = ConfigLoader()
    loader.logger = Mock()
    loader.path["task_root"] = cwd
    loader.path["jpath"] = cwd
    return loader


class CommandChainSemanticsTests(unittest.TestCase):
    def test_calculator_execution_cwd_inherits_after_cd(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dep_root = os.path.join(tmpdir, "deps")
            calc_exec_root = os.path.join(tmpdir, "calculators", "runtime", "program", "demo")
            loader = _build_loader(tmpdir)
            loader.config = {
                "LibDeps": {"path": dep_root, "Modules": []},
                "Calculators": {
                    "path": "&J/calculators/runtime/program",
                    "Modules": [
                        {
                            "name": "DemoCalc",
                            "path": "&J/calculators/runtime/program/demo",
                            "source": "&J/deps/demo-src",
                            "installation": [],
                            "initialization": [],
                            "execution": {
                                "path": "&J/calculators/runtime/program/demo",
                                "commands": [
                                    "cd ${LibDeps:path}",
                                    "cp ${LibDeps:path}/src.tar.gz ./",
                                ],
                                "input": [],
                                "output": [],
                            },
                        }
                    ],
                },
            }

            loader.analysis_calculator()
            commands = (
                loader.config["Calculators"]["Modules"][0]["execution"]["commands"]
            )
            self.assertEqual(commands[0]["cmd"], f"cd {dep_root}")
            self.assertEqual(commands[0]["cwd"], calc_exec_root)
            self.assertEqual(commands[1]["cmd"], f"cp {dep_root}/src.tar.gz ./")
            self.assertEqual(commands[1]["cwd"], dep_root)

    def test_unresolved_placeholder_fails_fast(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            loader = _build_loader(tmpdir)
            loader.config = {
                "LibDeps": {"path": os.path.join(tmpdir, "deps"), "Modules": []},
                "Calculators": {
                    "path": "&J/calculators/runtime/program",
                    "Modules": [
                        {
                            "name": "DemoCalc",
                            "path": "&J/calculators/runtime/program/demo",
                            "source": "&J/deps/demo-src",
                            "installation": ["cp ${MISSING:path} ./"],
                            "initialization": [],
                            "execution": {
                                "path": "&J/calculators/runtime/program/demo",
                                "commands": [],
                                "input": [],
                                "output": [],
                            },
                        }
                    ],
                },
            }

            with self.assertRaises(SystemExit) as cm:
                loader.analysis_calculator()
            self.assertEqual(cm.exception.code, 2)

    def test_cd_with_quoted_path_updates_cwd(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dep_root = os.path.join(tmpdir, "deps with space")
            start_cwd = os.path.join(tmpdir, "start")
            loader = _build_loader(tmpdir)
            loader.config = {
                "LibDeps": {"path": dep_root},
            }
            cmd, next_cwd = loader.decode_calc_command_via_config(
                'cd "${LibDeps:path}"',
                start_cwd,
                {},
            )
            self.assertEqual(cmd["cwd"], start_cwd)
            self.assertEqual(next_cwd, dep_root)

    def test_library_command_nonzero_exit_raises(self):
        module = LibraryModule(
            name="DemoLib",
            required_modules=[],
            installed=False,
            installation={"commands": []},
        )
        module.logger = _NoopLogger()
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaises(RuntimeError):
                loop = asyncio.new_event_loop()
                try:
                    loop.run_until_complete(
                        module.run_command(
                            {
                                "cmd": "sh -c 'exit 3'",
                                "cwd": tmpdir,
                            }
                        )
                    )
                finally:
                    loop.close()


if __name__ == "__main__":
    unittest.main()
