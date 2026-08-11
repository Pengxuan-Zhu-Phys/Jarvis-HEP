#!/usr/bin/env python3
"""D12.1 CalculatorSpec V1 YAML surface (string commands, selection, pools)."""

from __future__ import annotations

import os
import tempfile
import unittest

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.Module.calculator_spec import (
    CalculatorSpec,
    is_cwd_switch_command,
    next_cwd_from_command,
    normalize_command_list,
)
from jarvishep2.calculator_pools import resolve_calculator_pools
from jarvishep2.worker_config import build_worker_config

class CalculatorSpecV1YamlTests(unittest.TestCase):
    def test_string_installation_and_commands_are_normalized(self) -> None:
        spec = CalculatorSpec.from_config(
            "EggBox",
            {
                "name": "EggBox",
                "clone_shadow": True,
                "path": "/runtime/EggBox/@PackID",
                "source": "/assets/EggBox",
                "installation": ["cp -r ${source}/* ${path}"],
                "initialization": [
                    "cp -r ${source}/input.json ${path}/input.json",
                    "rm -f output.json",
                ],
                "execution": {
                    "path": "/runtime/EggBox/@PackID",
                    "commands": ["./eggbox.py", {"cmd": "echo ok", "cwd": "/tmp"}],
                    "input": [{"name": "inp", "type": "JSON"}],
                    "output": [{"name": "out", "type": "JSON"}],
                },
                "modes": False,
                "make_parallel": 16,
                "required_modules": [],
            },
        )

        self.assertEqual(len(spec.installation), 1)
        self.assertEqual(spec.installation[0]["cmd"], "cp -r ${source}/* ${path}")
        self.assertEqual(spec.installation[0]["cwd"], "/runtime/EggBox/@PackID")
        self.assertEqual(len(spec.initialization), 2)
        self.assertEqual(spec.initialization[1]["cmd"], "rm -f output.json")
        self.assertEqual(len(spec.commands), 2)
        self.assertEqual(
            spec.commands[0], {"cmd": "./eggbox.py", "cwd": "/runtime/EggBox/@PackID"}
        )
        self.assertEqual(spec.commands[1]["cmd"], "echo ok")
        self.assertEqual(spec.commands[1]["cwd"], "/tmp")
        self.assertEqual(len(spec.input_specs), 1)
        self.assertEqual(len(spec.output_specs), 1)
        # V1 unused keys tolerated and preserved in raw.
        self.assertIn("modes", spec.raw)
        self.assertEqual(spec.raw.get("make_parallel"), 16)

    def test_execution_path_is_the_default_command_cwd(self) -> None:
        spec = CalculatorSpec.from_config(
            "SeparateExecutionPath",
            {
                "path": "/runtime/calculator",
                "installation": ["make"],
                "initialization": ["prepare"],
                "execution": {
                    "path": "/runtime/calculator/bin",
                    "commands": ["./run", {"cmd": "echo explicit", "cwd": "/tmp"}],
                },
            },
        )

        self.assertEqual(spec.installation[0]["cwd"], "/runtime/calculator")
        self.assertEqual(spec.initialization[0]["cwd"], "/runtime/calculator")
        self.assertEqual(spec.commands[0]["cwd"], "/runtime/calculator/bin")
        self.assertEqual(spec.commands[1]["cwd"], "/tmp")

    def test_execution_commands_fall_back_to_calculator_path(self) -> None:
        spec = CalculatorSpec.from_config(
            "NoExecutionPath",
            {
                "path": "/runtime/calculator",
                "execution": {"commands": ["./run"]},
            },
        )

        self.assertEqual(spec.commands[0]["cwd"], "/runtime/calculator")

    def test_empty_and_non_command_entries_are_skipped(self) -> None:
        self.assertEqual(normalize_command_list(["", "  ", None], default_cwd="."), ())
        self.assertEqual(
            normalize_command_list([{"cwd": "."}], default_cwd="."),
            (),
        )

    def test_selection_is_captured(self) -> None:
        spec = CalculatorSpec.from_config(
            "Cut",
            {"name": "Cut", "selection": "x + y < 10", "execution": {"commands": []}},
        )
        self.assertEqual(spec.selection, "x + y < 10")
        self.assertIsNone(
            CalculatorSpec.from_config("NoCut", {"name": "NoCut"}).selection
        )

    def test_bare_string_command_is_one_element_list(self) -> None:
        self.assertEqual(
            normalize_command_list("./eggbox.py", default_cwd="/rt"),
            ({"cmd": "./eggbox.py", "cwd": "/rt"},),
        )

    def test_command_alias_is_accepted(self) -> None:
        entries = normalize_command_list(
            [{"command": "echo hi", "cwd": "/tmp"}], default_cwd="."
        )
        self.assertEqual(entries[0]["cmd"], "echo hi")

    def test_cd_updates_cwd_for_following_commands(self) -> None:
        """V1 chain: standalone ``cd path`` is inherited by later commands."""
        pack = "/runtime/microOMEGAs/@PackID"
        entries = normalize_command_list(
            [
                "make",
                "cd vector-iDM",
                "make main=main.cpp",
                'cd "/tmp/other pack"',
                "echo done",
            ],
            default_cwd=pack,
        )
        self.assertEqual(entries[0], {"cmd": "make", "cwd": pack})
        self.assertEqual(entries[1], {"cmd": "cd vector-iDM", "cwd": pack})
        self.assertEqual(
            entries[2],
            {
                "cmd": "make main=main.cpp",
                "cwd": f"{pack}/vector-iDM",
            },
        )
        self.assertEqual(entries[3]["cmd"], 'cd "/tmp/other pack"')
        self.assertEqual(entries[3]["cwd"], f"{pack}/vector-iDM")
        self.assertEqual(entries[4], {"cmd": "echo done", "cwd": "/tmp/other pack"})

    def test_cd_and_make_combined_does_not_chain(self) -> None:
        """Only pure ``cd <path>`` chains; ``cd x && make`` keeps default cwd."""
        root = "/calc/pack"
        entries = normalize_command_list(
            ["cd vector-iDM && make main=main.cpp", "ls"],
            default_cwd=root,
        )
        self.assertEqual(entries[0]["cwd"], root)
        self.assertEqual(entries[1]["cwd"], root)

    def test_next_cwd_from_command_helpers(self) -> None:
        self.assertTrue(is_cwd_switch_command("cd vector-iDM"))
        self.assertTrue(is_cwd_switch_command('cd "/tmp/a b"'))
        self.assertFalse(is_cwd_switch_command("cd vector-iDM && make"))
        self.assertFalse(is_cwd_switch_command("make main=main.cpp"))
        self.assertEqual(
            next_cwd_from_command("cd vector-iDM", "/pack/@PackID"),
            "/pack/@PackID/vector-iDM",
        )
        self.assertEqual(
            next_cwd_from_command("cd /abs/target", "/pack"),
            "/abs/target",
        )

    def test_explicit_mapping_cwd_resets_chain_base(self) -> None:
        entries = normalize_command_list(
            [
                {"cmd": "cd sub", "cwd": "/explicit/root"},
                "pwd",
            ],
            default_cwd="/default",
        )
        self.assertEqual(entries[0]["cwd"], "/explicit/root")
        self.assertEqual(entries[1], {"cmd": "pwd", "cwd": "/explicit/root/sub"})

    def test_installation_cd_chain_in_spec(self) -> None:
        spec = CalculatorSpec.from_config(
            "Micro",
            {
                "name": "Micro",
                "path": "/calc/micro/@PackID",
                "installation": [
                    "make",
                    "cd vector-iDM",
                    "make main=main.cpp",
                ],
                "execution": {"commands": ["./vector-iDM/main a b"]},
            },
        )
        self.assertEqual(spec.installation[1]["cmd"], "cd vector-iDM")
        self.assertEqual(spec.installation[1]["cwd"], "/calc/micro/@PackID")
        self.assertEqual(
            spec.installation[2]["cwd"],
            "/calc/micro/@PackID/vector-iDM",
        )


class CalculatorSelectionExecuteTests(unittest.TestCase):
    def test_selection_false_skips_execute_without_scheduler(self) -> None:
        module = CalculatorModule(
            "Cut",
            {
                "name": "Cut",
                "selection": "x > 5",
                "execution": {"commands": ["./never.py"]},
            },
        )
        result = module.execute({"params": {"x": 1.0}, "observables": {}})
        self.assertEqual(result, {})

    def test_selection_true_requires_scheduler(self) -> None:
        module = CalculatorModule(
            "Cut",
            {
                "name": "Cut",
                "selection": "x > 5",
                "execution": {"commands": ["./never.py"]},
            },
        )
        with self.assertRaisesRegex(RuntimeError, "scheduler is not attached"):
            module.execute({"params": {"x": 10.0}, "observables": {}})


class CalculatorMakeParallerConfigTests(unittest.TestCase):
    def test_top_level_make_parallel_feeds_pool_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            worker = build_worker_config(
                {
                    "project_root": tmpdir,
                    "Calculators": {
                        "make_parallel": 16,
                        "path": "&J/calculators/runtime/program",
                        "Modules": [
                            {
                                "name": "EggBox",
                                "path": ".",
                                "execution": {"commands": ["true"]},
                            }
                        ],
                    },
                },
                task_result_dir=tmpdir,
            )
            self.assertEqual(worker["calculator_make_parallel"], 16)
            self.assertEqual(
                worker["sample_config"].get("calculators_path"),
                os.path.join(tmpdir, "calculators", "runtime", "program"),
            )
            self.assertEqual(
                resolve_calculator_pools(worker),
                {"EggBox": 16},
            )

    def test_module_make_parallel_overrides_global(self) -> None:
        pools = resolve_calculator_pools(
            {
                "calculator_make_parallel": 16,
                "calculator_modules": [{"name": "EggBox", "make_parallel": 3}],
            }
        )
        self.assertEqual(pools, {"EggBox": 3})


if __name__ == "__main__":
    unittest.main()
