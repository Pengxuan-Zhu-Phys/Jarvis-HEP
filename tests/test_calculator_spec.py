#!/usr/bin/env python3
"""D12.1 CalculatorSpec V1 YAML surface (string commands, selection, pools)."""

from __future__ import annotations

import os
import tempfile
import unittest
from math import pi

import numpy as np

from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.Module.calculator_spec import (
    CalculatorSpec,
    is_cwd_switch_command,
    next_cwd_from_command,
    normalize_command_list,
)
from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig
from jarvishep2.calculator_pools import resolve_calculator_pools
from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.sample import Sample
from jarvishep2.worker_config import build_worker_config

TESTS_ROOT = os.path.dirname(__file__)
# Prefer the real Jarvis-Examples Eggbox asset when available; else local fixture.
_CANDIDATE_SOURCES = (
    os.path.abspath(
        os.path.join(
            TESTS_ROOT, "..", "..", "Jarvis-Examples", "Eggbox", "deps", "assets", "inertial", "EggBox"
        )
    ),
    os.path.abspath(
        os.path.join(TESTS_ROOT, "..", "Jarvis-Examples", "Eggbox", "deps", "assets", "inertial", "EggBox")
    ),
)
EGGBOX_SOURCE = next((p for p in _CANDIDATE_SOURCES if os.path.isdir(p)), None)

_TEST_SCHEDULERS: list[AsyncSubprocessScheduler] = []


def tearDownModule() -> None:
    for scheduler in _TEST_SCHEDULERS:
        try:
            scheduler.shutdown(wait=True)
        except Exception:
            pass
    _TEST_SCHEDULERS.clear()


def _attach_runtime(module: CalculatorModule, project_root: str) -> CalculatorModule:
    scheduler = AsyncSubprocessScheduler(
        SubprocessRuntimeConfig(max_concurrency=1, log_policy="quiet")
    )
    scheduler.start()
    _TEST_SCHEDULERS.append(scheduler)
    module.attach_scheduler(scheduler)
    module.attach_command_parser(CommandParser(project_root=project_root, scan_name="test"))
    return module


def _v1_eggbox_module(runtime_root: str, source: str) -> dict:
    """Mirror Example_Bridson_process.yaml Calculators.Modules[0] surface."""
    slot = os.path.join(runtime_root, "EggBox", "@PackID")
    return {
        "name": "EggBox",
        "required_modules": [],
        "clone_shadow": True,
        "path": slot,
        "source": source,
        "installation": [f"cp -r ${{source}}/* ${{path}}"],
        "initialization": [
            f"cp -r ${{source}}/input.json ${{path}}/input.json",
            "rm -f output.json",
        ],
        "modes": False,
        "make_paraller": 4,
        "execution": {
            "path": slot,
            "commands": ["./eggbox.py"],
            "input": [
                {
                    "name": "inpjson",
                    "path": os.path.join(slot, "input.json"),
                    "type": "JSON",
                    "save": False,
                    "actions": [
                        {
                            "type": "Dump",
                            "variables": [
                                {"name": "xx", "expression": "x * Pi"},
                                {"name": "yy", "expression": "y * Pi"},
                            ],
                        }
                    ],
                }
            ],
            "output": [
                {
                    "name": "oupjson",
                    "path": os.path.join(slot, "output.json"),
                    "type": "JSON",
                    "save": False,
                    "variables": [{"name": "z", "entry": "z"}],
                }
            ],
        },
    }


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
                "make_paraller": 16,
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
        self.assertEqual(spec.raw.get("make_paraller"), 16)

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
    def test_top_level_make_paraller_feeds_pool_resolution(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            worker = build_worker_config(
                {
                    "project_root": tmpdir,
                    "Calculators": {
                        "make_paraller": 16,
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
            self.assertEqual(worker["calculator_make_paraller"], 16)
            self.assertEqual(
                worker["sample_config"].get("calculators_path"),
                os.path.join(tmpdir, "calculators", "runtime", "program"),
            )
            self.assertEqual(
                resolve_calculator_pools(worker),
                {"EggBox": 16},
            )

    def test_module_make_paraller_overrides_global(self) -> None:
        pools = resolve_calculator_pools(
            {
                "calculator_make_paraller": 16,
                "calculator_modules": [{"name": "EggBox", "make_paraller": 3}],
            }
        )
        self.assertEqual(pools, {"EggBox": 3})


@unittest.skipUnless(EGGBOX_SOURCE, "Jarvis-Examples Eggbox assets not available")
class V1StringCloneShadowExecuteTests(unittest.TestCase):
    """Acceptance-shaped unit path: V1 string install/init/commands + clone_shadow."""

    def test_string_install_and_execute_produces_eggbox_z(self) -> None:
        assert EGGBOX_SOURCE is not None
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime_root = os.path.join(tmpdir, "runtime")
            raw_module = _v1_eggbox_module(runtime_root, EGGBOX_SOURCE)
            parser = CommandParser(project_root=tmpdir, scan_name="test")
            prepared = prepare_calculator_modules([raw_module], parser)[0]

            # Spec must keep install tokens until RuntimePreparer expands them.
            spec = CalculatorSpec.from_config("EggBox", prepared)
            self.assertIn("${source}", spec.installation[0]["cmd"])
            self.assertIn("${path}", spec.installation[0]["cmd"])
            self.assertEqual(spec.commands[0]["cmd"], "./eggbox.py")

            module = _attach_runtime(CalculatorModule("EggBox", prepared), tmpdir)
            module.preload_templates()
            module.acquire_pack_id("001")

            sample = Sample.from_params({"x": 0.25, "y": 0.25, "uuid": "v1-string"})
            sample.set_config(
                {
                    "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": False,
                }
            )
            sample.materialize()
            result = module.execute(sample.info)

            expected_z = float((np.sin(0.25 * pi) * np.cos(0.25 * pi) + 2) ** 5)
            self.assertAlmostEqual(float(result["z"]), expected_z, places=6)

            pack_dir = os.path.join(runtime_root, "EggBox", "001")
            self.assertTrue(os.path.isfile(os.path.join(pack_dir, "eggbox.py")))
            self.assertTrue(os.path.isfile(os.path.join(pack_dir, "output.json")))
            # V1: installation log lives under the pack shadow dir, not SAMPLE/.
            install_log = os.path.join(pack_dir, "Installation_EggBox-001.log")
            self.assertTrue(os.path.isfile(install_log), install_log)
            install_text = open(install_log, encoding="utf-8").read()
            self.assertIn("Start install EggBox-001", install_text)
            self.assertIn("Run installation command", install_text)
            self.assertIn("Command Summary -> [install#", install_text)

    def test_sample_running_log_matches_v1_shape(self) -> None:
        """Calculator + Portal + Likelihood must stream into Sample_running.log."""
        assert EGGBOX_SOURCE is not None
        with tempfile.TemporaryDirectory() as tmpdir:
            runtime_root = os.path.join(tmpdir, "runtime")
            raw_module = _v1_eggbox_module(runtime_root, EGGBOX_SOURCE)
            # Nested Dump entry like the Eggbox process card.
            raw_module["execution"]["input"][0]["actions"][0]["variables"].append(
                {"name": "cx", "expression": "(x + y) * Pi", "entry": "test.config.x"}
            )
            parser = CommandParser(project_root=tmpdir, scan_name="test")
            prepared = prepare_calculator_modules([raw_module], parser)[0]
            module = _attach_runtime(CalculatorModule("EggBox", prepared), tmpdir)
            module.acquire_pack_id("006")

            sample = Sample.from_params(
                {"x": 0.10026936915321401, "y": 3.557910680770874, "uuid": "log-shape"}
            )
            sample.set_config(
                {
                    "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
                    "task_result_dir": tmpdir,
                    "sample_artifacts": "always",
                    "workflow_has_calculator": True,
                    "workflow_references_sdir": False,
                }
            )
            sample.materialize()
            result = module.execute(sample.info)
            sample.merge_observables(result)
            sample.info["observables"] = dict(sample.observables)
            LogLikelihoodEvaluator(
                [{"name": "LogL_Z", "expression": "LogGauss(z, 100, 10)"}]
            ).calculate(sample.info)
            sample.observables = dict(sample.info["observables"])
            sample.close()

            run_log = str(sample.info.get("run_log") or "")
            self.assertTrue(os.path.isfile(run_log), run_log)
            text = open(run_log, encoding="utf-8").read()
            self.assertIn("Module load instance and logger is correctly set!", text)
            self.assertIn("(EggBox-No.006)", text)
            self.assertIn("Run initialize command", text)
            self.assertIn("Command Summary -> [initialize#", text)
            # Install stage must not pollute Sample_running.log (pack Installation_*.log).
            self.assertNotIn("Run installation command", text)
            self.assertNotIn("Command Summary -> [install#", text)
            self.assertIn("Adding the file inpjson as 'Portal:JSON' type", text)
            self.assertIn("Evaluating: expression", text)
            self.assertIn("Run execution command", text)
            self.assertIn("Command Summary -> [execution#", text)
            self.assertIn("Loading the file oupjson as 'Portal:JSON' type", text)
            self.assertIn("Evaluating   LogL_Z", text)
            self.assertIn("(Likelihood)", text)
            self.assertIn("Sample SUMMARY", text)
            self.assertIn("Sample closed", text)

            pack_install = os.path.join(
                runtime_root, "EggBox", "006", "Installation_EggBox-006.log"
            )
            self.assertTrue(os.path.isfile(pack_install), pack_install)
            self.assertIn(
                "Run installation command",
                open(pack_install, encoding="utf-8").read(),
            )


@unittest.skipUnless(EGGBOX_SOURCE, "Jarvis-Examples Eggbox assets not available")
class V1ProcessCardYamlLoadTests(unittest.TestCase):
    """Load the real Calculator card YAML and verify Spec surface after Phase-1."""

    def test_example_bridson_process_module_normalizes(self) -> None:
        import yaml

        # EGGBOX_SOURCE = <project>/deps/assets/inertial/EggBox
        eggbox_project = os.path.abspath(
            os.path.join(EGGBOX_SOURCE, "..", "..", "..", "..")  # type: ignore[arg-type]
        )
        card = os.path.join(eggbox_project, "bin", "Example_Bridson_process.yaml")
        self.assertTrue(os.path.isfile(card), card)
        with open(card, encoding="utf-8") as handle:
            config = yaml.safe_load(handle)
        modules = (config.get("Calculators") or {}).get("Modules") or []
        parser = CommandParser(
            project_root=eggbox_project, scan_name="EggBox_Bridson_process"
        )
        prepared = prepare_calculator_modules(modules, parser)
        self.assertEqual(len(prepared), 1)

        spec = CalculatorSpec.from_config("EggBox", prepared[0])
        self.assertTrue(spec.clone_shadow)
        self.assertEqual(len(spec.installation), 1)
        self.assertIn("${source}", spec.installation[0]["cmd"])
        self.assertEqual(spec.commands[0]["cmd"], "./eggbox.py")
        self.assertTrue(os.path.isabs(spec.source))
        self.assertTrue(os.path.isdir(spec.source))
        self.assertIn("@PackID", spec.basepath)
        # Top-level Calculators keys survive YAML load for worker_config.
        calc = config["Calculators"]
        self.assertEqual(calc.get("make_paraller"), 16)
        self.assertIn("path", calc)


if __name__ == "__main__":
    unittest.main()
