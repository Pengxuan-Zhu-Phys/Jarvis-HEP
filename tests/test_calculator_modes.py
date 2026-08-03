"""D20 shared-PackID calculator multi-mode tests."""

from __future__ import annotations

import os
import tempfile
import threading
import time
import unittest
from copy import deepcopy
from types import SimpleNamespace

import fakeredis
import numpy as np

from jarvishep2.calculator_modes import (
    expand_calculator_modes,
    validate_calculator_modes,
)
from jarvishep2.calculator_pools import register_calculator_pools, resolve_calculator_pools
from jarvishep2.core import Jarvis2Core
from jarvishep2.Module.calculator_spec import CalculatorSpec
from jarvishep2.Module.runtime_preparer import (
    RuntimePreparer,
    install_stamp_path,
    prepare_install_controls,
    read_install_control,
    read_install_stamp,
    write_install_control,
)
from jarvishep2.redis_queue import (
    RedisQueue,
    calc_shared_free_list_key,
    calc_shared_unassigned_list_key,
)
from jarvishep2.sample import ExecutionStep, Sample
from jarvishep2.task_validation import validate_task_config
from jarvishep2.worker import Worker
from jarvishep2.worker_config import build_worker_config
from jarvishep2.workflow import build_execution_plan, resolve_module_layers

from test_worker_calculator import _start_tcp_fakeredis, _worker_config


def _mode_card() -> dict:
    return {
        "Scan": {"name": "mode-test"},
        "Sampling": {
            "Method": "Random",
            "Bounds": {"Point number": 1},
            "Variables": [{
                "name": "x",
                "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
            }],
        },
        "Calculators": {"Modules": [
            {
                "name": "Prep",
                "path": "runtime/Prep/@PackID",
                "clone_shadow": True,
                "installation": ["parent-install"],
                "initialization": ["parent-init"],
                "modes": [
                    {
                        "name": "fast",
                        "installation": ["build-fast"],
                        "initialization": ["fast-init"],
                        "execution": {
                            "commands": ["fast"],
                            "output": [{"name": "prep-fast", "path": "out-fast.json", "type": "JSON", "variables": []}],
                        },
                    },
                    {
                        "name": "full",
                        "installation": ["build-full"],
                        "execution": {
                            "commands": ["full"],
                            "output": [{"name": "prep-full", "path": "out-full.json", "type": "JSON", "variables": []}],
                        },
                    },
                ],
            },
            {
                "name": "Analysis",
                "path": "runtime/Analysis/@PackID",
                "required_modules": ["Prep"],
                "execution": {"commands": ["analysis"], "output": []},
            },
            {
                "name": "Focused",
                "path": "runtime/Focused/@PackID",
                "required_modules": ["Prep.fast", "Prep.full"],
                "execution": {"commands": ["focused"], "output": []},
            },
        ]},
    }


class CalculatorModesTests(unittest.TestCase):
    def test_single_mode_module_is_left_unchanged(self) -> None:
        single = {
            "name": "Single",
            "path": "runtime/Single/@PackID",
            "make_paraller": 3,
            "execution": {"commands": ["single"], "output": []},
        }
        self.assertEqual(expand_calculator_modes([single]), [single])
        plan = build_execution_plan(calculator_modules=[single], include_likelihood=False)
        self.assertEqual([(step.name, step.layer) for step in plan], [("Single", 0)])
        self.assertEqual(resolve_calculator_pools({"calculator_modules": [single]}), {"Single": 3})

    def test_parent_modes_expand_to_logical_children_with_one_physical_path(self) -> None:
        card = _mode_card()
        report = validate_task_config(card)
        self.assertTrue(report.ok, [item.message for item in report.issues])

        modules = expand_calculator_modes(card["Calculators"]["Modules"])
        self.assertEqual([module["name"] for module in modules], ["Prep.fast", "Prep.full", "Analysis", "Focused"])
        fast, full, analysis, focused = modules
        self.assertEqual(fast["path"], "runtime/Prep/@PackID")
        self.assertEqual(full["path"], "runtime/Prep/@PackID")
        self.assertEqual(fast["_mode_parent"], "Prep")
        self.assertEqual(fast["_mode_name"], "fast")
        self.assertTrue(fast["_shared_mode_pack"])
        self.assertEqual(fast["installation"], ["parent-install", "build-fast"])
        self.assertEqual(fast["initialization"], ["parent-init", "fast-init"])
        self.assertEqual(full["initialization"], ["parent-init"])
        # A bare multi-mode parent is an all-modes dependency.
        self.assertEqual(analysis["required_modules"], ["Prep.fast", "Prep.full"])
        self.assertEqual(focused["required_modules"], ["Prep.fast", "Prep.full"])

    def test_expanded_modes_use_existing_layers_but_one_parent_pool(self) -> None:
        modules = expand_calculator_modes(_mode_card()["Calculators"]["Modules"])
        layers = resolve_module_layers(modules)
        self.assertEqual(layers["Prep.fast"], 0)
        self.assertEqual(layers["Prep.full"], 0)
        self.assertEqual(layers["Analysis"], 1)
        self.assertEqual(layers["Focused"], 1)
        plan = build_execution_plan(calculator_modules=_mode_card()["Calculators"]["Modules"], include_likelihood=False)
        self.assertEqual(
            [(step.name, step.layer) for step in plan],
            [("Prep.fast", 0), ("Prep.full", 0), ("Analysis", 1), ("Focused", 1)],
        )
        self.assertEqual(
            resolve_calculator_pools({"calculator_modules": modules}),
            {"Prep": 1, "Analysis": 1, "Focused": 1},
        )

    def test_worker_blueprint_expands_modes_but_preserves_parent_pool(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            card = _mode_card()
            card["task_result_dir"] = temporary
            card["project_root"] = temporary
            worker_config = build_worker_config(card, task_result_dir=temporary)
            self.assertEqual(
                [entry["name"] for entry in worker_config["calculator_modules"]],
                ["Prep.fast", "Prep.full", "Analysis", "Focused"],
            )
            self.assertEqual(
                resolve_calculator_pools(worker_config),
                {"Prep": 1, "Analysis": 1, "Focused": 1},
            )

    def test_unknown_qualified_mode_has_a_fixable_diagnostic(self) -> None:
        card = _mode_card()
        card["Calculators"]["Modules"][1]["required_modules"] = ["Prep.fats"]
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-MOD-005"]
        self.assertEqual(len(errors), 1)
        self.assertIn("Did you mean 'Prep.fast'?", errors[0].message)
        self.assertIn("Parent.mode", errors[0].hint or "")

    def test_bare_cross_block_dependencies_are_not_owned_by_mode_validation(self) -> None:
        card = _mode_card()
        card["Calculators"]["Modules"][1]["required_modules"] = [
            "LoopTools", "BuildBSMPTInput", "Parameters",
        ]
        errors = [
            item for item in validate_task_config(card).errors()
            if item.code == "JV2-MOD-005"
        ]
        self.assertEqual(errors, [])

    def test_bare_parent_dependency_excludes_the_expanded_child_itself(self) -> None:
        parent = deepcopy(_mode_card()["Calculators"]["Modules"][0])
        parent["modes"][1]["required_modules"] = ["Prep"]
        expanded = expand_calculator_modes([parent])
        self.assertEqual(expanded[1]["required_modules"], ["Prep.fast"])

    def test_parent_execution_and_duplicate_mode_output_are_rejected(self) -> None:
        card = _mode_card()
        parent = card["Calculators"]["Modules"][0]
        parent["execution"] = {"commands": ["wrong"]}
        parent["modes"][1]["execution"]["output"][0]["name"] = "prep-fast"
        codes = {item.code for item in validate_task_config(card).errors()}
        self.assertIn("JV2-MOD-002", codes)
        self.assertIn("JV2-MOD-004", codes)

    def test_mode_execution_receives_portal_io_schema_validation(self) -> None:
        card = _mode_card()
        card["Calculators"]["Modules"][0]["modes"][0]["execution"]["output"][0]["type"] = "NoSuchFormat"
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-SCH-002"]
        self.assertEqual(len(errors), 1)
        self.assertEqual(errors[0].path, "$.Calculators.Modules[0].modes[0].execution.output[0].type")

    def test_per_mode_directory_tokens_are_rejected(self) -> None:
        card = _mode_card()
        parent = card["Calculators"]["Modules"][0]
        parent["path"] = "runtime/Prep/@Mode/@PackID"
        parent["modes"][0]["installation"] = ["build ${mode_dir}"]
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-MOD-006"]
        self.assertEqual(len(errors), 2)
        self.assertTrue(any("@Mode" in item.message for item in errors))
        self.assertTrue(any("${mode_dir}" in item.message for item in errors))
        self.assertTrue(all("share" in (item.hint or "") for item in errors))

    def test_pools_name_the_parent_not_a_logical_mode(self) -> None:
        card = _mode_card()
        card["Calculators"]["Pools"] = {"Prep": 2}
        self.assertTrue(validate_task_config(card).ok)
        card["Calculators"]["Pools"] = {"Prep.fast": 2}
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-MOD-007"]
        self.assertEqual(len(errors), 1)
        self.assertIn("Calculators.Pools.Prep", errors[0].suggestion or "")

    def test_mode_execution_paths_are_allowed_but_physical_settings_are_not(self) -> None:
        card = _mode_card()
        mode = card["Calculators"]["Modules"][0]["modes"][0]
        mode["clone_shadow"] = False
        mode["source"] = "other-source"
        mode["path"] = "mode-runtime"
        mode["execution"]["path"] = "execution-runtime"
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-MOD-006"]
        self.assertEqual(
            {item.path for item in errors},
            {
                "Calculators.Modules[0].modes[0].clone_shadow",
                "Calculators.Modules[0].modes[0].source",
            },
        )

        expanded = expand_calculator_modes([card["Calculators"]["Modules"][0]])[0]
        self.assertTrue(expanded["clone_shadow"])
        self.assertNotIn("source", expanded)
        self.assertEqual(expanded["path"], "runtime/Prep/@PackID")
        self.assertEqual(expanded["execution"]["path"], "execution-runtime")

    def test_mode_execution_path_precedence(self) -> None:
        parent = _mode_card()["Calculators"]["Modules"][0]
        fast, full = parent["modes"]
        fast["path"] = "mode-fast/@PackID"
        fast["execution"]["path"] = "execution-fast/@PackID"
        full["path"] = "mode-full/@PackID"
        expanded_fast, expanded_full = expand_calculator_modes([parent])
        self.assertEqual(expanded_fast["path"], "runtime/Prep/@PackID")
        self.assertEqual(expanded_fast["execution"]["path"], "execution-fast/@PackID")
        self.assertEqual(expanded_full["execution"]["path"], "mode-full/@PackID")
        self.assertEqual(
            CalculatorSpec.from_config("Prep.fast", expanded_fast).commands[0]["cwd"],
            "execution-fast/@PackID",
        )
        self.assertEqual(
            CalculatorSpec.from_config("Prep.full", expanded_full).commands[0]["cwd"],
            "mode-full/@PackID",
        )

        full.pop("path")
        fallback = expand_calculator_modes([parent])[1]
        self.assertNotIn("path", fallback["execution"])
        self.assertEqual(
            CalculatorSpec.from_config("Prep.full", fallback).commands[0]["cwd"],
            "runtime/Prep/@PackID",
        )

    def test_unknown_pool_name_is_rejected_with_a_suggestion(self) -> None:
        card = _mode_card()
        card["Calculators"]["Pools"] = {"Prepp": 2}
        errors = [item for item in validate_task_config(card).errors() if item.code == "JV2-MOD-008"]
        self.assertEqual(len(errors), 1)
        self.assertIn("Did you mean 'Prep'?", errors[0].message)

    def test_small_shared_pool_emits_actionable_performance_warning(self) -> None:
        card = _mode_card()
        card["Calculators"]["Pools"] = {"Prep": 3}
        warnings = [
            item for item in validate_calculator_modes(
                card, resolved_config={"Runtime": {"workers": 3}},
            )
            if item.code == "JV2-MOD-009"
        ]
        self.assertEqual(len(warnings), 1)
        self.assertEqual(warnings[0].level, "warning")
        self.assertIn("approximately 5", warnings[0].message)
        self.assertIn("Calculators.Pools.Prep: 5", warnings[0].suggestion or "")

        card["Calculators"]["Pools"] = {"Prep": 5}
        self.assertFalse(any(
            item.code == "JV2-MOD-009"
            for item in validate_calculator_modes(
                card, resolved_config={"Runtime": {"workers": 3}},
            )
        ))

    def test_core_load_expands_only_after_validating_raw_card(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            task_path = os.path.join(temporary, "modes.yaml")
            with open(task_path, "w", encoding="ascii") as handle:
                handle.write(
                    "Scan:\n  name: modes\n"
                    "Sampling:\n  Method: Random\n  Point number: 1\n"
                    "  Variables:\n    - name: x\n      distribution:\n        type: Flat\n        parameters: {min: 0.0, max: 1.0}\n"
                    "Calculators:\n  Modules:\n    - name: Tool\n      path: runtime/Tool/@PackID\n      clone_shadow: true\n      modes:\n        - name: a\n          execution: {commands: [a], output: []}\n        - name: b\n          execution: {commands: [b], output: []}\n"
                )
            core = Jarvis2Core()
            core.load_task_yaml(task_path)
            self.assertEqual([entry["name"] for entry in core.config["Calculators"]["Modules"]], ["Tool.a", "Tool.b"])
            self.assertEqual(core.config["Calculators"]["Modules"][0]["path"], "runtime/Tool/@PackID")
            self.assertEqual(core.config.raw_task_card["Calculators"]["Modules"][0]["modes"][0]["name"], "a")

    def test_shared_pool_reuses_one_pack_and_records_mode_affinity(self) -> None:
        modules = expand_calculator_modes(_mode_card()["Calculators"]["Modules"])
        client = fakeredis.FakeRedis(decode_responses=True)
        queue = RedisQueue(client=client)
        register_calculator_pools(queue, {"calculator_modules": modules})
        self.assertEqual(client.llen(calc_shared_unassigned_list_key("Prep")), 1)
        self.assertEqual(client.exists("calc:free:Prep.fast"), 0)

        first = queue.acquire_shared_calc("Prep", "fast", modes=["fast", "full"], timeout=1)
        self.assertEqual(first, ("001", None))
        self.assertEqual(
            queue.shared_mode_busy_counts("Prep", ["fast", "full"]),
            {"fast": 1, "full": 0},
        )
        self.assertTrue(queue.release_shared_calc("Prep", "001", "fast"))
        self.assertEqual(client.lrange(calc_shared_free_list_key("Prep", "fast"), 0, -1), ["001"])

        warm = queue.acquire_shared_calc("Prep", "fast", modes=["fast", "full"], timeout=1)
        self.assertEqual(warm, ("001", "fast"))
        self.assertTrue(queue.release_shared_calc("Prep", "001", "fast"))
        borrowed = queue.acquire_shared_calc("Prep", "full", modes=["fast", "full"], timeout=1)
        self.assertEqual(borrowed, ("001", "fast"))
        self.assertTrue(queue.release_shared_calc("Prep", "001", "full"))
        self.assertEqual(client.hget("calc:packmode:Prep", "001"), "full")

    def test_shared_acquire_waits_for_busy_warm_pack_before_borrowing(self) -> None:
        client = fakeredis.FakeRedis(decode_responses=True)
        queue = RedisQueue(client=client)
        queue.register_shared_calc_pool("Prep", 2, modes=["fast", "full"])

        fast = queue.acquire_shared_calc(
            "Prep", "fast", modes=["fast", "full"], timeout=1,
        )
        self.assertEqual(fast, ("001", None))
        self.assertTrue(queue.release_shared_calc("Prep", "001", "fast"))
        full = queue.acquire_shared_calc(
            "Prep", "full", modes=["fast", "full"], timeout=1,
        )
        self.assertEqual(full, ("002", None))
        self.assertTrue(queue.release_shared_calc("Prep", "002", "full"))

        held = queue.acquire_shared_calc(
            "Prep", "fast", modes=["fast", "full"], timeout=1,
        )
        self.assertEqual(held, ("001", "fast"))

        releaser = threading.Thread(
            target=lambda: (
                time.sleep(0.05),
                queue.release_shared_calc("Prep", "001", "fast"),
            ),
            daemon=True,
        )
        releaser.start()
        started = time.monotonic()
        acquired = queue.acquire_shared_calc(
            "Prep", "fast", modes=["fast", "full"], timeout=1,
            affinity_wait_sec=0.5,
        )
        elapsed = time.monotonic() - started
        releaser.join(timeout=1)
        self.assertEqual(acquired, ("001", "fast"))
        self.assertGreaterEqual(elapsed, 0.03)
        self.assertTrue(queue.release_shared_calc("Prep", "001", "fast"))

    def test_shared_acquire_borrows_immediately_when_no_target_pack_is_busy(self) -> None:
        client = fakeredis.FakeRedis(decode_responses=True)
        queue = RedisQueue(client=client)
        queue.register_shared_calc_pool(
            "Prep", 1, modes=["fast", "full"], pack_modes={"001": "fast"},
        )
        started = time.monotonic()
        acquired = queue.acquire_shared_calc(
            "Prep", "full", modes=["fast", "full"], timeout=1,
            affinity_wait_sec=0.5,
        )
        self.assertEqual(acquired, ("001", "fast"))
        self.assertLess(time.monotonic() - started, 0.2)
        self.assertTrue(queue.release_shared_calc("Prep", "001", "full"))

    def test_three_workers_preserve_three_mode_affinity_under_contention(self) -> None:
        modes = ["a", "b", "c"]
        client = fakeredis.FakeRedis(decode_responses=True)
        queue = RedisQueue(client=client)
        queue.register_shared_calc_pool(
            "Tool",
            3,
            modes=modes,
            pack_modes={"001": "a", "002": "b", "003": "c"},
        )
        modules = {
            f"Tool.{mode}": SimpleNamespace(
                config={
                    "name": f"Tool.{mode}",
                    "_mode_parent": "Tool",
                    "_mode_name": mode,
                    "_shared_mode_pack": True,
                },
                should_run=lambda sample_info: True,
            )
            for mode in modes
        }
        steps = [
            ExecutionStep(type="calculator", name=f"Tool.{mode}", layer=0)
            for mode in modes
        ]
        barrier = threading.Barrier(3)
        rebuilds = 0
        rebuild_lock = threading.Lock()

        def run_worker(worker_id: int) -> None:
            nonlocal rebuilds
            worker = Worker(worker_id, {"host": "127.0.0.1", "port": 6379}, {})
            worker._redis = queue
            worker._calculators = modules

            def run_step(name, sample, selection_checked=False):
                nonlocal rebuilds
                target = name.split(".", 1)[1]
                acquired = queue.acquire_shared_calc(
                    "Tool", target, modes=modes, timeout=2,
                    affinity_wait_sec=0.25,
                )
                self.assertIsNotNone(acquired)
                pack_id, current_mode = acquired or ("", None)
                if current_mode != target:
                    with rebuild_lock:
                        rebuilds += 1
                time.sleep(0.02)
                self.assertTrue(queue.release_shared_calc("Tool", pack_id, target))

            worker._run_calculator_step = run_step  # type: ignore[method-assign]
            barrier.wait(timeout=2)
            for sample_index in range(6):
                worker._run_shared_mode_group(
                    steps,
                    Sample(
                        uuid=f"contention-{worker_id}-{sample_index}",
                        u_coords=np.array([0.0]),
                    ),
                )

        threads = [
            threading.Thread(target=run_worker, args=(index,), daemon=True)
            for index in range(3)
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join(timeout=5)
        self.assertTrue(all(not thread.is_alive() for thread in threads))
        self.assertEqual(rebuilds, 0)

    def test_mode_stamp_requires_rebuild_only_when_the_mode_changes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            card = _mode_card()
            parent = card["Calculators"]["Modules"][0]
            parent["path"] = os.path.join(temporary, "Prep", "@PackID")
            fast_config, full_config = expand_calculator_modes([parent])
            fast = RuntimePreparer(CalculatorSpec.from_config("Prep.fast", fast_config))
            full = RuntimePreparer(CalculatorSpec.from_config("Prep.full", full_config))
            fast.acquire_pack_id("001")
            full.acquire_pack_id("001")
            calls: list[tuple[str, list[dict[str, object]]]] = []

            def run_stage(commands, stage):
                calls.append((stage, [dict(command) for command in commands]))

            fast.ensure_shadow_installed(run_stage=run_stage)
            runtime = os.path.join(temporary, "Prep", "001")
            stamp = read_install_stamp(runtime)
            self.assertEqual(stamp["extra"]["mode"], "fast")
            self.assertEqual([command["cmd"] for command in calls[0][1]], ["parent-install"])
            self.assertEqual([command["cmd"] for command in calls[1][1]], ["build-fast"])
            fast.ensure_shadow_installed(run_stage=run_stage)
            self.assertEqual(len(calls), 2)

            stamp_missing_during_rebuild: list[bool] = []

            def rebuild(commands, stage):
                stamp_missing_during_rebuild.append(not os.path.exists(install_stamp_path(runtime)))
                calls.append((stage, [dict(command) for command in commands]))

            full.ensure_shadow_installed(run_stage=rebuild)
            stamp = read_install_stamp(runtime)
            self.assertEqual(stamp["extra"]["mode"], "full")
            self.assertEqual(len(calls), 3)
            self.assertEqual([command["cmd"] for command in calls[2][1]], ["build-full"])
            self.assertEqual(stamp_missing_during_rebuild, [True])

            forced = RuntimePreparer(
                CalculatorSpec.from_config("Prep.full", full_config),
                force_reinstall=True,
            )
            forced.acquire_pack_id("001")
            forced.ensure_shadow_installed(run_stage=run_stage)
            self.assertEqual([command["cmd"] for command in calls[3][1]], ["parent-install"])
            self.assertEqual([command["cmd"] for command in calls[4][1]], ["build-full"])

            # A resumed run restores affinity from the successful mode stamp.
            queue = RedisQueue(client=fakeredis.FakeRedis(decode_responses=True))
            register_calculator_pools(queue, {"calculator_modules": [fast_config, full_config]})
            self.assertEqual(
                queue.r.lrange(calc_shared_free_list_key("Prep", "full"), 0, -1),
                ["001"],
            )

    def test_modes_share_one_install_control_and_one_reinstall_epoch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            card = _mode_card()
            parent = card["Calculators"]["Modules"][0]
            parent["path"] = os.path.join(temporary, "Prep", "@PackID")
            modules = expand_calculator_modes([parent])
            first = prepare_install_controls(modules, env={})
            self.assertEqual([item["_install_epoch"] for item in first], [0, 0])
            spec = CalculatorSpec.from_config("Prep.fast", first[0])
            control = read_install_control(spec)
            self.assertEqual(control["module"], "Prep")
            control["reinstall"] = True
            write_install_control(spec, control)

            second = prepare_install_controls(modules, env={})
            self.assertEqual([item["_install_epoch"] for item in second], [1, 1])
            self.assertFalse(read_install_control(spec)["reinstall"])

    def test_worker_greedily_runs_the_pending_warm_mode_first(self) -> None:
        modules = expand_calculator_modes([_mode_card()["Calculators"]["Modules"][0]])

        class Counts:
            def shared_mode_free_counts(self, parent, modes):
                self.parent = parent
                self.modes = list(modes)
                return {"fast": 0, "full": 2}

        worker = Worker(0, {"host": "127.0.0.1", "port": 6379}, {})
        worker._redis = Counts()
        worker._calculators = {
            item["name"]: SimpleNamespace(
                config=item, should_run=lambda sample_info: True,
            )
            for item in modules
        }
        ran: list[str] = []
        worker._run_calculator_step = (  # type: ignore[method-assign]
            lambda name, sample, selection_checked=False: ran.append(name)
        )
        worker._run_shared_mode_group(
            [
                ExecutionStep(type="calculator", name="Prep.fast", layer=0),
                ExecutionStep(type="calculator", name="Prep.full", layer=0),
            ],
            Sample(uuid="greedy", u_coords=np.array([0.0])),
        )
        self.assertEqual(ran, ["Prep.full", "Prep.fast"])
        self.assertEqual(worker._redis.parent, "Prep")

    def test_worker_spreads_cold_builds_across_modes_in_use(self) -> None:
        modules = expand_calculator_modes([_mode_card()["Calculators"]["Modules"][0]])

        class ColdCounts:
            def shared_mode_free_counts(self, parent, modes):
                return {mode: 0 for mode in modes}

            def shared_mode_busy_counts(self, parent, modes):
                return {"fast": 1, "full": 0}

        worker = Worker(0, {"host": "127.0.0.1", "port": 6379}, {})
        worker._redis = ColdCounts()
        worker._calculators = {
            item["name"]: SimpleNamespace(
                config=item, should_run=lambda sample_info: True,
            )
            for item in modules
        }
        ran: list[str] = []
        worker._run_calculator_step = (  # type: ignore[method-assign]
            lambda name, sample, selection_checked=False: ran.append(name)
        )
        worker._run_shared_mode_group(
            [
                ExecutionStep(type="calculator", name="Prep.fast", layer=0),
                ExecutionStep(type="calculator", name="Prep.full", layer=0),
            ],
            Sample(uuid="cold-balance", u_coords=np.array([0.0])),
        )
        self.assertEqual(ran, ["Prep.full", "Prep.fast"])

    def test_selection_skips_mode_before_redis_acquire_and_greedy_ordering(self) -> None:
        modules = expand_calculator_modes([_mode_card()["Calculators"]["Modules"][0]])

        class Counts:
            def __init__(self) -> None:
                self.queries: list[list[str]] = []

            def shared_mode_free_counts(self, parent, modes):
                self.queries.append(list(modes))
                return {mode: 0 for mode in modes}

            def shared_mode_busy_counts(self, parent, modes):
                return {mode: 0 for mode in modes}

        worker = Worker(0, {"host": "127.0.0.1", "port": 6379}, {})
        worker._redis = Counts()
        worker._calculators = {
            "Prep.fast": SimpleNamespace(
                config=modules[0], should_run=lambda sample_info: False,
            ),
            "Prep.full": SimpleNamespace(
                config=modules[1], should_run=lambda sample_info: True,
            ),
        }
        ran: list[tuple[str, bool]] = []
        worker._run_calculator_step = (  # type: ignore[method-assign]
            lambda name, sample, selection_checked=False:
            ran.append((name, selection_checked))
        )
        worker._run_shared_mode_group(
            [
                ExecutionStep(type="calculator", name="Prep.fast", layer=0),
                ExecutionStep(type="calculator", name="Prep.full", layer=0),
            ],
            Sample(uuid="selection", u_coords=np.array([0.0])),
        )
        self.assertEqual(ran, [("Prep.full", True)])
        self.assertEqual(worker._redis.queries, [["full"]])

    def test_selection_false_step_never_calls_redis_acquire(self) -> None:
        class RedisMustNotBeUsed:
            def __getattr__(self, name):
                raise AssertionError(f"Redis method {name} must not be called")

        worker = Worker(0, {"host": "127.0.0.1", "port": 6379}, {})
        worker._redis = RedisMustNotBeUsed()
        worker._calculators = {
            "Prep.fast": SimpleNamespace(
                config={"name": "Prep.fast"},
                should_run=lambda sample_info: False,
            )
        }
        worker._run_calculator_step(
            "Prep.fast", Sample(uuid="skip", u_coords=np.array([0.0])),
        )

    def test_spawn_worker_switches_one_shared_runtime_and_prefers_warm_mode(self) -> None:
        """End-to-end: two logical modes share ``Slow/001`` and rebuild safely."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as temporary:
                base_script = os.path.join(temporary, "build_base.py")
                build_script = os.path.join(temporary, "build_mode.py")
                run_script = os.path.join(temporary, "run_mode.py")
                with open(base_script, "w", encoding="ascii") as handle:
                    handle.write(
                        "from pathlib import Path\nimport sys\n"
                        "runtime = Path(sys.argv[1])\n"
                        "runtime.mkdir(parents=True, exist_ok=True)\n"
                        "runtime.joinpath('base.txt').write_text('ready', encoding='ascii')\n"
                        "with runtime.parent.joinpath('build.log').open('a', encoding='ascii') as out: out.write('base\\n')\n"
                    )
                with open(build_script, "w", encoding="ascii") as handle:
                    handle.write(
                        "from pathlib import Path\nimport sys\n"
                        "runtime = Path(sys.argv[1])\nmode = sys.argv[2]\n"
                        "runtime.mkdir(parents=True, exist_ok=True)\n"
                        "assert runtime.joinpath('base.txt').read_text(encoding='ascii') == 'ready'\n"
                        "runtime.joinpath('mode.txt').write_text(mode, encoding='ascii')\n"
                        "with runtime.parent.joinpath('build.log').open('a', encoding='ascii') as out: out.write(mode + '\\n')\n"
                    )
                with open(run_script, "w", encoding="ascii") as handle:
                    handle.write(
                        "import json\nfrom pathlib import Path\nimport sys\n"
                        "output, key, value = sys.argv[1], sys.argv[2], int(sys.argv[3])\n"
                        "assert Path('mode.txt').read_text(encoding='ascii') == key\n"
                        "Path(output).parent.mkdir(parents=True, exist_ok=True)\n"
                        "Path(output).write_text(json.dumps({key: value}), encoding='ascii')\n"
                    )
                runtime_path = os.path.join(temporary, "runtime", "Slow", "@PackID")
                modes = []
                for name, value in (("a", 1), ("b", 2)):
                    output_path = f"@Sdir/mode_{name}.json"
                    modes.append({
                        "name": name,
                        "installation": [{"cmd": f"python3 {build_script} ${{path}} {name}", "cwd": runtime_path}],
                        "execution": {
                            "commands": [{"cmd": f"python3 {run_script} {output_path} {name} {value}", "cwd": runtime_path}],
                            "output": [{
                                "name": f"mode_{name}", "path": output_path,
                                "type": "JSON", "save": False,
                                "variables": [{"name": name, "entry": name}],
                            }],
                        },
                    })
                parent = {
                    "name": "Slow", "path": runtime_path, "clone_shadow": True,
                    "installation": [{"cmd": f"python3 {base_script} ${{path}}", "cwd": runtime_path}],
                    "modes": modes,
                }
                core = Jarvis2Core({"Runtime": {"mode": "redis", "workers": 1, "redis": redis_config}, "task_result_dir": temporary})
                core.init_redis()
                worker_config = _worker_config(temporary)
                worker_config.update({
                    "calculator_modules": [deepcopy(parent)],
                    "calculator_pools": {"Slow": 1},
                    "likelihood_expressions": [],
                    "handoff_to_staging": False,
                })
                core.init_factory(worker_config)
                assert core.redis is not None
                plan = build_execution_plan(calculator_modules=[parent], include_likelihood=False)
                self.assertEqual([step.name for step in plan], ["Slow.a", "Slow.b"])

                for uuid in ("mode-worker-1", "mode-worker-2"):
                    core.redis.push_task(Sample(uuid=uuid, u_coords=np.array([0.0, 0.0]), execution_plan=plan).to_task_dict())
                    deadline = time.monotonic() + 30.0
                    result = None
                    while time.monotonic() < deadline:
                        result = core.redis.pull_result(timeout=1)
                        if result is not None:
                            break
                    self.assertIsNotNone(result)
                    assert result is not None
                    self.assertEqual(result["status"], "Completed")
                    self.assertEqual(result["observables"]["a"], 1)
                    self.assertEqual(result["observables"]["b"], 2)
                core.shutdown()

                runtime = os.path.join(temporary, "runtime", "Slow", "001")
                self.assertTrue(os.path.isdir(runtime))
                self.assertFalse(os.path.exists(os.path.join(temporary, "runtime", "Slow", "a")))
                self.assertFalse(os.path.exists(os.path.join(temporary, "runtime", "Slow", "b")))
                with open(os.path.join(temporary, "runtime", "Slow", "build.log"), encoding="ascii") as handle:
                    self.assertEqual(handle.read().splitlines(), ["base", "a", "b", "a"])
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
