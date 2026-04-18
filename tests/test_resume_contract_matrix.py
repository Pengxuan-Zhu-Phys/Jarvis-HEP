#!/usr/bin/env python3
from __future__ import annotations

import os
import contextlib
import io
import tempfile
import unittest
from copy import deepcopy
from types import SimpleNamespace
from unittest import mock

from loguru import logger

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

from jarvishep.core import Core  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402
from jarvishep.moduleManager import ModuleManager  # noqa: E402
from jarvishep.Sampling.Source.MCMC.runtime_checkpoint import (  # noqa: E402
    StateSaver,
    build_state_payload,
    validate_state_payload,
)


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _MemorySink:
    def __init__(self) -> None:
        self.messages: list[str] = []

    def __call__(self, message):
        self.messages.append(message.record["message"])


class _FakePool:
    def __init__(self, name: str, pool_type: str) -> None:
        self.name = name
        self.type = pool_type
        self.restored: list[dict] = []

    def restore_blueprint(self, blueprint: dict) -> None:
        self.restored.append(deepcopy(blueprint))


class _TTYStringIO(io.StringIO):
    def isatty(self) -> bool:
        return True


def _minimal_config(method: str) -> dict:
    return {
        "Scan": {"sample_directory": {"limit": 4, "width": 2}},
        "Sampling": {
            "Method": method,
            "Variables": [
                {
                    "name": "x",
                    "description": "x",
                    "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                }
            ],
            "Bounds": {"num_chains": 1, "num_iters": 1, "proposal_scale": 0.1},
        },
    }


class TestResumeContractMatrix(unittest.TestCase):
    def test_distributor_resume_status_matrix_is_explicit(self):
        statuses = Distributor.list_resume_statuses()
        self.assertTrue(statuses)
        self.assertTrue(all(status in {"implemented", "intentionally unsupported"} for status in statuses.values()))
        self.assertTrue(all(status == "implemented" for status in statuses.values()))
        self.assertEqual(Distributor.get_resume_status("MCMC"), "implemented")
        self.assertEqual(Distributor.get_resume_status("Dynesty"), "implemented")
        self.assertEqual(Distributor.get_resume_status("MultiNest"), "implemented")

    def test_validate_state_payload_rejects_missing_required_sections(self):
        payload = {
            "format": "jarvis-hep.statesaver",
            "version": 1,
            "sampler_state": {},
        }
        ok, reason = validate_state_payload(payload)
        self.assertFalse(ok)
        self.assertIn("run_spec", reason)

        payload = {
            "format": "jarvis-hep.statesaver",
            "version": 1,
            "run_spec": {},
            "sampler_state": {},
        }
        ok, reason = validate_state_payload(payload)
        self.assertFalse(ok)
        self.assertIn("factory_blueprint", reason)

    def test_state_saver_rejects_corrupted_checkpoint_payload(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-corrupt-") as tmp:
            checkpoint = os.path.join(tmp, "state.pkl")
            with open(checkpoint, "wb") as handle:
                handle.write(b"not-a-pickle")

            saver = StateSaver(checkpoint, logger=_NoopLogger())
            with self.assertRaisesRegex(ValueError, "StateSaver checkpoint load failed"):
                saver.load()

    def test_core_resume_warns_on_yaml_drift_and_uses_frozen_checkpoint_spec(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-core-resume-") as tmp:
            scan_name = "resume-scan"
            checkpoint_root = os.path.join(tmp, "checkpoints", scan_name, "Grid")
            os.makedirs(checkpoint_root, exist_ok=True)

            frozen_cfg = _minimal_config("Grid")
            mutated_cfg = _minimal_config("Grid")
            mutated_cfg["Sampling"]["Bounds"]["num_iters"] = 99

            payload = build_state_payload(
                run_spec={
                    "raw_yaml_text": "frozen",
                    "normalized_config": deepcopy(frozen_cfg),
                    "scan_name": scan_name,
                    "task_root": tmp,
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "logs_dir": os.path.join(tmp, "LOGS"),
                    "images_dir": os.path.join(tmp, "IMAGES"),
                    "worker_parallel": 1,
                    "sampler_method": "Grid",
                    "workflow": {},
                    "workflow_layers": {},
                },
                factory_blueprint={"max_workers": 1, "workflow": {}, "calc_layer": {}, "module_pools": {}},
                sampler_state={
                    "sampler_signature": {
                        "method": "Grid",
                        "dimensions": 1,
                        "nchains": 1,
                        "niters": 1,
                        "variable_names": ["x"],
                        "task_result_dir": os.path.join(tmp, "RESULTS"),
                    },
                    "state_machine": {"state": "INIT", "ready_queue": [], "chains": [], "extras": {}},
                },
                reason="unit-test",
            )
            StateSaver(os.path.join(checkpoint_root, "state.pkl"), logger=_NoopLogger()).save(payload)

            core = Core()
            core.args = SimpleNamespace(debug=False, plot=False, skipFC=True, skiplibrary=True)
            core.path = {"task_root": tmp, "jpath": tmp, "logo": tmp, "args_info": tmp}
            core.info = {
                "scan_name": scan_name,
                "project_name": scan_name,
                "config_file": os.path.join(tmp, "current.yaml"),
                "sample": {
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "sample_dirs": os.path.join(tmp, "RESULTS", "SAMPLE"),
                    "archive_samples": False,
                },
                "logs_dir": os.path.join(tmp, "LOGS"),
                "images_dir": os.path.join(tmp, "IMAGES"),
            }
            core.yaml.config = deepcopy(mutated_cfg)
            core.yaml.get_sampling_method = lambda: "Grid"
            core.yaml.get_worker_parallel = lambda: 1
            core.sampler = SimpleNamespace(method="Grid")

            sink = _MemorySink()
            sink_id = logger.add(sink, level="WARNING")
            try:
                core._preload_resume_checkpoint()
            finally:
                logger.remove(sink_id)

            self.assertIsInstance(core._resume_checkpoint_payload, dict)
            self.assertEqual(core._resume_run_spec["normalized_config"], frozen_cfg)
            self.assertTrue(
                any("config drift" in message.lower() for message in sink.messages),
                msg=f"warning messages did not mention config drift: {sink.messages}",
            )

    def test_core_prompts_for_checkpoint_and_can_restart_fresh(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-core-fresh-") as tmp:
            scan_name = "fresh-scan"
            checkpoint_root = os.path.join(tmp, "checkpoints", scan_name, "Grid")
            os.makedirs(checkpoint_root, exist_ok=True)
            checkpoint_file = os.path.join(checkpoint_root, "state.pkl")

            payload = build_state_payload(
                run_spec={
                    "raw_yaml_text": "frozen",
                    "normalized_config": deepcopy(_minimal_config("Grid")),
                    "scan_name": scan_name,
                    "task_root": tmp,
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "logs_dir": os.path.join(tmp, "LOGS"),
                    "images_dir": os.path.join(tmp, "IMAGES"),
                    "worker_parallel": 1,
                    "sampler_method": "Grid",
                    "workflow": {},
                    "workflow_layers": {},
                },
                factory_blueprint={"max_workers": 1, "workflow": {}, "calc_layer": {}, "module_pools": {}},
                sampler_state={
                    "sampler_signature": {
                        "method": "Grid",
                        "dimensions": 1,
                        "nchains": 1,
                        "niters": 1,
                        "variable_names": ["x"],
                        "task_result_dir": os.path.join(tmp, "RESULTS"),
                    },
                    "state_machine": {"state": "INIT", "ready_queue": [], "chains": [], "extras": {}},
                },
                reason="unit-test",
            )
            StateSaver(checkpoint_file, logger=_NoopLogger()).save(payload)

            core = Core()
            core.args = SimpleNamespace(debug=False, plot=False, skipFC=True, skiplibrary=True, resume=False)
            core.path = {"task_root": tmp, "jpath": tmp, "logo": tmp, "args_info": tmp}
            core.info = {
                "scan_name": scan_name,
                "project_name": scan_name,
                "config_file": os.path.join(tmp, "current.yaml"),
                "sample": {
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "sample_dirs": os.path.join(tmp, "RESULTS", "SAMPLE"),
                    "archive_samples": False,
                },
                "logs_dir": os.path.join(tmp, "LOGS"),
                "images_dir": os.path.join(tmp, "IMAGES"),
            }
            core.yaml.config = deepcopy(_minimal_config("Grid"))
            core.yaml.get_sampling_method = lambda: "Grid"
            core.yaml.get_worker_parallel = lambda: 1
            core.sampler = SimpleNamespace(method="Grid")

            with mock.patch.object(core, "_prompt_resume_from_checkpoint", return_value=False) as prompt:
                core._preload_resume_checkpoint()

            prompt.assert_called_once()
            self.assertIsNone(core._resume_checkpoint_payload)
            self.assertFalse(os.path.exists(checkpoint_file))
            self.assertEqual(core._resume_checkpoint_policy, "fresh")

    def test_check_modules_mode_ignores_existing_checkpoint_by_default(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-core-opc-") as tmp:
            scan_name = "opc-scan"
            checkpoint_root = os.path.join(tmp, "checkpoints", scan_name, "Grid")
            os.makedirs(checkpoint_root, exist_ok=True)
            checkpoint_file = os.path.join(checkpoint_root, "state.pkl")

            payload = build_state_payload(
                run_spec={
                    "raw_yaml_text": "frozen",
                    "normalized_config": deepcopy(_minimal_config("Grid")),
                    "scan_name": scan_name,
                    "task_root": tmp,
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "logs_dir": os.path.join(tmp, "LOGS"),
                    "images_dir": os.path.join(tmp, "IMAGES"),
                    "worker_parallel": 1,
                    "sampler_method": "Grid",
                    "workflow": {},
                    "workflow_layers": {},
                },
                factory_blueprint={"max_workers": 1, "workflow": {}, "calc_layer": {}, "module_pools": {}},
                sampler_state={
                    "sampler_signature": {
                        "method": "Grid",
                        "dimensions": 1,
                        "nchains": 1,
                        "niters": 1,
                        "variable_names": ["x"],
                        "task_result_dir": os.path.join(tmp, "RESULTS"),
                    },
                    "state_machine": {"state": "INIT", "ready_queue": [], "chains": [], "extras": {}},
                },
                reason="unit-test",
            )
            StateSaver(checkpoint_file, logger=_NoopLogger()).save(payload)

            core = Core()
            core.args = SimpleNamespace(debug=False, plot=False, skipFC=True, skiplibrary=True, resume=False)
            core.mode = "1PC"
            core.path = {"task_root": tmp, "jpath": tmp, "logo": tmp, "args_info": tmp}
            core.info = {
                "scan_name": scan_name,
                "project_name": scan_name,
                "config_file": os.path.join(tmp, "current.yaml"),
                "sample": {
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "sample_dirs": os.path.join(tmp, "RESULTS", "SAMPLE"),
                    "archive_samples": False,
                },
                "logs_dir": os.path.join(tmp, "LOGS"),
                "images_dir": os.path.join(tmp, "IMAGES"),
            }
            core.yaml.config = deepcopy(_minimal_config("Grid"))
            core.yaml.get_sampling_method = lambda: "Grid"
            core.yaml.get_worker_parallel = lambda: 1
            core.sampler = SimpleNamespace(method="Grid")

            with mock.patch.object(core, "_prompt_resume_from_checkpoint", side_effect=AssertionError("prompt must not be called")):
                core._preload_resume_checkpoint()

            self.assertIsNone(core._resume_checkpoint_payload)
            self.assertIsNone(core._resume_run_spec)
            self.assertIsNone(core._resume_factory_blueprint)
            self.assertEqual(core._resume_checkpoint_policy, "fresh")
            self.assertTrue(os.path.exists(checkpoint_file))

    def test_core_resume_prompt_defaults_to_resume_on_blank_input(self):
        core = Core()
        with contextlib.redirect_stdout(io.StringIO()):
            with mock.patch("sys.stdin", _TTYStringIO("\n")):
                self.assertTrue(core._prompt_resume_from_checkpoint("/tmp/state.pkl", timeout_seconds=0.1))

    def test_core_resume_prompt_accepts_yes_for_fresh_run(self):
        core = Core()
        with contextlib.redirect_stdout(io.StringIO()):
            with mock.patch("sys.stdin", _TTYStringIO("y\n")):
                self.assertFalse(core._prompt_resume_from_checkpoint("/tmp/state.pkl", timeout_seconds=0.1))

    def test_core_resume_flag_skips_prompt_and_restores_checkpoint(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-core-resume-") as tmp:
            scan_name = "resume-scan-flag"
            checkpoint_root = os.path.join(tmp, "checkpoints", scan_name, "Grid")
            os.makedirs(checkpoint_root, exist_ok=True)
            checkpoint_file = os.path.join(checkpoint_root, "state.pkl")

            payload = build_state_payload(
                run_spec={
                    "raw_yaml_text": "frozen",
                    "normalized_config": deepcopy(_minimal_config("Grid")),
                    "scan_name": scan_name,
                    "task_root": tmp,
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "logs_dir": os.path.join(tmp, "LOGS"),
                    "images_dir": os.path.join(tmp, "IMAGES"),
                    "worker_parallel": 1,
                    "sampler_method": "Grid",
                    "workflow": {},
                    "workflow_layers": {},
                },
                factory_blueprint={"max_workers": 1, "workflow": {}, "calc_layer": {}, "module_pools": {}},
                sampler_state={
                    "sampler_signature": {
                        "method": "Grid",
                        "dimensions": 1,
                        "nchains": 1,
                        "niters": 1,
                        "variable_names": ["x"],
                        "task_result_dir": os.path.join(tmp, "RESULTS"),
                    },
                    "state_machine": {"state": "INIT", "ready_queue": [], "chains": [], "extras": {}},
                },
                reason="unit-test",
            )
            StateSaver(checkpoint_file, logger=_NoopLogger()).save(payload)

            core = Core()
            core.args = SimpleNamespace(debug=False, plot=False, skipFC=True, skiplibrary=True, resume=True)
            core.path = {"task_root": tmp, "jpath": tmp, "logo": tmp, "args_info": tmp}
            core.info = {
                "scan_name": scan_name,
                "project_name": scan_name,
                "config_file": os.path.join(tmp, "current.yaml"),
                "sample": {
                    "task_result_dir": os.path.join(tmp, "RESULTS"),
                    "sample_dirs": os.path.join(tmp, "RESULTS", "SAMPLE"),
                    "archive_samples": False,
                },
                "logs_dir": os.path.join(tmp, "LOGS"),
                "images_dir": os.path.join(tmp, "IMAGES"),
            }
            core.yaml.config = deepcopy(_minimal_config("Grid"))
            core.yaml.get_sampling_method = lambda: "Grid"
            core.yaml.get_worker_parallel = lambda: 1
            core.sampler = SimpleNamespace(method="Grid")

            with mock.patch.object(core, "_prompt_resume_from_checkpoint", side_effect=AssertionError("prompt must not be called")):
                core._preload_resume_checkpoint()

            self.assertIsInstance(core._resume_checkpoint_payload, dict)
            self.assertEqual(core._resume_checkpoint_policy, "resume")

    def test_init_statesaver_skips_auto_resume_when_policy_is_fresh(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-core-statesaver-") as tmp:
            class _FakeSampler:
                method = "Random"

                def __init__(self) -> None:
                    self.configure_calls = []
                    self.restore_calls = 0

                def supports_runtime_checkpointing(self):
                    return True

                def configure_runtime_checkpointing(self, checkpoint_root, **kwargs):
                    self.configure_calls.append((checkpoint_root, dict(kwargs)))

                def restore_runtime_checkpoint_if_available(self, _payload=None):
                    self.restore_calls += 1
                    return False

            core = Core()
            core._resume_checkpoint_policy = "fresh"
            core.path = {"task_root": tmp}
            core.info = {"scan_name": "fresh-opc", "project_name": "fresh-opc"}
            core.sampler = _FakeSampler()

            core.init_StateSaver()

            self.assertEqual(len(core.sampler.configure_calls), 1)
            _root, kwargs = core.sampler.configure_calls[0]
            self.assertFalse(kwargs.get("auto_resume", True))
            self.assertEqual(core.sampler.restore_calls, 0)

    def test_module_manager_restore_factory_blueprint_fails_fast(self):
        manager = ModuleManager()
        manager.module_pools = {
            "A": _FakePool("A", "Calculator"),
            "B": _FakePool("B", "Calculator"),
        }
        manager._max_worker = 1

        with self.assertRaises(ValueError):
            manager.restore_factory_blueprint(
                {
                    "max_workers": 1,
                    "workflow": {},
                    "module_pools": {
                        "A": {"name": "A", "type": "Calculator", "config": {}, "instances": {}},
                    },
                }
            )

        with self.assertRaises(ValueError):
            manager.restore_factory_blueprint(
                {
                    "max_workers": 1,
                    "workflow": {},
                    "module_pools": {
                        "A": {"name": "A", "type": "Calculator", "config": {}, "instances": {}},
                        "B": {"name": "B", "type": "Operas", "config": {}, "instances": {}},
                    },
                }
            )


if __name__ == "__main__":
    unittest.main()
