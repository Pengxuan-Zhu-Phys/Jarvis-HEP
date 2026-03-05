#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.core import Core  # noqa: E402
from jarvishep.factory import WorkerFactory  # noqa: E402
from jarvishep.moduleManager import ModuleManager  # noqa: E402
from jarvishep.workflow import Workflow  # noqa: E402


class _NoopSampleLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _InMemoryDatabase:
    def __init__(self):
        self.rows = []

    def add_data(self, row):
        self.rows.append(dict(row))


class _FakeYaml:
    def __init__(self):
        self.config = {"Sampling": {}}

    @staticmethod
    def get_worker_parallel():
        return 2


class _FakeSampler:
    def __init__(self):
        self._with_nuisance = False
        self.max_workers = None

    def set_max_workers(self, nworkers):
        self.max_workers = nworkers


class _FakeOperasModule:
    def __init__(self, name, inputs, outputs, callback):
        self.name = name
        self.type = "Operas"
        self.required_modules = []
        self.inputs = {item: None for item in inputs}
        self.outputs = {item: None for item in outputs}
        self._callback = callback
        self.calls = 0

    def set_funcs(self, funcs):
        self._funcs = funcs

    def set_logger(self, logger):
        self.logger = logger

    def execute(self, observables, sample_info):
        self.calls += 1
        return dict(self._callback(dict(observables), sample_info))


class ExecutionPathSmokeTests(unittest.TestCase):
    def setUp(self):
        WorkerFactory._instance = None
        ModuleManager._instance = None
        self.core = None

    def tearDown(self):
        if self.core is not None and getattr(self.core, "factory", None) is not None:
            if hasattr(self.core.factory, "shutdown"):
                self.core.factory.shutdown(wait=True, cancel_futures=True)
            else:
                executor = getattr(self.core.factory, "executor", None)
                if executor is not None:
                    executor.shutdown(wait=True, cancel_futures=True)
        WorkerFactory._instance = None
        ModuleManager._instance = None

    @staticmethod
    def _build_fake_workflow():
        workflow = Workflow()
        seed = _FakeOperasModule(
            "SeedX",
            inputs=[],
            outputs=["x"],
            callback=lambda obs, _sample: {"x": obs["x"]},
        )
        build_y = _FakeOperasModule(
            "BuildY",
            inputs=["x"],
            outputs=["y"],
            callback=lambda obs, _sample: {"y": obs["x"] + 1},
        )
        build_z = _FakeOperasModule(
            "BuildZ",
            inputs=["y"],
            outputs=["z"],
            callback=lambda obs, _sample: {"z": obs["y"] * 10},
        )
        workflow.add_module(seed)
        workflow.add_module(build_y)
        workflow.add_module(build_z)
        workflow.resolve_dependencies()
        workflow.workflow = {
            layer_id: list(layer_info["module"])
            for layer_id, layer_info in workflow.calc_layer.items()
            if layer_id > 1
        }
        return workflow

    def test_core_workflow_module_manager_smoke_path(self):
        workflow = self._build_fake_workflow()

        self.core = Core()
        self.core.yaml = _FakeYaml()
        self.core.sampler = _FakeSampler()
        self.core.workflow = workflow
        self.core.info = {"sample": {"task_result_dir": "/tmp"}}
        self.core._funcs = {}

        self.core.init_WorkerFactory()

        self.assertEqual(self.core.sampler.max_workers, 2)
        self.assertEqual(self.core.module_manager.workflow, {2: ["BuildY"], 3: ["BuildZ"]})
        self.assertIn("BuildY", self.core.module_manager.module_pools)
        self.assertIn("BuildZ", self.core.module_manager.module_pools)

        db = _InMemoryDatabase()
        self.core.module_manager._database = db

        sample_info = {
            "uuid": "smoke-001",
            "observables": {"x": 2},
            "logger": _NoopSampleLogger(),
        }
        result = self.core.module_manager.execute_workflow(sample_info)

        self.assertEqual(result, 1.0)
        self.assertEqual(len(db.rows), 1)
        self.assertEqual(db.rows[0]["x"], 2)
        self.assertEqual(db.rows[0]["y"], 3)
        self.assertEqual(db.rows[0]["z"], 30)
        self.assertEqual(sample_info["observables"]["z"], 30)


if __name__ == "__main__":
    unittest.main()
