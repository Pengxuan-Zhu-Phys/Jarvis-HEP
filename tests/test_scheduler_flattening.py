#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import threading
import time
import unittest
from unittest import mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.factory import WorkerFactory  # noqa: E402
from jarvishep.io_manager import IOManager  # noqa: E402
from jarvishep.runtime_config import workflow_has_calculator  # noqa: E402
from jarvishep.moduleManager import ModuleManager  # noqa: E402
from jarvishep.modulePool import ModulePool  # noqa: E402
from jarvishep.Module.calculator import CalculatorModule  # noqa: E402
from jarvishep.Module.module import Module  # noqa: E402


BENCHMARK_PROJECT = os.path.join(PROJECT_ROOT, "tests", "benchmark_project")
BENCHMARK_TASK = os.path.join(BENCHMARK_PROJECT, "bin", "benchmark_random_operas.yaml")
ENABLE_PERF_GATES = os.getenv("JARVIS_HEP_ENABLE_PERF_GATES") == "1"


def _thread_names() -> set[str]:
    return {thread.name for thread in threading.enumerate()}


def _count_threads(prefix: str) -> int:
    return sum(1 for name in _thread_names() if name.startswith(prefix))


class _FakeYaml:
    def __init__(self, *, with_calculator: bool = False):
        self.config = {"Sampling": {}}
        if with_calculator:
            self.config["Calculators"] = {"Modules": [{"name": "DemoCalc"}]}

    @staticmethod
    def get_worker_parallel():
        return 4

    def get_subprocess_runtime_options(self, worker_parallel=None):
        workers = max(1, int(worker_parallel or 4))
        return {
            "max_concurrency": workers,
            "max_pending": max(128, workers * 16),
            "per_task_timeout_sec": None,
            "progress_interval_sec": 5.0,
            "log_policy": "logger",
            "diagnostics_enabled": False,
            "diagnostics_interval_sec": 10.0,
            "terminate_grace_sec": 5.0,
        }


class _FakeSampler:
    def set_max_workers(self, nworkers):
        self.max_workers = nworkers


class _NoopLogger:
    def bind(self, **_kwargs):
        return self

    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _FakeModuleManager:
    def execute_workflow(self, sample_info):
        time.sleep(0.01)
        return 1.0

    def _module_failure_policy(self):
        return "fail-fast"


class SchedulerFlatteningUnitTests(unittest.TestCase):
    def setUp(self):
        WorkerFactory._instance = None
        ModuleManager._instance = None

    def tearDown(self):
        WorkerFactory._instance = None
        ModuleManager._instance = None

    def test_factory_uses_single_named_executor(self):
        factory = WorkerFactory()
        factory.configure(module_manager=_FakeModuleManager(), max_workers=3)
        self.assertFalse(hasattr(factory, "log_executor"))
        self.assertIsNotNone(factory.executor)
        factory.shutdown(wait=True, cancel_futures=True)

    def test_module_pool_has_no_installer_executor(self):
        module = CalculatorModule(
            "DemoCalc",
            {
                "modes": False,
                "required_modules": [],
                "clone_shadow": False,
                "installation": [],
                "initialization": [],
                "execution": {"commands": [], "input": [], "output": []},
                "path": "/tmp/demo",
            },
        )
        pool = ModulePool(module, max_workers=2)
        self.assertFalse(hasattr(pool, "executor"))

    def test_operas_workflow_skips_io_manager_pool(self):
        config = _FakeYaml(with_calculator=False).config
        self.assertFalse(workflow_has_calculator(config))
        io_manager = (
            IOManager(max_workers=4)
            if workflow_has_calculator(config)
            else None
        )
        self.assertIsNone(io_manager)

    def test_calculator_workflow_keeps_io_manager_pool(self):
        config = _FakeYaml(with_calculator=True).config
        self.assertTrue(workflow_has_calculator(config))
        io_manager = (
            IOManager(max_workers=4)
            if workflow_has_calculator(config)
            else None
        )
        self.assertIsNotNone(io_manager)
        self.assertEqual(io_manager.max_workers, 4)
        io_manager.shutdown(wait=True, cancel_futures=True)

    def test_module_pool_installs_on_factory_worker_thread(self):
        install_threads: list[str] = []

        def _record_install(instance):
            install_threads.append(threading.current_thread().name)
            instance.is_installed = True
            if instance.installation_event is not None:
                instance.installation_event.set()
            return instance

        module = CalculatorModule(
            "DemoCalc",
            {
                "modes": False,
                "required_modules": [],
                "clone_shadow": False,
                "installation": [],
                "initialization": [],
                "execution": {"commands": [], "input": [], "output": []},
                "path": "/tmp/demo",
            },
        )
        pool = ModulePool(module, max_workers=2)
        pool.set_logger()
        pool.set_funcs({})

        with mock.patch.object(ModulePool, "install_instance", side_effect=_record_install):
            with mock.patch.object(
                CalculatorModule,
                "execute",
                return_value={"x": 1.0},
            ):
                pool.execute({"x": 1.0}, {"uuid": "install-thread", "logger": _NoopLogger()})

        self.assertEqual(len(install_threads), 1)
        self.assertNotIn("ThreadPoolExecutor", install_threads[0])

    def test_operas_scan_has_no_io_manager_threads(self):
        factory = WorkerFactory()
        factory.configure(module_manager=_FakeModuleManager(), max_workers=4)
        factory.set_logger(_NoopLogger())

        barrier = threading.Event()
        observed_io_threads: list[int] = []

        class _SlowModuleManager:
            def execute_workflow(self, sample_info):
                observed_io_threads.append(_count_threads("jarvis-hep-io"))
                barrier.set()
                time.sleep(0.05)
                return 1.0

            def _module_failure_policy(self):
                return "fail-fast"

        factory.module_manager = _SlowModuleManager()
        futures = [
            factory.submit_task({"uuid": f"u-{idx}", "observables": {"x": idx}})
            for idx in range(8)
        ]
        self.assertTrue(barrier.wait(timeout=2.0))
        for fut in futures:
            fut.result(timeout=2.0)

        self.assertEqual(observed_io_threads[0], 0)
        factory_threads = _count_threads("jarvis-hep-factory")
        self.assertGreaterEqual(factory_threads, 1)
        self.assertLessEqual(factory_threads, 4)
        factory.shutdown(wait=True, cancel_futures=True)


class SchedulerFlatteningBenchmarkTests(unittest.TestCase):
    def setUp(self):
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "outputs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "logs"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "images"), ignore_errors=True)
        shutil.rmtree(os.path.join(BENCHMARK_PROJECT, "checkpoints"), ignore_errors=True)

    def tearDown(self):
        self.setUp()

    @unittest.skipUnless(
        ENABLE_PERF_GATES,
        "absolute throughput gates require JARVIS_HEP_ENABLE_PERF_GATES=1",
    )
    def test_benchmark_not_regressed_after_pool_flattening(self):
        env = dict(os.environ)
        env["PYTHONPATH"] = PROJECT_ROOT + os.pathsep + env.get("PYTHONPATH", "")
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "jarvishep",
                BENCHMARK_TASK,
                "--benchmark",
                "2",
                "--skip-draw-flowchart",
            ],
            cwd=BENCHMARK_PROJECT,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=60,
        )
        self.assertEqual(proc.returncode, 0, msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")

        benchmark_json = os.path.join(
            BENCHMARK_PROJECT,
            "outputs",
            "benchmark_random_operas",
            "benchmark.json",
        )
        with open(benchmark_json, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        self.assertGreaterEqual(payload["samples_per_sec"], 600.0)
        self.assertEqual(payload["total_failed"], 0)


if __name__ == "__main__":
    unittest.main()
