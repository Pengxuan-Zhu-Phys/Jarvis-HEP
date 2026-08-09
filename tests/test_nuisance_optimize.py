#!/usr/bin/env python3
"""D13.4 nuisance_optimize Worker step + Profile1D tests."""

from __future__ import annotations

import os
import tempfile
import threading
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Module.nuisance import (
    NuisanceExpressionRegistry,
    NuisancePassConditionRegistry,
    extract_nuisance_config,
    parse_nuisance_variable,
)
from jarvishep2.Module.profile1d import Profile1DProfiler, ProfileProbeResult
from jarvishep2.core import Jarvis2Core
from jarvishep2.expression import ExpressionContext
from jarvishep2.factory import TaskFactory
from jarvishep2.flowchart import build_flowchart_scene_from_config
from jarvishep2.sample import Sample
from jarvishep2.worker import Worker
from jarvishep2.workflow import build_execution_plan, execution_plan_template


def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0, "managed": False}


class RegistryUnitTests(unittest.TestCase):
    def test_logl_and_pass_compile_eval(self) -> None:
        ctx = ExpressionContext()
        logl = NuisanceExpressionRegistry(ctx)
        logl.set_config("L1", "(x - t)**2")
        pass_reg = NuisancePassConditionRegistry(ctx)
        pass_reg.set_config("ok", "x > 0")
        vals = {"x": 1.0, "t": 0.5}
        self.assertAlmostEqual(logl.total(vals), 0.25)
        self.assertTrue(pass_reg.all_pass(vals))
        self.assertFalse(pass_reg.all_pass({"x": -1.0}))

    def test_extract_and_parse(self) -> None:
        cfg = {
            "Sampling": {
                "Nuisance": {
                    "method": "Profile1D",
                    "max_attempt": 20,
                    "variables": [
                        {
                            "name": "ratio",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": -2.0, "max": 2.0},
                            },
                        }
                    ],
                    "log_likelihood": [{"name": "L", "expression": "ratio**2"}],
                    "target_mode": "min",
                }
            }
        }
        block = extract_nuisance_config(cfg)
        assert block is not None
        name, zmin, zmax = parse_nuisance_variable(block)
        self.assertEqual(name, "ratio")
        self.assertEqual(zmin, -2.0)
        self.assertEqual(zmax, 2.0)


class Profile1DUnitTests(unittest.TestCase):
    def test_re_run_physics_defaults_true_v1_parity(self) -> None:
        """D13.7c: default True matches V1 full-pipeline-per-probe cost."""
        cfg = {
            "Sampling": {
                "Nuisance": {
                    "method": "Profile1D",
                    "max_attempt": 10,
                    "variables": [
                        {
                            "name": "ratio",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": -1.0, "max": 1.0},
                            },
                        }
                    ],
                    "log_likelihood": [{"name": "L", "expression": "ratio**2"}],
                    "target_mode": "min",
                }
            }
        }
        profiler = Profile1DProfiler.from_config(cfg)
        assert profiler is not None
        self.assertTrue(profiler.re_run_physics)
        cfg["Sampling"]["Nuisance"]["rerun_physics"] = False
        profiler_off = Profile1DProfiler.from_config(cfg)
        assert profiler_off is not None
        self.assertFalse(profiler_off.re_run_physics)

    def test_golden_section_finds_minimum(self) -> None:
        profiler = Profile1DProfiler(
            var_name="z",
            zmin=-2.0,
            zmax=2.0,
            logl_registry=NuisanceExpressionRegistry(),
            pass_registry=NuisancePassConditionRegistry(),
            mode="min",
            max_attempt=40,
            re_run_physics=False,
        )

        def evaluate(z: float) -> ProfileProbeResult:
            # Minimum at z=0.3
            logl = (z - 0.3) ** 2
            return ProfileProbeResult(z=z, logl=logl, terms={"L": logl}, pass_ok=True)

        result = profiler.optimize(evaluate)
        self.assertAlmostEqual(result.best_z, 0.3, places=2)
        self.assertLess(result.best_logl, 1e-3)
        self.assertEqual(result.status, "Accept")

    def test_pass_condition_failure_status(self) -> None:
        ctx = ExpressionContext()
        logl = NuisanceExpressionRegistry(ctx)
        logl.set_config("L", "z**2")
        pass_reg = NuisancePassConditionRegistry(ctx)
        pass_reg.set_config("never", "z > 10")
        profiler = Profile1DProfiler(
            var_name="z",
            zmin=-1.0,
            zmax=1.0,
            logl_registry=logl,
            pass_registry=pass_reg,
            mode="min",
            max_attempt=15,
            re_run_physics=False,
        )

        def evaluate(z: float) -> ProfileProbeResult:
            payload = {"z": z}
            terms = logl.evaluate_all(payload)
            total = sum(terms.values())
            pass_terms = pass_reg.evaluate_all(payload)
            return ProfileProbeResult(
                z=z,
                logl=total,
                terms=terms,
                pass_ok=all(pass_terms.values()),
                pass_terms=pass_terms,
            )

        result = profiler.optimize(evaluate)
        self.assertFalse(result.pass_ok)
        self.assertEqual(result.status, "Failed")


class WorkflowPlanTests(unittest.TestCase):
    def test_plan_includes_nuisance_before_likelihood(self) -> None:
        steps = build_execution_plan(
            opera_modules=[{"name": "Egg"}],
            include_likelihood=True,
            include_nuisance=True,
        )
        types = [s.type for s in steps]
        self.assertIn("nuisance_optimize", types)
        self.assertIn("likelihood", types)
        self.assertLess(types.index("nuisance_optimize"), types.index("likelihood"))
        self.assertEqual(steps[types.index("nuisance_optimize")].name, "NuisanceOptimize")

    def test_template_json(self) -> None:
        plan = execution_plan_template(
            opera_modules=[{"name": "O"}],
            include_likelihood=True,
            include_nuisance=True,
        )
        self.assertTrue(any(s["type"] == "nuisance_optimize" for s in plan))


class FlowchartNuisanceTests(unittest.TestCase):
    def test_nuisance_node_on_flowchart(self) -> None:
        cfg = {
            "Operas": {
                "Modules": [
                    {
                        "name": "Egg",
                        "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                        "input": [{"name": "x", "expression": "x"}],
                        "output": [{"name": "z", "entry": "z"}],
                    }
                ]
            },
            "Sampling": {
                "Nuisance": {
                    "method": "Profile1D",
                    "variables": [
                        {
                            "name": "t",
                            "distribution": {
                                "type": "Flat",
                                "parameters": {"min": 0, "max": 1},
                            },
                        }
                    ],
                    "log_likelihood": [{"name": "L", "expression": "t**2"}],
                }
            },
        }
        scene = build_flowchart_scene_from_config(cfg, include_likelihood=True)
        names = [n.get("name") for n in scene.get("nodes") or []]
        # Scene may use modules list in different shapes — also check modules map.
        text = str(scene)
        self.assertIn("NuisanceOptimize", text)


class WorkerNuisanceInProcessTests(unittest.TestCase):
    def test_worker_profile_expression_only(self) -> None:
        """In-process Worker: expression-only profile (no physics re-run)."""
        from jarvishep2.redis_queue import make_fakeredis_queue

        queue = make_fakeredis_queue()
        nuisance = {
            "method": "Profile1D",
            "max_attempt": 30,
            "variables": [
                {
                    "name": "t",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": -1.0, "max": 1.0},
                    },
                }
            ],
            "log_likelihood": [{"name": "L", "expression": "(t - 0.25)**2"}],
            "target_mode": "min",
            "rerun_physics": False,
            "pass_condition": [{"name": "near", "expression": "Abs(t - 0.25) < 0.2"}],
        }
        worker = Worker(
            0,
            queue,
            {
                "opera_modules": [],
                "calculator_modules": [],
                "likelihood_expressions": [],
                "nuisance_config": nuisance,
                "mapper": {"type": "identity", "keys": ["x"]},
                "sample_config": {
                    "task_result_dir": tempfile.mkdtemp(),
                    "sample_artifacts": "never",
                },
                "file_operation_mode": "inline",
                "publish_feedback": False,
            },
        )
        worker._init_runtime()
        self.assertIsNotNone(worker._nuisance_profiler)

        sample = Sample(uuid="nuis-1", u_coords=np.array([0.5]), execution_plan=[])
        sample.params = {"x": 1.0}
        sample.observables = {"x": 1.0}
        sample.info = {"logger_name": "Sample@nuis-1"}
        sample.execution_plan = build_execution_plan(
            include_nuisance=True,
            include_likelihood=False,
        )
        # Force a bare nuisance step run.
        worker._run_nuisance_optimize(sample)
        self.assertAlmostEqual(float(sample.observables["t"]), 0.25, places=2)
        self.assertIn("NuisanceLogL", sample.observables)
        self.assertTrue(sample.info.get("nuisance", {}).get("pass"))
        self.assertEqual(sample.info["nuisance"]["status"], "Accept")


class CoreNuisanceIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_core_bridson_with_nuisance_expression(self) -> None:
        """End-to-end: Bridson + Opera + Profile1D expression path."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg: dict[str, Any] = {
                    "project_name": "nuisance_e2e",
                    "Scan": {"name": "nuis-scan"},
                    "task_result_dir": tmp,
                    "Runtime": {
                        "mode": "redis",
                        "workers": 1,
                        "batch_size": 1,
                        "redis": redis_config,
                        "sample_artifacts": "never",
                        "Watchdog": {"enabled": False},
                    },
                    "Sampling": {
                        "Method": "Bridson",
                        "Bounds": {"seed": 3, "radius": 0.4, "point_number": 4, "max_attempt": 30},
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {
                                        "min": 0,
                                        "max": 1,
                                        "length": 1.0,
                                    },
                                },
                            },
                            {
                                "name": "y",
                                "distribution": {
                                    "type": "Flat",
                                    "parameters": {
                                        "min": 0,
                                        "max": 1,
                                        "length": 1.0,
                                    },
                                },
                            },
                        ],
                        "LogLikelihood": [
                            {"name": "LogL_z", "expression": "LogGauss(z, 1, 1)"}
                        ],
                        "Nuisance": {
                            "method": "Profile1D",
                            "variables": [
                                {
                                    "name": "shift",
                                    "distribution": {
                                        "type": "Flat",
                                        "parameters": {"min": -0.5, "max": 0.5},
                                    },
                                }
                            ],
                            "log_likelihood": [
                                {
                                    "name": "LogLNP",
                                    "expression": "(shift - 0.1)**2",
                                }
                            ],
                            "target_mode": "min",
                            "max_attempt": 20,
                            "rerun_physics": False,
                            "pass_condition": [
                                {"name": "ok", "expression": "Abs(shift - 0.1) < 0.3"}
                            ],
                        },
                    },
                    "Operas": {
                        "Modules": [
                            {
                                "name": "Egg",
                                "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                                "call_mode": "call",
                                "input": [
                                    {"name": "x", "expression": "x"},
                                    {"name": "y", "expression": "y"},
                                ],
                                "output": [{"name": "z", "entry": "z"}],
                            }
                        ]
                    },
                }
                from jarvishep2.workflow import execution_plan_template

                plan = execution_plan_template(
                    opera_modules=cfg["Operas"]["Modules"],
                    include_likelihood=True,
                    include_nuisance=True,
                )
                self.assertIn("nuisance_optimize", [s["type"] for s in plan])
                core = Jarvis2Core(cfg)
                outcome = core.run()
                submitted = int(getattr(outcome, "submitted", 0) or 0)
                self.assertGreater(submitted, 0)
                # Flowchart should mention nuisance
                scene_path = os.path.join(tmp, "images", "nuis-scan", "flowchart.json")
                if not os.path.isfile(scene_path):
                    scene_path = os.path.join(tmp, "images", "scan", "flowchart.json")
                # Best-effort: graph export may land under scan_name
                found = False
                for root, _dirs, files in os.walk(tmp):
                    if "flowchart.json" in files:
                        text = open(
                            os.path.join(root, "flowchart.json"), encoding="utf-8"
                        ).read()
                        if "NuisanceOptimize" in text:
                            found = True
                            break
                self.assertTrue(found, "flowchart missing NuisanceOptimize node")
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
