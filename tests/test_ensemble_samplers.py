#!/usr/bin/env python3
"""D13.3 Ensemble / DEMCMC / PT sampler tests."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.slow

import tempfile
import threading
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.Source.MCMC.engine_demcmc import DEMCMCChain
from jarvishep2.Sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep2.Sampling.mcmc_sampler import MCMCSampler
from jarvishep2.core import Jarvis2Core
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue


def _start_tcp_fakeredis() -> tuple[TcpFakeServer, dict]:
    server = TcpFakeServer(("127.0.0.1", 0))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    host, port = server.server_address
    return server, {"host": host, "port": port, "db": 0, "managed": False}


def _flat_vars(names: list[str], lo: float = 0.0, hi: float = 1.0) -> list[dict[str, Any]]:
    return [
        {
            "name": name,
            "distribution": {
                "type": "Flat",
                "parameters": {"min": lo, "max": hi},
            },
        }
        for name in names
    ]


def _cfg(
    *,
    method: str,
    tmpdir: str,
    nchains: int = 4,
    niters: int = 10,
    workers: int = 2,
    seed: int = 7,
    redis_cfg: dict[str, Any] | None = None,
    extra_bounds: dict[str, Any] | None = None,
) -> dict[str, Any]:
    bounds: dict[str, Any] = {
        "num_chains": nchains,
        "num_iters": niters,
        "proposal_scale": 0.12,
    }
    if extra_bounds:
        bounds.update(extra_bounds)
    bounds.setdefault("seed", seed)
    runtime: dict[str, Any] = {
        "mode": "redis",
        "workers": workers,
        "batch_size": workers,
        "sample_artifacts": "never",
        "Watchdog": {"enabled": False},
    }
    if redis_cfg is not None:
        runtime["redis"] = redis_cfg
    return {
        "project_name": "ensemble_test",
        "Scan": {"name": f"{method.lower()}-scan"},
        "task_result_dir": tmpdir,
        "Runtime": runtime,
        "Sampling": {
            "Method": method,
            "Bounds": bounds,
            "Variables": _flat_vars(["x", "y"], 0.0, 5.0),
            "LogLikelihood": [
                {"name": "LogL_Z", "expression": "LogGauss(z, 100, 10)"}
            ],
        },
        "Operas": {
            "Modules": [
                {
                    "name": "EggBox",
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


class EngineUnitTests(unittest.TestCase):
    def test_ensemble_chain_seeded(self) -> None:
        pop = [np.array([0.2, 0.3]), np.array([0.7, 0.8])]

        def getter(_cid: int):
            return pop

        def run(seed: int) -> list[float]:
            rng = np.random.default_rng(seed)
            eng = EnsembleChain(
                rng.random(2), 0.1, 4, chain_id=0, population_getter=getter, rng=rng
            )
            vals = []
            for i in range(4):
                next(eng)
                eng.update(-0.1 * i)
                vals.append(float(eng.param[0]))
            return vals

        self.assertEqual(run(1), run(1))

    def test_demcmc_chain_seeded(self) -> None:
        pop = [np.array([0.1, 0.2]), np.array([0.8, 0.9]), np.array([0.4, 0.5])]

        def getter(_cid: int):
            return pop

        rng = np.random.default_rng(3)
        eng = DEMCMCChain(
            rng.random(2),
            0.1,
            3,
            chain_id=1,
            population_getter=getter,
            rng=rng,
        )
        next(eng)
        self.assertTrue(eng.update(0.0))


class DistributorEnsembleTests(unittest.TestCase):
    def test_methods_registered(self) -> None:
        for name in (
            "EnsembleMCMC",
            "Ensemble",
            "DEMCMC",
            "PTMCMC",
            "PTEnsemble",
        ):
            self.assertNotIn(name, STATELESS_METHODS)
            sampler = Distributor.set_method(name)
            self.assertIsInstance(sampler, MCMCSampler)
            self.assertEqual(Distributor.get_resume_status(name), "implemented")


class FeedbackLoopTests(unittest.TestCase):
    def _run_mock(
        self,
        method: str,
        *,
        nchains: int = 4,
        niters: int = 8,
        workers: int = 2,
        seed: int = 11,
        extra_bounds: dict[str, Any] | None = None,
    ) -> MCMCSampler:
        queue = make_fakeredis_queue()
        sampler = Distributor.set_method(method)
        assert isinstance(sampler, MCMCSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _cfg(
                method=method,
                tmpdir=tmp,
                nchains=nchains,
                niters=niters,
                workers=workers,
                seed=seed,
                extra_bounds=extra_bounds,
            )
            sampler.set_config(cfg)
            sampler.set_redis(queue)
            stop = threading.Event()

            def worker() -> None:
                while not stop.is_set():
                    task = queue.pull_task(timeout=1)
                    if task is None:
                        continue
                    u = np.asarray(task.get("u_coords") or [0.5, 0.5], dtype=float)
                    logl = -float(np.sum((u - 0.35) ** 2))
                    queue.publish_feedback(
                        {
                            "uuid": task["uuid"],
                            "status": "Completed",
                            "observables": {"LogL": logl},
                        }
                    )

            thread = threading.Thread(target=worker, daemon=True)
            thread.start()
            try:
                n = sampler.run_adaptive(timeout=60.0)
            finally:
                stop.set()
                thread.join(timeout=2.0)
            self.assertGreater(n, 0)
            self.assertTrue(sampler._all_finished())
            return sampler

    def test_ensemble_mcmc_loop(self) -> None:
        s = self._run_mock("EnsembleMCMC", nchains=4, niters=6)
        summary = s.summary()
        self.assertTrue(summary.get("half_ensemble"))
        self.assertEqual(summary["n_chains"], 4)

    def test_demcmc_loop(self) -> None:
        s = self._run_mock("DEMCMC", nchains=4, niters=6)
        self.assertEqual(s.method, "DEMCMC")
        self.assertGreater(s.summary()["total_proposed"], 0)

    def test_ptmcmc_exchange(self) -> None:
        s = self._run_mock(
            "PTMCMC",
            nchains=4,
            niters=8,
            extra_bounds={
                "temperature_ladder": [1.0, 2.0, 4.0, 8.0],
                "exchange_interval": 1,
            },
        )
        summary = s.summary()
        self.assertIn("swap_attempts", summary)
        self.assertGreaterEqual(summary["swap_attempts"], 1)
        self.assertEqual(summary.get("async_chains"), False)
        # Temperatures stay on ladder slots (params may have swapped).
        temps = [c.temperature for c in s._ensure_registry().all()]
        self.assertEqual(temps, [1.0, 2.0, 4.0, 8.0])
        self.assertFalse(s._uses_async_chain_pipeline())
        self.assertTrue(s._uses_pt())

    def test_pt_ensemble_loop(self) -> None:
        s = self._run_mock(
            "PTEnsemble",
            nchains=4,
            niters=6,
            extra_bounds={
                "temperature_ladder": [1.0, 1.5, 2.5, 4.0],
                "stretch_a": 2.0,
                "exchange_interval": 2,
            },
        )
        self.assertEqual(s.method, "PTEnsemble")
        self.assertTrue(s.summary().get("half_ensemble"))

    def test_worker_count_independent_ensemble(self) -> None:
        def finals(workers: int) -> list[list[float]]:
            queue = make_fakeredis_queue()
            sampler = Distributor.set_method("EnsembleMCMC")
            assert isinstance(sampler, MCMCSampler)
            with tempfile.TemporaryDirectory() as tmp:
                cfg = _cfg(
                    method="EnsembleMCMC",
                    tmpdir=tmp,
                    nchains=4,
                    niters=10,
                    workers=workers,
                    seed=42,
                )
                sampler.set_config(cfg)
                sampler.set_redis(queue)
                stop = threading.Event()

                def worker() -> None:
                    while not stop.is_set():
                        task = queue.pull_task(timeout=1)
                        if task is None:
                            continue
                        u = np.asarray(task.get("u_coords") or [0.5, 0.5], dtype=float)
                        logl = -float(np.sum((u - 0.4) ** 2))
                        queue.publish_feedback(
                            {
                                "uuid": task["uuid"],
                                "status": "Completed",
                                "observables": {"LogL": logl},
                            }
                        )

                thread = threading.Thread(target=worker, daemon=True)
                thread.start()
                try:
                    sampler.run_adaptive(timeout=60.0)
                finally:
                    stop.set()
                    thread.join(timeout=2.0)
                return [
                    np.asarray(c.engine.param, dtype=float).tolist()
                    for c in sampler._ensure_registry().all()
                ]

        a = finals(1)
        b = finals(3)
        for pa, pb in zip(a, b):
            np.testing.assert_allclose(pa, pb, atol=1e-12)


class CoreEnsembleIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_core_ensemble_eggbox_tiny(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = _cfg(
                    method="EnsembleMCMC",
                    tmpdir=tmp,
                    redis_cfg=redis_config,
                    nchains=4,
                    niters=3,
                    workers=1,
                    seed=5,
                    extra_bounds={"stretch_a": 2.0},
                )
                core = Jarvis2Core(cfg)
                outcome = core.run()
                submitted = int(getattr(outcome, "submitted", 0) or 0)
                self.assertGreater(submitted, 0)
                sampler = core.sampler
                assert isinstance(sampler, MCMCSampler)
                self.assertTrue(sampler.summary().get("half_ensemble"))
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
