#!/usr/bin/env python3
"""D13.2 MCMC / AM / DRAM FeedbackSampler tests."""

from __future__ import annotations

import tempfile
import threading
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.Source.MCMC.config_contract import parse_common_chain_counts
from jarvishep2.Sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep2.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep2.Sampling.mcmc_sampler import MCMCSampler, mcmc_sample_uuid
from jarvishep2.core import Jarvis2Core
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.testing.eggbox import eggbox2d_numpy


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


def _mcmc_config(
    *,
    method: str,
    tmpdir: str,
    redis_cfg: dict[str, Any] | None = None,
    nchains: int = 4,
    niters: int = 20,
    workers: int = 2,
    seed: int = 7,
    proposal_scale: float = 0.15,
    extra_bounds: dict[str, Any] | None = None,
) -> dict[str, Any]:
    bounds: dict[str, Any] = {
        "num_chains": nchains,
        "num_iters": niters,
        "proposal_scale": proposal_scale,
    }
    if extra_bounds:
        bounds.update(extra_bounds)
    bounds.setdefault("Seed", seed)
    runtime: dict[str, Any] = {
        "mode": "redis",
        "workers": workers,
        "batch_size": workers,
    }
    if redis_cfg is not None:
        runtime["redis"] = redis_cfg
    return {
        "project_name": "mcmc_test",
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


class ConfigContractTests(unittest.TestCase):
    def test_steps_alias(self) -> None:
        nchains, niters = parse_common_chain_counts({"chains": 3, "steps": 50})
        self.assertEqual(nchains, 3)
        self.assertEqual(niters, 50)


class EngineUnitTests(unittest.TestCase):
    def test_mcmc_chain_seeded_repro(self) -> None:
        def run(seed: int) -> list[float]:
            rng = np.random.default_rng(seed)
            eng = MCMCChain(rng.random(2), 0.2, 5, rng=rng)
            logls = []
            for i in range(5):
                next(eng)
                eng.update(-float(i))
                logls.append(float(eng.last_loglikelihood or 0.0))
            return logls

        self.assertEqual(run(1), run(1))
        self.assertNotEqual(run(1), run(2))

    def test_dram_stage_followup_path(self) -> None:
        rng = np.random.default_rng(0)
        eng = DRAMChain(np.array([0.5, 0.5]), 0.05, 10, dr_steps=2, rng=rng)
        # Bootstrap accept
        eng.propose_stage(0)
        out0 = eng.consume_stage_result(0, 0.0)
        self.assertTrue(out0["iteration_done"])
        # Force a reject then delayed rejection: propose far and give worse logl
        eng.propose_stage(0)
        # Very negative logl → likely reject
        out1 = eng.consume_stage_result(0, -1.0e6)
        if not out1["iteration_done"]:
            self.assertEqual(out1["next_stage"], 1)
            eng.propose_stage(1)
            out2 = eng.consume_stage_result(1, -1.0e6)
            self.assertTrue(out2["iteration_done"])


class DistributorMCMCTests(unittest.TestCase):
    def test_methods_registered_stateful(self) -> None:
        for name in ("MCMC", "AMMCMC", "AM", "DRAM"):
            self.assertNotIn(name, STATELESS_METHODS)
            sampler = Distributor.set_method(name)
            self.assertIsInstance(sampler, MCMCSampler)
            self.assertEqual(Distributor.get_resume_status(name), "implemented")

    def test_uuid_deterministic(self) -> None:
        a = mcmc_sample_uuid(method="DRAM", seed=1, chain_id=2, step=3, stage=0)
        b = mcmc_sample_uuid(method="DRAM", seed=1, chain_id=2, step=3, stage=0)
        c = mcmc_sample_uuid(method="DRAM", seed=1, chain_id=2, step=3, stage=1)
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)


class FeedbackLoopTests(unittest.TestCase):
    def _run_with_mock_worker(
        self,
        method: str,
        *,
        nchains: int = 3,
        niters: int = 8,
        workers: int = 2,
        seed: int = 11,
    ) -> MCMCSampler:
        queue = make_fakeredis_queue()
        sampler = Distributor.set_method(method)
        assert isinstance(sampler, MCMCSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method=method,
                tmpdir=tmp,
                nchains=nchains,
                niters=niters,
                workers=workers,
                seed=seed,
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
                    # Map unit cube → [0,5]^2 roughly as Flat
                    x, y = float(u[0]) * 5.0, float(u[1]) * 5.0
                    z = float(eggbox2d_numpy(x, y)["z"])
                    # Match LogGauss(z, 100, 10) roughly: -0.5*((z-100)/10)^2
                    logl = -0.5 * ((z - 100.0) / 10.0) ** 2
                    queue.publish_feedback(
                        {
                            "uuid": task["uuid"],
                            "status": "Completed",
                            "observables": {"z": z, "LogL": logl, "LogL_Z": logl},
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
            summary = sampler.summary()
            self.assertEqual(summary["n_chains"], nchains)
            self.assertEqual(summary["n_iters"], niters)
            self.assertGreater(summary["total_proposed"], 0)
            return sampler

    def test_mcmc_feedback_loop(self) -> None:
        s = self._run_with_mock_worker("MCMC", nchains=3, niters=6)
        self.assertGreaterEqual(s.summary()["accept_rate"], 0.0)

    def test_dram_feedback_loop(self) -> None:
        s = self._run_with_mock_worker("DRAM", nchains=3, niters=6)
        self.assertEqual(s.method, "DRAM")

    def test_worker_count_independent_trajectory(self) -> None:
        """Same seed → same chain final params for 1 vs 3 workers (mock LogL)."""

        def finals(workers: int) -> list[list[float]]:
            queue = make_fakeredis_queue()
            sampler = Distributor.set_method("MCMC")
            assert isinstance(sampler, MCMCSampler)
            with tempfile.TemporaryDirectory() as tmp:
                cfg = _mcmc_config(
                    method="MCMC",
                    tmpdir=tmp,
                    nchains=4,
                    niters=12,
                    workers=workers,
                    seed=99,
                    proposal_scale=0.12,
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
                        # Deterministic fake LogL from u only (no race).
                        logl = -float(np.sum((u - 0.3) ** 2))
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
                reg = sampler._ensure_registry()
                return [
                    np.asarray(c.engine.param, dtype=float).tolist() for c in reg.all()
                ]

        a = finals(1)
        b = finals(3)
        for pa, pb in zip(a, b):
            np.testing.assert_allclose(pa, pb, atol=1e-12)

    def test_checkpoint_roundtrip(self) -> None:
        queue = make_fakeredis_queue()
        sampler = Distributor.set_method("AMMCMC")
        assert isinstance(sampler, MCMCSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="AMMCMC",
                tmpdir=tmp,
                nchains=2,
                niters=30,
                seed=3,
                extra_bounds={
                    "adapt_enabled": True,
                    "adapt_start_iter": 5,
                    "adapt_window": 3,
                },
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
                    queue.publish_feedback(
                        {
                            "uuid": task["uuid"],
                            "status": "Completed",
                            "observables": {"LogL": -float(np.sum(u**2))},
                        }
                    )

            thread = threading.Thread(target=worker, daemon=True)
            thread.start()
            try:
                # Run a few steps by manually driving one generation then export.
                sampler._ensure_registry()
                props = sampler.propose_generation()
                assert props is not None
                sampler._submit_sample_batch(list(props))
                results = sampler.wait_for_generation(timeout=30.0)
                sampler.absorb_generation(results)
                blob = sampler.export_runtime_state()
                restored = Distributor.set_method("AMMCMC")
                assert isinstance(restored, MCMCSampler)
                restored.set_config(cfg)
                restored.set_redis(queue)
                restored.import_runtime_state(blob)
                self.assertEqual(
                    restored._ensure_registry().all()[0].engine.iterations,
                    sampler._ensure_registry().all()[0].engine.iterations,
                )
            finally:
                stop.set()
                thread.join(timeout=2.0)


class CoreDRAMIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_dram_config_and_bounds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="DRAM",
                tmpdir=tmp,
                nchains=2,
                niters=3,
                workers=1,
                seed=5,
                extra_bounds={
                    "adapt_start_iter": 2,
                    "adapt_window": 2,
                    "dr_steps": 2,
                    "dr_scale_factors": [1.0, 0.5],
                },
            )
            sampler = Distributor.set_method("DRAM")
            assert isinstance(sampler, MCMCSampler)
            sampler.set_config(cfg)
            self.assertEqual(sampler.method, "DRAM")
            self.assertEqual(sampler._nchains, 2)
            self.assertEqual(sampler._niters, 3)
            self.assertEqual(sampler._dr_steps, 2)

    def test_core_dram_eggbox_tiny(self) -> None:
        """Jarvis2Core end-to-end: DRAM + Operas eggbox + LogL (V1 card shape)."""
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmp:
                cfg = _mcmc_config(
                    method="DRAM",
                    tmpdir=tmp,
                    redis_cfg=redis_config,
                    nchains=2,
                    niters=4,
                    workers=1,
                    seed=5,
                    extra_bounds={
                        "adapt_start_iter": 2,
                        "adapt_window": 2,
                        "dr_steps": 2,
                        "dr_scale_factors": [1.0, 0.5],
                    },
                )
                cfg["Runtime"]["sample_artifacts"] = "never"
                cfg["Runtime"]["Watchdog"] = {"enabled": False}
                core = Jarvis2Core(cfg)
                outcome = core.run()
                submitted = int(getattr(outcome, "submitted", 0) or 0)
                self.assertGreater(submitted, 0)
                self.assertTrue(getattr(outcome, "ok", True))
                sampler = core.sampler
                assert isinstance(sampler, MCMCSampler)
                summary = sampler.summary()
                self.assertEqual(summary["n_iters"], 4)
                self.assertGreater(summary["total_proposed"], 0)
        finally:
            server.shutdown()
            server.server_close()


if __name__ == "__main__":
    unittest.main()
