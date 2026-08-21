#!/usr/bin/env python3
"""D13.2 MCMC / AM / DRAM FeedbackSampler tests."""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.slow

import tempfile
import threading
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.Sampling.Source.MCMC.config_contract import parse_common_chain_counts
from jarvishep2.Sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep2.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep2.Sampling.mcmc_base import MCMCBaseSampler
from jarvishep2.Sampling.mcmc_sampler import mcmc_sample_uuid
from jarvishep2.Sampling.ptmcmc import (
    PTMCMCSampler,
    normalize_temperature_ladder,
)
from jarvishep2.Sampling.toymcmc import ToyMCMCSampler
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
    bounds.setdefault("seed", seed)
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
    def test_methods_use_concrete_module_hierarchy(self) -> None:
        expected_modules = {
            "MCMC": "jarvishep2.sampling.mcmc",
            "ToyMCMC": "jarvishep2.sampling.toymcmc",
            "AMMCMC": "jarvishep2.sampling.ammcmc",
            "AM": "jarvishep2.sampling.am",
            "DRAM": "jarvishep2.sampling.dram",
            "EnsembleMCMC": "jarvishep2.sampling.ensemble_mcmc",
            "Ensemble": "jarvishep2.sampling.ensemble_mcmc",
            "DEMCMC": "jarvishep2.sampling.demcmc",
            "PTMCMC": "jarvishep2.sampling.ptmcmc",
            "PTEnsemble": "jarvishep2.sampling.ptensemble",
        }
        for method, module_name in expected_modules.items():
            with self.subTest(method=method):
                sampler = Distributor.set_method(method)
                self.assertIsInstance(sampler, MCMCBaseSampler)
                self.assertEqual(type(sampler).__module__, module_name)
                self.assertIsNot(type(sampler), MCMCBaseSampler)
                sampler.assert_checkpoint_attribute_contract()

    def test_methods_registered_stateful(self) -> None:
        for name in ("ToyMCMC", "MCMC", "AMMCMC", "AM", "DRAM", "PTMCMC"):
            self.assertNotIn(name, STATELESS_METHODS)
            sampler = Distributor.set_method(name)
            self.assertIsInstance(sampler, MCMCBaseSampler)
            self.assertEqual(Distributor.get_resume_status(name), "implemented")

    def test_ptmcmc_temperature_ladder_contract(self) -> None:
        ladder = normalize_temperature_ladder(
            [1.0, 1.5, 3.0, 6.0], nchains=4, sampler_method="PTMCMC"
        )
        self.assertEqual(ladder, [1.0, 1.5, 3.0, 6.0])
        with self.assertRaisesRegex(ValueError, "strictly increasing"):
            normalize_temperature_ladder(
                [1.0, 2.0, 2.0, 4.0], nchains=4, sampler_method="PTMCMC"
            )
        with self.assertRaisesRegex(ValueError, "cold temperature 1.0"):
            normalize_temperature_ladder(
                [1.1, 2.0, 3.0, 4.0], nchains=4, sampler_method="PTMCMC"
            )
        with self.assertRaisesRegex(ValueError, "num_chains"):
            normalize_temperature_ladder([1.0, 2.0], nchains=4, sampler_method="PTMCMC")

    def test_ptmcmc_config_and_adjacent_pairs(self) -> None:
        sampler = Distributor.set_method("PTMCMC")
        assert isinstance(sampler, PTMCMCSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="PTMCMC",
                tmpdir=tmp,
                nchains=4,
                niters=10,
                seed=3,
                extra_bounds={
                    "temperature_ladder": [1.0, 2.0, 4.0, 8.0],
                    "exchange_interval": 2,
                },
            )
            sampler.set_config(cfg)
            self.assertTrue(sampler._uses_pt())
            self.assertFalse(sampler._uses_async_chain_pipeline())
            self.assertEqual(sampler._exchange_interval, 2)
            self.assertEqual(sampler._temperature_ladder, [1.0, 2.0, 4.0, 8.0])
            registry = sampler._ensure_registry()
            self.assertEqual(
                [c.temperature for c in registry.all()],
                [1.0, 2.0, 4.0, 8.0],
            )
            self.assertTrue(registry.get(0).is_cold)
            # Even round: (0,1), (2,3); odd round: (1,2)
            sampler._exchange_offset = 0
            self.assertEqual(sampler._pair_chain_ids(), [(0, 1), (2, 3)])
            self.assertEqual(sampler._pair_chain_ids(), [(1, 2)])

    def test_ptmcmc_accepts_per_replica_scales_and_rejects_single_chain(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sampler = Distributor.set_method("PTMCMC")
            per_replica = _mcmc_config(
                method="PTMCMC",
                tmpdir=tmp,
                nchains=4,
                niters=5,
                extra_bounds={
                    "proposal_scale": [0.1, 0.2, 0.3, 0.4],
                    "temperature_ladder": [1.0, 2.0, 4.0, 8.0],
                },
            )
            sampler.set_config(per_replica)
            self.assertEqual(
                [chain.engine.proposal_scale for chain in sampler._ensure_registry().all()],
                [0.1, 0.2, 0.3, 0.4],
            )

            broadcast = Distributor.set_method("PTMCMC")
            scalar = _mcmc_config(
                method="PTMCMC",
                tmpdir=tmp,
                nchains=2,
                niters=5,
                proposal_scale=0.2,
                extra_bounds={"temperature_ladder": [1.0, 2.0]},
            )
            broadcast.set_config(scalar)
            self.assertEqual(
                [chain.engine.proposal_scale for chain in broadcast._ensure_registry().all()],
                [0.2, 0.2],
            )

            bad_length = _mcmc_config(
                method="PTMCMC",
                tmpdir=tmp,
                nchains=2,
                niters=5,
                extra_bounds={
                    "proposal_scale": [0.1, 0.2, 0.3],
                    "temperature_ladder": [1.0, 2.0],
                },
            )
            sampler2 = Distributor.set_method("PTMCMC")
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                sampler2.set_config(bad_length)

            sampler3 = Distributor.set_method("PTMCMC")
            bad_n = _mcmc_config(
                method="PTMCMC",
                tmpdir=tmp,
                nchains=1,
                niters=5,
                extra_bounds={"temperature_ladder": [1.0]},
            )
            with self.assertRaisesRegex(ValueError, "num_chains >= 2"):
                sampler3.set_config(bad_n)

    def test_toymcmc_native_pickle_checkpoint_roundtrip(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        assert isinstance(sampler, MCMCBaseSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="ToyMCMC",
                tmpdir=tmp,
                nchains=3,
                niters=4,
                workers=1,
            )
            sampler.set_config(cfg)
            sampler._ensure_registry()
            proposals = list(sampler.propose_generation() or [])
            sampler.absorb_generation(
                [
                    {
                        "uuid": sample.uuid,
                        "observables": {"LogL": -0.1},
                    }
                    for sample in proposals
                ]
            )
            expected_next = [
                np.asarray(chain.engine.param, dtype=float).copy()
                for chain in sampler._ensure_registry().all()
            ]
            state = sampler.export_runtime_state()
            self.assertEqual(state["native_mcmc_state_format"], "jarvis-hep.mcmc-native")
            self.assertIsInstance(state["native_mcmc_state_blob"], bytes)
            # Prove restore does not depend on rebuilding engines from the
            # explicit compatibility state.
            state.pop("chains")

            restored = Distributor.set_method("ToyMCMC")
            restored.set_config(cfg)
            restored.import_runtime_state(state)
            restored_chains = restored._ensure_registry().all()
            self.assertEqual(len(restored_chains), 3)
            for chain, expected in zip(restored_chains, expected_next):
                np.testing.assert_allclose(chain.engine.param, expected)
            self.assertEqual(len(restored_chains[0].history.all()), 1)

    def test_toymcmc_native_pickle_corruption_falls_back_to_explicit_state(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=3)
            sampler.set_config(cfg)
            sampler._ensure_registry()
            state = sampler.export_runtime_state()
            state["native_mcmc_state_blob"] = b"not-a-pickle"
            restored = Distributor.set_method("ToyMCMC")
            restored.set_config(cfg)
            restored.import_runtime_state(state)
            self.assertEqual(len(restored._ensure_registry().all()), 2)

    def test_checkpoint_rejects_num_chains_drift(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=4, niters=10, seed=1)
            sampler.set_config(cfg)
            sampler._ensure_registry()
            state = sampler.export_runtime_state()

            restored = Distributor.set_method("ToyMCMC")
            restored.set_config(
                _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=10, seed=1)
            )
            with self.assertRaisesRegex(ValueError, "num_chains"):
                restored.import_runtime_state(state)
            self.assertIsNone(restored._registry)
            self.assertFalse(restored._finished)

    def test_checkpoint_rejects_num_iters_and_seed_drift(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=20, seed=7)
            sampler.set_config(cfg)
            sampler._ensure_registry()
            state = sampler.export_runtime_state()

            restored_iters = Distributor.set_method("ToyMCMC")
            restored_iters.set_config(
                _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=5, seed=7)
            )
            with self.assertRaisesRegex(ValueError, "num_iters"):
                restored_iters.import_runtime_state(state)

            restored_seed = Distributor.set_method("ToyMCMC")
            restored_seed.set_config(
                _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=20, seed=99)
            )
            with self.assertRaisesRegex(ValueError, "seed"):
                restored_seed.import_runtime_state(state)

    def test_async_pipeline_stall_raises(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=1, niters=3)
            sampler.set_config(cfg)
            sampler.set_redis(make_fakeredis_queue())
            registry = sampler._ensure_registry()
            # Force an impossible state: unfinished engine but no ready proposals.
            chain = registry.all()[0]
            chain.engine.iterations = 0
            chain.engine.n_iterations = 3
            chain.open_stage = None
            sampler._pending_uuids.clear()
            # Monkey-patch ready set empty while _all_finished is false.
            sampler._ready_chain_ids = lambda: []  # type: ignore[method-assign]
            with self.assertRaisesRegex(RuntimeError, "async pipeline stalled"):
                sampler._run_async_independent_chains(timeout=1.0)
            self.assertFalse(sampler._finished)

    def test_mcmc_sample_carries_step_stage_on_wire(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=1, niters=2)
            sampler.set_config(cfg)
            samples = list(sampler.propose_generation() or [])
            self.assertEqual(len(samples), 1)
            task = samples[0].to_task_dict()
            self.assertEqual(task["chain_id"], 0)
            self.assertEqual(task["step"], 0)
            self.assertEqual(task["stage"], 0)
            self.assertNotIn("observables", task)

    def test_toymcmc_sampling_checkpoint_does_not_require_archive_ack(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=3)
            sampler.set_config(cfg)
            sampler._ensure_registry()
            sampler._register_pending("unrelated-archive-uuid")
            sampler._pending_uuids.clear()
            sampler._submitted_uuids.append("unrelated-archive-uuid")
            calls: list[str] = []
            sampler._save_checkpoint_callback = lambda **kwargs: calls.append(
                str(kwargs.get("reason"))
            ) or True
            self.assertTrue(
                sampler._checkpoint_at_sampling_barrier(
                    reason="ToyMCMC_sampling_step_1"
                )
            )
            self.assertEqual(calls, ["ToyMCMC_sampling_step_1"])

    def test_toymcmc_broadcasts_scalar_scale_to_all_chains(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        assert isinstance(sampler, MCMCBaseSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="ToyMCMC",
                tmpdir=tmp,
                nchains=4,
                niters=5,
                proposal_scale=0.17,
            )
            sampler.set_config(cfg)
            registry = sampler._ensure_registry()
            chains = registry.all()
            self.assertEqual(len(chains), 4)
            self.assertEqual(
                {float(chain.engine.proposal_scale) for chain in chains},
                {0.17},
            )
            self.assertEqual(
                {int(chain.engine.n_iterations) for chain in chains},
                {5},
            )
            self.assertEqual(len({id(chain.engine) for chain in chains}), 4)

    def test_toymcmc_accepts_per_chain_scales_and_rejects_bad_length(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        assert isinstance(sampler, MCMCBaseSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="ToyMCMC",
                tmpdir=tmp,
                nchains=2,
                niters=5,
                proposal_scale=0.17,
            )
            cfg["Sampling"]["Bounds"]["proposal_scale"] = [0.17, 0.2]
            sampler.set_config(cfg)
            scales = [
                float(chain.engine.proposal_scale)
                for chain in sampler._ensure_registry().all()
            ]
            self.assertEqual(scales, [0.17, 0.2])

        sampler = Distributor.set_method("ToyMCMC")
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(
                method="ToyMCMC",
                tmpdir=tmp,
                nchains=2,
                niters=5,
                proposal_scale=0.17,
            )
            cfg["Sampling"]["Bounds"]["proposal_scale"] = [0.17, 0.2, 0.3]
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                sampler.set_config(cfg)

    def test_toymcmc_progress_logs_permille_and_percent(self) -> None:
        sampler = Distributor.set_method("ToyMCMC")
        assert isinstance(sampler, ToyMCMCSampler)
        assert isinstance(sampler, MCMCBaseSampler)
        with tempfile.TemporaryDirectory() as tmp:
            sampler.set_config(
                _mcmc_config(
                    method="ToyMCMC",
                    tmpdir=tmp,
                    nchains=1,
                    niters=1000,
                    workers=1,
                )
            )
            info_messages: list[str] = []
            warning_messages: list[str] = []

            class _Capture:
                def info(self, msg, *args, **kwargs):  # noqa: ANN001
                    info_messages.append(str(msg) % args if args else str(msg))

                def warning(self, msg, *args, **kwargs):  # noqa: ANN001
                    warning_messages.append(str(msg) % args if args else str(msg))

            sampler._logger = _Capture()  # type: ignore[assignment]
            registry = sampler._ensure_registry()
            sampler._emit_progress()
            for iteration in range(1, 11):
                registry.get(0).engine.iterations = iteration
                sampler._emit_progress()

            self.assertTrue(
                any("Initializing the ToyMCMC Sampling" in line for line in warning_messages)
            )
            self.assertTrue(
                any("ToyMCMC Sampler configured" in line for line in info_messages)
            )
            self.assertTrue(
                any(
                    line.startswith("0‰ of 0/1000 ToyMCMC transitions completed")
                    for line in info_messages
                )
            )
            for permille in range(1, 10):
                self.assertTrue(
                    any(line.startswith(f"{permille}‰ of {permille}/1000") for line in info_messages),
                    f"missing ToyMCMC info heartbeat for {permille}‰",
                )
            self.assertTrue(any(line.startswith("10‰ of 10/1000") for line in warning_messages))
            self.assertFalse(any(line.startswith("10‰") for line in info_messages))
            registry.get(0).engine.iterations = 1000
            sampler._emit_progress()
            self.assertTrue(any(line.startswith("1000‰ of 1000/1000") for line in warning_messages))

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
    ) -> MCMCBaseSampler:
        queue = make_fakeredis_queue()
        sampler = Distributor.set_method(method)
        assert isinstance(sampler, MCMCBaseSampler)
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
                        },
                        queue=task.get("feedback_queue"),
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

    def test_toymcmc_feedback_loop(self) -> None:
        s = self._run_with_mock_worker("ToyMCMC", nchains=3, niters=6)
        self.assertEqual(s.method, "ToyMCMC")
        self.assertEqual(s.summary()["n_chains"], 3)
        self.assertTrue(s.summary().get("async_chains"))
        self.assertTrue(s.summary().get("sharded_feedback"))

    def test_toymcmc_independent_chain_pipeline(self) -> None:
        """ToyMCMC uses the base independent-chain async design."""
        s = self._run_with_mock_worker("ToyMCMC", nchains=4, niters=5)
        self.assertTrue(s._uses_async_chain_pipeline())
        self.assertTrue(s._uses_sharded_feedback())
        self.assertTrue(s._chains_are_independent())
        # Every chain completed its own length (not a forced global barrier).
        for chain in s._ensure_registry().all():
            self.assertEqual(int(chain.engine.iterations), 5)

    def test_toymcmc_task_carries_chain_feedback_route(self) -> None:
        queue = make_fakeredis_queue()
        sampler = Distributor.set_method("ToyMCMC")
        assert isinstance(sampler, ToyMCMCSampler)
        with tempfile.TemporaryDirectory() as tmp:
            cfg = _mcmc_config(method="ToyMCMC", tmpdir=tmp, nchains=2, niters=2)
            sampler.set_config(cfg)
            sampler.set_redis(queue)
            sampler._ensure_registry()
            props = list(sampler.propose_generation() or [])
            self.assertEqual(len(props), 2)
            for sample in props:
                task = sample.to_task_dict()
                self.assertIn("chain_id", task)
                self.assertEqual(
                    task["feedback_queue"],
                    f"hep:feedback:chain:{task['chain_id']}",
                )
            # Submit and complete via sharded feedback only.
            sampler._submit_sample_batch(props)

            def worker_once() -> None:
                for _ in range(2):
                    task = queue.pull_task(timeout=2)
                    self.assertIsNotNone(task)
                    assert task is not None
                    self.assertTrue(str(task["feedback_queue"]).startswith("hep:feedback:chain:"))
                    queue.publish_feedback(
                        {"uuid": task["uuid"], "logL": -1.0},
                        queue=task.get("feedback_queue"),
                    )

            thr = threading.Thread(target=worker_once, daemon=True)
            thr.start()
            try:
                results = sampler.wait_for_generation(
                    timeout=30.0, queues=sampler._pending_feedback_queues()
                )
                self.assertEqual(len(results), 2)
            finally:
                thr.join(timeout=2.0)

    def test_dram_feedback_loop(self) -> None:
        s = self._run_with_mock_worker("DRAM", nchains=3, niters=6)
        self.assertEqual(s.method, "DRAM")
        self.assertTrue(s.summary().get("async_chains"))

    def test_worker_count_independent_trajectory(self) -> None:
        """Same seed → same chain final params for 1 vs 3 workers (mock LogL)."""

        def finals(workers: int) -> list[list[float]]:
            queue = make_fakeredis_queue()
            sampler = Distributor.set_method("MCMC")
            assert isinstance(sampler, MCMCBaseSampler)
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
                            },
                            queue=task.get("feedback_queue"),
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
        assert isinstance(sampler, MCMCBaseSampler)
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
                        },
                        queue=task.get("feedback_queue"),
                    )

            thread = threading.Thread(target=worker, daemon=True)
            thread.start()
            try:
                # Run a few steps by manually driving one generation then export.
                sampler._ensure_registry()
                props = sampler.propose_generation()
                assert props is not None
                sampler._submit_sample_batch(list(props))
                results = sampler.wait_for_generation(
                    timeout=30.0, queues=sampler._pending_feedback_queues()
                )
                sampler.absorb_generation(results)
                blob = sampler.export_runtime_state()
                restored = Distributor.set_method("AMMCMC")
                assert isinstance(restored, MCMCBaseSampler)
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
            assert isinstance(sampler, MCMCBaseSampler)
            sampler.set_config(cfg)
            self.assertEqual(sampler.method, "DRAM")
            self.assertEqual(sampler._nchains, 2)
            self.assertEqual(sampler._niters, 3)
            self.assertEqual(sampler._dr_steps, 2)

if __name__ == "__main__":
    unittest.main()
