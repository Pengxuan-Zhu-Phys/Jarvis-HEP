#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.Source.MCMC.chain_history import ChainEvent, ChainHistory  # noqa: E402
from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime  # noqa: E402
from jarvishep.Sampling.Source.MCMC.controller import MCMCControlGuard, MCMCControlPatch  # noqa: E402
from jarvishep.Sampling.Source.MCMC.dram_chain import DRAMChain  # noqa: E402
from jarvishep.Sampling.ammcmc import AMMCMC  # noqa: E402
from jarvishep.Sampling.demcmc import DEMCMC  # noqa: E402
from jarvishep.Sampling.dream_lite import DREAMLite  # noqa: E402
from jarvishep.Sampling.dram import DRAM  # noqa: E402
from jarvishep.Sampling.ensemblemcmc import EnsembleMCMC  # noqa: E402
from jarvishep.Sampling.mcmc_standard import MCMC  # noqa: E402
from jarvishep.Sampling.pt_ensemble import PTEnsemble  # noqa: E402
from jarvishep.Sampling.robustam import RobustAM  # noqa: E402
from jarvishep.Sampling.tpmcmc import TPMCMC  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _FakeSample:
    close_calls = 0
    seq = 0

    @classmethod
    def reset(cls):
        cls.close_calls = 0
        cls.seq = 0

    def __init__(self, params):
        type(self).seq += 1
        self.uuid = f"mcmc-fake-{type(self).seq:04d}"
        self.params = dict(params)
        self.info = {
            "uuid": self.uuid,
            "params": dict(self.params),
            "observables": {**self.params, "uuid": self.uuid},
        }

    def set_config(self, config):
        if isinstance(config, dict):
            self.info.update(dict(config))
        self.info["uuid"] = self.uuid
        self.info["params"] = dict(self.params)
        obs = dict(self.info.get("observables", {}))
        obs.update(self.params)
        obs["uuid"] = self.uuid
        self.info["observables"] = obs

    def close(self):
        type(self).close_calls += 1


class _ImmediateFactory:
    def __init__(self):
        self.calls = 0

    def submit_task(self, sample_info):
        self.calls += 1
        fut = concurrent.futures.Future()
        params = sample_info.get("params", {})
        logl = float(sum(float(v) for v in params.values())) if isinstance(params, dict) else 0.0
        sample_info.setdefault("observables", {})
        sample_info["observables"]["LogL"] = logl
        fut.set_result(logl)
        return fut


class _BootstrapController:
    def __init__(self, ladder, exchange_interval):
        self._ladder = ladder
        self._exchange_interval = exchange_interval

    def on_run_start(self, _context):
        return MCMCControlPatch(
            temperature_ladder=self._ladder,
            exchange_interval=self._exchange_interval,
        )

    def on_pre_step(self, _snapshot):
        return None

    def on_post_step(self, _snapshot, _outcome):
        return None

    def on_pre_exchange(self, _snapshot):
        return None

    def on_post_exchange(self, _snapshot, _exchange_metrics):
        return None


class MCMCStateMachineTests(unittest.TestCase):
    def test_chain_history_tail(self):
        hist = ChainHistory()
        hist.append(
            ChainEvent(
                iter=1,
                state="accepted",
                proposal=[0.1],
                logl=1.0,
                accepted=True,
                temperature=1.0,
            )
        )
        hist.append(
            ChainEvent(
                iter=2,
                state="rejected",
                proposal=[0.2],
                logl=0.9,
                accepted=False,
                temperature=1.0,
            )
        )
        self.assertEqual(len(hist.all()), 2)
        self.assertEqual(len(hist.tail(1)), 1)
        self.assertEqual(hist.tail(1)[0].iter, 2)

    def test_chain_registry_cold_lookup(self):
        c0 = ChainRuntime(chain_id=0, engine=object(), is_cold=True)
        c1 = ChainRuntime(chain_id=1, engine=object(), is_cold=False)
        reg = ChainRegistry([c0, c1])
        self.assertTrue(reg.is_cold(0))
        self.assertFalse(reg.is_cold(1))
        reg.mark_cold(1, True)
        self.assertTrue(reg.is_cold(1))

    def test_control_guard_validation(self):
        ok, _ = MCMCControlGuard.validate_patch(
            MCMCControlPatch(temperature_ladder=[1.0, 2.0]), nchains=2
        )
        self.assertTrue(ok)

        ok_bad, _ = MCMCControlGuard.validate_patch(
            MCMCControlPatch(temperature_ladder=[1.0]), nchains=2
        )
        self.assertFalse(ok_bad)

    def test_mcmc_state_machine_history_api(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 2,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = MCMC()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertGreaterEqual(_FakeSample.close_calls, 1)
            self.assertTrue(sampler.is_cold_chain(0))
            self.assertEqual(len(sampler.history_all(0)), 4)
            self.assertEqual(len(sampler.history_tail(0, 2)), 2)
            metrics = sampler.metrics_latest()
            self.assertIsNotNone(metrics)
            self.assertIsNotNone(metrics.ess_proxy_mean)
            self.assertIsNotNone(metrics.autocorr_lag1_proxy_mean)

    def test_tpmcmc_controller_patch_applied(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 3,
                    "num_iters": 3,
                    "exchange_interval": 1,
                    "proposal_scales": [0.1, 0.2, 0.3],
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }
        ladder = [1.0, 3.0, 6.0]

        with tempfile.TemporaryDirectory() as td:
            sampler = TPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())
            sampler.set_controller(_BootstrapController(ladder=ladder, exchange_interval=1))

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            temps = [sampler.chain_snapshot(i)["temperature"] for i in range(3)]
            self.assertEqual(temps, ladder)
            self.assertTrue(sampler.is_cold_chain(0))
            self.assertGreaterEqual(_FakeSample.close_calls, 1)
            self.assertTrue(any(frame.event == "post_exchange" for frame in sampler.metrics_all()))

    def test_distributor_supports_toymcmc(self):
        sampler = Distributor.set_method("ToyMCMC")
        self.assertEqual(sampler.method, "ToyMCMC")

    def test_distributor_supports_ammcmc(self):
        sampler = Distributor.set_method("AMMCMC")
        self.assertEqual(sampler.method, "AMMCMC")

    def test_distributor_supports_robustam(self):
        sampler = Distributor.set_method("RobustAM")
        self.assertEqual(sampler.method, "RobustAM")

    def test_distributor_supports_dram(self):
        sampler = Distributor.set_method("DRAM")
        self.assertEqual(sampler.method, "DRAM")

    def test_distributor_supports_demcmc(self):
        sampler = Distributor.set_method("DEMCMC")
        self.assertEqual(sampler.method, "DEMCMC")

    def test_distributor_supports_dream_lite(self):
        sampler = Distributor.set_method("DREAMLite")
        self.assertEqual(sampler.method, "DREAMLite")

    def test_distributor_supports_ensemblemcmc(self):
        sampler = Distributor.set_method("EnsembleMCMC")
        self.assertEqual(sampler.method, "EnsembleMCMC")

    def test_distributor_supports_ptensemble(self):
        sampler = Distributor.set_method("PTEnsemble")
        self.assertEqual(sampler.method, "PTEnsemble")

    def test_ammcmc_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 2,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "adapt_enabled": True,
                    "adapt_start_iter": 2,
                    "adapt_window": 1,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = AMMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 8)
            self.assertEqual(_FakeSample.close_calls, 8)

    def test_demcmc_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 3,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "de_gamma": 0.0,
                    "de_noise": 1.0e-3,
                    "de_crossover": 0.9,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = DEMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 12)
            self.assertEqual(_FakeSample.close_calls, 12)

    def test_dream_lite_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 3,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "de_gamma": 0.0,
                    "de_noise": 1.0e-3,
                    "de_crossover": 0.9,
                    "dream_snooker_prob": 0.2,
                    "dream_archive_size": 64,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = DREAMLite()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 12)
            self.assertEqual(_FakeSample.close_calls, 12)

    def test_ensemblemcmc_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 4,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "stretch_a": 2.0,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = EnsembleMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 16)
            self.assertEqual(_FakeSample.close_calls, 16)

    def test_ptensemble_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 4,
                    "num_iters": 4,
                    "exchange_interval": 1,
                    "proposal_scales": [0.1, 0.1, 0.1, 0.1],
                    "temperature_ladder": [1.0, 1.6, 2.4, 3.6],
                    "stretch_a": 2.0,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = PTEnsemble()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 16)
            self.assertEqual(_FakeSample.close_calls, 16)

    def test_robustam_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 2,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "adapt_enabled": True,
                    "adapt_start_iter": 2,
                    "adapt_window": 1,
                    "global_jump_prob": 0.2,
                    "heavy_tail_df": 4.0,
                    "global_scale": 2.0,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = RobustAM()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertEqual(sampler.factory.calls, 8)
            self.assertEqual(_FakeSample.close_calls, 8)

    def test_dram_state_machine_smoke(self):
        cfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 2,
                    "num_iters": 4,
                    "proposal_scale": 0.1,
                    "adapt_enabled": True,
                    "adapt_start_iter": 2,
                    "adapt_window": 1,
                    "dr_scale_factors": [1.0, 0.5],
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            sampler = DRAM()
            sampler.set_logger(_NoopLogger())
            sampler.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(sampler.info["sample"]["sample_dirs"], exist_ok=True)
            sampler.set_config(cfg)
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch(
                "jarvishep.Sampling.Source.MCMC.state_machine_base.Sample",
                _FakeSample,
            ), patch(
                "jarvishep.Sampling.Source.MCMC.state_machine_multistage_base.Sample",
                _FakeSample,
            ):
                sampler.initialize()
                sampler.run_nested()

            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertGreaterEqual(sampler.factory.calls, 8)
            self.assertEqual(_FakeSample.close_calls, sampler.factory.calls)

    def test_dram_chain_uses_second_stage_after_first_rejection(self):
        chain = DRAMChain(
            initial_param=np.array([0.5], dtype=float),
            proposal_scale=0.2,
            n_iterations=4,
            adapt_enabled=False,
            dr_steps=2,
            dr_scale_factors=[1.0, 0.4],
        )

        chain.propose_stage(0)
        first = chain.consume_stage_result(0, 0.0, beta=1.0)
        self.assertTrue(first["iteration_done"])
        self.assertTrue(first["accepted"])

        chain.propose_stage(0)
        second = chain.consume_stage_result(0, -1e6, beta=1.0)
        self.assertFalse(second["iteration_done"])
        self.assertEqual(second["next_stage"], 1)

        chain.propose_stage(1)
        third = chain.consume_stage_result(1, -10.0, beta=1.0)
        self.assertTrue(third["iteration_done"])
        self.assertEqual(third["stage_attempts"], 2)

    def test_inflight_is_bounded_by_max_workers(self):
        mcfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 8,
                    "num_iters": 1,
                    "proposal_scale": 0.1,
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }
        tcfg = {
            "Sampling": {
                "Bounds": {
                    "num_chains": 4,
                    "num_iters": 1,
                    "exchange_interval": 1,
                    "proposal_scales": [0.1, 0.1, 0.1, 0.1],
                },
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
            },
            "Scan": {"sample_directory": {"limit": 20, "width": 4}},
        }

        with tempfile.TemporaryDirectory() as td:
            ms = MCMC()
            ms.set_logger(_NoopLogger())
            ms.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(ms.info["sample"]["sample_dirs"], exist_ok=True)
            ms.set_config(mcfg)
            ms.set_max_workers(3)
            ms.initialize()
            self.assertEqual(ms._max_inflight(), 3)

            ts = TPMCMC()
            ts.set_logger(_NoopLogger())
            ts.info["sample"] = {
                "task_result_dir": td,
                "sample_dirs": os.path.join(td, "SAMPLE"),
                "archive_samples": False,
            }
            os.makedirs(ts.info["sample"]["sample_dirs"], exist_ok=True)
            ts.set_config(tcfg)
            ts.set_max_workers(2)
            ts.initialize()
            self.assertEqual(ts._max_inflight(), 2)


if __name__ == "__main__":
    unittest.main()
