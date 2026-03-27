#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.config import ConfigValidator  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402
from jarvishep.Sampling.rltpmcmc import (  # noqa: E402
    RLTPActionMapper,
    RLTPMCMC,
    RLTPNumpyPPOTrainer,
    RLTPStateBuilder,
    RLTPTorchPPOTrainer,
)
from jarvishep.Sampling.Source.MCMC.state_machine_base import MCMCState  # noqa: E402


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
        self.uuid = f"rltpmcmc-fake-{type(self).seq:04d}"
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
        observables = dict(self.info.get("observables", {}))
        observables.update(self.params)
        observables["uuid"] = self.uuid
        self.info["observables"] = observables

    def close(self):
        type(self).close_calls += 1


class _ImmediateFactory:
    def __init__(self):
        self.calls = 0

    def submit_task(self, sample_info):
        self.calls += 1
        future = concurrent.futures.Future()
        params = sample_info.get("params", {})
        logl = float(sum(float(value) for value in params.values())) if isinstance(params, dict) else 0.0
        sample_info.setdefault("observables", {})
        sample_info["observables"]["LogL"] = logl
        future.set_result(logl)
        return future


class RLTPMCMCTests(unittest.TestCase):
    def _sample_info(self, task_root: str) -> dict:
        task_result_dir = os.path.join(task_root, "outputs", "rlscan")
        sample_dir = os.path.join(task_result_dir, "SAMPLE")
        os.makedirs(sample_dir, exist_ok=True)
        return {
            "task_result_dir": task_result_dir,
            "sample_dirs": sample_dir,
            "archive_samples": False,
        }

    def _config(self, *, controller="none", adaptation_mode="frozen", backend="torch") -> dict:
        return {
            "Scan": {
                "name": "rlscan",
                "save_dir": "&J/outputs",
                "sample_directory": {"limit": 20, "width": 4, "archive_samples": False},
            },
            "Sampling": {
                "Method": "RLTPMCMC",
                "Variables": [
                    {
                        "name": "x",
                        "description": "x",
                        "distribution": {"type": "Flat", "parameters": {"min": 0.0, "max": 1.0}},
                    }
                ],
                "LogLikelihood": [{"name": "LogL", "expression": "x"}],
                "Bounds": {
                    "num_chains": 3,
                    "num_iters": 3,
                    "exchange_interval": 1,
                    "proposal_scales": [0.1, 0.2, 0.3],
                    "temperature_ladder": [1.0, 2.0, 4.0],
                },
                "Control": {
                    "controller": controller,
                    "adaptation_mode": adaptation_mode,
                    "decision_interval": 1,
                    "smoothing": 0.5,
                    "max_action_magnitude": 0.15,
                    "freeze_after_burnin": False,
                    "burnin_iters": 0,
                    "deterministic_eval": False,
                    "min_beta_floor": 0.10,
                    "min_proposal_scale": 0.01,
                    "max_proposal_scale": 2.0
                },
                "Reward": {
                    "ess_weight": 1.0,
                    "round_trip_weight": 0.2,
                    "swap_weight": 0.5,
                    "discovery_weight": 1.0,
                    "swap_target": 0.25,
                    "discovery_metric": "best_loglike"
                },
                "PPO": {
                    "backend": backend,
                    "rollout_size": 1,
                    "epochs": 1,
                    "gamma": 0.95,
                    "gae_lambda": 0.90,
                    "clip_ratio": 0.2,
                    "learning_rate": 0.001,
                    "entropy_coef": 0.001,
                    "value_coef": 0.5,
                    "hidden_sizes": [8],
                    "seed": 7,
                    "device": "cpu",
                },
                "Diagnostics": {
                    "state_window": 8,
                    "reward_window": 8,
                    "log_interval": 1
                }
            },
        }

    def test_distributor_routes_rltpmcmc(self):
        sampler = Distributor.set_method("RLTPMCMC")
        self.assertEqual(sampler.method, "RLTPMCMC")
        self.assertTrue(sampler.schema.endswith("RLTPMCMC_schema.json"))

    def test_rltpmcmc_schema_accepts_rl_blocks(self):
        sampler = RLTPMCMC()
        validator = ConfigValidator()
        validator.logger = _NoopLogger()
        validator.set_schema(sampler.schema)
        validator.set_config(self._config(controller="ppo", adaptation_mode="train", backend="torch"))
        validator.validate_yaml()
        self.assertTrue(validator.passcheck)

    def test_rltpmcmc_vanilla_smoke_writes_rl_artifacts(self):
        with tempfile.TemporaryDirectory() as td:
            sampler = RLTPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.path["task_root"] = td
            sampler.path["jpath"] = td
            sampler.info["sample"] = self._sample_info(td)
            sampler.set_config(self._config(controller="none", adaptation_mode="frozen"))
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()
                sampler.finalize()

            rl_info = sampler.info["rl"]
            self.assertTrue(os.path.isdir(rl_info["outputs_root"]))
            self.assertTrue(os.path.isdir(rl_info["checkpoints_root"]))
            self.assertTrue(os.path.exists(os.path.join(rl_info["outputs_root"], "decision_metrics.jsonl")))
            self.assertTrue(os.path.exists(os.path.join(rl_info["outputs_root"], "manifest.json")))
            self.assertFalse(os.path.exists(os.path.join(rl_info["checkpoints_root"], "control_state.json")))
            self.assertFalse(os.path.exists(os.path.join(rl_info["checkpoints_root"], "policy.pt")))
            self.assertEqual(sampler.state.value, "TERMINATE")
            self.assertGreaterEqual(_FakeSample.close_calls, 1)

    def test_rltpmcmc_ppo_path_trains_and_records_metadata(self):
        with tempfile.TemporaryDirectory() as td:
            sampler = RLTPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.path["task_root"] = td
            sampler.path["jpath"] = td
            sampler.info["sample"] = self._sample_info(td)
            sampler.set_config(self._config(controller="ppo", adaptation_mode="train", backend="torch"))
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()
                sampler.finalize()

            rl_info = sampler.info["rl"]
            training_metrics_path = os.path.join(rl_info["outputs_root"], "training_metrics.jsonl")
            decision_metrics_path = os.path.join(rl_info["outputs_root"], "decision_metrics.jsonl")
            manifest_path = os.path.join(rl_info["outputs_root"], "manifest.json")

            self.assertTrue(os.path.exists(training_metrics_path))
            self.assertTrue(os.path.exists(decision_metrics_path))
            self.assertTrue(os.path.exists(manifest_path))
            self.assertFalse(os.path.exists(os.path.join(rl_info["checkpoints_root"], "policy.pt")))
            self.assertFalse(os.path.exists(os.path.join(rl_info["checkpoints_root"], "control_state.json")))

            with open(training_metrics_path, "r", encoding="utf-8") as handle:
                training_rows = [json.loads(line) for line in handle if line.strip()]
            with open(decision_metrics_path, "r", encoding="utf-8") as handle:
                decision_rows = [json.loads(line) for line in handle if line.strip()]
            with open(manifest_path, "r", encoding="utf-8") as handle:
                manifest = json.load(handle)

            self.assertGreaterEqual(len(training_rows), 1)
            self.assertGreaterEqual(len(decision_rows), 1)
            self.assertEqual(training_rows[0]["backend"], "torch")
            self.assertEqual(decision_rows[0]["ppo_backend"], "torch")
            self.assertEqual(manifest["ppo_backend"], "torch")
            self.assertEqual(manifest["control_state"]["ppo_backend"], "torch")
            self.assertIsNone(manifest["control_state"]["checkpoint"])
            for idx in range(3):
                snapshot = sampler.chain_snapshot(idx)
                self.assertGreater(float(snapshot["temperature"]), 0.0)
                self.assertGreater(float(snapshot["proposal_scale"]), 0.0)
            self.assertEqual(sampler.state.value, "TERMINATE")

    def test_rltpmcmc_numpy_reference_backend_trains_and_records_backend_identity(self):
        with tempfile.TemporaryDirectory() as td:
            sampler = RLTPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.path["task_root"] = td
            sampler.path["jpath"] = td
            sampler.info["sample"] = self._sample_info(td)
            sampler.set_config(self._config(controller="ppo", adaptation_mode="train", backend="numpy"))
            sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()
                sampler.finalize()

            rl_info = sampler.info["rl"]
            with open(os.path.join(rl_info["outputs_root"], "training_metrics.jsonl"), "r", encoding="utf-8") as handle:
                training_rows = [json.loads(line) for line in handle if line.strip()]
            with open(os.path.join(rl_info["outputs_root"], "decision_metrics.jsonl"), "r", encoding="utf-8") as handle:
                decision_rows = [json.loads(line) for line in handle if line.strip()]
            with open(os.path.join(rl_info["outputs_root"], "manifest.json"), "r", encoding="utf-8") as handle:
                manifest = json.load(handle)

            self.assertGreaterEqual(len(training_rows), 1)
            self.assertGreaterEqual(len(decision_rows), 1)
            self.assertEqual(training_rows[0]["backend"], "numpy")
            self.assertEqual(decision_rows[0]["ppo_backend"], "numpy")
            self.assertEqual(manifest["control_state"]["ppo_backend"], "numpy")
            self.assertEqual(sampler.state.value, "TERMINATE")

    def test_rltpmcmc_frozen_mode_loads_saved_checkpoint(self):
        with tempfile.TemporaryDirectory() as td:
            trainer_sampler = RLTPMCMC()
            trainer_sampler.set_logger(_NoopLogger())
            trainer_sampler.path["task_root"] = td
            trainer_sampler.path["jpath"] = td
            trainer_sampler.info["sample"] = self._sample_info(td)
            trainer_sampler.set_config(self._config(controller="ppo", adaptation_mode="train", backend="torch"))
            trainer_sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                trainer_sampler.initialize()
                trainer_sampler.run_nested()
                trainer_sampler.finalize()

            policy_path = trainer_sampler.save_rl_checkpoint("policy", trainer_sampler._trainer.export_state())
            self.assertTrue(os.path.exists(policy_path))

            frozen_cfg = self._config(controller="ppo", adaptation_mode="frozen", backend="torch")
            frozen_cfg["Sampling"]["PPO"]["resume_from"] = policy_path

            frozen_sampler = RLTPMCMC()
            frozen_sampler.set_logger(_NoopLogger())
            frozen_sampler.path["task_root"] = td
            frozen_sampler.path["jpath"] = td
            frozen_sampler.info["sample"] = self._sample_info(td)
            frozen_sampler.set_config(frozen_cfg)
            frozen_sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                frozen_sampler.initialize()
                frozen_sampler.run_nested()
                frozen_sampler.finalize()

            self.assertEqual(frozen_sampler.state.value, "TERMINATE")

    def test_rltpmcmc_resume_rejects_backend_mismatch(self):
        with tempfile.TemporaryDirectory() as td:
            trainer_sampler = RLTPMCMC()
            trainer_sampler.set_logger(_NoopLogger())
            trainer_sampler.path["task_root"] = td
            trainer_sampler.path["jpath"] = td
            trainer_sampler.info["sample"] = self._sample_info(td)
            trainer_sampler.set_config(self._config(controller="ppo", adaptation_mode="train", backend="numpy"))
            trainer_sampler.set_factory(_ImmediateFactory())

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                trainer_sampler.initialize()
                trainer_sampler.run_nested()
                trainer_sampler.finalize()

            policy_path = trainer_sampler.save_rl_checkpoint("policy", trainer_sampler._trainer.export_state())
            mismatch_cfg = self._config(controller="ppo", adaptation_mode="frozen", backend="torch")
            mismatch_cfg["Sampling"]["PPO"]["resume_from"] = policy_path

            sampler = RLTPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.path["task_root"] = td
            sampler.path["jpath"] = td
            sampler.info["sample"] = self._sample_info(td)
            sampler.set_config(mismatch_cfg)
            sampler.set_factory(_ImmediateFactory())

            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                with self.assertRaisesRegex(ValueError, "checkpoint backend mismatch"):
                    sampler.initialize()

    def test_rltpmcmc_runtime_checkpoint_roundtrip_restores_trainer_state(self):
        with tempfile.TemporaryDirectory() as td:
            checkpoint_root = os.path.join(td, "checkpoints", "rltpmcmc")

            sampler = RLTPMCMC()
            sampler.set_logger(_NoopLogger())
            sampler.path["task_root"] = td
            sampler.path["jpath"] = td
            sampler.info["sample"] = self._sample_info(td)
            sampler.set_config(self._config(controller="ppo", adaptation_mode="train", backend="numpy"))
            sampler.set_factory(_ImmediateFactory())
            sampler.configure_runtime_checkpointing(
                checkpoint_root,
                interval_seconds=3600.0,
                auto_resume=True,
                logger=_NoopLogger(),
            )

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                sampler.initialize()
                sampler.run_nested()
                sampler._state = MCMCState.UPDATE
                self.assertTrue(sampler.persist_runtime_checkpoint(force=True, reason="unit-test"))

            checkpoint_file = os.path.join(checkpoint_root, "state.pkl")
            self.assertTrue(os.path.exists(checkpoint_file))
            self.assertFalse(os.path.exists(os.path.join(checkpoint_root, "runtime_state.pkl")))
            original_update_count = sampler._trainer.update_count
            original_exchange_counter = sampler._exchange_counter
            original_replica_positions = list(sampler._replica_positions)

            restored = RLTPMCMC()
            restored.set_logger(_NoopLogger())
            restored.path["task_root"] = td
            restored.path["jpath"] = td
            restored.info["sample"] = self._sample_info(td)
            restored.set_config(self._config(controller="ppo", adaptation_mode="train", backend="numpy"))
            restored.set_factory(_ImmediateFactory())
            restored.set_runtime_checkpoint_resume_hint(True)
            restored.configure_runtime_checkpointing(
                checkpoint_root,
                interval_seconds=3600.0,
                auto_resume=True,
                logger=_NoopLogger(),
            )

            _FakeSample.reset()
            with patch("jarvishep.Sampling.Source.MCMC.state_machine_base.Sample", _FakeSample):
                restored.initialize()
                self.assertTrue(restored.restore_runtime_checkpoint_if_available())

            self.assertEqual(restored._exchange_counter, original_exchange_counter)
            self.assertEqual(restored._replica_positions, original_replica_positions)
            self.assertEqual(restored._trainer.update_count, original_update_count)
            self.assertEqual(restored._current_ppo_backend(), "numpy")
            self.assertEqual(restored._ppo_cfg["backend"], "numpy")

    def test_rltpmcmc_trainers_preserve_valid_control_patch_contract(self):
        summary = {
            "temperatures": [1.0, 2.0, 4.0],
            "proposal_scales": [0.1, 0.2, 0.3],
            "pair_acceptance_rates": [0.2, 0.3],
            "cold_ess_proxy": 8.0,
            "round_trip_proxy": 3.0,
            "best_improvement": 0.1,
            "best_logl": 0.5,
            "mean_acceptance": 0.4,
            "swap_acceptance_mean": 0.25,
            "training_enabled": True,
        }
        state_builder = RLTPStateBuilder(3, [0.1, 0.2, 0.3])
        state, _ = state_builder.build(summary)
        action_mapper = RLTPActionMapper(3, smoothing=0.5, max_action_magnitude=0.15)
        ppo_cfg = {"hidden_sizes": [8], "rollout_size": 1, "epochs": 1, "seed": 7, "device": "cpu"}

        trainers = [
            RLTPNumpyPPOTrainer(state_dim=state.size, action_dim=action_mapper.action_dim, ppo_cfg=dict(ppo_cfg)),
            RLTPTorchPPOTrainer(state_dim=state.size, action_dim=action_mapper.action_dim, ppo_cfg=dict(ppo_cfg)),
        ]
        for trainer in trainers:
            action, _, _ = trainer.select_action(state, deterministic=False)
            patch, _ = action_mapper.to_patch(summary, action)
            self.assertEqual(len(patch.temperature_ladder), 3)
            self.assertEqual(len(patch.proposal_scales), 3)
            self.assertAlmostEqual(float(patch.temperature_ladder[0]), 1.0, places=6)
            self.assertTrue(all(float(item) > 0.0 for item in patch.proposal_scales))
            self.assertTrue(
                all(
                    float(patch.temperature_ladder[idx]) < float(patch.temperature_ladder[idx + 1])
                    for idx in range(len(patch.temperature_ladder) - 1)
                )
            )


if __name__ == "__main__":
    unittest.main()
