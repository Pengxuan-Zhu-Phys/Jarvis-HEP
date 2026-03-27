#!/usr/bin/env python3
from __future__ import annotations

from collections import deque
import importlib
import math
from typing import Any, Dict, List, Sequence

import numpy as np
from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.Source.MCMC.controller import MCMCControlPatch
from jarvishep.Sampling.rl_sampler_base import RLSamplerBase, make_json_safe
from jarvishep.Sampling.tpmcmc import PTMCMC


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return float(default)


def _safe_mean(values: Sequence[float], default: float = 0.0) -> float:
    values = [float(v) for v in values if v is not None]
    if not values:
        return float(default)
    return float(sum(values) / len(values))


def _load_torch_module():
    try:
        return importlib.import_module("torch")
    except Exception as exc:
        raise RuntimeError(
            "RLTPMCMC torch backend initialization failed. "
            "Confirm that torch imports cleanly in the target Jarvis runtime context."
        ) from exc


class RLTPDiagnostics:
    def __init__(self, nchains: int, history_size: int = 64) -> None:
        self.nchains = int(max(2, nchains))
        self.history_size = int(max(4, history_size))
        self.exchange_history = deque(maxlen=self.history_size)
        self.current_window: List[Dict[str, Any]] = []
        self.best_logl: float | None = None
        self.window_start_best_logl: float | None = None
        self.window_best_improvement: float = 0.0

    def start_window(self) -> None:
        self.current_window = []
        self.window_start_best_logl = self.best_logl
        self.window_best_improvement = 0.0

    def update_best(self, logl: float) -> None:
        logl = float(logl)
        if self.best_logl is None:
            self.best_logl = logl
            if self.window_start_best_logl is None:
                self.window_start_best_logl = logl
            return
        if logl > self.best_logl:
            self.window_best_improvement += float(logl - self.best_logl)
            self.best_logl = logl

    def record_exchange(self, record: Dict[str, Any]) -> None:
        safe = make_json_safe(record)
        self.exchange_history.append(safe)
        self.current_window.append(safe)

    def _pair_rates_from(self, records: Sequence[Dict[str, Any]]) -> List[float]:
        attempted = [0 for _ in range(self.nchains - 1)]
        accepted = [0 for _ in range(self.nchains - 1)]
        for record in records:
            for pair_metric in record.get("pair_metrics", []):
                pair_index = int(pair_metric.get("pair_index", -1))
                if 0 <= pair_index < self.nchains - 1:
                    attempted[pair_index] += int(pair_metric.get("attempted", 0))
                    accepted[pair_index] += int(pair_metric.get("accepted", 0))
        return [
            float(accepted[ii]) / float(attempted[ii]) if attempted[ii] > 0 else 0.0
            for ii in range(self.nchains - 1)
        ]

    def collect_window_summary(self, sampler, *, training_enabled: bool) -> Dict[str, Any]:
        records = list(self.current_window) if self.current_window else list(self.exchange_history)
        pair_rates = self._pair_rates_from(records)
        round_trip_proxy = _safe_mean(
            [record.get("round_trip_proxy") for record in records],
            default=float(self.nchains),
        )
        swap_acceptance_mean = _safe_mean(pair_rates, default=0.0)

        registry = sampler._must_registry()
        chains = registry.all()
        temperatures = [float(chain.temperature) for chain in chains]
        proposal_scales = [
            _safe_float(getattr(chain.engine, "proposal_scale", None), default=0.1)
            for chain in chains
        ]
        acc_rates = []
        for chain in chains:
            denom = int(chain.accepted) + int(chain.rejected)
            acc_rates.append(float(chain.accepted) / float(denom) if denom > 0 else 0.0)

        cold_ids = registry.cold_chain_ids()
        cold_ess = 0.0
        if cold_ids:
            cold_chain = registry.get(cold_ids[0])
            cold_ess, _ = sampler._estimate_chain_proxies(cold_chain)
            cold_ess = 0.0 if cold_ess is None else float(cold_ess)

        return {
            "temperatures": temperatures,
            "proposal_scales": proposal_scales,
            "pair_acceptance_rates": pair_rates,
            "cold_ess_proxy": cold_ess,
            "round_trip_proxy": float(round_trip_proxy),
            "best_improvement": float(self.window_best_improvement),
            "best_logl": self.best_logl,
            "mean_acceptance": _safe_mean(acc_rates, default=0.0),
            "swap_acceptance_mean": float(swap_acceptance_mean),
            "replica_positions": list(getattr(sampler, "_replica_positions", [])),
            "training_enabled": bool(training_enabled),
        }


class RLTPStateBuilder:
    def __init__(self, nchains: int, initial_proposal_scales: Sequence[float]) -> None:
        self.nchains = int(max(2, nchains))
        self.initial_proposal_scales = np.asarray(initial_proposal_scales, dtype=np.float32)

    def build(self, summary: Dict[str, Any]) -> tuple[np.ndarray, Dict[str, Any]]:
        temperatures = np.asarray(summary["temperatures"], dtype=np.float32)
        temperatures = np.clip(temperatures, 1.0e-6, None)
        betas = 1.0 / temperatures
        gaps = np.clip(betas[:-1] - betas[1:], 1.0e-6, None)
        gaps = gaps / float(np.sum(gaps))

        pair_rates = np.asarray(summary["pair_acceptance_rates"], dtype=np.float32)
        if pair_rates.size != self.nchains - 1:
            pair_rates = np.zeros(self.nchains - 1, dtype=np.float32)

        proposal_scales = np.asarray(summary["proposal_scales"], dtype=np.float32)
        proposal_scales = np.clip(proposal_scales, 1.0e-6, None)
        reference = np.clip(self.initial_proposal_scales, 1.0e-6, None)
        proposal_features = np.log(proposal_scales / reference)

        cold_ess = _safe_float(summary.get("cold_ess_proxy"), 0.0)
        round_trip = _safe_float(summary.get("round_trip_proxy"), float(self.nchains))
        best_improvement = _safe_float(summary.get("best_improvement"), 0.0)
        best_logl = _safe_float(summary.get("best_logl"), 0.0)
        mean_accept = _safe_float(summary.get("mean_acceptance"), 0.0)
        swap_mean = _safe_float(summary.get("swap_acceptance_mean"), 0.0)
        training_enabled = 1.0 if bool(summary.get("training_enabled", False)) else 0.0

        summary_features = np.asarray(
            [
                np.tanh(cold_ess / 32.0),
                np.tanh(round_trip / float(max(1, self.nchains * 4))),
                np.tanh(best_improvement),
                np.tanh(best_logl / 50.0),
                np.clip(mean_accept, 0.0, 1.0),
                np.clip(swap_mean, 0.0, 1.0),
                training_enabled,
            ],
            dtype=np.float32,
        )

        state = np.concatenate([gaps, pair_rates, proposal_features, summary_features], axis=0)
        meta = {
            "beta_gaps": gaps.tolist(),
            "pair_acceptance_rates": pair_rates.tolist(),
            "proposal_features": proposal_features.tolist(),
            "summary_features": summary_features.tolist(),
        }
        return state.astype(np.float32), meta


class RLTPActionMapper:
    def __init__(
        self,
        nchains: int,
        *,
        smoothing: float = 1.0,
        max_action_magnitude: float = 0.20,
        min_beta_floor: float = 1.0e-3,
        min_proposal_scale: float = 1.0e-3,
        max_proposal_scale: float = 5.0,
    ) -> None:
        self.nchains = int(max(2, nchains))
        self.smoothing = float(np.clip(smoothing, 0.0, 1.0))
        self.max_action_magnitude = float(max(1.0e-6, max_action_magnitude))
        self.min_beta_floor = float(max(1.0e-6, min_beta_floor))
        self.min_proposal_scale = float(max(1.0e-6, min_proposal_scale))
        self.max_proposal_scale = float(max(self.min_proposal_scale, max_proposal_scale))

    @property
    def action_dim(self) -> int:
        return int((self.nchains - 1) + self.nchains)

    def to_patch(self, summary: Dict[str, Any], action: Sequence[float]) -> tuple[MCMCControlPatch, Dict[str, Any]]:
        action = np.asarray(action, dtype=np.float32).reshape(-1)
        if action.size != self.action_dim:
            raise ValueError(f"RLTPMCMC action size mismatch: expect {self.action_dim}, got {action.size}")

        ladder_delta = np.clip(action[: self.nchains - 1], -self.max_action_magnitude, self.max_action_magnitude)
        proposal_delta = np.clip(action[self.nchains - 1 :], -self.max_action_magnitude, self.max_action_magnitude)

        temperatures = np.asarray(summary["temperatures"], dtype=np.float32)
        temperatures = np.clip(temperatures, 1.0e-6, None)
        current_betas = 1.0 / temperatures
        current_gaps = np.clip(current_betas[:-1] - current_betas[1:], 1.0e-6, None)
        current_total_gap = float(np.sum(current_gaps))
        log_gaps = np.log(current_gaps)
        smoothed_log_gaps = log_gaps + (self.smoothing * ladder_delta)
        new_gaps = np.exp(smoothed_log_gaps)

        beta_floor = max(self.min_beta_floor, float(current_betas[-1]))
        target_total_gap = float(max(1.0e-6, 1.0 - beta_floor))
        new_gaps = new_gaps / float(np.sum(new_gaps)) * target_total_gap

        new_betas = [1.0]
        for gap in new_gaps:
            new_betas.append(new_betas[-1] - float(gap))
        new_betas[-1] = max(beta_floor, self.min_beta_floor)
        for ii in range(1, len(new_betas)):
            new_betas[ii] = max(self.min_beta_floor, min(new_betas[ii - 1] - 1.0e-6, new_betas[ii]))

        new_temperatures = [1.0 / float(beta) for beta in new_betas]

        current_scales = np.asarray(summary["proposal_scales"], dtype=np.float32)
        current_scales = np.clip(current_scales, self.min_proposal_scale, self.max_proposal_scale)
        new_scales = current_scales * np.exp(self.smoothing * proposal_delta)
        new_scales = np.clip(new_scales, self.min_proposal_scale, self.max_proposal_scale)

        patch = MCMCControlPatch(
            temperature_ladder=[float(value) for value in new_temperatures],
            proposal_scales=[float(value) for value in new_scales],
        )
        details = {
            "ladder_delta": ladder_delta.tolist(),
            "proposal_delta": proposal_delta.tolist(),
            "temperatures_before": temperatures.tolist(),
            "temperatures_after": [float(value) for value in new_temperatures],
            "proposal_scales_before": current_scales.tolist(),
            "proposal_scales_after": [float(value) for value in new_scales],
            "current_total_gap": float(current_total_gap),
        }
        return patch, details


class RLTPRewardBuilder:
    def __init__(self, reward_cfg: Dict[str, Any], nchains: int) -> None:
        self.nchains = int(max(2, nchains))
        self.ess_weight = _safe_float(reward_cfg.get("ess_weight", 1.0), 1.0)
        self.round_trip_weight = _safe_float(reward_cfg.get("round_trip_weight", 0.5), 0.5)
        self.swap_weight = _safe_float(reward_cfg.get("swap_weight", 0.5), 0.5)
        self.discovery_weight = _safe_float(reward_cfg.get("discovery_weight", 1.0), 1.0)
        self.swap_target = float(np.clip(_safe_float(reward_cfg.get("swap_target", 0.25), 0.25), 0.0, 1.0))
        self.discovery_metric = str(reward_cfg.get("discovery_metric", "best_loglike"))

    def compute(self, summary: Dict[str, Any]) -> tuple[float, Dict[str, float]]:
        cold_ess = _safe_float(summary.get("cold_ess_proxy"), 0.0)
        round_trip = max(1.0, _safe_float(summary.get("round_trip_proxy"), float(self.nchains)))
        pair_rates = [float(rate) for rate in summary.get("pair_acceptance_rates", [])]
        best_improvement = _safe_float(summary.get("best_improvement"), 0.0)

        ess_term = math.tanh(cold_ess / 32.0)
        round_trip_term = -math.tanh(round_trip / float(max(1, self.nchains * 4)))
        if pair_rates:
            swap_term = -float(sum(abs(rate - self.swap_target) for rate in pair_rates) / len(pair_rates))
        else:
            swap_term = 0.0
        discovery_term = math.tanh(best_improvement)

        total = (
            (self.ess_weight * ess_term)
            + (self.round_trip_weight * round_trip_term)
            + (self.swap_weight * swap_term)
            + (self.discovery_weight * discovery_term)
        )
        components = {
            "ess_term": float(ess_term),
            "round_trip_term": float(round_trip_term),
            "swap_term": float(swap_term),
            "discovery_term": float(discovery_term),
            "total_reward": float(total),
        }
        return float(total), components


class RLTPPPOTrainerBase:
    def __init__(self, state_dim: int, action_dim: int, ppo_cfg: Dict[str, Any], logger=None) -> None:
        self.backend = "base"
        self.state_dim = int(state_dim)
        self.action_dim = int(action_dim)
        self.logger = logger
        hidden_sizes = ppo_cfg.get("hidden_sizes", [64, 64])
        if not isinstance(hidden_sizes, list) or not hidden_sizes:
            hidden_sizes = [0]
        self.hidden_sizes = [int(max(0, value)) for value in hidden_sizes]
        self.hidden_dim = int(max(0, self.hidden_sizes[0])) if self.hidden_sizes else 0
        self.rollout_size = int(max(1, ppo_cfg.get("rollout_size", 8)))
        self.epochs = int(max(1, ppo_cfg.get("epochs", 4)))
        self.gamma = float(np.clip(_safe_float(ppo_cfg.get("gamma", 0.99), 0.99), 0.0, 0.9999))
        self.gae_lambda = float(np.clip(_safe_float(ppo_cfg.get("gae_lambda", 0.95), 0.95), 0.0, 1.0))
        self.clip_ratio = float(np.clip(_safe_float(ppo_cfg.get("clip_ratio", 0.2), 0.2), 1.0e-3, 1.0))
        self.learning_rate = float(max(1.0e-6, _safe_float(ppo_cfg.get("learning_rate", 3.0e-4), 3.0e-4)))
        self.entropy_coef = float(max(0.0, _safe_float(ppo_cfg.get("entropy_coef", 1.0e-3), 1.0e-3)))
        self.value_coef = float(max(0.0, _safe_float(ppo_cfg.get("value_coef", 0.5), 0.5)))
        self.rng_seed = int(ppo_cfg.get("seed", 12345))
        self.buffer: List[Dict[str, Any]] = []
        self.update_count = 0

    def select_action(self, state: np.ndarray, *, deterministic: bool = False) -> tuple[np.ndarray, float, float]:
        raise NotImplementedError

    def update(self) -> Dict[str, Any]:
        raise NotImplementedError

    def export_state(self) -> Dict[str, Any]:
        raise NotImplementedError

    def load_state(self, payload: Dict[str, Any]) -> None:
        raise NotImplementedError

    def add_transition(
        self,
        *,
        state: np.ndarray,
        action: Sequence[float],
        log_prob: float,
        value: float,
        reward: float,
        done: bool = False,
        train_enabled: bool = True,
    ) -> Dict[str, Any] | None:
        self.buffer.append(
            {
                "state": np.asarray(state, dtype=np.float32),
                "action": np.asarray(action, dtype=np.float32),
                "log_prob": float(log_prob),
                "value": float(value),
                "reward": float(reward),
                "done": bool(done),
            }
        )
        if not train_enabled or len(self.buffer) < self.rollout_size:
            return None
        return self.update()

    def finalize(self, *, train_enabled: bool) -> Dict[str, Any] | None:
        if not train_enabled or not self.buffer:
            self.buffer.clear()
            return None
        return self.update()

    def _compute_returns_and_advantages(self) -> tuple[np.ndarray, np.ndarray]:
        rewards = [float(item["reward"]) for item in self.buffer]
        values = [float(item["value"]) for item in self.buffer] + [0.0]
        dones = [bool(item["done"]) for item in self.buffer]

        advantages = [0.0 for _ in rewards]
        gae = 0.0
        for idx in reversed(range(len(rewards))):
            mask = 0.0 if dones[idx] else 1.0
            delta = rewards[idx] + (self.gamma * values[idx + 1] * mask) - values[idx]
            gae = delta + (self.gamma * self.gae_lambda * mask * gae)
            advantages[idx] = gae
        returns = np.asarray([advantages[idx] + values[idx] for idx in range(len(rewards))], dtype=np.float32)
        advantages_arr = np.asarray(advantages, dtype=np.float32)
        if advantages_arr.size > 1:
            advantages_arr = (advantages_arr - float(np.mean(advantages_arr))) / (float(np.std(advantages_arr)) + 1.0e-8)
        return returns, advantages_arr


class RLTPNumpyPPOTrainer(RLTPPPOTrainerBase):
    def __init__(self, state_dim: int, action_dim: int, ppo_cfg: Dict[str, Any], logger=None) -> None:
        super().__init__(state_dim=state_dim, action_dim=action_dim, ppo_cfg=ppo_cfg, logger=logger)
        self.backend = "numpy"
        self.rng = np.random.default_rng(seed=self.rng_seed)
        feature_dim = self.hidden_dim if self.hidden_dim > 0 else self.state_dim
        if self.hidden_dim > 0:
            self.feature_weight = self.rng.normal(0.0, 0.05, size=(self.hidden_dim, self.state_dim)).astype(np.float32)
            self.feature_bias = np.zeros(self.hidden_dim, dtype=np.float32)
        else:
            self.feature_weight = None
            self.feature_bias = None
        self.policy_weight = self.rng.normal(0.0, 0.05, size=(self.action_dim, feature_dim)).astype(np.float32)
        self.policy_bias = np.zeros(self.action_dim, dtype=np.float32)
        self.value_weight = self.rng.normal(0.0, 0.05, size=(feature_dim,)).astype(np.float32)
        self.value_bias = np.float32(0.0)
        self.log_std = np.zeros(self.action_dim, dtype=np.float32)

    def _forward(self, states: np.ndarray):
        states = np.asarray(states, dtype=np.float32)
        if self.hidden_dim > 0:
            pre_features = states @ self.feature_weight.T + self.feature_bias
            features = np.tanh(pre_features)
        else:
            features = states
        means = features @ self.policy_weight.T + self.policy_bias
        values = features @ self.value_weight + float(self.value_bias)
        return features, means, values

    def _log_prob(self, means: np.ndarray, actions: np.ndarray) -> np.ndarray:
        std = np.exp(self.log_std)
        var = std ** 2
        centered = actions - means
        return -0.5 * np.sum((centered ** 2) / var + (2.0 * self.log_std) + np.log(2.0 * np.pi), axis=1)

    def select_action(self, state: np.ndarray, *, deterministic: bool = False) -> tuple[np.ndarray, float, float]:
        state_batch = np.asarray(state, dtype=np.float32).reshape(1, -1)
        _, means, values = self._forward(state_batch)
        mean = means[0]
        if deterministic:
            action = mean
        else:
            action = mean + np.exp(self.log_std) * self.rng.standard_normal(self.action_dim).astype(np.float32)
        log_prob = float(self._log_prob(means, action.reshape(1, -1))[0])
        return action.astype(np.float32), log_prob, float(values[0])

    def update(self) -> Dict[str, Any]:
        if not self.buffer:
            return {}

        states = np.stack([item["state"] for item in self.buffer], axis=0).astype(np.float32)
        actions = np.stack([item["action"] for item in self.buffer], axis=0).astype(np.float32)
        old_log_probs = np.asarray([item["log_prob"] for item in self.buffer], dtype=np.float32)
        returns, advantages = self._compute_returns_and_advantages()
        sample_count = float(max(1, len(self.buffer)))

        metrics = {
            "policy_loss": 0.0,
            "value_loss": 0.0,
            "entropy": 0.0,
            "updates": self.update_count,
            "backend": self.backend,
        }
        for _ in range(self.epochs):
            features, means, values = self._forward(states)
            new_log_probs = self._log_prob(means, actions)
            ratio = np.exp(new_log_probs - old_log_probs)
            unclipped = ratio * advantages
            clipped_ratio = np.clip(ratio, 1.0 - self.clip_ratio, 1.0 + self.clip_ratio)
            clipped = clipped_ratio * advantages
            surrogate = np.minimum(unclipped, clipped)
            policy_loss = -float(np.mean(surrogate))
            value_residual = values - returns
            value_loss = float(np.mean(value_residual ** 2))
            entropy = float(np.sum(self.log_std + 0.5 * (1.0 + np.log(2.0 * np.pi))))

            std = np.exp(self.log_std)
            var = std ** 2
            grad_log_prob_mean = (actions - means) / var
            grad_log_prob_log_std = -1.0 + ((actions - means) ** 2) / var
            use_clipped = ((advantages >= 0.0) & (ratio > 1.0 + self.clip_ratio)) | (
                (advantages < 0.0) & (ratio < 1.0 - self.clip_ratio)
            )
            coeff = np.where(use_clipped, 0.0, -advantages * ratio) / sample_count

            grad_policy_output = coeff[:, None] * grad_log_prob_mean
            grad_policy_weight = grad_policy_output.T @ features
            grad_policy_bias = np.sum(grad_policy_output, axis=0)
            grad_log_std = np.sum(coeff[:, None] * grad_log_prob_log_std, axis=0) - self.entropy_coef

            grad_value = (2.0 * self.value_coef * value_residual) / sample_count
            grad_value_weight = features.T @ grad_value
            grad_value_bias = float(np.sum(grad_value))

            if self.hidden_dim > 0:
                feature_grad = grad_policy_output @ self.policy_weight
                feature_grad += grad_value[:, None] * self.value_weight[None, :]
                pre_grad = feature_grad * (1.0 - (features ** 2))
                grad_feature_weight = pre_grad.T @ states
                grad_feature_bias = np.sum(pre_grad, axis=0)
                self.feature_weight -= self.learning_rate * grad_feature_weight.astype(np.float32)
                self.feature_bias -= self.learning_rate * grad_feature_bias.astype(np.float32)

            self.policy_weight -= self.learning_rate * grad_policy_weight.astype(np.float32)
            self.policy_bias -= self.learning_rate * grad_policy_bias.astype(np.float32)
            self.value_weight -= self.learning_rate * grad_value_weight.astype(np.float32)
            self.value_bias -= np.float32(self.learning_rate * grad_value_bias)
            self.log_std -= self.learning_rate * grad_log_std.astype(np.float32)
            self.log_std = np.clip(self.log_std, -3.0, 1.0)

            metrics["policy_loss"] = float(policy_loss)
            metrics["value_loss"] = float(value_loss)
            metrics["entropy"] = float(entropy)

        self.buffer.clear()
        self.update_count += 1
        metrics["updates"] = int(self.update_count)
        if self.logger is not None:
            self.logger.info(
                "PPO update -> "
                f"backend={self.backend} "
                f"policy_loss={metrics['policy_loss']:.4f} "
                f"value_loss={metrics['value_loss']:.4f} "
                f"entropy={metrics['entropy']:.4f}"
            )
        return metrics

    def export_state(self) -> Dict[str, Any]:
        return {
            "backend": self.backend,
            "format_version": 1,
            "policy_weight": self.policy_weight,
            "policy_bias": self.policy_bias,
            "value_weight": self.value_weight,
            "value_bias": float(self.value_bias),
            "log_std": self.log_std,
            "feature_weight": self.feature_weight,
            "feature_bias": self.feature_bias,
            "trainer_state": {
                "state_dim": self.state_dim,
                "action_dim": self.action_dim,
                "hidden_dim": self.hidden_dim,
                "update_count": self.update_count,
                "rollout_size": self.rollout_size,
                "epochs": self.epochs,
                "gamma": self.gamma,
                "gae_lambda": self.gae_lambda,
                "clip_ratio": self.clip_ratio,
                "seed": self.rng_seed,
            },
        }

    def load_state(self, payload: Dict[str, Any]) -> None:
        if not payload:
            return
        payload_backend = str(payload.get("backend", "")).lower()
        if payload_backend and payload_backend != self.backend:
            raise ValueError(f"RLTPMCMC checkpoint backend mismatch: expect {self.backend}, got {payload_backend}")
        self.policy_weight = np.asarray(payload["policy_weight"], dtype=np.float32)
        self.policy_bias = np.asarray(payload["policy_bias"], dtype=np.float32)
        self.value_weight = np.asarray(payload["value_weight"], dtype=np.float32)
        self.value_bias = np.float32(payload["value_bias"])
        self.log_std = np.asarray(payload["log_std"], dtype=np.float32)
        self.feature_weight = (
            None if payload.get("feature_weight") is None else np.asarray(payload["feature_weight"], dtype=np.float32)
        )
        self.feature_bias = (
            None if payload.get("feature_bias") is None else np.asarray(payload["feature_bias"], dtype=np.float32)
        )
        trainer_state = payload.get("trainer_state", {})
        self.update_count = int(trainer_state.get("update_count", self.update_count))


class RLTPTorchPPOTrainer(RLTPPPOTrainerBase):
    def __init__(self, state_dim: int, action_dim: int, ppo_cfg: Dict[str, Any], logger=None) -> None:
        super().__init__(state_dim=state_dim, action_dim=action_dim, ppo_cfg=ppo_cfg, logger=logger)
        self.backend = "torch"
        self.torch = _load_torch_module()
        self.torch.manual_seed(self.rng_seed)
        self.device = self.torch.device(str(ppo_cfg.get("device", "cpu")))
        self._build_model()

    def _build_model(self) -> None:
        torch = self.torch
        nn = torch.nn
        state_dim = self.state_dim
        action_dim = self.action_dim
        hidden_dim = self.hidden_dim

        class _ActorCritic(nn.Module):
            def __init__(self, in_dim: int, out_dim: int, hidden: int) -> None:
                super().__init__()
                if hidden > 0:
                    self.feature = nn.Sequential(nn.Linear(in_dim, hidden), nn.Tanh())
                    feature_dim = hidden
                else:
                    self.feature = None
                    feature_dim = in_dim
                self.policy_head = nn.Linear(feature_dim, out_dim)
                self.value_head = nn.Linear(feature_dim, 1)
                self.log_std = nn.Parameter(torch.zeros(out_dim, dtype=torch.float32))

            def forward(self, states):
                features = states if self.feature is None else self.feature(states)
                means = self.policy_head(features)
                values = self.value_head(features).squeeze(-1)
                log_std = self.log_std.expand_as(means)
                return features, means, values, log_std

        self.model = _ActorCritic(state_dim, action_dim, hidden_dim).to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.learning_rate)

    def _log_prob(self, means, log_std, actions):
        centered = actions - means
        var = self.torch.exp(log_std * 2.0)
        return -0.5 * self.torch.sum((centered ** 2) / var + (2.0 * log_std) + math.log(2.0 * math.pi), dim=1)

    def select_action(self, state: np.ndarray, *, deterministic: bool = False) -> tuple[np.ndarray, float, float]:
        torch = self.torch
        state_tensor = torch.as_tensor(np.asarray(state, dtype=np.float32), dtype=torch.float32, device=self.device).unsqueeze(0)
        with torch.no_grad():
            _, means, values, log_std = self.model(state_tensor)
            if deterministic:
                action_tensor = means
            else:
                action_tensor = means + torch.exp(log_std) * torch.randn_like(means)
            log_prob = self._log_prob(means, log_std, action_tensor)
        action = action_tensor.squeeze(0).detach().cpu().numpy().astype(np.float32)
        return action, float(log_prob.item()), float(values.squeeze(0).item())

    def update(self) -> Dict[str, Any]:
        if not self.buffer:
            return {}

        torch = self.torch
        states = torch.as_tensor(np.stack([item["state"] for item in self.buffer], axis=0), dtype=torch.float32, device=self.device)
        actions = torch.as_tensor(np.stack([item["action"] for item in self.buffer], axis=0), dtype=torch.float32, device=self.device)
        old_log_probs = torch.as_tensor([item["log_prob"] for item in self.buffer], dtype=torch.float32, device=self.device)
        returns_np, advantages_np = self._compute_returns_and_advantages()
        returns = torch.as_tensor(returns_np, dtype=torch.float32, device=self.device)
        advantages = torch.as_tensor(advantages_np, dtype=torch.float32, device=self.device)

        metrics = {
            "policy_loss": 0.0,
            "value_loss": 0.0,
            "entropy": 0.0,
            "updates": self.update_count,
            "backend": self.backend,
        }

        for _ in range(self.epochs):
            _, means, values, log_std = self.model(states)
            new_log_probs = self._log_prob(means, log_std, actions)
            ratio = torch.exp(new_log_probs - old_log_probs)
            unclipped = ratio * advantages
            clipped = torch.clamp(ratio, 1.0 - self.clip_ratio, 1.0 + self.clip_ratio) * advantages
            policy_loss = -torch.mean(torch.minimum(unclipped, clipped))
            value_loss = torch.mean((values - returns) ** 2)
            entropy = torch.mean(torch.sum(log_std + 0.5 * (1.0 + math.log(2.0 * math.pi)), dim=1))
            loss = policy_loss + (self.value_coef * value_loss) - (self.entropy_coef * entropy)

            self.optimizer.zero_grad()
            loss.backward()
            self.optimizer.step()

            metrics["policy_loss"] = float(policy_loss.detach().cpu().item())
            metrics["value_loss"] = float(value_loss.detach().cpu().item())
            metrics["entropy"] = float(entropy.detach().cpu().item())

        self.buffer.clear()
        self.update_count += 1
        metrics["updates"] = int(self.update_count)
        if self.logger is not None:
            self.logger.info(
                "PPO update -> "
                f"backend={self.backend} "
                f"policy_loss={metrics['policy_loss']:.4f} "
                f"value_loss={metrics['value_loss']:.4f} "
                f"entropy={metrics['entropy']:.4f}"
            )
        return metrics

    def export_state(self) -> Dict[str, Any]:
        return {
            "backend": self.backend,
            "format_version": 1,
            "model_state": {key: value.detach().cpu() for key, value in self.model.state_dict().items()},
            "optimizer_state": self.optimizer.state_dict(),
            "trainer_state": {
                "state_dim": self.state_dim,
                "action_dim": self.action_dim,
                "hidden_dim": self.hidden_dim,
                "update_count": self.update_count,
                "rollout_size": self.rollout_size,
                "epochs": self.epochs,
                "gamma": self.gamma,
                "gae_lambda": self.gae_lambda,
                "clip_ratio": self.clip_ratio,
                "seed": self.rng_seed,
                "device": str(self.device),
                "torch_version": getattr(self.torch, "__version__", None),
            },
        }

    def load_state(self, payload: Dict[str, Any]) -> None:
        if not payload:
            return
        payload_backend = str(payload.get("backend", "")).lower()
        if payload_backend and payload_backend != self.backend:
            raise ValueError(f"RLTPMCMC checkpoint backend mismatch: expect {self.backend}, got {payload_backend}")
        self.model.load_state_dict(payload["model_state"])
        optimizer_state = payload.get("optimizer_state")
        if optimizer_state is not None:
            self.optimizer.load_state_dict(optimizer_state)
        trainer_state = payload.get("trainer_state", {})
        self.update_count = int(trainer_state.get("update_count", self.update_count))


class RLTPController:
    def __init__(
        self,
        sampler,
        diagnostics: RLTPDiagnostics,
        state_builder: RLTPStateBuilder,
        action_mapper: RLTPActionMapper,
        reward_builder: RLTPRewardBuilder,
        trainer: RLTPPPOTrainerBase | None,
    ) -> None:
        self.sampler = sampler
        self.diagnostics = diagnostics
        self.state_builder = state_builder
        self.action_mapper = action_mapper
        self.reward_builder = reward_builder
        self.trainer = trainer
        self.exchange_windows = 0
        self.pending_transition: Dict[str, Any] | None = None

    def on_run_start(self, context: Dict[str, Any]) -> MCMCControlPatch | None:
        _ = context
        self.diagnostics.start_window()
        self.sampler.rl_log("warning", "RLTPMCMC controller initialized")
        return None

    def on_pre_step(self, snapshot: Dict[str, Any]) -> MCMCControlPatch | None:
        _ = snapshot
        return None

    def on_post_step(
        self, snapshot: Dict[str, Any], outcome: Dict[str, Any]
    ) -> MCMCControlPatch | None:
        _ = snapshot
        _ = outcome
        return None

    def on_pre_exchange(self, snapshot: Dict[str, Any]) -> MCMCControlPatch | None:
        _ = snapshot
        return None

    def _control_enabled(self) -> bool:
        return str(self.sampler._control_cfg["controller"]).lower() == "ppo" and self.trainer is not None

    def _training_enabled(self, snapshot: Dict[str, Any]) -> bool:
        controller = str(self.sampler._control_cfg["controller"]).lower()
        if controller != "ppo" or self.trainer is None:
            return False

        adaptation_mode = str(self.sampler._control_cfg["adaptation_mode"]).lower()
        mean_iter = 0
        chains = snapshot.get("chains", [])
        if chains:
            mean_iter = int(sum(int(chain.get("iter", 0)) for chain in chains) / len(chains))

        freeze_after_burnin = bool(self.sampler._control_cfg.get("freeze_after_burnin", False))
        burnin_iters = int(max(0, self.sampler._control_cfg.get("burnin_iters", 0)))

        if adaptation_mode == "frozen":
            return False
        if adaptation_mode == "scheduled" and mean_iter >= burnin_iters:
            return False
        if freeze_after_burnin and mean_iter >= burnin_iters:
            return False
        return True

    def _deterministic_eval(self, training_enabled: bool) -> bool:
        if not training_enabled:
            return True
        return bool(self.sampler._control_cfg.get("deterministic_eval", False))

    def on_post_exchange(
        self, snapshot: Dict[str, Any], exchange_metrics: Dict[str, Any]
    ) -> MCMCControlPatch | None:
        _ = exchange_metrics
        self.exchange_windows += 1
        if self.exchange_windows % int(self.sampler._control_cfg["decision_interval"]) != 0:
            return None

        training_enabled = self._training_enabled(snapshot)
        summary = self.diagnostics.collect_window_summary(
            self.sampler,
            training_enabled=training_enabled,
        )
        reward, reward_components = self.reward_builder.compute(summary)

        training_metrics = None
        if self.pending_transition is not None and self.trainer is not None:
            training_metrics = self.trainer.add_transition(
                state=self.pending_transition["state"],
                action=self.pending_transition["action"],
                log_prob=self.pending_transition["log_prob"],
                value=self.pending_transition["value"],
                reward=reward,
                done=False,
                train_enabled=training_enabled,
            )
            if training_metrics:
                self.sampler.write_rl_record("training_metrics", training_metrics)

        state, state_meta = self.state_builder.build(summary)
        decision_record: Dict[str, Any] = {
            "decision_index": int(self.exchange_windows // int(self.sampler._control_cfg["decision_interval"])),
            "controller": str(self.sampler._control_cfg["controller"]),
            "adaptation_mode": str(self.sampler._control_cfg["adaptation_mode"]),
            "ppo_backend": self.sampler._current_ppo_backend(),
            "training_enabled": bool(training_enabled),
            "reward": float(reward),
            "reward_components": reward_components,
            "summary": summary,
            "state": state_meta,
        }

        patch = None
        if self._control_enabled():
            deterministic = self._deterministic_eval(training_enabled)
            action, log_prob, value = self.trainer.select_action(state, deterministic=deterministic)
            patch, action_details = self.action_mapper.to_patch(summary, action)
            decision_record["action"] = action.tolist()
            decision_record["action_details"] = action_details
            if training_enabled:
                self.pending_transition = {
                    "state": state,
                    "action": np.asarray(action, dtype=np.float32),
                    "log_prob": float(log_prob),
                    "value": float(value),
                }
            else:
                self.pending_transition = None
        else:
            self.pending_transition = None

        self.sampler.write_rl_record("decision_metrics", decision_record)
        self.sampler.write_rl_record(
            "state_action_reward",
            {
                "decision_index": decision_record["decision_index"],
                "ppo_backend": decision_record["ppo_backend"],
                "reward": float(reward),
                "reward_components": reward_components,
                "state": state_meta,
                "action": decision_record.get("action"),
            },
        )

        self.sampler.rl_log(
            "warning",
            format_two_column_log(
                "RLTPMCMC decision window",
                [
                    ("decision", decision_record["decision_index"]),
                    ("controller", self.sampler._control_cfg["controller"]),
                    ("reward", f"{reward:.4f}"),
                    ("swap_mean", f"{summary['swap_acceptance_mean']:.4f}"),
                    ("cold_ess", f"{summary['cold_ess_proxy']:.4f}"),
                ],
            ),
        )

        self.diagnostics.start_window()
        return patch

    def flush(self) -> Dict[str, Any] | None:
        if self.pending_transition is None or self.trainer is None:
            return None
        summary = self.diagnostics.collect_window_summary(
            self.sampler,
            training_enabled=True,
        )
        reward, _ = self.reward_builder.compute(summary)
        self.trainer.add_transition(
            state=self.pending_transition["state"],
            action=self.pending_transition["action"],
            log_prob=self.pending_transition["log_prob"],
            value=self.pending_transition["value"],
            reward=reward,
            done=True,
            train_enabled=True,
        )
        metrics = self.trainer.finalize(train_enabled=True)
        self.pending_transition = None
        return metrics


class RLTPMCMC(PTMCMC, RLSamplerBase):
    def __init__(self) -> None:
        super().__init__()
        self.method = "RLTPMCMC"
        self.load_schema_file()

        self._control_cfg: Dict[str, Any] = {
            "controller": "none",
            "adaptation_mode": "frozen",
            "decision_interval": 1,
            "smoothing": 0.50,
            "max_action_magnitude": 0.20,
            "freeze_after_burnin": False,
            "burnin_iters": 0,
            "deterministic_eval": False,
            "min_beta_floor": 1.0e-3,
            "min_proposal_scale": 1.0e-3,
            "max_proposal_scale": 5.0,
        }
        self._reward_cfg: Dict[str, Any] = {
            "ess_weight": 1.0,
            "round_trip_weight": 0.5,
            "swap_weight": 0.5,
            "discovery_weight": 1.0,
            "swap_target": 0.25,
            "discovery_metric": "best_loglike",
        }
        self._ppo_cfg: Dict[str, Any] = {
            "backend": "torch",
            "rollout_size": 8,
            "epochs": 4,
            "gamma": 0.99,
            "gae_lambda": 0.95,
            "clip_ratio": 0.20,
            "learning_rate": 3.0e-4,
            "entropy_coef": 1.0e-3,
            "value_coef": 0.5,
            "hidden_sizes": [64, 64],
            "seed": 12345,
            "device": "cpu",
            "resume_from": None,
        }
        self._diagnostics_cfg: Dict[str, Any] = {
            "state_window": 32,
            "reward_window": 32,
            "log_interval": 1,
        }

        self._diagnostics: RLTPDiagnostics | None = None
        self._state_builder: RLTPStateBuilder | None = None
        self._action_mapper: RLTPActionMapper | None = None
        self._reward_builder: RLTPRewardBuilder | None = None
        self._trainer: RLTPPPOTrainerBase | None = None
        self._rl_controller: RLTPController | None = None
        self._replica_positions: List[int] = []
        self._replica_hot_seen: Dict[int, int | None] = {}
        self._round_trip_intervals = deque(maxlen=128)
        self._exchange_counter = 0

    def supports_runtime_checkpointing(self) -> bool:
        return True

    def load_schema_file(self) -> None:
        self.schema = self.path["RLTPMCMCSchema"]

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("RLTPMCMC sampler logger initialized")

    def init_generator(self) -> None:
        super().init_generator()
        sampling_cfg = self.config.get("Sampling", {}) if isinstance(self.config, dict) else {}
        self._control_cfg.update(sampling_cfg.get("Control", {}) or {})
        self._reward_cfg.update(sampling_cfg.get("Reward", {}) or {})
        self._ppo_cfg.update(sampling_cfg.get("PPO", {}) or {})
        self._diagnostics_cfg.update(sampling_cfg.get("Diagnostics", {}) or {})
        self._ppo_cfg["backend"] = self._resolved_ppo_backend()
        self._selectionexp = sampling_cfg.get("selection")

    def _resolved_ppo_backend(self) -> str:
        backend = str(self._ppo_cfg.get("backend", "torch")).strip().lower() or "torch"
        if backend not in {"torch", "numpy"}:
            raise ValueError(f"RLTPMCMC unsupported PPO backend: {backend}")
        return backend

    def _export_runtime_extras(self) -> Dict[str, Any]:
        trainer_state = None
        if self._trainer is not None and hasattr(self._trainer, "export_state"):
            trainer_state = self._trainer.export_state()
        return {
            "control_cfg": dict(self._control_cfg),
            "reward_cfg": dict(self._reward_cfg),
            "ppo_cfg": dict(self._ppo_cfg),
            "diagnostics_cfg": dict(self._diagnostics_cfg),
            "replica_positions": list(self._replica_positions),
            "replica_hot_seen": {int(k): (None if v is None else int(v)) for k, v in self._replica_hot_seen.items()},
            "round_trip_intervals": list(self._round_trip_intervals),
            "exchange_counter": int(self._exchange_counter),
            "trainer_state": trainer_state,
            "ppo_backend": self._current_ppo_backend(),
        }

    def _import_runtime_extras(self, payload: Dict[str, Any]) -> None:
        if not isinstance(payload, dict):
            return
        self._control_cfg.update(payload.get("control_cfg", {}) or {})
        self._reward_cfg.update(payload.get("reward_cfg", {}) or {})
        self._ppo_cfg.update(payload.get("ppo_cfg", {}) or {})
        self._diagnostics_cfg.update(payload.get("diagnostics_cfg", {}) or {})
        self._replica_positions = list(payload.get("replica_positions", self._replica_positions))
        self._replica_hot_seen = {int(k): v for k, v in dict(payload.get("replica_hot_seen", {})).items()}
        self._round_trip_intervals = deque(payload.get("round_trip_intervals", []), maxlen=128)
        self._exchange_counter = int(payload.get("exchange_counter", self._exchange_counter))
        trainer_state = payload.get("trainer_state")
        if trainer_state is not None and self._trainer is not None and hasattr(self._trainer, "load_state"):
            try:
                self._trainer.load_state(trainer_state)
            except Exception as exc:
                self.rl_log("warning", f"Failed to restore PPO trainer state: {exc}")

    def _current_ppo_backend(self) -> str | None:
        if self._trainer is not None:
            return str(self._trainer.backend)
        controller = str(self._control_cfg.get("controller", "none")).lower()
        if controller != "ppo":
            return None
        return self._resolved_ppo_backend()

    def _build_trainer(self, state_dim: int, action_dim: int) -> RLTPPPOTrainerBase:
        backend = self._resolved_ppo_backend()
        trainer_cls = RLTPTorchPPOTrainer if backend == "torch" else RLTPNumpyPPOTrainer
        trainer = trainer_cls(
            state_dim=state_dim,
            action_dim=action_dim,
            ppo_cfg=dict(self._ppo_cfg),
            logger=self.rl_logger,
        )
        self._ppo_cfg["backend"] = backend
        return trainer

    def _create_chain_registry(self):
        registry = super()._create_chain_registry()
        self._replica_positions = list(range(len(registry)))
        self._replica_hot_seen = {int(replica_id): None for replica_id in self._replica_positions}
        self._round_trip_intervals = deque(maxlen=max(64, int(self._diagnostics_cfg.get("reward_window", 32))))
        self._exchange_counter = 0
        return registry

    def _pair_chain_ids(self):
        ids = self._must_registry().ids()
        start = int(self._exchange_offset % 2)
        pairs = [(ids[ii], ids[ii + 1]) for ii in range(start, len(ids) - 1, 2)]
        self._exchange_offset = 1 - start
        return pairs

    def _attempt_swap(self, cid1: int, cid2: int) -> bool:
        accepted = super()._attempt_swap(cid1, cid2)
        if accepted:
            self._replica_positions[cid1], self._replica_positions[cid2] = (
                self._replica_positions[cid2],
                self._replica_positions[cid1],
            )
        return accepted

    def _update_round_trip_proxy(self) -> float:
        self._exchange_counter += 1
        hot_chain_id = int(self._nchains - 1)
        for chain_id, replica_id in enumerate(self._replica_positions):
            if chain_id == hot_chain_id:
                self._replica_hot_seen[int(replica_id)] = int(self._exchange_counter)
            elif chain_id == 0:
                hot_seen = self._replica_hot_seen.get(int(replica_id))
                if hot_seen is not None:
                    interval = int(self._exchange_counter - hot_seen)
                    if interval > 0:
                        self._round_trip_intervals.append(interval)
                    self._replica_hot_seen[int(replica_id)] = None
        if self._round_trip_intervals:
            return float(sum(self._round_trip_intervals) / len(self._round_trip_intervals))
        return float(max(1, self._nchains * int(self._control_cfg.get("decision_interval", 1))))

    def _do_exchange(self):
        attempted = 0
        accepted = 0
        pair_metrics = []
        for cid1, cid2 in self._pair_chain_ids():
            attempted += 1
            is_accepted = self._attempt_swap(cid1, cid2)
            if is_accepted:
                accepted += 1
            pair_metrics.append(
                {
                    "pair": [int(cid1), int(cid2)],
                    "pair_index": int(min(cid1, cid2)),
                    "attempted": 1,
                    "accepted": 1 if is_accepted else 0,
                }
            )

        self._ready_queue.clear()
        self._ready_set.clear()
        for chain in self._must_registry().all():
            chain.window_iter = 0
            if chain.iter < int(self._niters):
                self._enqueue_chain(chain.chain_id)

        round_trip_proxy = self._update_round_trip_proxy()
        if self._diagnostics is not None:
            self._diagnostics.record_exchange(
                {
                    "attempted": attempted,
                    "accepted": accepted,
                    "pair_metrics": pair_metrics,
                    "round_trip_proxy": round_trip_proxy,
                    "temperatures": [chain.temperature for chain in self._must_registry().all()],
                    "proposal_scales": [getattr(chain.engine, "proposal_scale", None) for chain in self._must_registry().all()],
                    "replica_positions": list(self._replica_positions),
                }
            )

        self.logger.warning(
            format_two_column_log(
                "RLTPMCMC exchange summary",
                [
                    ("attempted", attempted),
                    ("accepted", accepted),
                    ("round_trip", f"{round_trip_proxy:.3f}"),
                ],
            )
        )
        return {
            "attempted": attempted,
            "accepted": accepted,
            "pair_metrics": pair_metrics,
            "round_trip_proxy": round_trip_proxy,
        }

    def _on_after_update(self, chain, sample_info: Dict[str, Any], accepted: bool) -> None:
        _ = chain
        _ = accepted
        if self._diagnostics is None:
            return
        try:
            logl = self._extract_logl(sample_info)
        except Exception:
            return
        self._diagnostics.update_best(logl)

    def initialize(self):
        self.init_rl_runtime("rltpmcmc")
        self.info.setdefault("rl", {})
        self.info["rl"]["ppo_backend"] = self._current_ppo_backend()

        history_size = max(
            int(self._diagnostics_cfg.get("state_window", 32)),
            int(self._diagnostics_cfg.get("reward_window", 32)),
        )
        self._diagnostics = RLTPDiagnostics(self._nchains, history_size=history_size)
        self._diagnostics.start_window()
        self._state_builder = RLTPStateBuilder(self._nchains, self._proposal_scales)
        self._action_mapper = RLTPActionMapper(
            self._nchains,
            smoothing=self._control_cfg["smoothing"],
            max_action_magnitude=self._control_cfg["max_action_magnitude"],
            min_beta_floor=self._control_cfg["min_beta_floor"],
            min_proposal_scale=self._control_cfg["min_proposal_scale"],
            max_proposal_scale=self._control_cfg["max_proposal_scale"],
        )
        self._reward_builder = RLTPRewardBuilder(self._reward_cfg, self._nchains)

        if str(self._control_cfg["controller"]).lower() == "ppo":
            state_dim = ((self._nchains - 1) * 2) + self._nchains + 7
            self._trainer = self._build_trainer(
                state_dim=state_dim,
                action_dim=self._action_mapper.action_dim,
            )
            resume_from = self._ppo_cfg.get("resume_from")
            if resume_from and not getattr(self, "_runtime_checkpoint_resume_hint", False):
                payload = self.rl_state_saver.load_checkpoint(str(resume_from), default=None)
                if payload:
                    self._trainer.load_state(payload)
        else:
            self._trainer = None

        self._rl_controller = RLTPController(
            sampler=self,
            diagnostics=self._diagnostics,
            state_builder=self._state_builder,
            action_mapper=self._action_mapper,
            reward_builder=self._reward_builder,
            trainer=self._trainer,
        )
        self.set_controller(self._rl_controller)
        super().initialize()

    def finalize(self):
        training_metrics = None
        if self._rl_controller is not None:
            training_metrics = self._rl_controller.flush()
        if training_metrics:
            self.write_rl_record("training_metrics", training_metrics)

        checkpoint_path = None
        control_state = {
            "method": self.method,
            "control": self._control_cfg,
            "reward": self._reward_cfg,
            "ppo": self._ppo_cfg,
            "diagnostics": self._diagnostics_cfg,
            "ppo_backend": self._current_ppo_backend(),
            "checkpoint": checkpoint_path,
        }

        super().finalize()
        self.close_rl_runtime(
            {
                "method": self.method,
                "controller": self._control_cfg["controller"],
                "ppo_backend": self._current_ppo_backend(),
                "checkpoint_path": checkpoint_path,
                "control_state": control_state,
            }
        )

    def set_factory(self, factory) -> None:
        super().set_factory(factory)
        self.rl_log(
            "warning",
            format_two_column_log(
                "RLTPMCMC runtime mode",
                [
                    ("controller", self._control_cfg["controller"]),
                    ("ppo_backend", self._current_ppo_backend() or "n/a"),
                    ("adaptation", self._control_cfg["adaptation_mode"]),
                    ("decision_interval", self._control_cfg["decision_interval"]),
                ],
            ),
        )
