#!/usr/bin/env python3
"""Adaptive Metropolis chain engine (V1 science, seeded RNG for V2)."""

from __future__ import annotations

from typing import Any

import numpy as np


class AMMCMCChain:
    """Adaptive Metropolis with bounded unit-cube proposals."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        adapt_enabled=True,
        adapt_start_iter=100,
        adapt_window=25,
        adapt_eps=1e-6,
        adapt_scale=2.38,
        *,
        rng: np.random.Generator | None = None,
    ) -> None:
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood: float | None = None
        self.proposed_param = None
        self._rng = rng if rng is not None else np.random.default_rng()

        self._dim = int(self.param.shape[0])
        self._adapt_enabled = bool(adapt_enabled)
        self._adapt_start_iter = max(1, int(adapt_start_iter))
        self._adapt_window = max(1, int(adapt_window))
        self._adapt_eps = float(adapt_eps)
        self._adapt_scale = float(adapt_scale)

        self._history: list[np.ndarray] = []
        self._cov = (self.proposal_scale**2) * np.eye(self._dim)
        self._accepted_since_adapt = 0

    def __iter__(self):
        return self

    def _draw_adaptive_step(self):
        if not self._adapt_enabled or self.iterations < self._adapt_start_iter:
            return self._rng.normal(0.0, self.proposal_scale, size=self._dim)

        scale = (self._adapt_scale**2) / max(1.0, float(self._dim))
        cov = self._cov * scale
        try:
            return self._rng.multivariate_normal(np.zeros(self._dim), cov)
        except Exception:
            return self._rng.normal(0.0, self.proposal_scale, size=self._dim)

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration

        if self.iterations == 0:
            proposed_param = self._rng.random(self._dim)
            self.proposed_param = proposed_param
            return proposed_param

        for _ in range(2048):
            step = self._draw_adaptive_step()
            proposed_param = self.param + step
            if np.all((proposed_param >= 0.0) & (proposed_param <= 1.0)):
                self.proposed_param = proposed_param
                return proposed_param

        proposed_param = np.clip(
            self.param + self._rng.normal(0.0, self.proposal_scale, size=self._dim),
            0.0,
            1.0,
        )
        self.proposed_param = proposed_param
        return proposed_param

    def _maybe_adapt_cov(self) -> None:
        if not self._adapt_enabled:
            return
        if self.iterations < self._adapt_start_iter:
            return
        if self._accepted_since_adapt < self._adapt_window:
            return
        if len(self._history) < 2:
            return

        arr = np.asarray(self._history, dtype=float)
        if arr.ndim != 2 or arr.shape[0] < 2:
            return
        cov = np.cov(arr, rowvar=False)
        if np.ndim(cov) == 0:
            cov = np.array([[float(cov)]], dtype=float)
        cov = np.asarray(cov, dtype=float)
        cov += self._adapt_eps * np.eye(cov.shape[0], dtype=float)
        self._cov = cov
        self._accepted_since_adapt = 0

    def update(self, new_loglikelihood, beta=1.0) -> bool:
        new_loglikelihood = float(new_loglikelihood)
        beta = float(beta)
        accepted = False

        if self.iterations == 0 or self.last_loglikelihood is None:
            accepted = True
        else:
            delta = (new_loglikelihood - float(self.last_loglikelihood)) * beta
            if delta >= 0.0:
                accepted = True
            else:
                accepted = bool(self._rng.random() < np.exp(delta))

        if accepted:
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood
            self._history.append(np.array(self.param, dtype=float))
            self._accepted_since_adapt += 1
            self._maybe_adapt_cov()

        self.iterations += 1
        return accepted

    def export_state(self) -> dict[str, Any]:
        return {
            "type": "AMMCMCChain",
            "param": np.asarray(self.param, dtype=float).tolist(),
            "proposal_scale": float(self.proposal_scale),
            "n_iterations": int(self.n_iterations),
            "iterations": int(self.iterations),
            "last_loglikelihood": self.last_loglikelihood,
            "proposed_param": (
                None
                if self.proposed_param is None
                else np.asarray(self.proposed_param, dtype=float).tolist()
            ),
            "adapt_enabled": self._adapt_enabled,
            "adapt_start_iter": self._adapt_start_iter,
            "adapt_window": self._adapt_window,
            "adapt_eps": self._adapt_eps,
            "adapt_scale": self._adapt_scale,
            "history": [np.asarray(h, dtype=float).tolist() for h in self._history],
            "cov": np.asarray(self._cov, dtype=float).tolist(),
            "accepted_since_adapt": int(self._accepted_since_adapt),
            "rng_state": self._rng.bit_generator.state,
        }

    def import_state(self, state: dict[str, Any]) -> None:
        self.param = np.asarray(state.get("param", self.param), dtype=float)
        self._dim = int(self.param.shape[0])
        self.proposal_scale = float(state.get("proposal_scale", self.proposal_scale))
        self.n_iterations = int(state.get("n_iterations", self.n_iterations))
        self.iterations = int(state.get("iterations", self.iterations))
        self.last_loglikelihood = state.get("last_loglikelihood")
        prop = state.get("proposed_param")
        self.proposed_param = None if prop is None else np.asarray(prop, dtype=float)
        self._adapt_enabled = bool(state.get("adapt_enabled", self._adapt_enabled))
        self._adapt_start_iter = int(state.get("adapt_start_iter", self._adapt_start_iter))
        self._adapt_window = int(state.get("adapt_window", self._adapt_window))
        self._adapt_eps = float(state.get("adapt_eps", self._adapt_eps))
        self._adapt_scale = float(state.get("adapt_scale", self._adapt_scale))
        self._history = [
            np.asarray(h, dtype=float) for h in (state.get("history") or [])
        ]
        cov = state.get("cov")
        if cov is not None:
            self._cov = np.asarray(cov, dtype=float)
        self._accepted_since_adapt = int(
            state.get("accepted_since_adapt", self._accepted_since_adapt)
        )
        rng_state = state.get("rng_state")
        if rng_state is not None:
            self._rng.bit_generator.state = rng_state


__all__ = ["AMMCMCChain"]
