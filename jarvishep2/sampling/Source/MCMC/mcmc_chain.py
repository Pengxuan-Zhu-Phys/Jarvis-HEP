#!/usr/bin/env python3
"""Standard Metropolis chain engine (V1 science, seeded RNG for V2)."""

from __future__ import annotations

from typing import Any

import numpy as np


class MCMCChain:
    """Random-walk Metropolis on the unit cube with Metropolis–Hastings accept."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
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

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposed_param = self._rng.random(self.param.shape[0])
            self.proposed_param = proposed_param
            return proposed_param
        while True:
            step = self._rng.normal(0.0, self.proposal_scale, size=self.param.shape)
            proposed_param = self.param + step
            if np.all((proposed_param >= 0.0) & (proposed_param <= 1.0)):
                self.proposed_param = proposed_param
                return proposed_param

    def update(self, new_loglikelihood, beta=1.0) -> bool:
        accepted = False
        new_loglikelihood = float(new_loglikelihood)
        beta = float(beta)

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

        self.iterations += 1
        return accepted

    def export_state(self) -> dict[str, Any]:
        return {
            "type": "MCMCChain",
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
            "rng_state": self._rng.bit_generator.state,
        }

    def import_state(self, state: dict[str, Any]) -> None:
        self.param = np.asarray(state.get("param", self.param), dtype=float)
        self.proposal_scale = float(state.get("proposal_scale", self.proposal_scale))
        self.n_iterations = int(state.get("n_iterations", self.n_iterations))
        self.iterations = int(state.get("iterations", self.iterations))
        self.last_loglikelihood = state.get("last_loglikelihood")
        prop = state.get("proposed_param")
        self.proposed_param = None if prop is None else np.asarray(prop, dtype=float)
        rng_state = state.get("rng_state")
        if rng_state is not None:
            self._rng.bit_generator.state = rng_state


__all__ = ["MCMCChain"]
