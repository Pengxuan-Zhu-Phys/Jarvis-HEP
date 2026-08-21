#!/usr/bin/env python3
"""Differential Evolution MCMC chain engine (V1 science, seeded RNG for V2)."""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

PopulationGetter = Callable[[int], list[np.ndarray]]


class DEMCMCChain:
    """DE-MCMC engine for one chain."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        chain_id,
        population_getter: PopulationGetter,
        de_gamma: float = 0.0,
        de_noise: float = 1.0e-3,
        de_crossover: float = 1.0,
        *,
        rng: np.random.Generator | None = None,
    ) -> None:
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood: float | None = None
        self.proposed_param = None

        self.chain_id = int(chain_id)
        self._population_getter = population_getter
        self._dim = int(self.param.shape[0])
        if float(de_gamma) > 0:
            self._gamma = float(de_gamma)
        else:
            self._gamma = float(2.38 / np.sqrt(max(1.0, 2.0 * self._dim)))
        self._noise = max(1.0e-9, float(de_noise))
        self._crossover = min(1.0, max(1.0e-6, float(de_crossover)))
        self._rng = rng if rng is not None else np.random.default_rng()

    def __iter__(self):
        return self

    def _other_population(self) -> list[np.ndarray]:
        try:
            population = self._population_getter(self.chain_id)
        except TypeError:
            population = self._population_getter()  # type: ignore[call-arg]
        if population is None:
            return []
        out: list[np.ndarray] = []
        for point in population:
            arr = np.asarray(point, dtype=float)
            if arr.shape == self.param.shape:
                out.append(arr)
        return out

    def _draw_demove(self) -> np.ndarray:
        others = self._other_population()
        if len(others) < 2:
            return self._rng.normal(0.0, self.proposal_scale, size=self._dim)

        ids = self._rng.choice(len(others), size=2, replace=False)
        xr1 = others[int(ids[0])]
        xr2 = others[int(ids[1])]
        diff = xr1 - xr2
        step = self._gamma * diff + self._noise * self._rng.normal(
            0.0, 1.0, size=self._dim
        )
        if self._crossover < 1.0:
            mask = self._rng.random(self._dim) < self._crossover
            if not np.any(mask):
                mask[int(self._rng.integers(0, self._dim))] = True
            step = step * mask
        return step

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = self._rng.random(self._dim)
            self.proposed_param = proposal
            return proposal

        proposal = None
        for _ in range(2048):
            step = self._draw_demove()
            cand = self.param + step
            if np.all((cand >= 0.0) & (cand <= 1.0)):
                proposal = cand
                break
        if proposal is None:
            proposal = np.clip(
                self.param + self._rng.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
        self.proposed_param = np.asarray(proposal, dtype=float)
        return self.proposed_param

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
                accepted = bool(
                    self._rng.random() < np.exp(np.clip(delta, -700.0, 0.0))
                )

        if accepted:
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood

        self.iterations += 1
        return accepted

    def export_state(self) -> dict[str, Any]:
        return {
            "type": "DEMCMCChain",
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
            "chain_id": int(self.chain_id),
            "gamma": float(self._gamma),
            "noise": float(self._noise),
            "crossover": float(self._crossover),
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
        self.chain_id = int(state.get("chain_id", self.chain_id))
        self._gamma = float(state.get("gamma", self._gamma))
        self._noise = float(state.get("noise", self._noise))
        self._crossover = float(state.get("crossover", self._crossover))
        rng_state = state.get("rng_state")
        if rng_state is not None:
            self._rng.bit_generator.state = rng_state


__all__ = ["DEMCMCChain"]
