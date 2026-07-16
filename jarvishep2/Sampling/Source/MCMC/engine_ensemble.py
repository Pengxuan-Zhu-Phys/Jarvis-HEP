#!/usr/bin/env python3
"""Goodman–Weare stretch-move ensemble engine (V1 science, seeded RNG for V2)."""

from __future__ import annotations

from typing import Any, Callable

import numpy as np

PopulationGetter = Callable[[int], list[np.ndarray]]


class EnsembleChain:
    """Stretch-move engine for one walker."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        chain_id,
        population_getter: PopulationGetter,
        stretch_a: float = 2.0,
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
        self._stretch_a = max(1.01, float(stretch_a))
        self._last_log_jacobian = 0.0
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

    def _draw_stretch_factor(self) -> float:
        u = float(self._rng.random())
        a = self._stretch_a
        return (((a - 1.0) * u + 1.0) ** 2) / a

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration

        if self.iterations == 0:
            proposal = self._rng.random(self._dim)
            self.proposed_param = proposal
            self._last_log_jacobian = 0.0
            return proposal

        others = self._other_population()
        if not others:
            proposal = np.clip(
                self.param + self._rng.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
            self.proposed_param = proposal
            self._last_log_jacobian = 0.0
            return proposal

        proposal = None
        log_jacobian = 0.0
        for _ in range(2048):
            partner = np.asarray(
                others[int(self._rng.integers(0, len(others)))], dtype=float
            )
            z = float(self._draw_stretch_factor())
            cand = partner + z * (self.param - partner)
            if np.all((cand >= 0.0) & (cand <= 1.0)):
                proposal = cand
                log_jacobian = (self._dim - 1.0) * np.log(max(1.0e-12, z))
                break
        if proposal is None:
            proposal = np.clip(
                self.param + self._rng.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
            log_jacobian = 0.0

        self.proposed_param = np.asarray(proposal, dtype=float)
        self._last_log_jacobian = float(log_jacobian)
        return self.proposed_param

    def update(self, new_loglikelihood, beta=1.0) -> bool:
        new_loglikelihood = float(new_loglikelihood)
        beta = float(beta)
        accepted = False

        if self.iterations == 0 or self.last_loglikelihood is None:
            accepted = True
        else:
            delta = (new_loglikelihood - float(self.last_loglikelihood)) * beta
            delta += float(self._last_log_jacobian)
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
            "type": "EnsembleChain",
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
            "stretch_a": float(self._stretch_a),
            "last_log_jacobian": float(self._last_log_jacobian),
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
        self._stretch_a = max(1.01, float(state.get("stretch_a", self._stretch_a)))
        self._last_log_jacobian = float(
            state.get("last_log_jacobian", self._last_log_jacobian)
        )
        rng_state = state.get("rng_state")
        if rng_state is not None:
            self._rng.bit_generator.state = rng_state


__all__ = ["EnsembleChain"]
