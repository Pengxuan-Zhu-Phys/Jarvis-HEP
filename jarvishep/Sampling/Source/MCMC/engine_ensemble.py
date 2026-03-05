#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class EnsembleChain:
    """Goodman-Weare stretch-move engine for one walker."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        chain_id,
        population_getter,
        stretch_a=2.0,
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self.chain_id = int(chain_id)
        self._population_getter = population_getter
        self._dim = int(self.param.shape[0])
        self._stretch_a = max(1.01, float(stretch_a))
        self._last_log_jacobian = 0.0

    def __iter__(self):
        return self

    def _other_population(self):
        try:
            population = self._population_getter(self.chain_id)
        except TypeError:
            population = self._population_getter()
        if population is None:
            return []
        out = []
        for point in population:
            arr = np.asarray(point, dtype=float)
            if arr.shape == self.param.shape:
                out.append(arr)
        return out

    def _draw_stretch_factor(self):
        u = np.random.rand()
        a = self._stretch_a
        return (((a - 1.0) * u + 1.0) ** 2) / a

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration

        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
            self.proposed_param = proposal
            self._last_log_jacobian = 0.0
            return proposal

        others = self._other_population()
        if not others:
            proposal = np.clip(
                self.param + np.random.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
            self.proposed_param = proposal
            self._last_log_jacobian = 0.0
            return proposal

        proposal = None
        log_jacobian = 0.0
        for _ in range(2048):
            partner = np.asarray(others[np.random.randint(0, len(others))], dtype=float)
            z = float(self._draw_stretch_factor())
            cand = partner + z * (self.param - partner)
            if np.all((cand >= 0.0) & (cand <= 1.0)):
                proposal = cand
                log_jacobian = (self._dim - 1.0) * np.log(max(1.0e-12, z))
                break
        if proposal is None:
            proposal = np.clip(
                self.param + np.random.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
            log_jacobian = 0.0

        self.proposed_param = np.asarray(proposal, dtype=float)
        self._last_log_jacobian = float(log_jacobian)
        return self.proposed_param

    def update(self, new_loglikelihood, beta=1.0):
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
                accepted = bool(np.random.rand() < np.exp(np.clip(delta, -700.0, 0.0)))

        if accepted:
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood

        self.iterations += 1
        return accepted
