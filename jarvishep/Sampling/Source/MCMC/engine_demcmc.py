#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class DEMCMCChain:
    """Differential Evolution MCMC chain engine."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        chain_id,
        population_getter,
        de_gamma=0.0,
        de_noise=1.0e-3,
        de_crossover=1.0,
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
        if float(de_gamma) > 0:
            self._gamma = float(de_gamma)
        else:
            self._gamma = float(2.38 / np.sqrt(max(1.0, 2.0 * self._dim)))
        self._noise = max(1.0e-9, float(de_noise))
        self._crossover = min(1.0, max(1.0e-6, float(de_crossover)))

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

    def _draw_demove(self):
        others = self._other_population()
        if len(others) < 2:
            return np.random.normal(0.0, self.proposal_scale, size=self._dim)

        ids = np.random.choice(len(others), size=2, replace=False)
        xr1 = others[int(ids[0])]
        xr2 = others[int(ids[1])]
        diff = xr1 - xr2
        step = self._gamma * diff + self._noise * np.random.normal(0.0, 1.0, size=self._dim)
        if self._crossover < 1.0:
            mask = np.random.rand(self._dim) < self._crossover
            if not np.any(mask):
                mask[np.random.randint(0, self._dim)] = True
            step = step * mask
        return step

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
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
                self.param + np.random.normal(0.0, self.proposal_scale, size=self._dim),
                0.0,
                1.0,
            )
        self.proposed_param = np.asarray(proposal, dtype=float)
        return self.proposed_param

    def update(self, new_loglikelihood, beta=1.0):
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
                accepted = bool(np.random.rand() < np.exp(np.clip(delta, -700.0, 0.0)))

        if accepted:
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood

        self.iterations += 1
        return accepted
