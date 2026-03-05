#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class NUTSChain:
    """No-U-Turn style placeholder chain with stochastic tree-depth expansion."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        step_size=0.05,
        max_depth=6,
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self._dim = int(self.param.shape[0])
        self._step_size = max(1.0e-4, float(step_size))
        self._max_depth = max(1, int(max_depth))

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
            self.proposed_param = proposal
            return proposal

        proposal = None
        for _ in range(256):
            depth = np.random.randint(1, self._max_depth + 1)
            step = np.zeros(self._dim, dtype=float)
            scale = self._step_size
            for _level in range(depth):
                direction = 1.0 if np.random.rand() < 0.5 else -1.0
                step += direction * np.random.normal(0.0, scale, size=self._dim)
                scale *= 0.7
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
