#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class MALAChain:
    """Gradient-free MALA placeholder chain.

    Uses a local drift estimated from accepted trajectory and Gaussian noise.
    """

    def __init__(self, initial_param, proposal_scale, n_iterations, step_size=0.1):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self._dim = int(self.param.shape[0])
        self._step_size = max(1.0e-4, float(step_size))
        self._prev_accepted_param = None

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
            self.proposed_param = proposal
            return proposal

        drift = np.zeros(self._dim, dtype=float)
        if self._prev_accepted_param is not None:
            drift = np.asarray(self.param - self._prev_accepted_param, dtype=float)

        proposal = None
        for _ in range(2048):
            noise = np.random.normal(0.0, self._step_size, size=self._dim)
            cand = self.param + 0.5 * self._step_size * drift + noise
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
            self._prev_accepted_param = np.asarray(self.param, dtype=float)
            self.param = np.asarray(self.proposed_param, dtype=float)
            self.last_loglikelihood = new_loglikelihood

        self.iterations += 1
        return accepted
