#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class HMCChain:
    """Gradient-free HMC placeholder chain with leapfrog-style random dynamics."""

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        step_size=0.05,
        leapfrog_steps=8,
        gradient_provider=None,
        gradient_clip_norm=None,
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None
        self.gradient_provider = gradient_provider
        self.gradient_clip_norm = gradient_clip_norm
        self.gradient_contract_level = "placeholder"

        self._dim = int(self.param.shape[0])
        self._step_size = max(1.0e-4, float(step_size))
        self._leapfrog_steps = max(1, int(leapfrog_steps))

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
            pos = np.asarray(self.param, dtype=float)
            mom = np.random.normal(0.0, 1.0, size=self._dim)
            for _k in range(self._leapfrog_steps):
                mom = mom + np.random.normal(0.0, 0.05, size=self._dim)
                pos = pos + self._step_size * mom
            if np.all((pos >= 0.0) & (pos <= 1.0)):
                proposal = pos
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
