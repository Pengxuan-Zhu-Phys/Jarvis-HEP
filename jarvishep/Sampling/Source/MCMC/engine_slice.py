#!/usr/bin/env python3
from __future__ import annotations

import numpy as np


class SliceChain:
    """Slice-inspired directional proposal chain.

    Notes:
    - Designed for external expensive logL where each proposal is evaluated by
      WorkerFactory.
    - Uses one proposal per iteration with adaptive directional width.
    """

    def __init__(
        self,
        initial_param,
        proposal_scale,
        n_iterations,
        mode="random_direction",
        width=0.2,
        max_steps_out=16,
        max_shrink=32,
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self._dim = int(self.param.shape[0])
        self._mode = str(mode)
        self._width = max(1.0e-4, float(width))
        self._max_steps_out = max(1, int(max_steps_out))
        self._max_shrink = max(1, int(max_shrink))

    def __iter__(self):
        return self

    def _draw_direction(self):
        if self._mode in ("univariate", "coordinate"):
            direction = np.zeros(self._dim, dtype=float)
            idx = np.random.randint(0, self._dim)
            direction[idx] = 1.0 if np.random.rand() < 0.5 else -1.0
            return direction
        vec = np.random.normal(0.0, 1.0, size=self._dim)
        norm = float(np.linalg.norm(vec))
        if norm <= 1.0e-12:
            vec = np.zeros(self._dim, dtype=float)
            vec[np.random.randint(0, self._dim)] = 1.0
            return vec
        return vec / norm

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration
        if self.iterations == 0:
            proposal = np.random.rand(self._dim)
            self.proposed_param = proposal
            return proposal

        proposal = None
        for _ in range(max(self._max_steps_out, self._max_shrink)):
            direction = self._draw_direction()
            step = np.random.uniform(-self._width, self._width)
            cand = self.param + step * direction
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
            self._width = min(0.8, self._width * 1.01)
        else:
            self._width = max(1.0e-4, self._width * 0.95)

        self.iterations += 1
        return accepted
