#!/usr/bin/env python3

import numpy as np


class AMMCMCChain:
    """Adaptive Metropolis chain with bounded unit-cube proposals."""

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
    ):
        self.param = np.asarray(initial_param, dtype=float)
        self.proposal_scale = float(proposal_scale)
        self.n_iterations = int(n_iterations)
        self.iterations = 0
        self.last_loglikelihood = None
        self.proposed_param = None

        self._dim = int(self.param.shape[0])
        self._adapt_enabled = bool(adapt_enabled)
        self._adapt_start_iter = max(1, int(adapt_start_iter))
        self._adapt_window = max(1, int(adapt_window))
        self._adapt_eps = float(adapt_eps)
        self._adapt_scale = float(adapt_scale)

        self._history = []
        self._cov = (self.proposal_scale ** 2) * np.eye(self._dim)
        self._accepted_since_adapt = 0

    def __iter__(self):
        return self

    def _draw_adaptive_step(self):
        if not self._adapt_enabled or self.iterations < self._adapt_start_iter:
            return np.random.normal(0.0, self.proposal_scale, size=self._dim)

        scale = (self._adapt_scale ** 2) / max(1.0, float(self._dim))
        cov = self._cov * scale
        try:
            return np.random.multivariate_normal(np.zeros(self._dim), cov)
        except Exception:
            return np.random.normal(0.0, self.proposal_scale, size=self._dim)

    def __next__(self):
        if self.iterations >= self.n_iterations:
            raise StopIteration

        if self.iterations == 0:
            proposed_param = np.random.rand(self._dim)
            self.proposed_param = proposed_param
            return proposed_param

        for _ in range(2048):
            step = self._draw_adaptive_step()
            proposed_param = self.param + step
            if np.all((proposed_param >= 0.0) & (proposed_param <= 1.0)):
                self.proposed_param = proposed_param
                return proposed_param

        # Fallback when repeated proposals leave the domain.
        proposed_param = np.clip(self.param + np.random.normal(0.0, self.proposal_scale, size=self._dim), 0.0, 1.0)
        self.proposed_param = proposed_param
        return proposed_param

    def _maybe_adapt_cov(self):
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
                accepted = bool(np.random.rand() < np.exp(delta))

        if accepted:
            self.param = self.proposed_param
            self.last_loglikelihood = new_loglikelihood
            self._history.append(np.array(self.param, dtype=float))
            self._accepted_since_adapt += 1
            self._maybe_adapt_cov()

        self.iterations += 1
        return accepted
