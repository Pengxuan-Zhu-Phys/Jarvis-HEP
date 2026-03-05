#!/usr/bin/env python3

import numpy as np

from .ammcmc_chain import AMMCMCChain


class RobustAMChain(AMMCMCChain):
    """Adaptive Metropolis with occasional heavy-tail/global-jump proposals."""

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
        global_jump_prob=0.08,
        heavy_tail_df=5.0,
        global_scale=1.8,
    ):
        super().__init__(
            initial_param=initial_param,
            proposal_scale=proposal_scale,
            n_iterations=n_iterations,
            adapt_enabled=adapt_enabled,
            adapt_start_iter=adapt_start_iter,
            adapt_window=adapt_window,
            adapt_eps=adapt_eps,
            adapt_scale=adapt_scale,
        )
        self._global_jump_prob = float(global_jump_prob)
        self._heavy_tail_df = max(1e-6, float(heavy_tail_df))
        self._global_scale = max(1e-6, float(global_scale))

    def _draw_global_jump(self):
        base_cov = np.asarray(self._cov, dtype=float) * self._global_scale
        if np.ndim(base_cov) == 0:
            base_cov = np.array([[float(base_cov)]], dtype=float)
        if base_cov.shape != (self._dim, self._dim):
            base_cov = (self.proposal_scale ** 2) * np.eye(self._dim)

        z = np.random.multivariate_normal(np.zeros(self._dim), base_cov)
        chi = np.random.chisquare(self._heavy_tail_df)
        scale = np.sqrt(max(1e-12, chi / self._heavy_tail_df))
        return z / scale

    def _draw_adaptive_step(self):
        if np.random.rand() < self._global_jump_prob:
            return self._draw_global_jump()
        return super()._draw_adaptive_step()
