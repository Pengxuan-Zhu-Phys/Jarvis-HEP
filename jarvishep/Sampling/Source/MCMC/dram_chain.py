#!/usr/bin/env python3

import numpy as np

from .ammcmc_chain import AMMCMCChain


class DRAMChain(AMMCMCChain):
    """DRAM-lite chain: adaptive MH + rejection-triggered scale reduction."""

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
        dr_scale_factors=None,
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
        if dr_scale_factors is None:
            dr_scale_factors = [1.0, 0.5]
        self._dr_scale_factors = [max(1e-6, float(x)) for x in dr_scale_factors]
        if not self._dr_scale_factors:
            self._dr_scale_factors = [1.0, 0.5]
        self._reject_streak = 0

    def _draw_adaptive_step(self):
        step = super()._draw_adaptive_step()
        stage = min(self._reject_streak, len(self._dr_scale_factors) - 1)
        return step * self._dr_scale_factors[stage]

    def update(self, new_loglikelihood, beta=1.0):
        accepted = super().update(new_loglikelihood, beta=beta)
        if accepted:
            self._reject_streak = 0
        else:
            self._reject_streak = min(self._reject_streak + 1, len(self._dr_scale_factors) - 1)
        return accepted
