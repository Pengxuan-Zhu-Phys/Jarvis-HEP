#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.robustam_chain import RobustAMChain
from jarvishep.Sampling.ammcmc import AMMCMC


class RobustAM(AMMCMC):
    """Adaptive Metropolis with heavy-tail/global-jump mixture proposals."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "RobustAM"
        self._global_jump_prob = 0.08
        self._heavy_tail_df = 5.0
        self._global_scale = 1.8

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._global_jump_prob = float(smp.get("global_jump_prob", 0.08))
        self._heavy_tail_df = float(smp.get("heavy_tail_df", 5.0))
        self._global_scale = float(smp.get("global_scale", 1.8))

    def _create_chain_registry(self) -> ChainRegistry:
        chains = []
        for ii in range(self._nchains):
            engine = RobustAMChain(
                np.random.rand(self._dimensions),
                self._proposal_scales,
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
                global_jump_prob=self._global_jump_prob,
                heavy_tail_df=self._heavy_tail_df,
                global_scale=self._global_scale,
            )
            chains.append(
                ChainRuntime(
                    chain_id=ii,
                    engine=engine,
                    temperature=1.0,
                    is_cold=(ii == 0),
                )
            )
        return ChainRegistry(chains)
