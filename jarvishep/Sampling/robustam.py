#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get_float
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
        self._global_jump_prob = bounds_get_float(
            smp,
            "global_jump_prob",
            aliases=("robust.global_jump_prob",),
            default=0.08,
            minimum=0.0,
        )
        self._heavy_tail_df = bounds_get_float(
            smp,
            "heavy_tail_df",
            aliases=("robust.heavy_tail_df",),
            default=5.0,
            minimum=1e-9,
        )
        self._global_scale = bounds_get_float(
            smp,
            "global_scale",
            aliases=("robust.global_scale",),
            default=1.8,
            minimum=1e-9,
        )

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = []
        for ii in range(self._nchains):
            engine = RobustAMChain(
                np.random.rand(self._dimensions),
                proposal_scales[ii],
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
