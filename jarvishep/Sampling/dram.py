#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep.Sampling.ammcmc import AMMCMC


class DRAM(AMMCMC):
    """DRAM-lite sampler (adaptive MH + delayed-rejection style scaling)."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "DRAM"
        self._dr_scale_factors = [1.0, 0.5]

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        factors = smp.get("dr_scale_factors", [1.0, 0.5])
        if isinstance(factors, (int, float)):
            factors = [1.0, float(factors)]
        self._dr_scale_factors = [float(x) for x in factors]

    def _create_chain_registry(self) -> ChainRegistry:
        chains = []
        for ii in range(self._nchains):
            engine = DRAMChain(
                np.random.rand(self._dimensions),
                self._proposal_scales,
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
                dr_scale_factors=self._dr_scale_factors,
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
