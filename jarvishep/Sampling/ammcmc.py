#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.ammcmc_chain import AMMCMCChain
from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.mcmc_standard import MCMC


class AMMCMC(MCMC):
    """Adaptive Metropolis sampler on top of MCMC state-machine runtime."""

    def __init__(self) -> None:
        super().__init__(method_name="AMMCMC")
        self._adapt_enabled = True
        self._adapt_start_iter = 100
        self._adapt_window = 25
        self._adapt_eps = 1e-6
        self._adapt_scale = 2.38

    def load_schema_file(self):
        # Reuse MCMC schema; adaptive keys are optional bounds extras.
        self.schema = self.path["MCMCSchema"]

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._adapt_enabled = bool(smp.get("adapt_enabled", True))
        self._adapt_start_iter = int(smp.get("adapt_start_iter", 100))
        self._adapt_window = int(smp.get("adapt_window", 25))
        self._adapt_eps = float(smp.get("adapt_eps", 1e-6))
        self._adapt_scale = float(smp.get("adapt_scale", 2.38))

    def _create_chain_registry(self) -> ChainRegistry:
        chains = []
        for ii in range(self._nchains):
            engine = AMMCMCChain(
                np.random.rand(self._dimensions),
                self._proposal_scales,
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
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
