#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get_float, normalize_proposal_scales
from jarvishep.Sampling.Source.MCMC.engine_demcmc import DEMCMCChain
from jarvishep.Sampling.mcmc_standard import MCMC


class DEMCMC(MCMC):
    """Differential Evolution MCMC on top of shared state-machine runtime."""

    def __init__(self) -> None:
        super().__init__(method_name="DEMCMC")
        self._de_gamma = 0.0
        self._de_noise = 1.0e-3
        self._de_crossover = 1.0

    def load_schema_file(self):
        self.schema = self.path.get("DEMCMCSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._de_gamma = bounds_get_float(
            smp,
            "de_gamma",
            aliases=("de.gamma",),
            default=0.0,
            minimum=0.0,
        )
        self._de_noise = bounds_get_float(
            smp,
            "de_noise",
            aliases=("de.noise",),
            default=1.0e-3,
            minimum=0.0,
        )
        self._de_crossover = bounds_get_float(
            smp,
            "de_crossover",
            aliases=("de.crossover",),
            default=1.0,
            minimum=0.0,
        )

    def _normalize_proposal_scales(self):
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = [
            ChainRuntime(
                chain_id=ii,
                engine=None,
                temperature=1.0,
                is_cold=(ii == 0),
            )
            for ii in range(self._nchains)
        ]
        registry = ChainRegistry(chains)

        def _population_getter(cid):
            return [
                np.asarray(registry.get(jj).engine.param, dtype=float)
                for jj in registry.ids()
                if jj != int(cid) and registry.get(jj).engine is not None
            ]

        for ii in registry.ids():
            registry.get(ii).engine = DEMCMCChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                chain_id=ii,
                population_getter=_population_getter,
                de_gamma=self._de_gamma,
                de_noise=self._de_noise,
                de_crossover=self._de_crossover,
            )
        return registry
