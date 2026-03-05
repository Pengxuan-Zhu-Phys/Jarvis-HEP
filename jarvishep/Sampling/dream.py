#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.engine_dream import DREAMChain
from jarvishep.Sampling.demcmc import DEMCMC


class DREAM(DEMCMC):
    """DREAM sampler (DE + snooker + archive + adaptive crossover)."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "DREAM"
        self._dream_snooker_prob = 0.1
        self._dream_archive_size = 512
        self._dream_crossover_values = [0.2, 0.5, 0.9]
        self._dream_crossover_adapt_interval = 64
        self._dream_scale_jitter = 0.1

    def load_schema_file(self):
        self.schema = self.path.get("DREAMSchema", self.path.get("DREAMLiteSchema", self.path.get("DEMCMCSchema", self.path["MCMCSchema"])))

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._dream_snooker_prob = float(smp.get("dream_snooker_prob", 0.1))
        self._dream_archive_size = int(smp.get("dream_archive_size", 512))

        values = smp.get("dream_crossover_values", [0.2, 0.5, 0.9])
        if isinstance(values, (int, float)):
            values = [float(values)]
        self._dream_crossover_values = [float(x) for x in values]
        if not self._dream_crossover_values:
            self._dream_crossover_values = [0.9]

        self._dream_crossover_adapt_interval = int(smp.get("dream_crossover_adapt_interval", 64))
        self._dream_scale_jitter = float(smp.get("dream_scale_jitter", 0.1))

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
            registry.get(ii).engine = DREAMChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                chain_id=ii,
                population_getter=_population_getter,
                de_gamma=self._de_gamma,
                de_noise=self._de_noise,
                de_crossover=self._de_crossover,
                dream_snooker_prob=self._dream_snooker_prob,
                dream_archive_size=self._dream_archive_size,
                dream_crossover_values=self._dream_crossover_values,
                dream_crossover_adapt_interval=self._dream_crossover_adapt_interval,
                dream_scale_jitter=self._dream_scale_jitter,
            )
        return registry
