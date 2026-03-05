#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.engine_dream_lite import DREAMLiteChain
from jarvishep.Sampling.demcmc import DEMCMC


class DREAMLite(DEMCMC):
    """DREAM-lite sampler (DE + snooker/global archive moves)."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "DREAMLite"
        self._dream_snooker_prob = 0.1
        self._dream_archive_size = 256

    def load_schema_file(self):
        self.schema = self.path.get("DREAMLiteSchema", self.path.get("DEMCMCSchema", self.path["MCMCSchema"]))

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._dream_snooker_prob = float(smp.get("dream_snooker_prob", 0.1))
        self._dream_archive_size = int(smp.get("dream_archive_size", 256))

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
            registry.get(ii).engine = DREAMLiteChain(
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
            )
        return registry
