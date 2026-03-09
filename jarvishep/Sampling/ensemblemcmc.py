#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get_float, normalize_proposal_scales
from jarvishep.Sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep.Sampling.mcmc_standard import MCMC


class EnsembleMCMC(MCMC):
    """Ensemble stretch-move MCMC sampler."""

    def __init__(self) -> None:
        super().__init__(method_name="EnsembleMCMC")
        self._stretch_a = 2.0

    def load_schema_file(self):
        self.schema = self.path.get("EnsembleSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._stretch_a = bounds_get_float(
            smp,
            "stretch_a",
            aliases=("ensemble.stretch_a",),
            default=2.0,
            minimum=1.01,
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
            registry.get(ii).engine = EnsembleChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                chain_id=ii,
                population_getter=_population_getter,
                stretch_a=self._stretch_a,
            )
        return registry
