#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get_float
from jarvishep.Sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep.Sampling.tpmcmc import PTMCMC


class PTEnsemble(PTMCMC):
    """Parallel-tempered EnsembleMCMC sampler."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "PTEnsemble"
        self._stretch_a = 2.0

    def load_schema_file(self):
        self.schema = self.path.get("PTEnsembleSchema", self.path["TPMCMCSchema"])

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

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        if len(self._temperature_ladder) != self._nchains:
            raise ValueError(
                f"temperature_ladder size mismatch for PTEnsemble: expect {self._nchains}, got {len(self._temperature_ladder)}"
            )

        cold_idx = 0
        cold_temp = min(self._temperature_ladder) if self._temperature_ladder else 1.0
        for ii, temp in enumerate(self._temperature_ladder):
            if float(temp) == float(cold_temp):
                cold_idx = ii
                break

        chains = [
            ChainRuntime(
                chain_id=ii,
                engine=None,
                temperature=float(self._temperature_ladder[ii]),
                is_cold=(ii == cold_idx),
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
