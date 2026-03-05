#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.engine_hmc import HMCChain
from jarvishep.Sampling.mcmc_standard import MCMC


class HMC(MCMC):
    """HMC placeholder sampler with unified external sampler interface."""

    def __init__(self) -> None:
        super().__init__(method_name="HMC")
        self._hmc_step_size = 0.05
        self._hmc_leapfrog_steps = 8

    def load_schema_file(self):
        self.schema = self.path.get("HMCSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._hmc_step_size = float(smp.get("hmc_step_size", smp.get("step_size", 0.05)))
        self._hmc_leapfrog_steps = int(
            smp.get("hmc_leapfrog_steps", smp.get("leapfrog_steps", 8))
        )

    def _normalize_proposal_scales(self):
        scales = self._proposal_scales
        if isinstance(scales, (int, float)):
            return [float(scales) for _ in range(self._nchains)]
        values = [float(x) for x in scales]
        if len(values) != self._nchains:
            raise ValueError(
                f"proposal_scale size mismatch for HMC: expect {self._nchains}, got {len(values)}"
            )
        return values

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = []
        for ii in range(self._nchains):
            engine = HMCChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                step_size=self._hmc_step_size,
                leapfrog_steps=self._hmc_leapfrog_steps,
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
