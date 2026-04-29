#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get, normalize_proposal_scales
from jarvishep.Sampling.Source.MCMC.engine_ess import ESSChain
from jarvishep.Sampling.mcmc_standard import MCMC


class ESS(MCMC):
    """ESS-style sampler with prior covariance aware elliptical proposals."""

    def __init__(self) -> None:
        super().__init__(method_name="ESS")
        self._ess_prior_cov = None

    def load_schema_file(self):
        self.schema = self.path.get("ESSSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._ess_prior_cov = bounds_get(
            smp,
            "ess_prior_cov",
            aliases=("ess.prior_cov",),
            default=(
                smp.get("ess", {}).get("prior_cov", None)
                if isinstance(smp.get("ess", {}), dict)
                else None
            ),
        )

    def _normalize_proposal_scales(self):
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = []
        for ii in range(self._nchains):
            engine = ESSChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                prior_cov=self._ess_prior_cov,
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
