#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
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
        ess_cfg = smp.get("ess", {})
        if not isinstance(ess_cfg, dict):
            ess_cfg = {}
        self._ess_prior_cov = smp.get("ess_prior_cov", ess_cfg.get("prior_cov", None))

    def _normalize_proposal_scales(self):
        scales = self._proposal_scales
        if isinstance(scales, (int, float)):
            return [float(scales) for _ in range(self._nchains)]
        values = [float(x) for x in scales]
        if len(values) != self._nchains:
            raise ValueError(
                f"proposal_scale size mismatch for ESS: expect {self._nchains}, got {len(values)}"
            )
        return values

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
