#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import (
    bounds_get,
    bounds_get_float,
    bounds_get_int,
    normalize_proposal_scales,
)
from jarvishep.Sampling.Source.MCMC.engine_nuts import NUTSChain
from jarvishep.Sampling.Source.MCMC.gradient_provider import validate_gradient_provider
from jarvishep.Sampling.mcmc_standard import MCMC


class NUTS(MCMC):
    """NUTS placeholder sampler with adaptive-tree proposal depth control."""

    def __init__(self) -> None:
        super().__init__(method_name="NUTS")
        self._nuts_step_size = 0.05
        self._nuts_max_depth = 6
        self._gradient_provider = None
        self._gradient_provider_spec = None
        self._gradient_clip_norm = None
        self.gradient_contract_level = "placeholder-reference"
        self.gradient_contract_target = "v1.7.0"

    def load_schema_file(self):
        self.schema = self.path.get("NUTSSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._nuts_step_size = bounds_get_float(
            smp,
            "nuts_step_size",
            aliases=("nuts.step_size", "step_size"),
            default=0.05,
            minimum=1e-6,
        )
        self._nuts_max_depth = bounds_get_int(
            smp,
            "nuts_max_depth",
            aliases=("nuts.max_depth", "max_depth"),
            default=6,
            minimum=1,
        )
        provider = bounds_get(
            smp,
            "gradient_provider",
            aliases=("grad.provider", "gradient.provider"),
            default=None,
        )
        self._gradient_provider_spec = provider
        self._gradient_provider = None
        if provider is not None and not isinstance(provider, str):
            validate_gradient_provider(provider, method_name=self.method)
            self._gradient_provider = provider
        self._gradient_clip_norm = bounds_get(
            smp,
            "gradient_clip_norm",
            aliases=("grad.clip_norm", "gradient.clip_norm"),
            default=None,
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
            engine = NUTSChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                step_size=self._nuts_step_size,
                max_depth=self._nuts_max_depth,
                gradient_provider=self._gradient_provider,
                gradient_clip_norm=self._gradient_clip_norm,
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
