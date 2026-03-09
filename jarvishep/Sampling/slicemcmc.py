#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import (
    bounds_get_float,
    bounds_get_int,
    normalize_proposal_scales,
)
from jarvishep.Sampling.Source.MCMC.engine_slice import SliceChain
from jarvishep.Sampling.mcmc_standard import MCMC


class SliceMCMC(MCMC):
    """Slice-inspired MCMC sampler with directional proposals."""

    def __init__(self) -> None:
        super().__init__(method_name="SliceMCMC")
        self._slice_mode = "random_direction"
        self._slice_width = 0.2
        self._slice_max_steps_out = 16
        self._slice_max_shrink = 32

    def load_schema_file(self):
        self.schema = self.path.get("SliceSchema", self.path["MCMCSchema"])

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._slice_mode = str(
            smp.get(
                "slice_mode",
                smp.get("slice", {}).get("mode", "random_direction")
                if isinstance(smp.get("slice", {}), dict)
                else "random_direction",
            )
        )
        self._slice_width = bounds_get_float(
            smp,
            "slice_width",
            aliases=("slice.width",),
            default=(
                smp.get("slice", {}).get("width", 0.2)
                if isinstance(smp.get("slice", {}), dict)
                else 0.2
            ),
            minimum=1e-6,
        )
        self._slice_max_steps_out = bounds_get_int(
            smp,
            "slice_max_steps_out",
            aliases=("slice.max_steps_out",),
            default=(
                smp.get("slice", {}).get("max_steps_out", 16)
                if isinstance(smp.get("slice", {}), dict)
                else 16
            ),
            minimum=1,
        )
        self._slice_max_shrink = bounds_get_int(
            smp,
            "slice_max_shrink",
            aliases=("slice.max_shrink",),
            default=(
                smp.get("slice", {}).get("max_shrink", 32)
                if isinstance(smp.get("slice", {}), dict)
                else 32
            ),
            minimum=1,
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
            engine = SliceChain(
                initial_param=np.random.rand(self._dimensions),
                proposal_scale=proposal_scales[ii],
                n_iterations=self._niters,
                mode=self._slice_mode,
                width=self._slice_width,
                max_steps_out=self._slice_max_steps_out,
                max_shrink=self._slice_max_shrink,
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
