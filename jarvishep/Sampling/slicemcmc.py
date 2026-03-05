#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
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
        slice_cfg = smp.get("slice", {})
        if not isinstance(slice_cfg, dict):
            slice_cfg = {}
        self._slice_mode = str(smp.get("slice_mode", slice_cfg.get("mode", "random_direction")))
        self._slice_width = float(smp.get("slice_width", slice_cfg.get("width", 0.2)))
        self._slice_max_steps_out = int(
            smp.get("slice_max_steps_out", slice_cfg.get("max_steps_out", 16))
        )
        self._slice_max_shrink = int(smp.get("slice_max_shrink", slice_cfg.get("max_shrink", 32)))

    def _normalize_proposal_scales(self):
        scales = self._proposal_scales
        if isinstance(scales, (int, float)):
            return [float(scales) for _ in range(self._nchains)]
        values = [float(x) for x in scales]
        if len(values) != self._nchains:
            raise ValueError(
                f"proposal_scale size mismatch for SliceMCMC: expect {self._nchains}, got {len(values)}"
            )
        return values

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
