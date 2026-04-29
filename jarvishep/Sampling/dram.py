#!/usr/bin/env python3
from __future__ import annotations

import numpy as np

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import bounds_get_int, bounds_get_list
from jarvishep.Sampling.Source.MCMC.dram_chain import DRAMChain
from jarvishep.Sampling.Source.MCMC.state_machine_multistage_base import (
    MCMCMultiStageStateMachineBase,
)
from jarvishep.Sampling.ammcmc import AMMCMC


class DRAM(AMMCMC):
    """Full DRAM sampler with delayed rejection and adaptive covariance."""

    def __init__(self) -> None:
        super().__init__()
        self.method = "DRAM"
        self._dr_steps = 2
        self._dr_scale_factors = [1.0, 0.5]
        self._future_to_stage_idx = {}
        self._chain_stage_idx = {}

    def init_generator(self):
        super().init_generator()
        smp = self.config["Sampling"]["Bounds"]
        self._dr_steps = bounds_get_int(
            smp,
            "dr_steps",
            aliases=("dr.steps",),
            default=2,
            minimum=1,
        )
        factors = bounds_get_list(
            smp,
            "dr_scale_factors",
            aliases=("dr.scale_factors",),
            default=[1.0, 0.5],
        )
        if isinstance(factors, (int, float)):
            factors = [1.0, float(factors)]
        self._dr_scale_factors = [max(1e-6, float(x)) for x in factors]

    def _initialize_state_machine(self) -> None:
        MCMCMultiStageStateMachineBase._initialize_state_machine(self)

    def _export_runtime_extras(self):
        return MCMCMultiStageStateMachineBase._export_runtime_extras(self)

    def _import_runtime_extras(self, payload):
        MCMCMultiStageStateMachineBase._import_runtime_extras(self, payload)

    def _initial_stage_index(self, chain: ChainRuntime) -> int:
        _ = chain
        return 0

    def _proposal_unit_point_for_stage(self, chain: ChainRuntime, stage_index: int):
        return chain.engine.propose_stage(stage_index)

    def _propose_params_for_stage(self, chain: ChainRuntime, stage_index: int):
        while True:
            proposal = self._proposal_unit_point_for_stage(chain, stage_index)
            params = self.map_point_into_distribution(proposal)
            if self._selectionexp:
                if self.evaluate_selection(self._selectionexp, params):
                    return params
                continue
            return params

    def _consume_stage_result(self, chain: ChainRuntime, sample_info, stage_index):
        logl_new = self._extract_logl(sample_info)
        beta = 1.0
        if chain.temperature > 0:
            beta = 1.0 / float(chain.temperature)
        outcome = chain.engine.consume_stage_result(stage_index, logl_new, beta=beta)
        outcome.setdefault("logl", getattr(chain.engine, "last_loglikelihood", logl_new))
        return outcome

    def run_nested(self):
        return MCMCMultiStageStateMachineBase.run_nested(self)

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = []
        for ii in range(self._nchains):
            engine = DRAMChain(
                np.random.rand(self._dimensions),
                proposal_scales[ii],
                self._niters,
                adapt_enabled=self._adapt_enabled,
                adapt_start_iter=self._adapt_start_iter,
                adapt_window=self._adapt_window,
                adapt_eps=self._adapt_eps,
                adapt_scale=self._adapt_scale,
                dr_steps=self._dr_steps,
                dr_scale_factors=self._dr_scale_factors,
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
