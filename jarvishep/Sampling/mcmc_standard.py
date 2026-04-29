#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import time

from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import (
    normalize_proposal_scales,
    parse_common_chain_counts,
    parse_proposal_scale_value,
)
from jarvishep.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep.Sampling.Source.MCMC.state_machine_base import MCMCStateMachineBase
from jarvishep.sample import Sample  # Backward-compatible patch point for tests/tools.


class MCMC(MCMCStateMachineBase):
    def __init__(self, method_name: str = "MCMC") -> None:
        super().__init__()
        self.load_schema_file()
        self.method = str(method_name)
        self._selectionexp = None
        self._nchains = 1
        self._niters = 1
        self._proposal_scales = 0.1

    def load_schema_file(self):
        self.schema = self.path["MCMCSchema"]

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def __iter__(self):
        return self

    def __next__(self):
        if not self._state_machine_ready:
            self._initialize_state_machine()
        cid = self._pop_ready_chain()
        if cid is None:
            raise StopIteration
        chain = self._must_registry().get(cid)
        params = self._propose_params_for_chain(chain)
        self._enqueue_chain(cid)
        return params

    def next_sample(self):
        return self.__next__()

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._dimensions = len(self.vars)
        smp = self.config["Sampling"]["Bounds"]
        self._nchains, self._niters = parse_common_chain_counts(smp)
        self._proposal_scales = parse_proposal_scale_value(smp, default=0.1)
        self._selectionexp = self.config["Sampling"].get("selection")

    def initialize(self):
        self.logger.warning(f"Initializing the {self.method} Sampling")
        t0 = time.time()
        self._initialize_state_machine()
        self.info["t0"] = time.time() - t0
        self.logger.info("{} Sampler initialized in {:.2f} sec".format(self.method, self.info["t0"]))

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        chains = []
        for ii in range(self._nchains):
            engine = MCMCChain(
                np.random.rand(self._dimensions),
                proposal_scales[ii],
                self._niters,
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

    def _max_inflight(self) -> int:
        return self._bounded_inflight(self._nchains)

    def _normalize_proposal_scales(self):
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _next_save_dir(self) -> str | None:
        return self._next_bucket_dir_for_sample()

    def _apply_proposal_scales(self, proposal_scales):
        if isinstance(proposal_scales, (int, float)):
            scale = float(proposal_scales)
            for chain in self._must_registry().all():
                chain.engine.proposal_scale = scale
            self._proposal_scales = scale
            return

        values = normalize_proposal_scales(
            proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )
        for ii, chain in enumerate(self._must_registry().all()):
            chain.engine.proposal_scale = values[ii]
        self._proposal_scales = values

    def finalize(self):
        self._emit_chain_summary(self.method)

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning(f"WorkerFactory is ready for {self.method} sampler")
