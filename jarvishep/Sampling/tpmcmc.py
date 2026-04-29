#!/usr/bin/env python3
from __future__ import annotations

import numpy as np
import time

from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.Source.MCMC.chain_runtime import ChainRegistry, ChainRuntime
from jarvishep.Sampling.Source.MCMC.config_contract import (
    bounds_get,
    bounds_get_int,
    normalize_proposal_scales,
    parse_common_chain_counts,
    parse_proposal_scale_value,
)
from jarvishep.Sampling.Source.MCMC.mcmc_chain import MCMCChain
from jarvishep.Sampling.Source.MCMC.state_machine_base import MCMCStateMachineBase
from jarvishep.sample import Sample  # Backward-compatible patch point for tests/tools.


class PTMCMC(MCMCStateMachineBase):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "PTMCMC"
        self._selectionexp = None
        self._nchains = 2
        self._niters = 1
        self._exchange_interval = 1
        self._proposal_scales = [0.1, 0.1]
        self._temperature_ladder = [1.0, 2.0]
        self._exchange_offset = 0
        self._swap_pairing_mode = "round_robin"
        self._exchange_rule = "metropolis"

    def load_schema_file(self):
        self.schema = self.path["TPMCMCSchema"]

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
        self._exchange_interval = bounds_get_int(
            smp,
            "exchange_interval",
            aliases=("exchange.interval",),
            default=1,
            minimum=1,
        )
        proposal_value = parse_proposal_scale_value(smp, default=0.1)
        self._proposal_scales = normalize_proposal_scales(
            proposal_value,
            nchains=self._nchains,
            sampler_method=self.method,
        )
        self._exchange_rule = str(
            bounds_get(
                smp,
                "exchange_rule",
                aliases=("exchange.rule",),
                default=self._exchange_rule,
            )
        )
        self._swap_pairing_mode = str(
            bounds_get(
                smp,
                "swap_pairing_mode",
                aliases=("swap.pairing_mode",),
                default=self._swap_pairing_mode,
            )
        )
        self._selectionexp = self.config["Sampling"].get("selection")

        ladder = bounds_get(
            smp,
            "temperature_ladder",
            aliases=("temperature.ladder",),
            default=None,
        )
        if ladder is None:
            # Backward compatible default: chain-0 is cold chain.
            ladder = [1.0 + float(ii) for ii in range(self._nchains)]
        self._temperature_ladder = [float(x) for x in ladder]

    def initialize(self):
        self.logger.warning("Initializing the temporal parallel MCMC (PTMCMC) Sampling")
        t0 = time.time()
        self._initialize_state_machine()
        self.info["t0"] = time.time() - t0
        self.logger.info("PTMCMC Sampler initialized in {:.2f} sec".format(self.info["t0"]))

    def _normalize_proposal_scales(self):
        return normalize_proposal_scales(
            self._proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _create_chain_registry(self) -> ChainRegistry:
        proposal_scales = self._normalize_proposal_scales()
        if len(self._temperature_ladder) != self._nchains:
            raise ValueError(
                f"temperature_ladder size mismatch for PTMCMC: expect {self._nchains}, got {len(self._temperature_ladder)}"
            )

        cold_idx = 0
        cold_temp = min(self._temperature_ladder) if self._temperature_ladder else 1.0
        for ii, temp in enumerate(self._temperature_ladder):
            if float(temp) == float(cold_temp):
                cold_idx = ii
                break

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
                    temperature=float(self._temperature_ladder[ii]),
                    is_cold=(ii == cold_idx),
                )
            )
        return ChainRegistry(chains)

    def _max_inflight(self) -> int:
        return self._bounded_inflight(self._nchains)

    def _next_save_dir(self) -> str | None:
        return self._next_bucket_dir_for_sample()

    def _apply_proposal_scales(self, proposal_scales):
        values = normalize_proposal_scales(
            proposal_scales,
            nchains=self._nchains,
            sampler_method=self.method,
        )
        for ii, chain in enumerate(self._must_registry().all()):
            chain.engine.proposal_scale = values[ii]
        self._proposal_scales = values

    def _should_exchange(self) -> bool:
        if int(self._exchange_interval) <= 0:
            return False
        if self._pending_futures:
            return False
        for chain in self._must_registry().all():
            if chain.iter >= int(self._niters):
                continue
            if chain.window_iter < int(self._exchange_interval):
                return False
            if chain.last_logl is None:
                return False
        return True

    def _pair_chain_ids(self):
        ids = self._must_registry().ids()
        if len(ids) < 2:
            return []
        offset = self._exchange_offset % len(ids)
        rotated = ids[offset:] + ids[:offset]
        pairs = []
        for ii in range(0, len(rotated) - 1, 2):
            pairs.append((rotated[ii], rotated[ii + 1]))
        self._exchange_offset = (self._exchange_offset + 1) % max(len(ids), 1)
        return pairs

    @staticmethod
    def _swap_acceptance(logl1: float, temp1: float, logl2: float, temp2: float) -> bool:
        beta1 = 1.0 / float(temp1) if float(temp1) > 0 else 1.0
        beta2 = 1.0 / float(temp2) if float(temp2) > 0 else 1.0
        delta = (beta1 - beta2) * (float(logl2) - float(logl1))
        if delta >= 0.0:
            return True
        return bool(np.random.rand() < np.exp(delta))

    def _attempt_swap(self, cid1: int, cid2: int) -> bool:
        chain1 = self._must_registry().get(cid1)
        chain2 = self._must_registry().get(cid2)
        if chain1.last_logl is None or chain2.last_logl is None:
            return False

        accepted = self._swap_acceptance(
            chain1.last_logl,
            chain1.temperature,
            chain2.last_logl,
            chain2.temperature,
        )
        if not accepted:
            return False

        chain1.engine.param, chain2.engine.param = chain2.engine.param, chain1.engine.param
        chain1.engine.last_loglikelihood, chain2.engine.last_loglikelihood = (
            chain2.engine.last_loglikelihood,
            chain1.engine.last_loglikelihood,
        )
        chain1.last_logl, chain2.last_logl = chain2.last_logl, chain1.last_logl
        return True

    def _do_exchange(self):
        attempted = 0
        accepted = 0
        for cid1, cid2 in self._pair_chain_ids():
            attempted += 1
            if self._attempt_swap(cid1, cid2):
                accepted += 1

        self._ready_queue.clear()
        self._ready_set.clear()
        for chain in self._must_registry().all():
            chain.window_iter = 0
            if chain.iter < int(self._niters):
                self._enqueue_chain(chain.chain_id)

        self.logger.warning(
            format_two_column_log(
                "PTMCMC exchange summary",
                [("attempted", attempted), ("accepted", accepted)],
            )
        )
        return {"attempted": attempted, "accepted": accepted}

    def finalize(self):
        self._emit_chain_summary("PTMCMC")

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for PTMCMC sampler")
