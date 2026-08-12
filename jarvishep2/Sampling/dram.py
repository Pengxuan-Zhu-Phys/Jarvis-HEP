"""Delayed-rejection adaptive Metropolis sampler."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.adaptive_mcmc import AdaptiveMCMCBase
from jarvishep2.Sampling.Source.MCMC.config_contract import bounds_get_int, bounds_get_list
from jarvishep2.Sampling.Source.MCMC.dram_chain import DRAMChain


class DRAMSampler(AdaptiveMCMCBase):
    method = "DRAM"
    _checkpoint_excluded_attributes = frozenset(
        {"_dr_steps", "_dr_scale_factors"}
    )

    def __init__(self) -> None:
        super().__init__()
        self._dr_steps = 2
        self._dr_scale_factors = [1.0, 0.5]

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        super()._configure_method(bounds)
        self._dr_steps = bounds_get_int(
            bounds, "dr_steps", aliases=("dr.steps",), default=2, minimum=1
        )
        factors = bounds_get_list(
            bounds, "dr_scale_factors", aliases=("dr.scale_factors",), default=[1.0, 0.5]
        )
        if isinstance(factors, (int, float)):
            factors = [1.0, float(factors)]
        self._dr_scale_factors = [max(1e-6, float(x)) for x in factors]

    def _make_engine(self, chain_id, scale, rng):  # noqa: ANN001
        initial = rng.random(self._dim)
        return DRAMChain(
            initial,
            scale,
            self._niters,
            adapt_enabled=self._adapt_enabled,
            adapt_start_iter=self._adapt_start_iter,
            adapt_window=self._adapt_window,
            adapt_eps=self._adapt_eps,
            adapt_scale=self._adapt_scale,
            dr_steps=self._dr_steps,
            dr_scale_factors=self._dr_scale_factors,
            rng=rng,
        )


def create_dram() -> DRAMSampler:
    return DRAMSampler()


__all__ = ["DRAMSampler", "create_dram"]
