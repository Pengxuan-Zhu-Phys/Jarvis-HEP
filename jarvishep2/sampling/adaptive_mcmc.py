"""Shared base for adaptive covariance MCMC samplers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.sampling.Source.MCMC.config_contract import (
    bounds_get_bool,
    bounds_get_float,
    bounds_get_int,
)


class AdaptiveMCMCBase(MCMCBaseSampler):
    """Common configuration for AM, AMMCMC, and DRAM."""

    _checkpoint_excluded_attributes = frozenset(
        {
            "_adapt_enabled",
            "_adapt_start_iter",
            "_adapt_window",
            "_adapt_eps",
            "_adapt_scale",
        }
    )

    def __init__(self) -> None:
        super().__init__()
        self._adapt_enabled = True
        self._adapt_start_iter = 100
        self._adapt_window = 25
        self._adapt_eps = 1e-6
        self._adapt_scale = 2.38

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._adapt_enabled = bounds_get_bool(
            bounds, "adapt_enabled", aliases=("adapt.enabled",), default=True
        )
        self._adapt_start_iter = bounds_get_int(
            bounds,
            "adapt_start_iter",
            aliases=("adapt.start_iter",),
            default=100,
            minimum=1,
        )
        self._adapt_window = bounds_get_int(
            bounds,
            "adapt_window",
            aliases=("adapt.window",),
            default=25,
            minimum=1,
        )
        self._adapt_eps = bounds_get_float(
            bounds,
            "adapt_eps",
            aliases=("adapt.eps",),
            default=1e-6,
            minimum=0.0,
        )
        self._adapt_scale = bounds_get_float(
            bounds,
            "adapt_scale",
            aliases=("adapt.scale",),
            default=2.38,
            minimum=1e-9,
        )


__all__ = ["AdaptiveMCMCBase"]
