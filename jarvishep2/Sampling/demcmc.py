"""Differential-evolution MCMC sampler."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.Sampling.Source.MCMC.config_contract import bounds_get_float
from jarvishep2.Sampling.Source.MCMC.engine_demcmc import DEMCMCChain


class DEMCMCSampler(MCMCBaseSampler):
    method = "DEMCMC"
    _checkpoint_excluded_attributes = frozenset(
        {"_de_gamma", "_de_noise", "_de_crossover"}
    )

    def __init__(self) -> None:
        super().__init__()
        self._de_gamma = 0.0
        self._de_noise = 1.0e-3
        self._de_crossover = 1.0

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._de_gamma = bounds_get_float(
            bounds, "de_gamma", aliases=("de.gamma",), default=0.0, minimum=0.0
        )
        self._de_noise = bounds_get_float(
            bounds, "de_noise", aliases=("de.noise",), default=1.0e-3, minimum=0.0
        )
        self._de_crossover = bounds_get_float(
            bounds,
            "de_crossover",
            aliases=("de.crossover",),
            default=1.0,
            minimum=0.0,
        )
        if self._nchains < 2:
            raise ValueError(f"{self.method} requires num_chains >= 2")

    def _uses_half_ensemble(self) -> bool:
        return True

    def _make_engine(self, chain_id, scale, rng):  # noqa: ANN001
        initial = rng.random(self._dim)
        return DEMCMCChain(
            initial,
            scale,
            self._niters,
            chain_id=chain_id,
            population_getter=self._population_for,
            de_gamma=self._de_gamma,
            de_noise=self._de_noise,
            de_crossover=self._de_crossover,
            rng=rng,
        )


def create_demcmc() -> DEMCMCSampler:
    return DEMCMCSampler()


__all__ = ["DEMCMCSampler", "create_demcmc"]
