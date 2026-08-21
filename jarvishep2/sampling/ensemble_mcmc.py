"""Stretch-move ensemble MCMC sampler."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.sampling.Source.MCMC.config_contract import bounds_get_float
from jarvishep2.sampling.Source.MCMC.engine_ensemble import EnsembleChain


class EnsembleMCMCBase(MCMCBaseSampler):
    """Common stretch-move configuration and engine construction."""

    _checkpoint_excluded_attributes = frozenset({"_stretch_a"})

    def __init__(self) -> None:
        super().__init__()
        self._stretch_a = 2.0

    def _configure_ensemble(self, bounds: Mapping[str, Any]) -> None:
        self._stretch_a = bounds_get_float(
            bounds,
            "stretch_a",
            aliases=("ensemble.stretch_a",),
            default=2.0,
            minimum=1.01,
        )
        if self._nchains < 2:
            raise ValueError(f"{self.method} requires num_chains >= 2")

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._configure_ensemble(bounds)

    def _uses_half_ensemble(self) -> bool:
        return True

    def _make_engine(self, chain_id, scale, rng):  # noqa: ANN001
        initial = rng.random(self._dim)
        return EnsembleChain(
            initial,
            scale,
            self._niters,
            chain_id=chain_id,
            population_getter=self._population_for,
            stretch_a=self._stretch_a,
            rng=rng,
        )


class EnsembleMCMCSampler(EnsembleMCMCBase):
    method = "EnsembleMCMC"


class EnsembleSampler(EnsembleMCMCSampler):
    method = "Ensemble"


def create_ensemble() -> EnsembleMCMCSampler:
    return EnsembleMCMCSampler()


def create_ensemble_alias() -> EnsembleSampler:
    return EnsembleSampler()


__all__ = [
    "EnsembleMCMCBase",
    "EnsembleMCMCSampler",
    "EnsembleSampler",
    "create_ensemble",
    "create_ensemble_alias",
]
