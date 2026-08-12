"""Parallel-tempering ensemble MCMC sampler."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.Source.MCMC.engine_ensemble import EnsembleChain
from jarvishep2.Sampling.ensemble_mcmc import EnsembleMCMCBase
from jarvishep2.Sampling.ptmcmc import PTMCMCBase


class PTEnsembleSampler(PTMCMCBase, EnsembleMCMCBase):
    method = "PTEnsemble"

    def __init__(self) -> None:
        # PTMCMCBase is first in the MRO, so initialize the ensemble-only
        # setting explicitly instead of relying on EnsembleMCMCBase.__init__.
        super().__init__()
        self._stretch_a = 2.0

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._configure_pt(bounds)
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


def create_pt_ensemble() -> PTEnsembleSampler:
    return PTEnsembleSampler()


__all__ = ["PTEnsembleSampler", "create_pt_ensemble"]
