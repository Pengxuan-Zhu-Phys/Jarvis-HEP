"""Parallel-tempering random-walk MCMC sampler."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.Sampling.Source.MCMC.config_contract import bounds_get, bounds_get_int


class PTMCMCBase(MCMCBaseSampler):
    """Common temperature ladder, exchange, and checkpoint state."""

    _checkpoint_excluded_attributes = frozenset(
        {
            "_exchange_interval",
            "_exchange_offset",
            "_swap_attempts",
            "_swap_accepts",
            "_exchange_rng",
            "_temperature_ladder",
        }
    )

    def __init__(self) -> None:
        super().__init__()
        self._exchange_interval = 1
        self._exchange_offset = 0
        self._swap_attempts = 0
        self._swap_accepts = 0
        self._exchange_rng = None
        self._temperature_ladder: list[float] = [1.0]

    def _configure_pt(self, bounds: Mapping[str, Any]) -> None:
        self._exchange_interval = bounds_get_int(
            bounds,
            "exchange_interval",
            aliases=("exchange.interval",),
            default=1,
            minimum=1,
        )
        ladder = bounds_get(
            bounds,
            "temperature_ladder",
            aliases=("temperature.ladder",),
            default=None,
        )
        if ladder is None:
            ladder = [1.0 + float(ii) for ii in range(self._nchains)]
        self._temperature_ladder = [float(x) for x in ladder]
        if len(self._temperature_ladder) != self._nchains:
            raise ValueError(
                f"temperature_ladder size mismatch for {self.method}: "
                f"expect {self._nchains}, got {len(self._temperature_ladder)}"
            )

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._configure_pt(bounds)

    def _uses_pt(self) -> bool:
        return True

    def _export_method_state(self) -> dict[str, Any]:
        return {
            "swap_attempts": int(self._swap_attempts),
            "swap_accepts": int(self._swap_accepts),
            "exchange_offset": int(self._exchange_offset),
            "temperature_ladder": list(self._temperature_ladder),
        }

    def _import_method_state(self, state: Mapping[str, Any]) -> None:
        self._swap_attempts = int(state.get("swap_attempts", 0) or 0)
        self._swap_accepts = int(state.get("swap_accepts", 0) or 0)
        self._exchange_offset = int(state.get("exchange_offset", 0) or 0)
        ladder = state.get("temperature_ladder")
        if ladder is not None:
            self._temperature_ladder = [float(x) for x in ladder]


class PTMCMCSampler(PTMCMCBase):
    method = "PTMCMC"


class PTSampler(PTMCMCSampler):
    method = "PT"


def create_ptmcmc() -> PTMCMCSampler:
    return PTMCMCSampler()


def create_pt() -> PTSampler:
    return PTSampler()


__all__ = ["PTMCMCBase", "PTMCMCSampler", "PTSampler", "create_pt", "create_ptmcmc"]
