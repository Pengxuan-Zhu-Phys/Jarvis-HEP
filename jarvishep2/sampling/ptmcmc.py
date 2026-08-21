"""PTMCMC — parallel-tempering random-walk Metropolis (ToyMCMC-style profile).

Open/closed layering:

* Transport / MH bookkeeping live in
  :class:`~jarvishep2.sampling.mcmc_sampler.MCMCBaseSampler`
  (generation barrier schedule + control-side temperature exchange).
* This module only adds PT science knobs: temperature ladder, exchange
  interval, adjacent-swap pairing, and progress observers.

PT is **not** independent-chain async: replicas must hit a generation barrier
so exchange can run with consistent logL.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.sampling.Source.MCMC.config_contract import bounds_get, bounds_get_int
from jarvishep2.sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.log_kv import PermilleProgress


def normalize_temperature_ladder(
    ladder: Any,
    *,
    nchains: int,
    sampler_method: str,
) -> list[float]:
    """Validate / coerce a PT temperature ladder.

    Design contract (DESIGN_MCMC_ARCHITECTURE §8.5):

    * length equals ``num_chains``
    * all temperatures positive
    * first entry is the cold temperature ``1.0``
    * strictly increasing (hotter replicas at higher indices)
    """
    if ladder is None:
        raise ValueError(
            f"{sampler_method} requires Sampling.Bounds.temperature_ladder "
            f"(length {int(nchains)})"
        )
    if isinstance(ladder, (list, tuple)):
        values = [float(x) for x in ladder]
    else:
        raise ValueError(
            f"{sampler_method} temperature_ladder must be a list of positive numbers"
        )
    if len(values) != int(nchains):
        raise ValueError(
            f"{sampler_method} temperature_ladder length must equal num_chains "
            f"({int(nchains)}), got {len(values)}"
        )
    if any(t <= 0.0 for t in values):
        raise ValueError(
            f"{sampler_method} temperature_ladder values must all be positive"
        )
    if abs(float(values[0]) - 1.0) > 1e-12:
        raise ValueError(
            f"{sampler_method} temperature_ladder[0] must be the cold temperature 1.0, "
            f"got {values[0]!r}"
        )
    for lo, hi in zip(values, values[1:]):
        if not (float(hi) > float(lo)):
            raise ValueError(
                f"{sampler_method} temperature_ladder must be strictly increasing, "
                f"got {values!r}"
            )
    return values


class PTMCMCBase(MCMCBaseSampler):
    """Temperature ladder + exchange configuration shared by PT profiles."""

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
        if int(self._nchains) < 2:
            raise ValueError(f"{self.method} requires num_chains >= 2")
        self._exchange_interval = bounds_get_int(
            bounds,
            "exchange_interval",
            aliases=("exchange.interval",),
            default=1,
            minimum=1,
        )
        raw_ladder = bounds_get(
            bounds,
            "temperature_ladder",
            aliases=("temperature.ladder",),
            default=None,
        )
        if raw_ladder is None:
            # Sensible default when omitted: cold T=1, then 2, 3, …
            raw_ladder = [1.0 + float(ii) for ii in range(self._nchains)]
        self._temperature_ladder = normalize_temperature_ladder(
            raw_ladder,
            nchains=self._nchains,
            sampler_method=self.method,
        )

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        self._configure_pt(bounds)

    def _uses_pt(self) -> bool:
        return True

    def _pair_chain_ids(self) -> list[tuple[int, int]]:
        """Adjacent ladder swaps; alternate even/odd edges each exchange round.

        Chain ids are assigned in ladder order (cid → temperature_ladder[cid]),
        so adjacent ids are adjacent temperatures. Non-adjacent pairings (e.g.
        hottest↔coldest) are never proposed.
        """
        ids = self._ensure_registry().ids()
        if len(ids) < 2:
            return []
        start = int(self._exchange_offset) % 2
        pairs: list[tuple[int, int]] = []
        for ii in range(start, len(ids) - 1, 2):
            pairs.append((ids[ii], ids[ii + 1]))
        self._exchange_offset = (int(self._exchange_offset) + 1) % 2
        return pairs

    def _export_method_state(self) -> dict[str, Any]:
        return {
            "swap_attempts": int(self._swap_attempts),
            "swap_accepts": int(self._swap_accepts),
            "exchange_offset": int(self._exchange_offset),
            "exchange_interval": int(self._exchange_interval),
            "temperature_ladder": list(self._temperature_ladder),
        }

    def _import_method_state(self, state: Mapping[str, Any]) -> None:
        self._swap_attempts = int(state.get("swap_attempts", 0) or 0)
        self._swap_accepts = int(state.get("swap_accepts", 0) or 0)
        self._exchange_offset = int(state.get("exchange_offset", 0) or 0)
        if state.get("exchange_interval") is not None:
            self._exchange_interval = max(1, int(state.get("exchange_interval") or 1))
        ladder = state.get("temperature_ladder")
        if ladder is not None:
            # Resume may carry the ladder; live Bounds still win at set_config time.
            # Only restore when lengths match the configured chain count.
            try:
                restored = [float(x) for x in ladder]
            except (TypeError, ValueError):
                restored = []
            if len(restored) == int(self._nchains):
                self._temperature_ladder = restored


class PTMCMCSampler(PTMCMCBase):
    """Parallel-tempering random-walk MH with scalar or per-replica scales.

    Science: Metropolis on each temperature replica (β = 1/T), adjacent
    temperature swaps at generation barriers every ``exchange_interval`` draws.
    Runtime: inherits barrier-coupled MCMC transport (not async independent).
    """

    method = "PTMCMC"
    _checkpoint_excluded_attributes = frozenset(
        set(PTMCMCBase._checkpoint_excluded_attributes) | {"_submit_progress"}
    )

    def __init__(self) -> None:
        super().__init__()
        self._submit_progress: PermilleProgress | None = None

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        # Validate scalar/list scale cardinality here; scalar values and one-item
        # lists broadcast, while a full list maps one scale to each replica.
        self._normalize_scales()
        self._configure_pt(bounds)
        self._submit_progress = None

    def _progress_done(self) -> int:
        registry = self._ensure_registry()
        return sum(max(0, int(chain.engine.iterations)) for chain in registry.all())

    def _progress_extra(self) -> str:
        registry = self._ensure_registry()
        max_iteration = max(
            (int(chain.engine.iterations) for chain in registry.all()),
            default=0,
        )
        swaps = int(self._swap_accepts)
        attempts = int(self._swap_attempts)
        return (
            f"replicas={self._nchains} iterations={max_iteration}/{self._niters} "
            f"swaps={swaps}/{attempts}"
        )

    def _ensure_progress(self) -> None:
        total = max(1, int(self._nchains) * int(self._niters))
        if self._submit_progress is not None and self._submit_progress.total == total:
            return
        self._submit_progress = PermilleProgress(
            self._logger,
            total=total,
            label="PTMCMC transitions completed",
        )
        self._logger.warning("Initializing the PTMCMC Sampling")
        self._logger.info(
            "PTMCMC Sampler configured: replicas=%d iterations=%d "
            "total_transitions=%d proposal_scale=%s ladder=%s exchange_interval=%d",
            self._nchains,
            self._niters,
            total,
            self._proposal_scales,
            self._temperature_ladder,
            self._exchange_interval,
        )
        self._submit_progress.update(
            self._progress_done(),
            extra=self._progress_extra(),
            force=True,
        )

    def _emit_progress(self) -> None:
        self._ensure_progress()
        assert self._submit_progress is not None
        self._submit_progress.update(
            self._progress_done(),
            extra=self._progress_extra(),
        )

    def _on_run_started(self) -> None:
        self._ensure_progress()

    def _on_generation_completed(self) -> None:
        self._emit_progress()

    def _on_runtime_state_imported(self) -> None:
        self._submit_progress = None


def create_ptmcmc() -> PTMCMCSampler:
    return PTMCMCSampler()


__all__ = [
    "PTMCMCBase",
    "PTMCMCSampler",
    "create_ptmcmc",
    "normalize_temperature_ladder",
]
