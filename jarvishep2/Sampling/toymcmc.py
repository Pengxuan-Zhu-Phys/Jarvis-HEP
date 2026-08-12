"""ToyMCMC — reference independent-chain MCMC profile.

This module is intentionally thin. All high-throughput transport lives in
:class:`~jarvishep2.Sampling.mcmc_sampler.MCMCBaseSampler`:

* shared Redis task queue;
* per-chain feedback shards (``hep:feedback:chain:{id}``);
* async 1-inflight-per-chain event loop (no cross-chain generation barrier).

New independent MCMC methods should copy this pattern: set ``method``,
override ``_make_engine`` / Bounds hooks, and optionally attach progress
observers. Do not re-implement the Redis driver.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.Sampling.mcmc_sampler import MCMCBaseSampler
from jarvishep2.log_kv import PermilleProgress


class ToyMCMCSampler(MCMCBaseSampler):
    """Baseline random-walk MCMC with scalar or per-chain settings.

    Science: independent chains with plain MH and either a shared or
    per-chain ``proposal_scale``.
    Runtime: inherits the independent-chain async multi-chain design from
    ``MCMCBaseSampler`` (shared task queue + per-chain feedback shards).
    """

    method = "ToyMCMC"
    _checkpoint_excluded_attributes = frozenset({"_submit_progress"})

    def __init__(self) -> None:
        super().__init__()
        self._submit_progress: PermilleProgress | None = None

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        # Validate scalar/list cardinality here; scalar values and one-item
        # lists broadcast, while a full list maps one scale to each chain.
        self._normalize_scales()
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
        return f"chains={self._nchains} iterations={max_iteration}/{self._niters}"

    def _ensure_progress(self) -> None:
        total = max(1, int(self._nchains) * int(self._niters))
        if self._submit_progress is not None and self._submit_progress.total == total:
            return
        self._submit_progress = PermilleProgress(
            self._logger,
            total=total,
            label="ToyMCMC transitions completed",
        )
        self._logger.warning(
            "Initializing the ToyMCMC Sampling",
        )
        self._logger.info(
            "ToyMCMC Sampler configured: chains=%d iterations=%d "
            "total_transitions=%d proposal_scale=%s",
            self._nchains,
            self._niters,
            total,
            self._proposal_scales,
        )
        self._submit_progress.update(
            self._progress_done(),
            extra=self._progress_extra(),
            force=True,
        )

    def _emit_progress(self) -> None:
        """Emit one INFO line per ‰ and WARNING at every whole percent."""
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


def create_toymcmc() -> ToyMCMCSampler:
    return ToyMCMCSampler()


__all__ = ["ToyMCMCSampler", "create_toymcmc"]
