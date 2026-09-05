"""ToyMCMC — reference independent-chain MCMC profile.

This module is intentionally thin. All high-throughput transport lives in
:class:`~jarvishep2.sampling.mcmc_sampler.MCMCBaseSampler`:

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

from jarvishep2.sampling.mcmc_sampler import MCMCBaseSampler


class ToyMCMCSampler(MCMCBaseSampler):
    """Baseline random-walk MCMC with scalar or per-chain settings.

    Science: independent chains with plain MH and either a shared or
    per-chain ``proposal_scale``.
    Runtime: inherits the independent-chain async multi-chain design from
    ``MCMCBaseSampler`` (shared task queue + per-chain feedback shards).
    Progress: Jarvis ``PermilleProgress`` from the MCMC base logger.
    """

    method = "ToyMCMC"

    def _configure_method(self, bounds: Mapping[str, Any]) -> None:
        # Validate scalar/list cardinality here; scalar values and one-item
        # lists broadcast, while a full list maps one scale to each chain.
        self._normalize_scales()


def create_toymcmc() -> ToyMCMCSampler:
    return ToyMCMCSampler()


__all__ = ["ToyMCMCSampler", "create_toymcmc"]
