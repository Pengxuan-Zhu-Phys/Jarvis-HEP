#!/usr/bin/env python3
"""Chain registry / runtime records for V2 MCMC (ported from V1)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from jarvishep2.Sampling.Source.MCMC.chain_history import ChainHistory


@dataclass
class ChainRuntime:
    chain_id: int
    engine: Any
    temperature: float = 1.0
    is_cold: bool = False
    iter: int = 0
    accepted: int = 0
    rejected: int = 0
    window_iter: int = 0
    last_logl: float | None = None
    history: ChainHistory = field(default_factory=ChainHistory)
    meta: dict[str, Any] = field(default_factory=dict)
    # Open delayed-rejection stage (None when idle between iterations).
    open_stage: int | None = None

    def snapshot(self) -> dict[str, Any]:
        return {
            "chain_id": int(self.chain_id),
            "temperature": float(self.temperature),
            "is_cold": bool(self.is_cold),
            "iter": int(self.iter),
            "accepted": int(self.accepted),
            "rejected": int(self.rejected),
            "window_iter": int(self.window_iter),
            "last_logl": self.last_logl,
            "history_size": len(self.history),
            "proposal_scale": getattr(self.engine, "proposal_scale", None),
            "open_stage": self.open_stage,
            "meta": dict(self.meta),
        }


class ChainRegistry:
    """O(1) chain lookup."""

    def __init__(self, chains: list[ChainRuntime] | None = None) -> None:
        self._chains: dict[int, ChainRuntime] = {}
        if chains is not None:
            for chain in chains:
                self.add(chain)

    def add(self, chain: ChainRuntime) -> None:
        self._chains[int(chain.chain_id)] = chain

    def get(self, chain_id: int) -> ChainRuntime:
        return self._chains[int(chain_id)]

    def ids(self) -> list[int]:
        return sorted(self._chains.keys())

    def all(self) -> list[ChainRuntime]:
        return [self._chains[cid] for cid in self.ids()]

    def __len__(self) -> int:
        return len(self._chains)


__all__ = ["ChainRegistry", "ChainRuntime"]
