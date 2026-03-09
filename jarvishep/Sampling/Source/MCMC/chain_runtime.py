#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List

from .chain_history import ChainHistory
from .engine_contract import MCMCEngineProtocol


@dataclass
class ChainRuntime:
    chain_id: int
    engine: MCMCEngineProtocol | Any
    temperature: float = 1.0
    is_cold: bool = False
    iter: int = 0
    accepted: int = 0
    rejected: int = 0
    window_iter: int = 0
    last_logl: float | None = None
    history: ChainHistory = field(default_factory=ChainHistory)
    meta: Dict[str, Any] = field(default_factory=dict)

    def snapshot(self) -> Dict[str, Any]:
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
        }


class ChainRegistry:
    """O(1) chain lookup + O(1) cold-chain checks."""

    def __init__(self, chains: Iterable[ChainRuntime] | None = None) -> None:
        self._chains: Dict[int, ChainRuntime] = {}
        self._cold_chain_ids: set[int] = set()
        if chains is not None:
            for chain in chains:
                self.add(chain)

    def add(self, chain: ChainRuntime) -> None:
        cid = int(chain.chain_id)
        self._chains[cid] = chain
        if bool(chain.is_cold):
            self._cold_chain_ids.add(cid)
        else:
            self._cold_chain_ids.discard(cid)

    def get(self, chain_id: int) -> ChainRuntime:
        return self._chains[int(chain_id)]

    def ids(self) -> List[int]:
        return sorted(self._chains.keys())

    def all(self) -> List[ChainRuntime]:
        return [self._chains[cid] for cid in self.ids()]

    def is_cold(self, chain_id: int) -> bool:
        return int(chain_id) in self._cold_chain_ids

    def mark_cold(self, chain_id: int, flag: bool = True) -> None:
        cid = int(chain_id)
        chain = self._chains[cid]
        chain.is_cold = bool(flag)
        if chain.is_cold:
            self._cold_chain_ids.add(cid)
        else:
            self._cold_chain_ids.discard(cid)

    def set_temperature(self, chain_id: int, temperature: float) -> None:
        cid = int(chain_id)
        self._chains[cid].temperature = float(temperature)

    def cold_chain_ids(self) -> List[int]:
        return sorted(self._cold_chain_ids)

    def __len__(self) -> int:
        return len(self._chains)
