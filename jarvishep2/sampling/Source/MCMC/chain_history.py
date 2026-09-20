#!/usr/bin/env python3
from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Deque, Dict, Sequence

# RAM tail only — full per-iteration rows stream to DATABASE/chain_history.csv.
# A million-iteration ensemble cannot keep every ChainEvent in the sampler
# process or in state.pkl (D21.10: checkpoints must stay O(1) in scan size).
CHAIN_HISTORY_RAM_TAIL = 256


@dataclass
class ChainEvent:
    iter: int
    state: str
    proposal: Any
    logl: float | None
    accepted: bool
    temperature: float
    timestamp_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    meta: Dict[str, Any] = field(default_factory=dict)


class ChainHistory:
    """Bounded in-memory tail with O(1) append.

    ``all()`` / ``__len__`` report the RAM tail, not the lifetime total.
    ``total_appended`` is the lifetime count (survives pickle).
    """

    def __init__(self, maxlen: int | None = CHAIN_HISTORY_RAM_TAIL) -> None:
        self._maxlen = None if maxlen is None else max(1, int(maxlen))
        self._events: Deque[ChainEvent] = deque(maxlen=self._maxlen)
        self._total = 0

    def append(self, event: ChainEvent) -> None:
        self._events.append(event)
        self._total += 1

    def append_from_values(
        self,
        *,
        iter: int,
        state: str,
        proposal: Any,
        logl: float | None,
        accepted: bool,
        temperature: float,
        meta: Dict[str, Any] | None = None,
    ) -> None:
        self.append(
            ChainEvent(
                iter=int(iter),
                state=str(state),
                proposal=proposal,
                logl=None if logl is None else float(logl),
                accepted=bool(accepted),
                temperature=float(temperature),
                meta=dict(meta or {}),
            )
        )

    def all(self) -> Sequence[ChainEvent]:
        return tuple(self._events)

    def tail(self, n: int) -> Sequence[ChainEvent]:
        n = int(n)
        if n <= 0:
            return tuple()
        if n >= len(self._events):
            return tuple(self._events)
        return tuple(list(self._events)[-n:])

    def last(self) -> ChainEvent | None:
        if not self._events:
            return None
        return self._events[-1]

    @property
    def total_appended(self) -> int:
        return int(self._total)

    def __len__(self) -> int:
        return len(self._events)

    def __getstate__(self) -> dict[str, Any]:
        return {
            "_maxlen": self._maxlen,
            "_events": list(self._events),
            "_total": int(self._total),
        }

    def __setstate__(self, state: dict[str, Any]) -> None:
        events = list(state.get("_events") or [])
        maxlen = state.get("_maxlen", CHAIN_HISTORY_RAM_TAIL)
        if maxlen is not None:
            maxlen = max(1, int(maxlen))
        self._maxlen = maxlen
        self._events = deque(events, maxlen=maxlen)
        total = state.get("_total")
        self._total = int(total) if total is not None else len(self._events)
