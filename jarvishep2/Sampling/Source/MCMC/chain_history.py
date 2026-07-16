#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Sequence


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
    """Append-only chain history with O(1) append and fast tail slicing."""

    def __init__(self) -> None:
        self._events: List[ChainEvent] = []

    def append(self, event: ChainEvent) -> None:
        self._events.append(event)

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
        self._events.append(
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
        return tuple(self._events[-n:])

    def last(self) -> ChainEvent | None:
        if not self._events:
            return None
        return self._events[-1]

    def __len__(self) -> int:
        return len(self._events)

