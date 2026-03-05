#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Sequence


@dataclass
class MCMCMetricsFrame:
    step: int
    event: str
    state: str
    pending_futures: int
    ready_queue: int
    acceptance_rate_mean: float
    ess_proxy_mean: float | None = None
    autocorr_lag1_proxy_mean: float | None = None
    swap_acceptance_rate: float | None = None
    queue_depth_hint: int = 0
    meta: Dict[str, Any] = field(default_factory=dict)
    timestamp_utc: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())


class MCMCMetricsBus:
    """Lightweight in-process metrics bus for controller/analysis consumption."""

    def __init__(self, max_frames: int = 4096) -> None:
        self.max_frames = int(max_frames)
        self._frames: List[MCMCMetricsFrame] = []
        self._dropped = 0

    def publish(self, frame: MCMCMetricsFrame) -> None:
        if len(self._frames) >= self.max_frames:
            self._frames.pop(0)
            self._dropped += 1
        self._frames.append(frame)

    def all(self) -> Sequence[MCMCMetricsFrame]:
        return tuple(self._frames)

    def tail(self, n: int) -> Sequence[MCMCMetricsFrame]:
        n = int(n)
        if n <= 0:
            return tuple()
        if n >= len(self._frames):
            return tuple(self._frames)
        return tuple(self._frames[-n:])

    def latest(self) -> MCMCMetricsFrame | None:
        if not self._frames:
            return None
        return self._frames[-1]

    @property
    def dropped(self) -> int:
        return int(self._dropped)
