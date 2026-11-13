"""Feedback wait budgets exclude only measured archive admission delay."""

from __future__ import annotations

import time
from collections.abc import Collection


class DispatchWait:
    """A host-side wait window shared by completion, feedback and nested waits."""

    def __init__(self, redis, timeout: float, pending: Collection[str]):
        self.redis = redis
        self.timeout = timeout
        self.restart(pending)

    def _credits(self, pending: Collection[str]) -> dict[str, float]:
        if getattr(self.redis, "_dispatch_has_paused", False) is True:
            return self.redis.dispatch_wait_credits(tuple(pending))
        return {}

    def restart(self, pending: Collection[str]) -> None:
        self.started = time.monotonic()
        self._baseline = self._credits(pending)

    def remaining(self, pending: Collection[str]) -> float:
        remaining = self.timeout - (time.monotonic() - self.started)
        if remaining > 0 or not pending:
            return remaining
        credits = self._credits(pending)
        # A running task's credit is frozen at claim. A queued peer therefore
        # cannot extend the running task's feedback budget.
        extra = min(max(0.0, credits.get(uuid, 0.0) - self._baseline.get(uuid, 0.0))
                    for uuid in pending)
        return remaining + extra
