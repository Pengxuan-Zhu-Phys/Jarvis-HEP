#!/usr/bin/env python3
"""Archive receipts: which submitted samples are durable in DATABASE."""

from __future__ import annotations

import threading
from collections.abc import Iterable
from typing import Any

# One SISMEMBER pipeline per chunk keeps each Redis round trip short.
_CONFIRM_CHUNK = 10_000
# ``add`` prunes confirmed receipts once the pending set reaches this size and
# again each time it doubles, so memory stays bounded without checkpoints.
_PRUNE_MIN = 20_000


class ArchiveReceipts:
    """Track submitted samples until the Archiver posts their receipts.

    The Archiver SADDs each UUID to ``hep:archived:{scan}`` after its HDF5
    batch is fsynced.  :meth:`confirm` checks only the pending UUIDs, removes
    the receipts it consumed, and returns the lowest sample index still
    pending: every submitted index below it is in DATABASE.
    """

    def __init__(self, redis: Any, scan_name: str) -> None:
        self._redis = redis
        self._scan = str(scan_name or "scan")
        self._pending: dict[int, str] = {}
        self._next_index: int | None = None
        self._prune_at = _PRUNE_MIN
        self._lock = threading.Lock()

    def add(self, samples: Iterable[Any]) -> None:
        """Record samples just before their tasks are pushed."""
        with self._lock:
            for sample in samples:
                index = getattr(sample, "sample_index", None)
                uuid = str(getattr(sample, "uuid", "") or "").strip()
                if index is None or not uuid:
                    continue
                index = int(index)
                self._pending[index] = uuid
                if self._next_index is None or index >= self._next_index:
                    self._next_index = index + 1
            prune = len(self._pending) >= self._prune_at
        if prune:
            self.confirm()
            with self._lock:
                self._prune_at = max(_PRUNE_MIN, 2 * len(self._pending))

    def pending_count(self) -> int:
        with self._lock:
            return len(self._pending)

    def confirm(self) -> int | None:
        """Return the durable index cursor, or None before the first ``add``."""
        with self._lock:
            if self._next_index is None:
                return None
            pending = list(self._pending.items())
        confirmed: list[tuple[int, str]] = []
        for start in range(0, len(pending), _CONFIRM_CHUNK):
            chunk = pending[start : start + _CONFIRM_CHUNK]
            missing = set(
                self._redis.missing_archived_uuids(self._scan, [uuid for _, uuid in chunk])
            )
            confirmed.extend(item for item in chunk if item[1] not in missing)
        if confirmed:
            self._redis.remove_archived_uuids(self._scan, [uuid for _, uuid in confirmed])
        with self._lock:
            for index, uuid in confirmed:
                if self._pending.get(index) == uuid:
                    del self._pending[index]
            return min(self._pending) if self._pending else self._next_index


__all__ = ["ArchiveReceipts"]
