#!/usr/bin/env python3
"""Checkpoint/resume mixin for distributed stateless samplers."""

from __future__ import annotations

import threading
from collections import deque
from abc import ABC, abstractmethod
from copy import deepcopy
from typing import Any, Mapping

from jarvishep2.Sampling.runtime_checkpoint import (
    CHECKPOINT_HEARTBEAT_SEC,
    CheckpointHeartbeat,
)
from jarvishep2.Sampling.sampler import SamplingVirtial


class CheckpointedSampler(SamplingVirtial, ABC):
    """Adds lightweight checkpoint capture and safe resume support to samplers."""

    # Fixed/CSV enumerators override this: they use the DATABASE's contiguous
    # ``sample_index`` prefix rather than a full UUID acknowledgement set.
    uses_indexed_resume = False
    # Every instance attribute is explicitly saved by a subclass export or
    # deliberately excluded.  ``assert_checkpoint_attribute_contract`` makes
    # a newly added mutable field fail tests instead of silently turning a
    # checkpoint into an incomplete or O(N) snapshot (D21.10).
    _checkpoint_saved_attributes = frozenset()
    _checkpoint_excluded_attributes = frozenset(
        {
            "config", "info", "redis", "runtime_mode", "_execution_plan_template",
            "_opera_modules", "_calculator_modules", "_sample_artifacts",
            "_expression_context", "_runtime_checkpoint_save_lock", "_submitted_uuids",
            "_persisted_uuids", "_persisted_index_prefix", "_last_safe_state",
            "_checkpoint_snapshots", "_checkpoint_requested", "_checkpoint_file",
            "_checkpoint_heartbeat", "_save_checkpoint_callback", "_repropose_after_resume",
            # D22: rebuilt from task card on every process start (not checkpointed).
            "_mapper_pipeline",
        }
    )

    def __init__(self) -> None:
        super().__init__()
        self._runtime_checkpoint_save_lock = threading.RLock()
        # UUID acknowledgements remain for feedback samplers which genuinely
        # need an acknowledgement barrier.  Stateless samplers resume from an
        # index prefix instead; never grow this list on their hot path.
        self._submitted_uuids: list[str] = []
        self._persisted_uuids: set[str] = set()
        self._persisted_index_prefix = 0
        self._last_safe_state: dict[str, Any] | None = None
        self._checkpoint_snapshots: deque[tuple[int, dict[str, Any]]] = deque(maxlen=8)
        self._checkpoint_requested = False
        self._checkpoint_file = ""
        self._checkpoint_heartbeat: CheckpointHeartbeat | None = None
        self._save_checkpoint_callback = None
        self._repropose_after_resume = False

    def configure_checkpoint(
        self,
        *,
        checkpoint_file: str,
        save_callback,
    ) -> None:
        self._checkpoint_file = str(checkpoint_file)
        self._save_checkpoint_callback = save_callback
        # Initial state is always a valid recovery point: no work is in flight.
        initial = deepcopy(self.export_runtime_state())
        initial_index = self._state_index(initial)
        initial["checkpoint_index"] = initial_index
        self._last_safe_state = initial
        self._checkpoint_snapshots.clear()
        self._checkpoint_snapshots.append((initial_index, initial))
        if self._checkpoint_heartbeat is not None:
            self._checkpoint_heartbeat.stop()
        self._checkpoint_heartbeat = CheckpointHeartbeat(
            interval_sec=float(getattr(self, "_checkpoint_heartbeat_sec", CHECKPOINT_HEARTBEAT_SEC)),
            # The heartbeat must not snapshot sampler state, acquire the
            # runtime-checkpoint lock, or touch Redis.  It only requests a
            # capture at the next natural batch boundary (D21.9).
            save_callback=lambda reason="checkpoint_heartbeat": self.request_checkpoint_capture(),
        )
        self._checkpoint_heartbeat.start()

    def shutdown_checkpointing(self) -> None:
        if self._checkpoint_heartbeat is not None:
            self._checkpoint_heartbeat.stop()
            self._checkpoint_heartbeat = None

    @property
    def submitted_uuids(self) -> frozenset[str]:
        return frozenset(self._submitted_uuids)

    def set_persisted_uuids(self, uuids: set[str] | frozenset[str]) -> None:
        """Install DATABASE-derived completion truth for resume and dedup."""
        self._persisted_uuids = {str(item) for item in uuids if str(item)}

    def set_persisted_prefix(self, prefix: int) -> None:
        """Advance the DATABASE-proven contiguous sample-index prefix.

        A snapshot at index ``j`` is safe iff ``j <= prefix``.  Keeping a
        bounded snapshot window makes this check O(1) per archived batch and
        avoids keeping either UUID maps or every previous sampler state.
        """
        self._persisted_index_prefix = max(self._persisted_index_prefix, int(prefix or 0))
        for index, state in self._checkpoint_snapshots:
            if index <= self._persisted_index_prefix:
                self._last_safe_state = deepcopy(state)

    def assert_checkpoint_attribute_contract(self) -> None:
        """Raise when a sampler field lacks an explicit checkpoint decision."""
        saved: set[str] = set()
        excluded: set[str] = set()
        for klass in type(self).__mro__:
            saved.update(getattr(klass, "_checkpoint_saved_attributes", ()))
            excluded.update(getattr(klass, "_checkpoint_excluded_attributes", ()))
        overlap = saved & excluded
        if overlap:
            raise AssertionError(f"checkpoint attributes both saved and excluded: {sorted(overlap)}")
        unknown = set(vars(self)) - saved - excluded
        if unknown:
            raise AssertionError(
                f"checkpoint attributes need an explicit saved/excluded decision: {sorted(unknown)}"
            )

    def should_submit_uuid(self, uuid: str) -> bool:
        return str(uuid) not in self._persisted_uuids

    def record_submitted_batch(self, uuids: list[str]) -> None:
        self._submitted_uuids.extend(str(item) for item in uuids if str(item))

    def record_persisted_proposals(self) -> None:
        self._last_safe_state = deepcopy(self.export_runtime_state())

    @staticmethod
    def _state_index(state: Mapping[str, Any]) -> int:
        """Return the next logical sample index represented by ``state``."""
        raw = state.get("checkpoint_index", state.get("accepted_index", state.get("sample_index", 0)))
        try:
            return max(0, int(raw or 0))
        except (TypeError, ValueError):
            return 0

    def request_checkpoint_capture(self) -> None:
        """O(1) heartbeat-side request; no state is inspected here."""
        # A plain assignment is sufficient under CPython's GIL.  If the
        # batch-boundary reader races a heartbeat, at worst one request folds
        # into the next heartbeat; no recovery state is advanced unsafely.
        self._checkpoint_requested = True

    def capture_checkpoint_at_batch_boundary(self, *, reason: str = "checkpoint_heartbeat") -> bool:
        """Capture one small state only after a complete submission batch.

        This is called by enumerating samplers, never by the heartbeat thread.
        The callback writes the newest DATABASE-safe snapshot selected through
        ``set_persisted_prefix``.
        """
        if not self._checkpoint_requested:
            return False
        self._checkpoint_requested = False
        state = deepcopy(self.export_runtime_state())
        index = self._state_index(state)
        state["checkpoint_index"] = index
        self._checkpoint_snapshots.append((index, state))
        if index <= self._persisted_index_prefix:
            self._last_safe_state = deepcopy(state)
        if self._save_checkpoint_callback is not None:
            with self._runtime_checkpoint_save_lock:
                self._save_checkpoint_callback(reason=reason)
        return True

    def checkpoint_runtime_state(self, *, safe: bool) -> dict[str, Any]:
        """Return the newest state whose submitted work is fully durable."""
        if safe:
            self._last_safe_state = deepcopy(self.export_runtime_state())
        if self._last_safe_state is None:
            self._last_safe_state = deepcopy(self.export_runtime_state())
        return deepcopy(self._last_safe_state)

    def set_resume_repropose_hint(self, enabled: bool = True) -> None:
        self._repropose_after_resume = bool(enabled)

    def _checkpoint_control_state(self) -> dict[str, Any]:
        return {"repropose_after_resume": bool(self._repropose_after_resume)}

    def _import_checkpoint_control_state(self, state: Mapping[str, Any]) -> None:
        control = dict(state.get("control_state") or {})
        self._repropose_after_resume = bool(control.get("repropose_after_resume", False))
        self._submitted_uuids = [str(item) for item in state.get("submitted_uuids") or []]
        self._last_safe_state = deepcopy(dict(state))

    def persist_runtime_checkpoint(
        self,
        *,
        force: bool = False,
        reason: str = "",
        archiver_persistence: Mapping[str, Any] | None = None,
    ) -> bool:
        if self._save_checkpoint_callback is None:
            return False
        if not force:
            from jarvishep2.Sampling.runtime_checkpoint import safe_barrier_ready

            if not safe_barrier_ready(
                sampler_at_barrier=self.at_safe_barrier(),
                submitted_uuids=self.submitted_uuids,
                archiver_persistence=archiver_persistence,
            ):
                return False
        with self._runtime_checkpoint_save_lock:
            return bool(self._save_checkpoint_callback(reason=reason))

    @abstractmethod
    def at_safe_barrier(self) -> bool:
        """Return True when the sampler may safely write a resume checkpoint."""


__all__ = ["CheckpointedSampler"]
