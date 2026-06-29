#!/usr/bin/env python3
"""Checkpoint/resume mixin for distributed stateless samplers."""

from __future__ import annotations

import threading
from typing import Any, Mapping

from jarvishep2.Sampling.runtime_checkpoint import (
    CHECKPOINT_HEARTBEAT_SEC,
    CheckpointHeartbeat,
)
from jarvishep2.Sampling.sampler import SamplingVirtial


class CheckpointedSampler(SamplingVirtial):
    """Adds WP-D6.2 checkpoint heartbeat + resume repropose to samplers."""

    def __init__(self) -> None:
        super().__init__()
        self._runtime_checkpoint_save_lock = threading.RLock()
        self._submitted_uuids: list[str] = []
        self._completed_uuids: set[str] = set()
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
        if self._checkpoint_heartbeat is not None:
            self._checkpoint_heartbeat.stop()
        self._checkpoint_heartbeat = CheckpointHeartbeat(
            interval_sec=CHECKPOINT_HEARTBEAT_SEC,
            save_callback=lambda reason="checkpoint_heartbeat": self.persist_runtime_checkpoint(
                force=True,
                reason=reason,
            ),
        )
        self._checkpoint_heartbeat.start()

    def shutdown_checkpointing(self) -> None:
        if self._checkpoint_heartbeat is not None:
            self._checkpoint_heartbeat.stop()
            self._checkpoint_heartbeat = None

    @property
    def submitted_uuids(self) -> frozenset[str]:
        return frozenset(self._submitted_uuids)

    def mark_completed(self, uuid: str) -> None:
        self._completed_uuids.add(str(uuid))

    def set_resume_repropose_hint(self, enabled: bool = True) -> None:
        self._repropose_after_resume = bool(enabled)

    def _checkpoint_control_state(self) -> dict[str, Any]:
        return {"repropose_after_resume": bool(self._repropose_after_resume)}

    def _import_checkpoint_control_state(self, state: Mapping[str, Any]) -> None:
        control = dict(state.get("control_state") or {})
        self._repropose_after_resume = bool(control.get("repropose_after_resume", False))
        self._submitted_uuids = [str(item) for item in state.get("submitted_uuids") or []]
        self._completed_uuids = {str(item) for item in state.get("completed_uuids") or []}

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

    def at_safe_barrier(self) -> bool:
        raise NotImplementedError


__all__ = ["CheckpointedSampler"]