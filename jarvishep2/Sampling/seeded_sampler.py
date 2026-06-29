#!/usr/bin/env python3
"""SeedSequence-driven sampler for distributed reproducibility (WP-D6.2)."""

from __future__ import annotations

import hashlib
import threading
from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.runtime_checkpoint import (
    CHECKPOINT_HEARTBEAT_SEC,
    CheckpointHeartbeat,
    derive_sample_seed,
    deserialize_seed_sequence,
    serialize_seed_sequence,
)
from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.sample import Sample


def deterministic_uuid(*, master: np.random.SeedSequence, sample_index: int) -> str:
    digest = hashlib.sha256(
        f"{int(master.entropy)}:{int(sample_index)}".encode("utf-8")
    ).hexdigest()
    return (
        f"{digest[0:8]}-{digest[8:12]}-{digest[12:16]}-{digest[16:20]}-{digest[20:32]}"
    )


class SeededOperaSampler(SamplingVirtial):
    """Deterministic opera sampler for checkpoint/resume acceptance tests."""

    def __init__(
        self,
        *,
        seed: int = 0,
        total_points: int = 10,
        dimensions: int = 3,
    ) -> None:
        super().__init__()
        self._runtime_checkpoint_save_lock = threading.RLock()
        self._master_seq = np.random.SeedSequence(int(seed))
        self._next_sample_index = 0
        self._total_points = max(0, int(total_points))
        self._dimensions = max(1, int(dimensions))
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

    def derive_u_coords(self, sample_index: int) -> np.ndarray:
        child = derive_sample_seed(self._master_seq, sample_index)
        rng = np.random.default_rng(child)
        return rng.uniform(0.0, 1.0, size=self._dimensions).astype(np.float64)

    def propose_next(self) -> Sample | None:
        if self._next_sample_index >= self._total_points:
            return None
        index = self._next_sample_index
        sample = self._build_sample(self.derive_u_coords(index))
        sample.uuid = deterministic_uuid(master=self._master_seq, sample_index=index)
        self._next_sample_index += 1
        return sample

    def propose_remaining(self) -> list[Sample]:
        samples: list[Sample] = []
        while True:
            sample = self.propose_next()
            if sample is None:
                break
            samples.append(sample)
        return samples

    def submit_next(self) -> str | None:
        sample = self.propose_next()
        if sample is None:
            return None
        self._submit(sample)
        self._submitted_uuids.append(sample.uuid)
        return sample.uuid

    def submit_all_remaining(self) -> list[str]:
        uuids: list[str] = []
        while True:
            submitted = self.submit_next()
            if submitted is None:
                break
            uuids.append(submitted)
        return uuids

    @property
    def submitted_uuids(self) -> frozenset[str]:
        return frozenset(self._submitted_uuids)

    def mark_completed(self, uuid: str) -> None:
        self._completed_uuids.add(str(uuid))

    def at_safe_barrier(self) -> bool:
        if self._next_sample_index < self._total_points:
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def export_runtime_state(self) -> dict[str, Any]:
        return {
            "seed_sequence": serialize_seed_sequence(self._master_seq),
            "next_sample_index": int(self._next_sample_index),
            "total_points": int(self._total_points),
            "dimensions": int(self._dimensions),
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "chains": [],
            "ready_queue": [],
            "control_state": {"repropose_after_resume": bool(self._repropose_after_resume)},
            "numpy_random_state": None,
        }

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        seed_payload = state.get("seed_sequence")
        if isinstance(seed_payload, Mapping):
            self._master_seq = deserialize_seed_sequence(seed_payload)
        self._next_sample_index = int(state.get("next_sample_index", 0) or 0)
        self._total_points = int(state.get("total_points", self._total_points) or self._total_points)
        self._dimensions = int(state.get("dimensions", self._dimensions) or self._dimensions)
        self._submitted_uuids = [str(item) for item in state.get("submitted_uuids") or []]
        self._completed_uuids = {str(item) for item in state.get("completed_uuids") or []}
        control = dict(state.get("control_state") or {})
        self._repropose_after_resume = bool(control.get("repropose_after_resume", False))

    def set_resume_repropose_hint(self, enabled: bool = True) -> None:
        self._repropose_after_resume = bool(enabled)

    def repropose_unfinished(self) -> list[str]:
        """Re-submit Samples not yet archived after resume (design §10)."""
        if not self._repropose_after_resume:
            return []
        pending = [
            uuid
            for uuid in self._submitted_uuids
            if uuid not in self._completed_uuids
        ]
        requeued: list[str] = []
        for uuid in pending:
            try:
                index = self._submitted_uuids.index(uuid)
            except ValueError:
                continue
            sample = self._build_sample(self.derive_u_coords(index))
            sample.uuid = uuid
            self._submit(sample)
            requeued.append(uuid)
        return requeued

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


__all__ = ["SeededOperaSampler", "deterministic_uuid"]