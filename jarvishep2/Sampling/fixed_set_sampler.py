#!/usr/bin/env python3
"""Shared base for finite, stateless batch samplers (D9.7)."""

from __future__ import annotations

import time
from typing import Any, Mapping

import numpy as np

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.Sampling.variables import Variable, load_variables
from jarvishep2.log_kv import format_duration
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample

# Poll interval while waiting for in-flight slots (backpressure).
_INFLIGHT_POLL_SEC = 0.05
# How often to log a backpressure wait line (seconds).
_BACKPRESSURE_LOG_INTERVAL_SEC = 2.0


class FixedSetSampler(CheckpointedSampler):
    """Template base: index/accepted maps + batch flush + common set_config fields.

    Subclasses implement generation via ``propose_next`` and any sampler-specific
    initialize/export state. Shared batch submission lives here so modules no
    longer poke ``_submit`` private members from the outside.
    """

    uuid_prefix: str = "fixed"

    def __init__(self) -> None:
        super().__init__()
        self.vars: list[Variable] = []
        self._index = 0
        self._accepted_index = 0
        self._selectionexp: str | None = None
        self._seed = 0
        self._batch_size = 16
        # High-water mark for Redis in-flight tasks (submitted - finished).
        self._max_inflight = 16
        self._uuid_by_accepted_index: dict[int, str] = {}
        self._u_by_uuid: dict[str, np.ndarray] = {}
        self._logger = get_jarvis_logger("sampler.fixed")
        self._backpressure_last_log = 0.0

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        self.vars = load_variables(self.config)
        self._selectionexp = sampling.get("selection")
        self._seed = int(sampling.get("Seed", sampling.get("seed", 0)) or 0)
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        # Default backpressure: keep the pipeline full but bounded.
        # Subclasses (e.g. Bridson) may override with Sampling.MaxWorker.
        explicit = runtime.get("max_inflight", sampling.get("MaxWorker"))
        if explicit is not None:
            self._max_inflight = max(1, int(explicit or 1))
        else:
            self._max_inflight = max(self._batch_size, workers * 4, 1)

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        if accepted_index in self._uuid_by_accepted_index:
            return self._uuid_by_accepted_index[accepted_index]
        uuid = deterministic_sampler_uuid(
            prefix=self.uuid_prefix,
            seed=self._seed,
            sample_index=accepted_index,
        )
        self._uuid_by_accepted_index[accepted_index] = uuid
        return uuid

    def flush_batch(self, batch: list[Sample]) -> int:
        if not batch:
            return 0
        if len(batch) == 1:
            self._submit(batch[0])
        else:
            self._submit_group(batch)
        self._submitted_uuids.extend(sample.uuid for sample in batch)
        return len(batch)

    def _finished_count(self) -> int:
        """Completed + failed samples observed in Redis (worker results)."""
        if self.redis is None:
            return 0
        try:
            stats = self.redis.fetch_sample_stats()
        except Exception:
            return 0
        return int(stats.get("completed", 0) or 0) + int(stats.get("failed", 0) or 0)

    def _queue_depth(self) -> int:
        if self.redis is None:
            return 0
        try:
            lengths = self.redis.get_queue_lengths()
        except Exception:
            return 0
        return int(lengths.get("task_queue_length", 0) or 0)

    def _inflight(self, pushed: int, *, base_finished: int) -> int:
        finished = max(0, self._finished_count() - base_finished)
        return max(0, int(pushed) - finished)

    def _wait_for_inflight_room(
        self,
        pushed: int,
        *,
        base_finished: int,
        max_inflight: int,
    ) -> None:
        """Block until fewer than ``max_inflight`` samples are outstanding."""
        while self._inflight(pushed, base_finished=base_finished) >= max_inflight:
            now = time.time()
            if now - self._backpressure_last_log >= _BACKPRESSURE_LOG_INTERVAL_SEC:
                self._backpressure_last_log = now
                finished = max(0, self._finished_count() - base_finished)
                inflight = self._inflight(pushed, base_finished=base_finished)
                self._logger.info(
                    "backpressure wait: inflight=%d/%d submitted=%d finished=%d queued=%d",
                    inflight,
                    max_inflight,
                    pushed,
                    finished,
                    self._queue_depth(),
                )
            time.sleep(_INFLIGHT_POLL_SEC)

    def run_distributed(self) -> int:
        """Propose → Redis with in-flight backpressure (V1-like pipeline bound)."""
        pushed = 0
        batch: list[Sample] = []
        max_inflight = max(1, int(self._max_inflight or self._batch_size or 1))
        base_finished = self._finished_count()
        self._backpressure_last_log = 0.0
        self._logger.info(
            "distributed submit start: max_inflight=%d batch_size=%d",
            max_inflight,
            self._batch_size,
        )
        t0 = time.time()

        while True:
            self._wait_for_inflight_room(
                pushed, base_finished=base_finished, max_inflight=max_inflight
            )
            room = max_inflight - self._inflight(pushed, base_finished=base_finished)
            if room <= 0:
                continue
            target = min(self._batch_size, room)

            while len(batch) < target:
                sample = self.propose_next()
                if sample is None:
                    if batch:
                        pushed += self.flush_batch(batch)
                        batch = []
                    self._logger.info(
                        "distributed submit done: pushed=%d in %s (max_inflight=%d)",
                        pushed,
                        format_duration(time.time() - t0),
                        max_inflight,
                    )
                    return pushed
                batch.append(sample)

            pushed += self.flush_batch(batch)
            batch = []

    def propose_next(self) -> Sample | None:
        raise NotImplementedError

    def _stream_exhausted(self) -> bool:
        """Return True when the sampler has no more candidates to propose."""
        raise NotImplementedError

    def at_safe_barrier(self) -> bool:
        if not self._stream_exhausted():
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def _common_export_fields(self) -> dict[str, Any]:
        return {
            "index": int(self._index),
            "accepted_index": int(self._accepted_index),
            "selectionexp": self._selectionexp,
            "seed": int(self._seed),
            "uuid_by_accepted_index": dict(self._uuid_by_accepted_index),
            "u_by_uuid": {key: value.tolist() for key, value in self._u_by_uuid.items()},
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "chains": [],
            "ready_queue": [],
            "control_state": self._checkpoint_control_state(),
        }

    def _import_common_fields(self, state: Mapping[str, Any]) -> None:
        self._index = int(state.get("index", self._index) or 0)
        self._accepted_index = int(state.get("accepted_index", self._accepted_index) or 0)
        self._selectionexp = state.get("selectionexp", self._selectionexp)
        self._seed = int(state.get("seed", self._seed) or self._seed)
        raw_uuid_map = state.get("uuid_by_accepted_index") or {}
        self._uuid_by_accepted_index = {int(k): str(v) for k, v in raw_uuid_map.items()}
        raw_u_by_uuid = state.get("u_by_uuid") or {}
        self._u_by_uuid = {
            str(key): np.asarray(value, dtype=np.float64) for key, value in raw_u_by_uuid.items()
        }
        self._import_checkpoint_control_state(state)


__all__ = ["FixedSetSampler"]
