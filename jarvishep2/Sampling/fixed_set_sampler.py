#!/usr/bin/env python3
"""Shared base for finite, stateless batch samplers (D9.7)."""

from __future__ import annotations

import time
from typing import Any, Mapping

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
    uses_indexed_resume = True
    _checkpoint_saved_attributes = frozenset(
        {"_index", "_accepted_index", "_selectionexp", "_seed"}
    )
    _checkpoint_excluded_attributes = frozenset(
        {
            "vars", "_batch_size", "_max_inflight", "_logger",
            "_backpressure_last_log", "_checkpoint_heartbeat_sec",
        }
    )

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
        self._checkpoint_heartbeat_sec = float(runtime.get("checkpoint_heartbeat_sec", 30.0) or 30.0)
        # Default backpressure: keep the pipeline full but bounded.
        # Subclasses (e.g. Bridson) may override with Sampling.MaxWorker.
        explicit = runtime.get("max_inflight", sampling.get("MaxWorker"))
        if explicit is not None:
            self._max_inflight = max(1, int(explicit or 1))
        else:
            self._max_inflight = max(self._batch_size, workers * 4, 1)

    def _uuid_for_accepted_index(self, accepted_index: int) -> str:
        return deterministic_sampler_uuid(
            prefix=self.uuid_prefix,
            seed=self._seed,
            sample_index=accepted_index,
        )

    def flush_batch(self, batch: list[Sample]) -> int:
        if not batch:
            return 0
        if len(batch) == 1:
            self._submit(batch[0])
        else:
            self._submit_group(batch)
        # At most one boolean check per batch in the normal submission path.
        self.capture_checkpoint_at_batch_boundary()
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
        return self._stream_exhausted()

    def advance_to_persisted_prefix(self, prefix: int) -> int:
        """Replay cheap proposals locally through the durable prefix.

        No worker task is submitted here.  This is a resume-only optimisation:
        DATABASE has proved ``[0, prefix)`` durable, so rebuilding the
        deterministic generator state locally avoids recalculating those
        physics points when the latest checkpoint is older than the watermark
        (or absent altogether).

        Bridson-style submit ‰ progress is suppressed for the replay so a
        resume does not look like it is re-submitting durable points (D21.14).
        """
        target = max(0, int(prefix or 0))
        previous = bool(getattr(self, "_suppress_submit_progress", False))
        self._suppress_submit_progress = True
        try:
            while self._accepted_index < target:
                if self.propose_next() is None:
                    break
        finally:
            self._suppress_submit_progress = previous
            # Fresh bar after local replay so the first real submit is 0‰ again.
            if hasattr(self, "barinfo"):
                self.barinfo = {}
        return int(self._accepted_index)

    def _common_export_fields(self) -> dict[str, Any]:
        return {
            "index": int(self._index),
            "accepted_index": int(self._accepted_index),
            "selectionexp": self._selectionexp,
            "seed": int(self._seed),
            "chains": [],
            "ready_queue": [],
            "control_state": self._checkpoint_control_state(),
        }

    def _import_common_fields(self, state: Mapping[str, Any]) -> None:
        self._index = int(state.get("index", self._index) or 0)
        self._accepted_index = int(state.get("accepted_index", self._accepted_index) or 0)
        self._selectionexp = state.get("selectionexp", self._selectionexp)
        self._seed = int(state.get("seed", self._seed) or self._seed)
        self._import_checkpoint_control_state(state)


__all__ = ["FixedSetSampler"]
