#!/usr/bin/env python3
"""Feedback-driven sampler base (D13.1): propose → Redis → hep:feedback → absorb.

Extracted from AdaptiveBridson so MCMC / nested / ensemble methods share one
barrier-synchronized generation loop. Method-specific science lives in
``propose_generation`` / ``absorb_generation``; transport, pending-uuid
bookkeeping, SeedSequence spawning, and checkpoint barriers live here.

Design: ``Jarvis-Books/Jarvis-HEP V2/DESIGN_SAMPLERS_2.0.md`` §3.2 and
``components/feedback_sampler.md``.
"""

from __future__ import annotations

import time
from abc import ABC, abstractmethod
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import PROC_BOARD_TTL_SEC, RedisQueue
from jarvishep2.sample import Sample

_FEEDBACK_BLPOP_CAP_SEC = 1
_BLOCKING_SCAN_MODES = frozenset({"paused", "stopping", "draining"})


class FeedbackSampler(CheckpointedSampler, ABC):
    """Barrier-synchronized feedback sampler.

    Subclasses implement:

    * ``propose_generation()`` → samples for the next generation, or ``None``
      when the run is finished.
    * ``absorb_generation(results)`` → fold feedback records into method state.

    The default ``run_adaptive`` loop is:

    .. code-block:: text

        while True:
            proposals = propose_generation()
            if proposals is None: break
            publish batch (pending uuids tracked)
            results = wait_for_generation()   # drain hep:feedback
            absorb_generation(results)
            checkpoint_at_barrier()

    Methods with richer control flow (AdaptiveBridson: converge / refine /
    max_points) may override ``run_adaptive`` while still reusing
    ``wait_for_generation``, ``_submit_sample_batch``, and seed helpers.
    """

    method = ""

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.feedback")
        self._pending_uuids: set[str] = set()
        self._seed = 0
        self._seed_seq: np.random.SeedSequence | None = None
        self._batch_size = 16
        self._generation = 0
        # Physics-safe default for failed samples at a barrier (D13 §3.3).
        # ALS treats Failed as f=None (crossing-conservative); MCMC will reject.
        self._on_failure = "reject"  # reject | halt

    # ----------------------------------------------------------------- seeds
    def _init_seed_sequence(self, seed: int) -> None:
        """Install the root SeedSequence used by ``_generation_rng``."""
        self._seed = int(seed)
        self._seed_seq = np.random.SeedSequence(self._seed)

    def _ensure_seed_sequence(self) -> np.random.SeedSequence:
        if self._seed_seq is None:
            self._seed_seq = np.random.SeedSequence(self._seed)
        return self._seed_seq

    def _generation_rng(self, generation: int) -> np.random.Generator:
        """Deterministic per-generation child of the root SeedSequence (D6.2)."""
        root = self._ensure_seed_sequence()
        gen = max(0, int(generation))
        child = root.spawn(gen + 1)[gen]
        return np.random.default_rng(child)

    # ---------------------------------------------------------------- redis
    def _require_redis(self, what: str) -> RedisQueue:
        if self.redis is None:
            raise RuntimeError(f"{what} requires redis")
        return self.redis

    # ---------------------------------------------------------- pending/submit
    def _register_pending(self, uuid: str) -> None:
        self._pending_uuids.add(str(uuid))

    def _clear_pending(self, uuid: str) -> None:
        self._pending_uuids.discard(str(uuid))

    def _submit_sample_batch(self, samples: Sequence[Sample]) -> int:
        """Push samples in groups of ``_batch_size``; track pending + submitted.

        Caller is responsible for assigning deterministic uuids and any
        method-specific bookkeeping *before* calling this method. Each sample
        uuid is registered as pending.
        """
        if not samples:
            return 0
        self._require_redis(f"{type(self).__name__}._submit_sample_batch")
        batch: list[Sample] = []
        n = 0
        for sample in samples:
            uuid = str(sample.uuid)
            self._register_pending(uuid)
            batch.append(sample)
            n += 1
            if len(batch) >= self._batch_size:
                self._submit_group(list(batch))
                self._submitted_uuids.extend(s.uuid for s in batch)
                batch = []
        if batch:
            self._submit_group(list(batch))
            self._submitted_uuids.extend(s.uuid for s in batch)
        return n

    # --------------------------------------------------------- feedback drain
    def wait_for_any_feedback(
        self,
        *,
        timeout: float,
        queues: Sequence[str] | None = None,
    ) -> dict[str, Any]:
        """Block until one pending uuid returns on the feedback channel(s).

        Used by independent-chain MCMC: advance a single chain as soon as its
        evaluation completes, without waiting for a generation barrier.

        ``queues`` selects one or more Redis feedback lists (default
        ``hep:feedback``). Per-chain shards use multi-key BLPOP so the control
        process wakes on whichever chain finishes first. Same 1s BLPOP cap and
        ``scan_mode`` checks as ``wait_for_generation``.
        """
        if not self._pending_uuids:
            raise RuntimeError(
                f"{type(self).__name__}.wait_for_any_feedback called with no pending uuids"
            )
        redis = self._require_redis(f"{type(self).__name__}.wait_for_any_feedback")
        deadline = time.monotonic() + max(1.0, float(timeout))
        last_scan_mode: str | None = None
        missing_since: float | None = None
        board_ttl = self._core_board_ttl_sec()
        what = f"feedback wait (generation={self._generation})"
        while self._pending_uuids:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(
                    f"{type(self).__name__} feedback wait timed out "
                    f"with {len(self._pending_uuids)} pending sample(s) "
                    f"(generation={self._generation})"
                )
            record, last_scan_mode, missing_since = self._await_feedback_round(
                redis,
                queues=queues,
                last_scan_mode=last_scan_mode,
                missing_since=missing_since,
                board_ttl=board_ttl,
                what=what,
            )
            if record is None:
                continue
            uuid = str(record.get("uuid", ""))
            if uuid not in self._pending_uuids:
                self._logger.warning(
                    "dropping unmatched hep:feedback uuid=%s "
                    "(pending=%d generation=%d logL=%s)",
                    uuid or "<empty>",
                    len(self._pending_uuids),
                    self._generation,
                    record.get("logL", ""),
                )
                continue
            self._clear_pending(uuid)
            self._on_feedback_record(record)
            return dict(record)
        raise RuntimeError(
            f"{type(self).__name__}.wait_for_any_feedback: pending set emptied mid-wait"
        )

    def wait_for_generation(
        self,
        *,
        timeout: float,
        queues: Sequence[str] | None = None,
    ) -> list[dict[str, Any]]:
        """Drain feedback until every pending uuid has a record.

        Returns the matched feedback records (order = arrival order). Unrelated
        uuids on the channel are ignored (watchdog re-queues stay transparent:
        the barrier waits on uuids, not workers).

        ``queues`` may list per-chain feedback shards; when omitted the default
        ``hep:feedback`` list is used.

        Feedback BLPOP is capped at 1s so ``scan_mode`` can be re-read each
        round. A core board that was never published is treated as running.
        """
        redis = self._require_redis(f"{type(self).__name__}.wait_for_generation")
        deadline = time.monotonic() + max(1.0, float(timeout))
        matched: list[dict[str, Any]] = []
        last_scan_mode: str | None = None
        missing_since: float | None = None
        board_ttl = self._core_board_ttl_sec()
        what = f"generation {self._generation}"
        while self._pending_uuids:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(
                    f"{type(self).__name__} generation {self._generation} timed out "
                    f"with {len(self._pending_uuids)} pending sample(s)"
                )
            record, last_scan_mode, missing_since = self._await_feedback_round(
                redis,
                queues=queues,
                last_scan_mode=last_scan_mode,
                missing_since=missing_since,
                board_ttl=board_ttl,
                what=what,
            )
            if record is None:
                continue
            uuid = str(record.get("uuid", ""))
            if uuid not in self._pending_uuids:
                self._logger.warning(
                    "dropping unmatched hep:feedback uuid=%s "
                    "(pending=%d generation=%d logL=%s)",
                    uuid or "<empty>",
                    len(self._pending_uuids),
                    self._generation,
                    record.get("logL", ""),
                )
                continue
            self._clear_pending(uuid)
            matched.append(dict(record))
            self._on_feedback_record(record)
        return matched

    def _core_board_ttl_sec(self) -> float:
        override = getattr(self, "_board_ttl_sec", None)
        if override is not None:
            return float(override)
        stale_sec = 30.0
        try:
            from jarvishep2.runtime_config import get_watchdog_config

            stale_sec = float(
                get_watchdog_config(getattr(self, "config", None)).get("stale_sec", 30.0)
            )
        except Exception:
            stale_sec = 30.0
        return max(float(PROC_BOARD_TTL_SEC), 2.0 * stale_sec)

    def _generation_interrupt_requested(self) -> bool:
        if bool(getattr(self, "_interrupt_requested", False)):
            return True
        core = getattr(self, "_core", None)
        return bool(core is not None and getattr(core, "_interrupt_requested", False))

    def _fail_if_generation_blocked(self, *, what: str) -> None:
        if not self._generation_interrupt_requested():
            return
        raise RuntimeError(
            f"{type(self).__name__} {what} aborted: "
            "scan_mode=stopping pause_reason=interrupt"
        )

    def _await_feedback_round(
        self,
        redis: RedisQueue,
        *,
        queues: Sequence[str] | None,
        last_scan_mode: str | None,
        missing_since: float | None,
        board_ttl: float,
        what: str,
    ) -> tuple[dict[str, Any] | None, str | None, float | None]:
        """One 1s BLPOP plus scan_mode / interrupt checks (KD-15)."""
        self._fail_if_generation_blocked(what=what)
        mode, pause_reason, seen = self._read_core_scan_mode(redis)
        if seen:
            last_scan_mode = mode
            missing_since = None
            if mode in _BLOCKING_SCAN_MODES:
                raise RuntimeError(
                    f"{type(self).__name__} {what} aborted: "
                    f"scan_mode={mode} pause_reason={pause_reason or mode}"
                )
        elif last_scan_mode is not None:
            if missing_since is None:
                missing_since = time.monotonic()
            if (time.monotonic() - missing_since) >= board_ttl:
                raise RuntimeError(
                    f"{type(self).__name__} {what} aborted: "
                    "core proc board missing (scan_mode unknown)"
                )
        try:
            record = redis.pull_feedback(
                timeout=_FEEDBACK_BLPOP_CAP_SEC, queues=queues
            )
        except Exception as exc:
            if self._generation_interrupt_requested() or _redis_closed(exc):
                raise RuntimeError(
                    f"{type(self).__name__} {what} aborted: "
                    "scan_mode=stopping pause_reason=interrupt"
                ) from exc
            raise
        return record, last_scan_mode, missing_since

    def _read_core_scan_mode(self, redis: RedisQueue) -> tuple[str, str, bool]:
        reader = getattr(redis, "read_proc_board", None)
        if not callable(reader):
            return "", "", False
        try:
            board = reader("core") or {}
        except Exception as exc:
            if _redis_closed(exc):
                return "stopping", "interrupt", True
            return "", "", False
        if not isinstance(board, dict) or not board:
            return "", "", False
        mode = str(board.get("scan_mode") or "").strip().lower()
        if not mode:
            return "", "", False
        pause_reason = str(board.get("pause_reason") or "").strip()
        return mode, pause_reason, True

    def _on_feedback_record(self, record: Mapping[str, Any]) -> None:
        """Optional streaming hook fired once per matched feedback record.

        Default is a no-op. Prefer ``absorb_generation`` for bulk processing
        after the barrier; override this only when state must update mid-drain
        (rare — ALS historically filled ``f`` here; it now uses absorb).
        """
        return None

    def _failure_policy_halt(self, record: Mapping[str, Any]) -> bool:
        """Return True if an unusable feedback record should abort the run.

        D13.8: no sample ``status`` on feedback — unusable means non-finite
        ``logL`` (−∞ / nan).
        """
        from jarvishep2.feedback_return import (
            WIRE_LOGL_KEY,
            feedback_logl,
            is_unusable_logl,
            normalize_feedback_record,
        )

        flat = normalize_feedback_record(record)
        if WIRE_LOGL_KEY not in flat:
            return False
        if not is_unusable_logl(feedback_logl(flat)):
            return False
        return str(self._on_failure).lower() == "halt"

    # -------------------------------------------------------------- barrier
    def at_safe_barrier(self) -> bool:
        """Safe to checkpoint when no generation is in flight."""
        return not self._pending_uuids

    def checkpoint_at_barrier(self, *, reason: str = "generation_barrier") -> bool:
        """Checkpoint only after the completed generation is durable.

        Feedback arrives before the process Archiver necessarily fsyncs its
        HDF5 batch.  Wait briefly for Redis's post-fsync UUID set; if it does
        not catch up, the checkpoint callback records the previous safe state
        with ``safe_barrier_confirmed=False`` rather than advancing past disk.
        """
        redis = self.redis
        scan_block = self.config.get("Scan")
        scan_mapping = scan_block if isinstance(scan_block, Mapping) else {}
        scan_name = str(
            self.config.get("scan_name") or scan_mapping.get("name") or "scan"
        )
        if redis is not None and self._submitted_uuids:
            deadline = time.monotonic() + 5.0
            while time.monotonic() < deadline:
                persisted = redis.get_archived_uuids(scan_name)
                self.set_persisted_uuids(persisted)
                if set(self._submitted_uuids) <= persisted:
                    break
                time.sleep(0.05)
        return self.persist_runtime_checkpoint(force=True, reason=reason)

    # ---------------------------------------------------------- method hooks
    @abstractmethod
    def propose_generation(self) -> Sequence[Sample] | None:
        """Build the next generation of Samples, or ``None`` when finished.

        Return an empty sequence to absorb an empty generation (no Redis work)
        without terminating — useful for bookkeeping-only steps. Return
        ``None`` to exit ``run_adaptive``.
        """

    @abstractmethod
    def absorb_generation(self, results: Sequence[Mapping[str, Any]]) -> None:
        """Fold barrier results into method state (accept/reject, fill f, …)."""

    # ----------------------------------------------------------------- driver
    def run_adaptive(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Default generation loop with a per-generation timeout.

        Returns the number of samples submitted across all generations.
        Subclasses with non-linear control flow (e.g. AdaptiveBridson) override
        this method and call ``wait_for_generation`` / ``_submit_sample_batch``
        directly.
        """
        if timeout is not None:
            generation_timeout = timeout
        self._require_redis(f"{type(self).__name__}.run_adaptive")
        self._ensure_seed_sequence()
        total_submitted = 0
        while True:
            proposals = self.propose_generation()
            if proposals is None:
                break
            batch = list(proposals)
            if batch:
                total_submitted += self._submit_sample_batch(batch)
                results = self.wait_for_generation(timeout=generation_timeout)
            else:
                results = []
            self.absorb_generation(results)
            self._generation += 1
            self.checkpoint_at_barrier()
        return total_submitted

    def run_distributed(self) -> int:
        """Alias so generic drivers can call run_distributed when registered."""
        return self.run_adaptive()

    # ----------------------------------------------------------- checkpoint IO
    def _feedback_export_state(self) -> dict[str, Any]:
        """Shared feedback fields for ``export_runtime_state`` subclasses."""
        return {
            "generation": self._generation,
            "pending_uuids": sorted(self._pending_uuids),
            "seed": self._seed,
            "submitted_uuids": list(self._submitted_uuids),
            "control_state": self._checkpoint_control_state(),
            "on_failure": self._on_failure,
        }

    def _feedback_import_state(self, state: Mapping[str, Any]) -> None:
        """Restore shared feedback fields from a runtime checkpoint."""
        self._generation = int(state.get("generation", 0) or 0)
        self._pending_uuids = {str(u) for u in (state.get("pending_uuids") or [])}
        self._seed = int(state.get("seed", self._seed) or self._seed)
        self._seed_seq = np.random.SeedSequence(self._seed)
        self._on_failure = str(state.get("on_failure") or self._on_failure or "reject")
        self._import_checkpoint_control_state(state)


def _redis_closed(exc: BaseException) -> bool:
    text = str(exc).lower()
    return "not connected" in text or "connection closed" in text


__all__ = ["FeedbackSampler"]
