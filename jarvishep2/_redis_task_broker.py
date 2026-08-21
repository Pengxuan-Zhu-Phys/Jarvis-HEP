#!/usr/bin/env python3
"""Task / archive / feedback Redis lists (D25.4)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.redis_queue import (
    ARCHIVE_QUEUE,
    ARCHIVED_INDEX_PREFIX,
    ARCHIVED_UUIDS,
    CHAIN_FEEDBACK_QUEUE_PATTERN,
    CodecError,
    FEEDBACK_QUEUE,
    OP_COUNT,
    RESULTS,
    SAMPLE_STATS,
    TASK_QUEUE,
    TaskValidationError,
    _normalize_feedback_queue_name,
    _redis_text,
    _validate_result_payload,
    _validate_task_payload,
    decode_payload,
    encode_payload,
)


class _TaskBroker:
    """Private RedisQueue mixin (D25.4)."""

    def push_task(self, task: Mapping[str, Any]) -> None:
        self._require_client()
        _validate_task_payload(task)
        encoded = encode_payload(dict(task), codec=self._codec)
        pipe = self.r.pipeline(transaction=True)
        pipe.rpush(TASK_QUEUE, encoded)
        pipe.incr(OP_COUNT.format(kind="task"))
        pipe.execute()

    def pull_task(self, timeout: int = 5) -> dict[str, Any] | None:
        self._require_client()
        raw = self._blpop(TASK_QUEUE, timeout=timeout)
        if raw is None:
            return None
        _, payload = raw
        decoded = decode_payload(payload, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("task payload must decode to a dict")
        # Sample left the queue and is now in-flight on a Worker. Only the
        # stats hash moves here: the D1.1 contract is a single sample
        # op_count increment per sample, in submit_result. Backpressure
        # reads fetch_sample_stats() directly, so it needs no op bump.
        self.r.hincrby(SAMPLE_STATS, "running", 1)
        return decoded

    def drain_task_queue(self) -> int:
        """Discard all queued tasks (resume path — WP-D6.2)."""
        self._require_client()
        drained = 0
        while int(self.r.llen(TASK_QUEUE)) > 0:
            # Bypass pull_task so we do not bump SAMPLE_STATS.running for discarded work.
            raw = self._blpop(TASK_QUEUE, timeout=1)
            if raw is None:
                break
            drained += 1
        return drained

    def add_archived_uuids(self, scan_name: str, uuids: Sequence[str]) -> int:
        """Publish UUIDs only after their HDF5 batch is durable."""
        self._require_client()
        values = [str(item).strip() for item in uuids if str(item).strip()]
        if not values:
            return 0
        key = ARCHIVED_UUIDS.format(scan=str(scan_name or "scan").strip() or "scan")
        return int(self.r.sadd(key, *values))

    def get_archived_uuids(self, scan_name: str) -> set[str]:
        """Read the process-visible durable-archive acknowledgement cache."""
        self._require_client()
        key = ARCHIVED_UUIDS.format(scan=str(scan_name or "scan").strip() or "scan")
        return {_redis_text(item) for item in (self.r.smembers(key) or set())}

    def clear_archived_uuids(self, scan_name: str) -> int:
        """Clear the Redis cache for a fresh run; DATABASE remains authoritative."""
        self._require_client()
        key = ARCHIVED_UUIDS.format(scan=str(scan_name or "scan").strip() or "scan")
        return int(self.r.delete(key))

    def set_archived_index_prefix(self, scan_name: str, prefix: int) -> None:
        """Publish the HDF5-durable contiguous sample-index prefix."""
        self._require_client()
        key = ARCHIVED_INDEX_PREFIX.format(scan=str(scan_name or "scan").strip() or "scan")
        self.r.set(key, max(0, int(prefix)))

    def get_archived_index_prefix(self, scan_name: str) -> int:
        """Read the durable prefix in O(1); Redis remains only a cache."""
        self._require_client()
        key = ARCHIVED_INDEX_PREFIX.format(scan=str(scan_name or "scan").strip() or "scan")
        raw = self.r.get(key)
        try:
            return max(0, int(_redis_text(raw))) if raw is not None else 0
        except (TypeError, ValueError):
            return 0

    def clear_archived_index_prefix(self, scan_name: str) -> int:
        self._require_client()
        key = ARCHIVED_INDEX_PREFIX.format(scan=str(scan_name or "scan").strip() or "scan")
        return int(self.r.delete(key))

    def push_many_tasks(self, tasks: list[Mapping[str, Any]]) -> None:
        if not tasks:
            return
        self._require_client()
        for task in tasks:
            _validate_task_payload(task)
        encoded = [encode_payload(dict(task), codec=self._codec) for task in tasks]
        pipe = self.r.pipeline(transaction=True)
        for item in encoded:
            pipe.rpush(TASK_QUEUE, item)
        pipe.incrby(OP_COUNT.format(kind="task"), len(tasks))
        pipe.execute()

    def submit_result(self, info: Mapping[str, Any]) -> None:
        self._require_client()
        _validate_result_payload(info)
        encoded = encode_payload(dict(info), codec=self._codec)
        status = str(info.get("status", "Completed"))
        pipe = self.r.pipeline(transaction=True)
        pipe.rpush(ARCHIVE_QUEUE, encoded)
        if status == "Failed":
            pipe.hincrby(SAMPLE_STATS, "failed", 1)
            pipe.hincrby(SAMPLE_STATS, "running", -1)
        else:
            pipe.hincrby(SAMPLE_STATS, "completed", 1)
            pipe.hincrby(SAMPLE_STATS, "running", -1)
        pipe.incr(OP_COUNT.format(kind="sample"))
        pipe.execute()

    def publish_feedback(
        self,
        info: Mapping[str, Any],
        *,
        queue: str | None = None,
    ) -> None:
        """Push a flat sampler barrier record (D13.8).

        Expected shape::

            {"uuid": str, "logL": float|str, ...optional flat fields}

        Nested ``observables`` / sample ``status`` are not part of the contract;
        if a legacy caller still passes them they are dropped (not re-nested).

        ``queue`` selects the Redis list (default ``hep:feedback``). Independent
        MCMC chains may use ``hep:feedback:chain:{id}`` while still sharing the
        single task queue.
        """
        self._require_client()
        from jarvishep2.feedback_return import (
            WIRE_LOGL_KEY,
            encode_feedback_float,
            extract_logl_total,
            normalize_feedback_record,
        )

        raw = dict(info or {})
        # Accept already-flat records or legacy nested bags from older tests.
        if "observables" in raw or (
            WIRE_LOGL_KEY not in raw and "LogL" in dict(raw.get("observables") or {})
        ):
            payload = normalize_feedback_record(raw)
        else:
            payload = {
                key: value
                for key, value in raw.items()
                if key not in {"status", "observables"}
            }
            payload["uuid"] = str(raw.get("uuid") or payload.get("uuid") or "")
            if WIRE_LOGL_KEY not in payload and "LogL" in raw:
                payload[WIRE_LOGL_KEY] = raw["LogL"]
            if WIRE_LOGL_KEY not in payload:
                nested = raw.get("observables")
                if isinstance(nested, Mapping):
                    logl = extract_logl_total(nested)
                    if logl is not None:
                        payload[WIRE_LOGL_KEY] = logl
        if not payload.get("uuid"):
            raise TaskValidationError("feedback payload requires uuid")
        # Encode non-finite floats for JSON codecs.
        if WIRE_LOGL_KEY in payload:
            try:
                payload[WIRE_LOGL_KEY] = encode_feedback_float(
                    float(payload[WIRE_LOGL_KEY])
                )
            except (TypeError, ValueError):
                from jarvishep2.feedback_return import decode_feedback_float

                payload[WIRE_LOGL_KEY] = encode_feedback_float(
                    decode_feedback_float(payload[WIRE_LOGL_KEY])
                )
        encoded = encode_payload(payload, codec=self._codec)
        self.r.rpush(_normalize_feedback_queue_name(queue), encoded)

    def pull_feedback(
        self,
        timeout: int = 1,
        *,
        queues: Sequence[str] | None = None,
    ) -> dict[str, Any] | None:
        """Blocking pop from one or more feedback lists; normalizes flat floats.

        When ``queues`` is omitted, pops from the default ``hep:feedback``.
        When several keys are given, Redis returns the first non-empty list
        (MCMC async control uses this to wait on any ready chain shard).
        """
        self._require_client()
        keys = [
            _normalize_feedback_queue_name(key)
            for key in (list(queues) if queues is not None else [FEEDBACK_QUEUE])
        ]
        # Preserve order, drop duplicates (BLPOP rejects empty key lists).
        seen: set[str] = set()
        unique_keys: list[str] = []
        for key in keys:
            if key not in seen:
                seen.add(key)
                unique_keys.append(key)
        if not unique_keys:
            unique_keys = [FEEDBACK_QUEUE]
        if len(unique_keys) == 1:
            raw = self._blpop(unique_keys[0], timeout=timeout)
        else:
            raw = self._blpop_many(unique_keys, timeout=float(timeout))
        if raw is None:
            return None
        _, payload = raw
        decoded = decode_payload(payload, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("feedback payload must decode to a dict")
        from jarvishep2.feedback_return import normalize_feedback_record

        return normalize_feedback_record(decoded)

    def list_feedback_queue_keys(self) -> list[str]:
        """Default feedback list plus any live per-chain shards."""
        self._require_client()
        keys: list[str] = [FEEDBACK_QUEUE]
        try:
            for key in self.r.scan_iter(match=CHAIN_FEEDBACK_QUEUE_PATTERN):
                text = key.decode() if isinstance(key, (bytes, bytearray)) else str(key)
                if text and text not in keys:
                    keys.append(text)
        except Exception:
            pass
        return keys

    def drain_feedback_queue(self) -> int:
        """Discard all queued feedback records (resume path).

        Drains the default list and any ``hep:feedback:chain:*`` shards so a
        crash-era independent-chain run cannot leak stale feedback after resume.
        """
        self._require_client()
        drained = 0
        for key in self.list_feedback_queue_keys():
            while int(self.r.llen(key) or 0) > 0:
                item = self.pull_feedback(timeout=1, queues=[key])
                if item is None:
                    break
                drained += 1
        return drained

    def pull_result(self, timeout: int = 1) -> dict[str, Any] | None:
        self._require_client()
        raw = self._blpop(ARCHIVE_QUEUE, timeout=timeout)
        if raw is None:
            return None
        _, payload = raw
        decoded = decode_payload(payload, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("result payload must decode to a dict")
        return decoded

    def store_result_hash(self, uuid: str, info: Mapping[str, Any]) -> None:
        """Optional per-uuid result hash for debug/sync handoff."""
        self._require_client()
        _validate_result_payload(info)
        key = RESULTS.format(uuid=uuid)
        encoded = encode_payload(dict(info), codec=self._codec)
        self.r.hset(key, mapping={"payload": encoded})
