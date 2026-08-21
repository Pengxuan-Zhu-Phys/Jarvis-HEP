#!/usr/bin/env python3
"""SAMPLE bucket allocator keyspace (D25.4)."""

from __future__ import annotations

import os
import time
from typing import Any

from jarvishep2.queue.redis_queue import (
    BUCKET_LOCK,
    BUCKET_META,
    BUCKET_READY_QUEUE,
    BUCKET_STATE,
    _redis_text,
    encode_payload,
)


class _SampleBuckets:
    """Private RedisQueue mixin (D25.4)."""

    def init_sample_buckets(
        self,
        *,
        sample_root: str,
        limit: int = 200,
        width: int = 6,
        start_bucket: int = 1,
        pack: bool = True,
        enabled: bool = True,
    ) -> dict[str, Any]:
        """Initialize Redis SAMPLE bucket meta for this run."""
        self._require_client()
        root = os.path.abspath(str(sample_root))
        os.makedirs(root, exist_ok=True)
        mapping = {
            "sample_root": root,
            "limit": int(max(1, limit)),
            "width": int(max(1, width)),
            "current": int(max(1, start_bucket)),
            "count": 0,
            "pack": 1 if pack else 0,
            "enabled": 1 if enabled else 0,
        }
        self.r.delete(BUCKET_META, BUCKET_READY_QUEUE)
        self.r.hset(BUCKET_META, mapping=mapping)
        return dict(mapping)

    def _acquire_bucket_lock(
        self,
        *,
        timeout_sec: float = 15.0,
        sleep_sec: float = 0.005,
        hold_ms: int = 5000,
        # Legacy kwargs kept for callers/tests.
        retries: int | None = None,
        sleep_sec_legacy: float | None = None,
    ) -> bool:
        """Acquire SAMPLE bucket meta lock (fakeredis-safe; no Lua required).

        Dynesty / multi-worker scans call allocate+finish on every sample. A
        ~100ms spin (old default) dies under 10+ workers. Wait up to
        *timeout_sec* with a short sleep; SET NX + TTL recovers dead holders.
        """
        if sleep_sec_legacy is not None:
            sleep_sec = float(sleep_sec_legacy)
        if retries is not None:
            # Historical API: retries * sleep ≈ total wait.
            timeout_sec = max(float(timeout_sec), float(retries) * float(sleep_sec))
        deadline = time.monotonic() + max(0.05, float(timeout_sec))
        hold_ms = max(500, int(hold_ms))
        while time.monotonic() < deadline:
            # SET NX with TTL keeps a crashed holder from blocking forever.
            try:
                ok = self.r.set(BUCKET_LOCK, "1", nx=True, px=hold_ms)
            except TypeError:
                # Older redis-py / fakeredis spellings.
                ok = self.r.set(
                    BUCKET_LOCK, "1", nx=True, ex=max(1, int(hold_ms / 1000))
                )
            if ok:
                return True
            time.sleep(max(0.001, float(sleep_sec)))
        return False

    def _release_bucket_lock(self) -> None:
        try:
            self.r.delete(BUCKET_LOCK)
        except Exception:
            pass

    def _maybe_enqueue_bucket_pack_locked(self, bucket_id: int) -> bool:
        """If sealed, idle, fully archived, and not packing/packed → ready queue.

        Pack must wait until Archiver has written DATABASE rows for every sample
        in the bucket. Worker-side ``active==0`` alone is insufficient: packing
        early prunes SAMPLE/<bucket>/ and drops remaining archive payloads.
        """
        state_key = BUCKET_STATE.format(id=int(bucket_id))
        state = self.r.hgetall(state_key) or {}
        try:
            active = int(state.get("active") or 0)
            sealed = int(state.get("sealed") or 0)
            packed = int(state.get("packed") or 0)
            packing = int(state.get("packing") or 0)
            assigned = int(state.get("assigned") or 0)
            archived = int(state.get("archived") or 0)
        except (TypeError, ValueError):
            return False
        if (
            sealed != 1
            or active != 0
            or packed != 0
            or packing != 0
            or assigned <= 0
            or archived < assigned
        ):
            return False
        self.r.hset(state_key, mapping={"packing": 1})
        self.r.rpush(BUCKET_READY_QUEUE, str(int(bucket_id)))
        return True

    def allocate_sample_bucket(self) -> dict[str, Any] | None:
        """Assign the next SAMPLE bucket slot to a pulled sample.

        Returns ``None`` when bucket layout is disabled. Otherwise::

            {bucket_id, bucket_name, bucket_dir, count, width, sealed_id}
        """
        self._require_client()
        meta = self.r.hgetall(BUCKET_META) or {}
        if str(meta.get("enabled", "1")) in {"0", "false", "False"}:
            return None
        sample_root = str(meta.get("sample_root") or "").strip()
        if not sample_root:
            return None
        try:
            limit = max(1, int(meta.get("limit") or 200))
            width = max(1, int(meta.get("width") or 6))
        except (TypeError, ValueError):
            limit, width = 200, 6

        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            meta = self.r.hgetall(BUCKET_META) or {}
            try:
                current = max(1, int(meta.get("current") or 1))
                count = max(0, int(meta.get("count") or 0))
            except (TypeError, ValueError):
                current, count = 1, 0
            sealed_id = 0
            if count >= limit:
                sealed_id = current
                self.r.hset(BUCKET_STATE.format(id=sealed_id), mapping={"sealed": 1})
                current += 1
                count = 0
                self.r.hset(BUCKET_META, mapping={"current": current, "count": 0})
            count += 1
            bucket_id = current
            self.r.hset(
                BUCKET_META,
                mapping={
                    "current": bucket_id,
                    "count": count,
                    "limit": limit,
                    "width": width,
                },
            )
            state_key = BUCKET_STATE.format(id=bucket_id)
            self.r.hincrby(state_key, "active", 1)
            self.r.hincrby(state_key, "assigned", 1)
            # If we sealed the previous bucket and it is already idle, enqueue pack.
            if sealed_id:
                self._maybe_enqueue_bucket_pack_locked(int(sealed_id))
        finally:
            self._release_bucket_lock()

        from jarvishep2.sample_bucket import bucket_dir_path, format_bucket_name

        path = bucket_dir_path(sample_root, bucket_id, width=width)
        os.makedirs(path, exist_ok=True)
        return {
            "bucket_id": bucket_id,
            "bucket_name": format_bucket_name(bucket_id, width=width),
            "bucket_dir": path,
            "count": count,
            "width": width,
            "sealed_id": sealed_id or None,
            "sample_root": sample_root,
        }

    def finish_sample_bucket(self, bucket_id: int | str | None) -> bool:
        """Mark one Worker-finished sample (active--). Pack only after Archiver catches up."""
        self._require_client()
        if bucket_id is None or str(bucket_id).strip() == "":
            return False
        bid = int(bucket_id)
        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            state_key = BUCKET_STATE.format(id=bid)
            active = int(self.r.hincrby(state_key, "active", -1))
            if active < 0:
                self.r.hset(state_key, mapping={"active": 0})
            self.r.hincrby(state_key, "completed", 1)
            # Usually still archived < assigned here; pack waits for note_sample_archived.
            return self._maybe_enqueue_bucket_pack_locked(bid)
        finally:
            self._release_bucket_lock()

    def note_sample_archived(self, bucket_id: int | str | None) -> bool:
        """Archiver wrote one sample's DATABASE row. May enqueue pack when bucket is ready."""
        self._require_client()
        if bucket_id is None or str(bucket_id).strip() == "":
            return False
        bid = int(bucket_id)
        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            state_key = BUCKET_STATE.format(id=bid)
            self.r.hincrby(state_key, "archived", 1)
            return self._maybe_enqueue_bucket_pack_locked(bid)
        finally:
            self._release_bucket_lock()

    def seal_current_sample_bucket(self) -> bool:
        """Seal the open bucket at end-of-run; enqueue pack if fully archived and idle."""
        self._require_client()
        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            meta = self.r.hgetall(BUCKET_META) or {}
            try:
                current = int(meta.get("current") or 0)
            except (TypeError, ValueError):
                current = 0
            if current <= 0:
                return False
            self.r.hset(BUCKET_STATE.format(id=current), mapping={"sealed": 1})
            return self._maybe_enqueue_bucket_pack_locked(current)
        finally:
            self._release_bucket_lock()

    def pull_ready_bucket(self, timeout: int = 0) -> dict[str, Any] | None:
        """Pop a sealed+idle bucket ready for tar packing."""
        self._require_client()
        if timeout and timeout > 0:
            raw = self._blpop(BUCKET_READY_QUEUE, timeout=timeout)
            if raw is None:
                return None
            _, bucket_token = raw
        else:
            bucket_token = self.r.lpop(BUCKET_READY_QUEUE)
            if bucket_token is None:
                return None
        try:
            bucket_id = int(bucket_token)
        except (TypeError, ValueError):
            return None
        meta = self.r.hgetall(BUCKET_META) or {}
        sample_root = str(meta.get("sample_root") or "").strip()
        try:
            width = max(1, int(meta.get("width") or 6))
        except (TypeError, ValueError):
            width = 6
        from jarvishep2.sample_bucket import bucket_dir_path, format_bucket_name

        return {
            "bucket_id": bucket_id,
            "bucket_name": format_bucket_name(bucket_id, width=width),
            "bucket_dir": bucket_dir_path(sample_root, bucket_id, width=width) if sample_root else "",
            "sample_root": sample_root,
            "pack": str(meta.get("pack", "1")) not in {"0", "false", "False"},
        }

    def mark_bucket_packed(self, bucket_id: int | str) -> None:
        self._require_client()
        key = BUCKET_STATE.format(id=int(bucket_id))
        self.r.hset(key, mapping={"packed": 1, "packing": 0})

    def get_bucket_state(self, bucket_id: int | str) -> dict[str, Any]:
        self._require_client()
        raw = self.r.hgetall(BUCKET_STATE.format(id=int(bucket_id))) or {}
        out: dict[str, Any] = {}
        for key, value in raw.items():
            try:
                out[str(key)] = int(value)
            except (TypeError, ValueError):
                out[str(key)] = value
        return out
