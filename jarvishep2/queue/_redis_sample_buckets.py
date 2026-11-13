#!/usr/bin/env python3
"""SAMPLE bucket allocator keyspace (D25.4)."""

from __future__ import annotations

import os
import time
from typing import Any

from jarvishep2.queue._redis_calc_pool import _lua_eval_unavailable
from jarvishep2.queue.redis_queue import (
    BUCKET_LOCK,
    BUCKET_META,
    BUCKET_READY_QUEUE,
    BUCKET_STATE,
    _redis_text,
)

BUCKET_PACK_ACTIVE = "hep:sample:bucket:pack_active"


_UPDATE_BUCKET = """
local meta, ready = KEYS[1], KEYS[2]
local operation, state_prefix = ARGV[1], ARGV[2]

local function enqueue(bucket, state)
    if redis.call('HGET', meta, 'pack') == '0' then return 0 end
    local v = redis.call('HMGET', state, 'sealed', 'active', 'packing',
                         'packed', 'assigned', 'archived')
    local assigned = tonumber(v[5]) or 0
    if v[1] ~= '1' or (tonumber(v[2]) or 0) ~= 0 or v[3] == '1' or
       v[4] == '1' or assigned <= 0 or (tonumber(v[6]) or 0) < assigned then
        return 0
    end
    redis.call('HSET', state, 'packing', 1)
    redis.call('RPUSH', ready, bucket)
    if redis.call('HGET', meta, 'observe') == '1' then
        redis.call('HINCRBY', meta, 'ready_buckets', 1)
        redis.call('HINCRBY', meta, 'outstanding_points', assigned)
    end
    return 1
end

if operation == 'allocate' then
    local v = redis.call('HMGET', meta, 'enabled', 'sample_root', 'limit',
                         'width', 'current', 'count')
    if v[1] == '0' or not v[2] or v[2] == '' then return {} end
    local limit = tonumber(v[3]) or 200
    local width = tonumber(v[4]) or 6
    local bucket = v[5] or '1'
    local count = tonumber(v[6]) or 0
    local sealed = ''
    if count >= limit then
        sealed = bucket
        redis.call('HINCRBY', meta, 'current', 1)
        bucket = redis.call('HGET', meta, 'current')
        count = 0
        redis.call('HSET', state_prefix .. sealed, 'sealed', 1)
        enqueue(sealed, state_prefix .. sealed)
    end
    count = count + 1
    redis.call('HSET', meta, 'count', count)
    local state = state_prefix .. bucket
    redis.call('HINCRBY', state, 'active', 1)
    redis.call('HINCRBY', state, 'assigned', 1)
    return {bucket, count, width, sealed, v[2]}
end

local bucket = ARGV[3]
if operation == 'seal' then
    bucket = redis.call('HGET', meta, 'current')
    if not bucket then return 0 end
end
local state = state_prefix .. bucket
if operation == 'finish' then
    if redis.call('HINCRBY', state, 'active', -1) < 0 then
        redis.call('HSET', state, 'active', 0)
    end
    redis.call('HINCRBY', state, 'completed', 1)
elseif operation == 'archived' then
    redis.call('HINCRBY', state, 'archived', ARGV[4])
elseif operation == 'seal' then
    redis.call('HSET', state, 'sealed', 1)
end
return enqueue(bucket, state)
"""

_PACK_TRANSITION = """
local meta, ready, active = KEYS[1], KEYS[2], KEYS[3]
local op, owner, bucket, arg = ARGV[1], ARGV[2], ARGV[3], ARGV[4]
if redis.call('HGET', meta, 'observe') ~= '1' then return false end
local now = redis.call('TIME')
local function progress()
    redis.call('HSET', meta, 'pack_progress', now[1] .. '.' .. string.format('%06d', now[2]))
end
if op == 'start' then
    if redis.call('HGET', meta, 'pack_owner') == owner then return 1 end
    local abandoned = redis.call('HVALS', active)
    for _, count in ipairs(abandoned) do
        redis.call('HINCRBY', meta, 'outstanding_points', -tonumber(count))
    end
    redis.call('HINCRBY', meta, 'failed_buckets', #abandoned)
    redis.call('DEL', active)
    redis.call('HSET', meta, 'pack_owner', owner, 'packer_slots', arg,
               'active_buckets', 0, 'completed_buckets', 0, 'completed_points', 0)
    progress()
    return 1
end
if redis.call('HGET', meta, 'pack_owner') ~= owner then return false end
if op == 'claim' then
    bucket = redis.call('LPOP', ready)
    if not bucket then return false end
    local assigned = redis.call('HGET', ARGV[5] .. bucket, 'assigned') or '0'
    redis.call('HSET', active, bucket, assigned)
    redis.call('HINCRBY', meta, 'ready_buckets', -1)
    redis.call('HINCRBY', meta, 'active_buckets', 1)
    progress()
    return bucket
end
if op == 'progress' then
    if redis.call('HEXISTS', active, bucket) == 1 then progress() end
    return 1
end
local assigned = redis.call('HGET', active, bucket)
if not assigned then return 0 end
redis.call('HDEL', active, bucket)
redis.call('HINCRBY', meta, 'active_buckets', -1)
redis.call('HINCRBY', meta, 'outstanding_points', -tonumber(assigned))
if arg == '1' then
    redis.call('HSET', ARGV[5] .. bucket, 'packed', 1, 'packing', 0)
    redis.call('HINCRBY', meta, 'completed_buckets', 1)
    redis.call('HINCRBY', meta, 'completed_points', assigned)
else
    redis.call('HINCRBY', meta, 'failed_buckets', 1)
end
progress()
return 1
"""


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
        observe: bool = False,
    ) -> dict[str, Any]:
        """Initialize Redis SAMPLE bucket meta for this run."""
        self._require_client()
        root = os.path.abspath(str(sample_root))
        if enabled:
            os.makedirs(root, exist_ok=True)
        mapping = {
            "sample_root": root,
            "limit": int(max(1, limit)),
            "width": int(max(1, width)),
            "current": int(max(1, start_bucket)),
            "count": 0,
            "pack": 1 if pack else 0,
            "enabled": 1 if enabled else 0,
            "observe": int(observe),
        }
        self.r.delete(BUCKET_META, BUCKET_READY_QUEUE, BUCKET_PACK_ACTIVE)
        self.r.hset(BUCKET_META, mapping=mapping)
        return dict(mapping)

    def _update_bucket(self, operation: str, *args: Any) -> Any:
        """Apply one bucket transition atomically on Redis."""
        if getattr(self, "_bucket_eval_unavailable", False):
            return None
        try:
            return self.r.eval(
                _UPDATE_BUCKET, 2, BUCKET_META, BUCKET_READY_QUEUE,
                operation, BUCKET_STATE.format(id=""), *args,
            )
        except Exception as exc:
            if not _lua_eval_unavailable(exc):
                raise
            # Offline clients such as fakeredis without Lua use the lock path.
            self._bucket_eval_unavailable = True
            return None

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

    def _bucket_meta(self) -> dict[str, str]:
        return {
            _redis_text(key): _redis_text(value)
            for key, value in (self.r.hgetall(BUCKET_META) or {}).items()
        }

    def _maybe_enqueue_bucket_pack_locked(self, bucket_id: int) -> bool:
        """If sealed, idle, fully archived, and not packing/packed → ready queue.

        Pack must wait until Archiver has written DATABASE rows for every sample
        in the bucket. Worker-side ``active==0`` alone is insufficient: packing
        early prunes SAMPLE/<bucket>/ and drops remaining archive payloads.
        """
        if _redis_text(self.r.hget(BUCKET_META, "pack")) == "0":
            return False
        state_key = BUCKET_STATE.format(id=int(bucket_id))
        state = {
            _redis_text(key): value
            for key, value in (self.r.hgetall(state_key) or {}).items()
        }
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
        allocation = self._update_bucket("allocate")
        if allocation is None:
            allocation = self._allocate_bucket_with_lock()
        if not allocation:
            return None
        bucket, count, width, sealed, root = allocation
        bucket_id = int(bucket)
        width = int(width)
        sample_root = _redis_text(root)

        from jarvishep2.sample_bucket import bucket_dir_path, format_bucket_name

        path = bucket_dir_path(sample_root, bucket_id, width=width)
        os.makedirs(path, exist_ok=True)
        return {
            "bucket_id": bucket_id,
            "bucket_name": format_bucket_name(bucket_id, width=width),
            "bucket_dir": path,
            "count": int(count),
            "width": width,
            "sealed_id": int(sealed) if sealed else None,
            "sample_root": sample_root,
        }

    def _allocate_bucket_with_lock(self) -> tuple[int, int, int, int, str] | None:
        meta = self._bucket_meta()
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
            meta = self._bucket_meta()
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

        return bucket_id, count, width, sealed_id, sample_root

    def finish_sample_bucket(self, bucket_id: int | str | None) -> bool:
        """Mark one Worker-finished sample (active--). Pack only after Archiver catches up."""
        self._require_client()
        if bucket_id is None or str(bucket_id).strip() == "":
            return False
        bid = int(bucket_id)
        result = self._update_bucket("finish", bid)
        if result is not None:
            return bool(result)
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

    def note_sample_archived(self, bucket_id: int | str | None, *, count: int = 1) -> bool:
        """Record durable DATABASE rows in one bucket and enqueue it when ready."""
        self._require_client()
        if bucket_id is None or str(bucket_id).strip() == "" or count <= 0:
            return False
        bid = int(bucket_id)
        result = self._update_bucket("archived", bid, count)
        if result is not None:
            return bool(result)
        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            state_key = BUCKET_STATE.format(id=bid)
            self.r.hincrby(state_key, "archived", count)
            return self._maybe_enqueue_bucket_pack_locked(bid)
        finally:
            self._release_bucket_lock()

    def seal_current_sample_bucket(self) -> bool:
        """Seal the open bucket at end-of-run; enqueue pack if fully archived and idle."""
        self._require_client()
        result = self._update_bucket("seal")
        if result is not None:
            return bool(result)
        if not self._acquire_bucket_lock():
            raise RuntimeError("timed out acquiring SAMPLE bucket lock")
        try:
            meta = self._bucket_meta()
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

    def _pack_transition(self, op: str, owner: str, bucket: str = "", arg: int = 0) -> Any:
        return self._eval_ctrl(
            _PACK_TRANSITION, 3, BUCKET_META, BUCKET_READY_QUEUE, BUCKET_PACK_ACTIVE,
            op, owner, bucket, arg, BUCKET_STATE.format(id=""),
        )

    def pull_ready_bucket(
        self, timeout: float = 0, *, owner: str | None = None
    ) -> dict[str, Any] | None:
        """Pop a sealed+idle bucket ready for tar packing."""
        self._require_client()
        if owner is not None:
            bucket_token = self._pack_transition("claim", owner)
            if bucket_token is None:
                return None
        elif timeout and timeout > 0:
            raw = self._blpop_many([BUCKET_READY_QUEUE], timeout=timeout)
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
        meta = self._bucket_meta()
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

    def sample_pack_pressure(self) -> dict[str, str]:
        return {
            _redis_text(key): _redis_text(value)
            for key, value in self._ctrl().hgetall(BUCKET_META).items()
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
                out[_redis_text(key)] = int(value)
            except (TypeError, ValueError):
                out[_redis_text(key)] = value
        return out
