#!/usr/bin/env python3
"""Thin Redis access layer for Jarvis-HEP V2 distributed runtime."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
import os
import time
from typing import Any, Mapping
from uuid import UUID, uuid4

import numpy as np

# --- Key namespace (DESIGN §7, redis_queue.md §2) ---

TASK_QUEUE = "hep:task_queue"
CALC_FREE = "calc:free:{name}"
CALC_BUSY_PACKS = "calc:busy:{name}"
RESULTS = "hep:results:{uuid}"
ARCHIVE_QUEUE = "hep:archive_queue"
FEEDBACK_QUEUE = "hep:feedback"
WORKER_STATUS = "hep:worker:status:{id}"
CALC_STATUS = "hep:calculator:status"
SAMPLE_STATS = "hep:sample:stats"
OP_COUNT = "hep:{kind}:op_count"
# SAMPLE bucket allocator (limit samples per bucket; seal → ready for tar pack).
BUCKET_META = "hep:sample:bucket:meta"
BUCKET_STATE = "hep:sample:bucket:state:{id}"
BUCKET_READY_QUEUE = "hep:sample:bucket:ready"
BUCKET_LOCK = "hep:sample:bucket:lock"
# Single control-process lease for one Redis DB (hard guard against stacked Jarvis2 runs).
CONTROL_LOCK = "hep:control:lock"
CONTROL_LOCK_TTL_SEC = 120
_VALID_OP_KINDS = frozenset({"worker", "calculator", "sample", "task"})
_VALID_SAMPLE_ARTIFACTS = frozenset({"auto", "always", "never"})
_VALID_RESULT_STATUSES = frozenset({"Created", "Init", "Running", "Completed", "Failed"})
_VALID_EXECUTION_STEP_TYPES = frozenset(
    {"calculator", "opera", "likelihood", "nuisance_optimize"}
)
# The distributed broker is a V2 implementation detail, not task-YAML input.
INTERNAL_REDIS_CONFIG = {"host": "127.0.0.1", "port": 6379, "db": 0}


class CodecError(RuntimeError):
    """Raised when payload encoding or decoding fails."""


class TaskValidationError(ValueError):
    """Raised when a task or result payload fails schema validation."""


def calc_free_list_key(name: str) -> str:
    """Redis list key for calculator free slots."""
    return CALC_FREE.format(name=name)


def calc_busy_packs_key(name: str) -> str:
    """Redis hash key tracking active pack_id owners per calculator."""
    return CALC_BUSY_PACKS.format(name=name)


def format_calc_pack_id(index: int, *, slots: int) -> str:
    """V1-style zero-padded PackID (``001`` …) for a 1-based slot index.

    Width is at least 3, and grows with ``slots`` so pools larger than 999 stay
    sortable (``0001`` … ``1000``).
    """
    n = max(1, int(index))
    width = max(3, len(str(max(1, int(slots)))))
    return f"{n:0{width}d}"


def is_stable_calc_pack_id(pack_id: str | None) -> bool:
    """Return True for Redis-registered numeric PackIDs (``001``, ``16``, …).

    Rejects legacy free-list junk such as the pre-migration token ``ready``,
    empty strings, and UUID-shaped ids minted by older ``acquire_calc``.
    """
    text = str(pack_id or "").strip()
    if not text or not text.isdigit():
        return False
    return int(text) >= 1


def calc_status_free_field(name: str) -> str:
    """Hash field inside CALC_STATUS for free-slot count."""
    return f"{name}:free"


def calc_status_busy_field(name: str) -> str:
    """Hash field inside CALC_STATUS for busy-slot count."""
    return f"{name}:busy"


def _json_default(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, UUID):
        return str(obj)
    if is_dataclass(obj) and not isinstance(obj, type):
        return asdict(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def encode_payload(payload: Any, *, codec: str = "json") -> str | bytes:
    """Serialize a payload for Redis storage."""
    if codec == "json":
        return json.dumps(payload, default=_json_default, separators=(",", ":"))
    if codec == "msgpack":
        try:
            import msgpack
        except ImportError as exc:
            raise CodecError("msgpack codec requires the msgpack package") from exc
        return msgpack.packb(payload, default=_json_default, use_bin_type=True)
    raise CodecError(f"unsupported codec: {codec}")


def decode_payload(raw: str | bytes | None, *, codec: str = "json") -> Any:
    """Deserialize a payload from Redis."""
    if raw is None:
        return None
    if codec == "json":
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        return json.loads(raw)
    if codec == "msgpack":
        try:
            import msgpack
        except ImportError as exc:
            raise CodecError("msgpack codec requires the msgpack package") from exc
        return msgpack.unpackb(raw, raw=False)
    raise CodecError(f"unsupported codec: {codec}")


def _encode_heartbeat_value(value: Any) -> str | int | float:
    if isinstance(value, (str, int, float)):
        return value
    if isinstance(value, bool):
        return int(value)
    return json.dumps(value, default=_json_default, separators=(",", ":"))


def _validate_u_coords(value: Any) -> None:
    if value is None:
        return
    if isinstance(value, np.ndarray):
        return
    if isinstance(value, (list, tuple)):
        return
    raise TaskValidationError("u_coords must be a list, tuple, or numpy array")


def _validate_execution_plan(plan: Any) -> None:
    if plan is None:
        return
    if not isinstance(plan, list):
        raise TaskValidationError("execution_plan must be a list")
    for index, step in enumerate(plan):
        if not isinstance(step, Mapping):
            raise TaskValidationError(f"execution_plan[{index}] must be a mapping")
        step_type = str(step.get("type", "")).strip()
        if step_type not in _VALID_EXECUTION_STEP_TYPES:
            allowed = ", ".join(sorted(_VALID_EXECUTION_STEP_TYPES))
            raise TaskValidationError(
                f"execution_plan[{index}].type '{step_type}' is invalid; allowed: {allowed}"
            )
        if "layer" not in step:
            raise TaskValidationError(f"execution_plan[{index}] requires 'layer'")


def _validate_task_payload(task: Mapping[str, Any]) -> None:
    """Validate a lightweight task dict before it enters hep:task_queue."""
    if not isinstance(task, Mapping):
        raise TaskValidationError("task payload must be a mapping")
    uuid = task.get("uuid")
    if not uuid or not str(uuid).strip():
        raise TaskValidationError("task payload requires non-empty 'uuid'")
    has_coords = "u_coords" in task and task.get("u_coords") is not None
    has_plan = bool(task.get("execution_plan"))
    if not has_coords and not has_plan:
        raise TaskValidationError("task payload requires 'u_coords' and/or 'execution_plan'")
    if has_coords:
        _validate_u_coords(task.get("u_coords"))
    if has_plan:
        _validate_execution_plan(task.get("execution_plan"))
    if "sample_artifacts" in task:
        mode = str(task["sample_artifacts"]).strip().lower()
        if mode not in _VALID_SAMPLE_ARTIFACTS:
            raise TaskValidationError(
                f"sample_artifacts must be one of {sorted(_VALID_SAMPLE_ARTIFACTS)}"
            )


def _validate_result_payload(info: Mapping[str, Any]) -> None:
    """Validate a result/info dict before it enters the archive queue."""
    if not isinstance(info, Mapping):
        raise TaskValidationError("result payload must be a mapping")
    uuid = info.get("uuid")
    if not uuid or not str(uuid).strip():
        raise TaskValidationError("result payload requires non-empty 'uuid'")
    status = info.get("status")
    if status is not None and str(status) not in _VALID_RESULT_STATUSES:
        raise TaskValidationError(
            f"result status '{status}' is invalid; allowed: {sorted(_VALID_RESULT_STATUSES)}"
        )
    observables = info.get("observables")
    if observables is not None and not isinstance(observables, Mapping):
        raise TaskValidationError("result observables must be a mapping when provided")


class RedisQueue:
    """Redis broker for tasks, calculator pools, results, and monitor counters."""

    def __init__(self, config: Mapping[str, Any] | None = None, *, client: Any = None) -> None:
        self.config = dict(config or {})
        self._codec = str(self.config.get("codec", "json")).strip().lower()
        self.r = client

    def connect(self) -> None:
        """Build a redis client from config when one was not injected."""
        if self.r is not None:
            return

        url = self.config.get("url")
        if url:
            import redis

            self.r = redis.Redis.from_url(str(url), decode_responses=self._codec == "json")
            return

        import redis

        self.r = redis.Redis(
            host=str(self.config.get("host", "localhost")),
            port=int(self.config.get("port", 6379)),
            db=int(self.config.get("db", 0)),
            decode_responses=self._codec == "json",
        )

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
        raw = self.r.blpop(TASK_QUEUE, timeout=timeout)
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
            raw = self.r.blpop(TASK_QUEUE, timeout=1)
            if raw is None:
                break
            drained += 1
        return drained

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

    def register_calc_pool(self, name: str, n: int) -> None:
        """Register *n* exclusive PackID slots owned by Redis.

        Lifecycle (one sample / one worker per slot at a time)::

            free list  calc:free:<name>  = [001, 002, …, N]
            acquire    BLPOP free → busy[pack]=running   (exclusive)
            release    HDEL busy  → RPUSH free           (reuse same PackID)

        Always rebuilds the pool from scratch so leftover ``ready`` / UUID
        tokens from older builds cannot re-enter the free list.
        """
        self._require_client()
        if n <= 0:
            return
        slots = int(n)
        pool_key = calc_free_list_key(name)
        busy_key = calc_busy_packs_key(name)
        pack_ids = [format_calc_pack_id(i, slots=slots) for i in range(1, slots + 1)]
        pipe = self.r.pipeline(transaction=True)
        # Hard reset: wipe free + busy so no UUID/ready residue survives.
        pipe.delete(pool_key)
        pipe.delete(busy_key)
        for pack_id in pack_ids:
            pipe.rpush(pool_key, pack_id)
        pipe.hset(
            CALC_STATUS,
            mapping={
                calc_status_free_field(name): slots,
                calc_status_busy_field(name): 0,
            },
        )
        pipe.execute()

    def _discard_stale_free_token(self, name: str, token: str) -> None:
        """Drop a non-PackID free-list entry (e.g. legacy ``ready``/UUID) permanently."""
        free_field = calc_status_free_field(name)
        try:
            current = int(self.r.hget(CALC_STATUS, free_field) or 0)
        except (TypeError, ValueError):
            current = 0
        if current > 0:
            self.r.hincrby(CALC_STATUS, free_field, -1)
        # Never RPUSH: stale tokens must not re-enter the pool.

    def acquire_calc(self, name: str, timeout: int = 30) -> str | None:
        """Claim one free PackID exclusively (blocks until a slot or timeout).

        Only stable numeric PackIDs (``001`` …) are returned. Junk free-list
        values (``ready``, UUIDs, empty) are discarded and the wait continues.
        Returns ``None`` only when *timeout* elapses with no valid free slot.
        """
        self._require_client()
        pool_key = calc_free_list_key(name)
        busy_key = calc_busy_packs_key(name)
        deadline = time.monotonic() + max(0.0, float(timeout))

        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return None
            # redis-py BLPOP timeout is whole seconds; keep waiting in chunks.
            wait_sec = max(1, int(remaining) if remaining >= 1 else 1)
            raw = self.r.blpop(pool_key, timeout=wait_sec)
            if raw is None:
                if time.monotonic() >= deadline:
                    return None
                continue

            pack_id = str(raw[1]).strip()
            if not is_stable_calc_pack_id(pack_id):
                # Legacy free-list pollution — drop and wait for a real slot.
                self._discard_stale_free_token(name, pack_id)
                continue

            # free → running (exclusive ownership for this worker/sample).
            pipe = self.r.pipeline(transaction=True)
            pipe.hset(busy_key, pack_id, "running")
            pipe.hincrby(CALC_STATUS, calc_status_free_field(name), -1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), 1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return pack_id

    def release_calc(self, name: str, pack_id: str) -> None:
        """Return a PackID to free after the sample finishes (running → free)."""
        if not pack_id or not str(pack_id).strip():
            raise ValueError("pack_id is required for release_calc")
        pack_id = str(pack_id).strip()
        self._require_client()

        # Never put non-stable ids back into the free list (closes the UUID loop).
        if not is_stable_calc_pack_id(pack_id):
            busy_key = calc_busy_packs_key(name)
            # Best-effort cleanup if an old worker still held a UUID slot.
            self.r.hdel(busy_key, pack_id)
            return

        busy_key = calc_busy_packs_key(name)
        # Atomic guard: HDEL's return value decides the single winner when a
        # Worker's finally-release races the watchdog sweep for the same slot.
        if not self.r.hdel(busy_key, pack_id):
            raise ValueError(f"unknown pack_id '{pack_id}' for calculator '{name}'")

        pool_key = calc_free_list_key(name)
        pipe = self.r.pipeline(transaction=True)
        # running → free: same PackID re-enters the pool for the next sample.
        pipe.rpush(pool_key, pack_id)
        pipe.hincrby(CALC_STATUS, calc_status_free_field(name), 1)
        pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), -1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
        pipe.execute()

    def force_release_calc(self, name: str, pack_id: str) -> bool:
        """Best-effort PackID return for failure/cleanup paths (never raises).

        Returns True when this call transitioned the slot from busy → free.
        Safe under double-release races (second caller gets False).
        """
        try:
            if not pack_id or not str(pack_id).strip():
                return False
            pack_id = str(pack_id).strip()
            name = str(name or "").strip()
            if not name:
                return False
            self._require_client()
            busy_key = calc_busy_packs_key(name)
            if not is_stable_calc_pack_id(pack_id):
                self.r.hdel(busy_key, pack_id)
                return True
            if not self.r.hdel(busy_key, pack_id):
                return False
            pool_key = calc_free_list_key(name)
            pipe = self.r.pipeline(transaction=True)
            pipe.rpush(pool_key, pack_id)
            pipe.hincrby(CALC_STATUS, calc_status_free_field(name), 1)
            pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), -1)
            pipe.incr(OP_COUNT.format(kind="calculator"))
            pipe.execute()
            return True
        except Exception:
            return False

    def claim_control_lock(
        self,
        owner: str,
        *,
        ttl_sec: int = CONTROL_LOCK_TTL_SEC,
    ) -> bool:
        """Acquire exclusive control lease (SET NX EX). False if another owner holds it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        if not owner_text:
            raise ValueError("control lock owner is required")
        ttl = max(5, int(ttl_sec))
        return bool(self.r.set(CONTROL_LOCK, owner_text, nx=True, ex=ttl))

    def refresh_control_lock(
        self,
        owner: str,
        *,
        ttl_sec: int = CONTROL_LOCK_TTL_SEC,
    ) -> bool:
        """Extend the lease only if *owner* still holds it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        current = self.r.get(CONTROL_LOCK)
        if current is None:
            return bool(self.r.set(CONTROL_LOCK, owner_text, nx=True, ex=max(5, int(ttl_sec))))
        if str(current) != owner_text:
            return False
        self.r.expire(CONTROL_LOCK, max(5, int(ttl_sec)))
        return True

    def release_control_lock(self, owner: str) -> None:
        """Drop the control lease if we still own it."""
        self._require_client()
        owner_text = str(owner or "").strip()
        current = self.r.get(CONTROL_LOCK)
        if current is not None and str(current) == owner_text:
            self.r.delete(CONTROL_LOCK)

    def get_control_lock_owner(self) -> str | None:
        self._require_client()
        value = self.r.get(CONTROL_LOCK)
        return None if value is None else str(value)

    def reset_run_ephemeral_keys(
        self,
        *,
        calculator_names: list[str] | None = None,
        worker_ids: list[str] | int | None = None,
    ) -> dict[str, int]:
        """Hard-reset queues/stats/calc pools for a fresh run (not a full FLUSHDB).

        Calculator free/busy keys are deleted here; the caller must
        ``register_calc_pool`` again before Workers start.
        """
        self._require_client()
        keys: list[str] = [
            TASK_QUEUE,
            ARCHIVE_QUEUE,
            FEEDBACK_QUEUE,
            SAMPLE_STATS,
            CALC_STATUS,
            BUCKET_META,
            BUCKET_READY_QUEUE,
            BUCKET_LOCK,
        ]
        for kind in sorted(_VALID_OP_KINDS):
            keys.append(OP_COUNT.format(kind=kind))
        for name in calculator_names or []:
            text = str(name or "").strip()
            if not text:
                continue
            keys.append(calc_free_list_key(text))
            keys.append(calc_busy_packs_key(text))
        if worker_ids is None:
            # Best-effort: clear a reasonable worker status range.
            worker_ids = list(range(0, 64))
        elif isinstance(worker_ids, int):
            worker_ids = list(range(0, max(0, int(worker_ids))))
        for worker_id in worker_ids:
            keys.append(WORKER_STATUS.format(id=str(worker_id)))
        # Drop known SAMPLE bucket state hashes (current + a generous lookback window).
        try:
            meta = self.r.hgetall(BUCKET_META) or {}
            current = int(meta.get("current") or 0)
        except (TypeError, ValueError):
            current = 0
        for bucket_id in range(1, max(current, 0) + 64):
            keys.append(BUCKET_STATE.format(id=bucket_id))
        # Deduplicate while preserving order.
        seen: set[str] = set()
        unique_keys = []
        for key in keys:
            if key not in seen:
                seen.add(key)
                unique_keys.append(key)
        deleted = 0
        if unique_keys:
            deleted = int(self.r.delete(*unique_keys) or 0)
        # Explicit zeroed sample stats so readers never see missing keys as stale.
        self.r.hset(
            SAMPLE_STATS,
            mapping={"completed": 0, "failed": 0, "running": 0},
        )
        return {"deleted_keys": deleted, "sample_stats_reset": 1}

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
            raw = self.r.blpop(BUCKET_READY_QUEUE, timeout=timeout)
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

    def publish_feedback(self, info: Mapping[str, Any]) -> None:
        """Push a flat sampler barrier record (D13.8).

        Expected shape::

            {"uuid": str, "logL": float|str, ...optional flat fields}

        Nested ``observables`` / sample ``status`` are not part of the contract;
        if a legacy caller still passes them they are dropped (not re-nested).
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
        self.r.rpush(FEEDBACK_QUEUE, encoded)

    def pull_feedback(self, timeout: int = 1) -> dict[str, Any] | None:
        """Blocking pop from hep:feedback; normalizes flat floats."""
        self._require_client()
        raw = self.r.blpop(FEEDBACK_QUEUE, timeout=timeout)
        if raw is None:
            return None
        _, payload = raw
        decoded = decode_payload(payload, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("feedback payload must decode to a dict")
        from jarvishep2.feedback_return import normalize_feedback_record

        return normalize_feedback_record(decoded)

    def drain_feedback_queue(self) -> int:
        """Discard all queued feedback records (resume path)."""
        self._require_client()
        drained = 0
        while int(self.r.llen(FEEDBACK_QUEUE)) > 0:
            item = self.pull_feedback(timeout=1)
            if item is None:
                break
            drained += 1
        return drained

    def pull_result(self, timeout: int = 1) -> dict[str, Any] | None:
        self._require_client()
        raw = self.r.blpop(ARCHIVE_QUEUE, timeout=timeout)
        if raw is None:
            return None
        _, payload = raw
        decoded = decode_payload(payload, codec=self._codec)
        if not isinstance(decoded, dict):
            raise CodecError("result payload must decode to a dict")
        return decoded

    def heartbeat(self, worker_id: str, **fields: Any) -> None:
        self._require_client()
        if "last_heartbeat" not in fields and "ts" in fields:
            fields["last_heartbeat"] = fields["ts"]
        key = WORKER_STATUS.format(id=worker_id)
        mapping = {k: _encode_heartbeat_value(v) for k, v in fields.items()}
        pipe = self.r.pipeline(transaction=True)
        if mapping:
            pipe.hset(key, mapping=mapping)
        pipe.incr(OP_COUNT.format(kind="worker"))
        pipe.execute()

    def get_op_count(self, kind: str) -> int:
        self._require_client()
        if kind not in _VALID_OP_KINDS:
            raise ValueError(f"invalid op_count kind: {kind}")
        value = self.r.get(OP_COUNT.format(kind=kind))
        return int(value or 0)

    def get_all_op_counts(self) -> dict[str, int]:
        """Return all subsystem op_count values in one pipeline round-trip."""
        self._require_client()
        kinds = sorted(_VALID_OP_KINDS)
        pipe = self.r.pipeline(transaction=False)
        for kind in kinds:
            pipe.get(OP_COUNT.format(kind=kind))
        values = pipe.execute()
        return {kind: int(value or 0) for kind, value in zip(kinds, values)}

    def get_queue_lengths(self) -> dict[str, int]:
        """Return task, archive, and feedback queue lengths in one pipeline."""
        self._require_client()
        pipe = self.r.pipeline(transaction=False)
        pipe.llen(TASK_QUEUE)
        pipe.llen(ARCHIVE_QUEUE)
        pipe.llen(FEEDBACK_QUEUE)
        task_len, archive_len, feedback_len = pipe.execute()
        return {
            "task_queue_length": int(task_len),
            "archive_queue_length": int(archive_len),
            "feedback_queue_length": int(feedback_len),
        }

    def fetch_calculator_status(self) -> dict[str, int | float | str]:
        """Read the calculator status hash (monitor subsystem fetch)."""
        self._require_client()
        calc_status = self.r.hgetall(CALC_STATUS) or {}
        return {key: _coerce_numeric(value) for key, value in calc_status.items()}

    def fetch_sample_stats(self) -> dict[str, int | float | str]:
        """Read the sample stats hash (monitor subsystem fetch)."""
        self._require_client()
        sample_stats = self.r.hgetall(SAMPLE_STATS) or {}
        result = {key: _coerce_numeric(value) for key, value in sample_stats.items()}
        # Clamp: submit_result may outpace pull_task in unit tests / crash paths.
        if "running" in result:
            try:
                result["running"] = max(0, int(result["running"] or 0))
            except (TypeError, ValueError):
                result["running"] = 0
        return result

    def fetch_worker_status(self, worker_ids: list[str]) -> dict[str, dict[str, Any]]:
        """Read heartbeat hashes for the given Worker ids."""
        self._require_client()
        if not worker_ids:
            return {}
        pipe = self.r.pipeline(transaction=False)
        for worker_id in worker_ids:
            pipe.hgetall(WORKER_STATUS.format(id=worker_id))
        rows = pipe.execute()
        return {
            str(worker_id): dict(row or {})
            for worker_id, row in zip(worker_ids, rows)
        }

    def sweep_held_calc_slots(self, held_packs: Mapping[str, Any]) -> int:
        """Release calculator slots recorded for a dead Worker (WP-D6.1)."""
        released = 0
        for name, pack_id in dict(held_packs or {}).items():
            if not pack_id or not str(pack_id).strip():
                continue
            try:
                self.release_calc(str(name), str(pack_id))
            except ValueError:
                continue
            released += 1
        return released

    def encode_task_for_heartbeat(self, task: Mapping[str, Any]) -> str:
        """Serialize an in-flight task for the Worker heartbeat hash."""
        return encode_payload(dict(task), codec=self._codec)

    def decode_heartbeat_task(self, heartbeat: Mapping[str, Any]) -> dict[str, Any] | None:
        """Decode a Worker heartbeat's serialized in-flight task payload."""
        raw = heartbeat.get("current_task")
        if raw is None or raw == "":
            return None
        if isinstance(raw, Mapping):
            return dict(raw)
        decoded = decode_payload(str(raw), codec=self._codec)
        return dict(decoded) if isinstance(decoded, dict) else None

    def decode_heartbeat_subprocess_pids(self, heartbeat: Mapping[str, Any]) -> list[int]:
        """Decode active calculator subprocess PIDs from a Worker heartbeat."""
        raw = heartbeat.get("active_subprocess_pids")
        if raw is None or raw == "":
            return []
        if isinstance(raw, (list, tuple)):
            decoded: Any = list(raw)
        else:
            try:
                decoded = json.loads(str(raw))
            except (TypeError, ValueError, json.JSONDecodeError):
                return []
        if not isinstance(decoded, list):
            return []
        pids: list[int] = []
        for item in decoded:
            try:
                pid = int(item)
            except (TypeError, ValueError):
                continue
            if pid > 0:
                pids.append(pid)
        return pids

    def decode_heartbeat_held_packs(self, heartbeat: Mapping[str, Any]) -> dict[str, str]:
        """Decode held calculator slots from a Worker heartbeat hash."""
        raw = heartbeat.get("held_calc_packs")
        if raw is None or raw == "":
            return {}
        if isinstance(raw, Mapping):
            return {str(key): str(value) for key, value in raw.items() if value}
        try:
            decoded = json.loads(str(raw))
        except (TypeError, ValueError, json.JSONDecodeError):
            return {}
        if not isinstance(decoded, dict):
            return {}
        return {str(key): str(value) for key, value in decoded.items() if value}

    def incr_op(self, kind: str, amount: int = 1) -> int:
        self._require_client()
        if kind not in _VALID_OP_KINDS:
            raise ValueError(f"invalid op_count kind: {kind}")
        return int(self.r.incrby(OP_COUNT.format(kind=kind), amount))

    def snapshot_raw(self) -> dict[str, Any]:
        self._require_client()
        calc_status = self.r.hgetall(CALC_STATUS) or {}
        sample_stats = self.r.hgetall(SAMPLE_STATS) or {}
        return {
            "task_queue_length": int(self.r.llen(TASK_QUEUE)),
            "archive_queue_length": int(self.r.llen(ARCHIVE_QUEUE)),
            "calculator_status": {k: _coerce_numeric(v) for k, v in calc_status.items()},
            "sample_stats": {k: _coerce_numeric(v) for k, v in sample_stats.items()},
            "op_counts": {kind: self.get_op_count(kind) for kind in sorted(_VALID_OP_KINDS)},
        }

    def store_result_hash(self, uuid: str, info: Mapping[str, Any]) -> None:
        """Optional per-uuid result hash for debug/sync handoff."""
        self._require_client()
        _validate_result_payload(info)
        key = RESULTS.format(uuid=uuid)
        encoded = encode_payload(dict(info), codec=self._codec)
        self.r.hset(key, mapping={"payload": encoded})

    def connection_config(self) -> dict[str, Any]:
        """Return picklable connection settings for spawn child processes."""
        return dict(self.config)

    def ping(self) -> bool:
        """Verify that the configured Redis server is reachable."""
        self._require_client()
        return bool(self.r.ping())

    @staticmethod
    def extract_connection_config(
        source: RedisQueue | Mapping[str, Any],
    ) -> dict[str, Any]:
        """Normalize a live queue or mapping into picklable connection settings."""
        if isinstance(source, RedisQueue):
            return source.connection_config()
        return dict(source)

    def close(self) -> None:
        """Close the underlying Redis client and release the handle."""
        client = self.r
        self.r = None
        if client is None:
            return
        closer = getattr(client, "close", None)
        if callable(closer):
            try:
                closer()
            except Exception:
                pass

    def _require_client(self) -> None:
        if self.r is None:
            raise RuntimeError("Redis client is not connected; call connect() or inject client")


def _coerce_numeric(value: Any) -> int | float | str:
    if isinstance(value, (int, float)):
        return value
    text = str(value)
    try:
        if "." in text:
            return float(text)
        return int(text)
    except ValueError:
        return text


__all__ = [
    "ARCHIVE_QUEUE",
    "BUCKET_META",
    "BUCKET_READY_QUEUE",
    "BUCKET_STATE",
    "CALC_BUSY_PACKS",
    "CALC_FREE",
    "CALC_STATUS",
    "CONTROL_LOCK",
    "CONTROL_LOCK_TTL_SEC",
    "CodecError",
    "OP_COUNT",
    "RESULTS",
    "RedisQueue",
    "SAMPLE_STATS",
    "TASK_QUEUE",
    "TaskValidationError",
    "WORKER_STATUS",
    "calc_busy_packs_key",
    "calc_free_list_key",
    "calc_status_busy_field",
    "calc_status_free_field",
    "decode_payload",
    "encode_payload",
    "format_calc_pack_id",
    "is_stable_calc_pack_id",
    "INTERNAL_REDIS_CONFIG",
    "make_fakeredis_queue",
]


def make_fakeredis_queue(**config: Any) -> RedisQueue:
    """Construct a RedisQueue backed by fakeredis (for unit tests)."""
    import fakeredis

    codec = str(config.get("codec", "json")).strip().lower()
    client = fakeredis.FakeStrictRedis(decode_responses=codec == "json")
    return RedisQueue(config, client=client)
