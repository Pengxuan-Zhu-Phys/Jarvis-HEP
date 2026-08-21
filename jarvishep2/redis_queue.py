#!/usr/bin/env python3
"""Thin Redis access layer for Jarvis-HEP V2 distributed runtime.

D25.4: ``RedisQueue`` is the public connection façade. Keyspaces live on
private mixins (``_TaskBroker``, ``_CalcPool``, ``_SampleBuckets``,
``_ControlAndHeartbeat``). Callers keep importing ``RedisQueue``.
"""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from typing import Any, Mapping, Sequence
from uuid import UUID

import numpy as np

# --- Key namespace (DESIGN §7, redis_queue.md §2) ---

TASK_QUEUE = "hep:task_queue"
CALC_FREE = "calc:free:{name}"
CALC_BUSY_PACKS = "calc:busy:{name}"
CALC_SHARED_FREE = "calc:free:{name}:mode:{mode}"
CALC_SHARED_UNASSIGNED = "calc:free:{name}:unassigned"
CALC_SHARED_PACK_MODE = "calc:packmode:{name}"
SHARED_HELD_PREFIX = "@shared:"
RESULTS = "hep:results:{uuid}"
ARCHIVE_QUEUE = "hep:archive_queue"
ARCHIVED_UUIDS = "hep:archived:{scan}"
# Largest durable contiguous logical sample prefix [0, P).  Unlike
# ARCHIVED_UUIDS this is O(1) to read from the control heartbeat.
ARCHIVED_INDEX_PREFIX = "hep:archived-prefix:{scan}"
FEEDBACK_QUEUE = "hep:feedback"
# Per-chain feedback shards for independent MCMC pipelines. Task submit stays
# on the single ``hep:task_queue``; workers route by task.feedback_queue.
CHAIN_FEEDBACK_QUEUE = "hep:feedback:chain:{chain_id}"
CHAIN_FEEDBACK_QUEUE_PATTERN = "hep:feedback:chain:*"
WORKER_STATUS = "hep:worker:status:{id}"
CALC_STATUS = "hep:calculator:status"
SAMPLE_STATS = "hep:sample:stats"
OP_COUNT = "hep:{kind}:op_count"
# SAMPLE bucket allocator (limit samples per bucket; seal → ready for tar pack).
BUCKET_META = "hep:sample:bucket:meta"
BUCKET_STATE = "hep:sample:bucket:state:{id}"
BUCKET_READY_QUEUE = "hep:sample:bucket:ready"
BUCKET_LOCK = "hep:sample:bucket:lock"
# Single control-process lease for one Redis DB (hard guard against stacked Jarvis runs).
CONTROL_LOCK = "hep:control:lock"
CONTROL_LOCK_TTL_SEC = 30
_VALID_OP_KINDS = frozenset({"worker", "calculator", "sample", "task"})
_VALID_SAMPLE_ARTIFACTS = frozenset({"auto", "always", "never"})
_VALID_RESULT_STATUSES = frozenset({"Created", "Init", "Running", "Completed", "Failed"})
_VALID_EXECUTION_STEP_TYPES = frozenset(
    {"calculator", "opera", "likelihood", "nuisance_optimize"}
)

# The distributed broker is a V2 implementation detail, not task-YAML input.
INTERNAL_REDIS_CONFIG = {"host": "127.0.0.1", "port": 6379, "db": 0}
RUNTIME_METADATA_PATH = "hep:runtime:metadata_path"


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


def calc_shared_free_list_key(name: str, mode: str) -> str:
    """Redis free list for PackIDs currently built as one shared-pack mode."""
    return CALC_SHARED_FREE.format(name=name, mode=mode)


def calc_shared_unassigned_list_key(name: str) -> str:
    """Redis free list for shared packs with no trustworthy installed mode."""
    return CALC_SHARED_UNASSIGNED.format(name=name)


def calc_shared_pack_mode_key(name: str) -> str:
    """Redis hash of PackID -> successfully installed shared-pack mode."""
    return CALC_SHARED_PACK_MODE.format(name=name)


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


def _redis_text(value: Any) -> str:
    """Normalize redis-py decode-responses and byte-client return values."""
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


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


def chain_feedback_queue(chain_id: int) -> str:
    """Redis list key for one independent MCMC chain's feedback stream."""
    return CHAIN_FEEDBACK_QUEUE.format(chain_id=int(chain_id))


def _normalize_feedback_queue_name(queue: str | None) -> str:
    text = str(queue or "").strip()
    return text if text else FEEDBACK_QUEUE


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
    if "feedback_queue" in task and task.get("feedback_queue") is not None:
        fbq = str(task.get("feedback_queue") or "").strip()
        if not fbq:
            raise TaskValidationError("task feedback_queue must be non-empty when set")
    if "chain_id" in task and task.get("chain_id") is not None:
        try:
            int(task.get("chain_id"))
        except (TypeError, ValueError) as exc:
            raise TaskValidationError("task chain_id must be an integer") from exc
    for key in ("step", "stage"):
        if key in task and task.get(key) is not None:
            try:
                int(task.get(key))
            except (TypeError, ValueError) as exc:
                raise TaskValidationError(f"task {key} must be an integer") from exc


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


from jarvishep2._redis_calc_pool import _CalcPool
from jarvishep2._redis_control import _ControlAndHeartbeat
from jarvishep2._redis_sample_buckets import _SampleBuckets
from jarvishep2._redis_task_broker import _TaskBroker


class RedisQueue(_TaskBroker, _CalcPool, _SampleBuckets, _ControlAndHeartbeat):
    """Redis broker for tasks, calculator pools, results, and monitor counters.

    Internals (D25.4): task/archive/feedback, calculator pools, SAMPLE buckets,
    and the control lock + worker heartbeats are mixin keyspaces. This class
    owns the connection, codec, and cross-keyspace reset/monitor helpers.
    """

    # redis-py 5+/8 default socket_timeout=5 races with BLPOP(timeout<=5) used
    # throughout the control/worker loops (feedback barrier, task pull, archive
    # drain). Blocking commands need socket_timeout=None so server-side BLPOP
    # timeouts return None cleanly instead of raising TimeoutError mid-wait.
    _SOCKET_CONNECT_TIMEOUT_SEC = 5.0

    def __init__(self, config: Mapping[str, Any] | None = None, *, client: Any = None) -> None:
        self.config = dict(config or {})
        self._codec = str(self.config.get("codec", "json")).strip().lower()
        self.r = client

    def _client_kwargs(self) -> dict[str, Any]:
        """Shared redis-py kwargs for long-lived blocking BLPOP clients."""
        return {
            "decode_responses": self._codec == "json",
            # Override redis-py 8 default (5s) — must exceed any BLPOP wait.
            "socket_timeout": None,
            "socket_connect_timeout": float(
                self.config.get(
                    "socket_connect_timeout", self._SOCKET_CONNECT_TIMEOUT_SEC
                )
            ),
        }

    def connect(self) -> None:
        """Build a redis client from config when one was not injected."""
        if self.r is not None:
            return

        import redis

        kwargs = self._client_kwargs()
        url = self.config.get("url")
        if url:
            self.r = redis.Redis.from_url(str(url), **kwargs)
            return

        self.r = redis.Redis(
            host=str(self.config.get("host", "localhost")),
            port=int(self.config.get("port", 6379)),
            db=int(self.config.get("db", 0)),
            **kwargs,
        )

    def _blpop(self, key: str, *, timeout: int = 1) -> Any | None:
        """BLPOP wrapper that treats client socket timeouts as empty pops.

        With ``socket_timeout=None`` this should not fire; kept as a defensive
        guard so a future redis-py default change cannot crash generation
        barriers (AdaptiveBridson / FeedbackSampler wait loops).
        """
        self._require_client()
        try:
            return self.r.blpop(key, timeout=max(0, int(timeout)))
        except Exception as exc:
            # redis.exceptions.TimeoutError (and socket.timeout wrappers).
            name = type(exc).__name__
            if name in {"TimeoutError", "Timeout"} or "timeout" in str(exc).lower():
                return None
            raise

    def _blpop_many(self, keys: list[str], *, timeout: float = 1) -> Any | None:
        """Blocking pop from several candidate lists, preserving Redis affinity order."""
        if not keys:
            return None
        try:
            return self.r.blpop(keys, timeout=max(0.0, float(timeout)))
        except Exception as exc:
            name = type(exc).__name__
            if name in {"TimeoutError", "Timeout"} or "timeout" in str(exc).lower():
                return None
            raise

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
        keys.extend(self.r.scan_iter(match=CHAIN_FEEDBACK_QUEUE_PATTERN))
        for kind in sorted(_VALID_OP_KINDS):
            keys.append(OP_COUNT.format(kind=kind))
        for name in calculator_names or []:
            text = str(name or "").strip()
            if not text:
                continue
            keys.append(calc_free_list_key(text))
            keys.append(calc_busy_packs_key(text))
            keys.append(calc_shared_pack_mode_key(text))
            keys.extend(self.r.scan_iter(match=f"calc:free:{text}:*"))
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

    def reconcile_resume_ephemeral(
        self,
        *,
        completed: int,
        failed: int = 0,
    ) -> dict[str, int]:
        """Discard crash-era transport state and seed counters from DATABASE."""
        self._require_client()
        keys: list[str] = [TASK_QUEUE, ARCHIVE_QUEUE, FEEDBACK_QUEUE]
        keys.extend(self.r.scan_iter(match=CHAIN_FEEDBACK_QUEUE_PATTERN))
        keys.extend(self.r.scan_iter(match="hep:worker:status:*"))
        deleted = int(self.r.delete(*keys) or 0) if keys else 0
        self.r.hset(
            SAMPLE_STATS,
            mapping={
                "completed": max(0, int(completed)),
                "failed": max(0, int(failed)),
                "running": 0,
            },
        )
        return {
            "deleted_keys": deleted,
            "completed": max(0, int(completed)),
            "failed": max(0, int(failed)),
        }

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

    def set_runtime_metadata_path(self, path: str) -> None:
        self._require_client()
        self.r.set(RUNTIME_METADATA_PATH, str(path))

    def get_runtime_metadata_path(self) -> str | None:
        self._require_client()
        value = self.r.get(RUNTIME_METADATA_PATH)
        return str(value) if value else None

    def clear_runtime_metadata_path(self) -> None:
        self._require_client()
        self.r.delete(RUNTIME_METADATA_PATH)

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
    "ARCHIVED_UUIDS",
    "BUCKET_META",
    "BUCKET_READY_QUEUE",
    "BUCKET_STATE",
    "CALC_BUSY_PACKS",
    "CALC_FREE",
    "CALC_SHARED_FREE",
    "CALC_SHARED_PACK_MODE",
    "CALC_SHARED_UNASSIGNED",
    "CALC_STATUS",
    "CHAIN_FEEDBACK_QUEUE",
    "CHAIN_FEEDBACK_QUEUE_PATTERN",
    "CONTROL_LOCK",
    "CONTROL_LOCK_TTL_SEC",
    "CodecError",
    "FEEDBACK_QUEUE",
    "OP_COUNT",
    "RESULTS",
    "RedisQueue",
    "SAMPLE_STATS",
    "SHARED_HELD_PREFIX",
    "TASK_QUEUE",
    "TaskValidationError",
    "WORKER_STATUS",
    "calc_busy_packs_key",
    "calc_free_list_key",
    "calc_shared_free_list_key",
    "calc_shared_pack_mode_key",
    "calc_shared_unassigned_list_key",
    "calc_status_busy_field",
    "calc_status_free_field",
    "chain_feedback_queue",
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
