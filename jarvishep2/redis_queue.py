#!/usr/bin/env python3
"""Thin Redis access layer for Jarvis-HEP V2 distributed runtime."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from typing import Any, Mapping
from uuid import UUID, uuid4

import numpy as np

# --- Key namespace (DESIGN §7, redis_queue.md §2) ---

TASK_QUEUE = "hep:task_queue"
CALC_FREE = "calc:free:{name}"
CALC_BUSY_PACKS = "calc:busy:{name}"
RESULTS = "hep:results:{uuid}"
ARCHIVE_QUEUE = "hep:archive_queue"
WORKER_STATUS = "hep:worker:status:{id}"
CALC_STATUS = "hep:calculator:status"
SAMPLE_STATS = "hep:sample:stats"
OP_COUNT = "hep:{kind}:op_count"

_CALC_SLOT_TOKEN = "ready"
_VALID_OP_KINDS = frozenset({"worker", "calculator", "sample", "task"})
_VALID_SAMPLE_ARTIFACTS = frozenset({"auto", "always", "never"})
_VALID_RESULT_STATUSES = frozenset({"Created", "Init", "Running", "Completed", "Failed"})
_VALID_EXECUTION_STEP_TYPES = frozenset(
    {"calculator", "opera", "likelihood", "nuisance_optimize"}
)


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
        return decoded

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
        self._require_client()
        if n <= 0:
            return
        pool_key = calc_free_list_key(name)
        busy_key = calc_busy_packs_key(name)
        pipe = self.r.pipeline(transaction=True)
        pipe.delete(pool_key)
        pipe.delete(busy_key)
        for _ in range(n):
            pipe.rpush(pool_key, _CALC_SLOT_TOKEN)
        pipe.hset(
            CALC_STATUS,
            mapping={
                calc_status_free_field(name): n,
                calc_status_busy_field(name): 0,
            },
        )
        pipe.execute()

    def acquire_calc(self, name: str, timeout: int = 30) -> str | None:
        self._require_client()
        pool_key = calc_free_list_key(name)
        raw = self.r.blpop(pool_key, timeout=timeout)
        if raw is None:
            return None

        pack_id = str(uuid4())
        busy_key = calc_busy_packs_key(name)
        pipe = self.r.pipeline(transaction=True)
        pipe.hset(busy_key, pack_id, "active")
        pipe.hincrby(CALC_STATUS, calc_status_free_field(name), -1)
        pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), 1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
        pipe.execute()
        return pack_id

    def release_calc(self, name: str, pack_id: str) -> None:
        if not pack_id or not str(pack_id).strip():
            raise ValueError("pack_id is required for release_calc")
        self._require_client()
        busy_key = calc_busy_packs_key(name)
        if not self.r.hexists(busy_key, pack_id):
            raise ValueError(f"unknown pack_id '{pack_id}' for calculator '{name}'")

        pool_key = calc_free_list_key(name)
        pipe = self.r.pipeline(transaction=True)
        pipe.hdel(busy_key, pack_id)
        pipe.rpush(pool_key, _CALC_SLOT_TOKEN)
        pipe.hincrby(CALC_STATUS, calc_status_free_field(name), 1)
        pipe.hincrby(CALC_STATUS, calc_status_busy_field(name), -1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
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
        """Return task and archive queue lengths in one pipeline round-trip."""
        self._require_client()
        pipe = self.r.pipeline(transaction=False)
        pipe.llen(TASK_QUEUE)
        pipe.llen(ARCHIVE_QUEUE)
        task_len, archive_len = pipe.execute()
        return {
            "task_queue_length": int(task_len),
            "archive_queue_length": int(archive_len),
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
        return {key: _coerce_numeric(value) for key, value in sample_stats.items()}

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
    "CALC_BUSY_PACKS",
    "CALC_FREE",
    "CALC_STATUS",
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
    "make_fakeredis_queue",
]


def make_fakeredis_queue(**config: Any) -> RedisQueue:
    """Construct a RedisQueue backed by fakeredis (for unit tests)."""
    import fakeredis

    codec = str(config.get("codec", "json")).strip().lower()
    client = fakeredis.FakeStrictRedis(decode_responses=codec == "json")
    return RedisQueue(config, client=client)