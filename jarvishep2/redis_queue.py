#!/usr/bin/env python3
"""Thin Redis access layer for Jarvis-HEP V2 distributed runtime."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
import json
from typing import Any, Mapping
from uuid import UUID, uuid4

import numpy as np

# --- Key namespace (DESIGN §7) ---

TASK_QUEUE = "hep:task_queue"
CALC_FREE = "calc:free:{name}"
RESULTS = "hep:results:{uuid}"
ARCHIVE_QUEUE = "hep:archive_queue"
WORKER_STATUS = "hep:worker:status:{id}"
CALC_STATUS = "hep:calculator:status"
SAMPLE_STATS = "hep:sample:stats"
OP_COUNT = "hep:{kind}:op_count"

_CALC_SLOT_TOKEN = "ready"
_VALID_OP_KINDS = frozenset({"worker", "calculator", "sample", "task"})


class CodecError(RuntimeError):
    """Raised when payload encoding or decoding fails."""


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
        key = CALC_FREE.format(name=name)
        pipe = self.r.pipeline(transaction=True)
        pipe.delete(key)
        for _ in range(n):
            pipe.rpush(key, _CALC_SLOT_TOKEN)
        self._set_calc_counts(pipe, name, free=n, busy=0)
        pipe.execute()

    def acquire_calc(self, name: str, timeout: int = 30) -> str | None:
        self._require_client()
        key = CALC_FREE.format(name=name)
        raw = self.r.blpop(key, timeout=timeout)
        if raw is None:
            return None

        pack_id = str(uuid4())
        pipe = self.r.pipeline(transaction=True)
        self._adjust_calc_counts(pipe, name, free_delta=-1, busy_delta=1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
        pipe.execute()
        return pack_id

    def release_calc(self, name: str, pack_id: str | None = None) -> None:
        del pack_id  # pack_id is traceability-only; slots are fungible
        self._require_client()
        key = CALC_FREE.format(name=name)
        pipe = self.r.pipeline(transaction=True)
        pipe.rpush(key, _CALC_SLOT_TOKEN)
        self._adjust_calc_counts(pipe, name, free_delta=1, busy_delta=-1)
        pipe.incr(OP_COUNT.format(kind="calculator"))
        pipe.execute()

    def submit_result(self, info: Mapping[str, Any]) -> None:
        self._require_client()
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
        mapping = {k: str(v) for k, v in fields.items()}
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
            "calculator_status": dict(calc_status),
            "sample_stats": dict(sample_stats),
            "op_counts": {kind: self.get_op_count(kind) for kind in sorted(_VALID_OP_KINDS)},
        }

    def store_result_hash(self, uuid: str, info: Mapping[str, Any]) -> None:
        """Optional per-uuid result hash for debug/sync handoff."""
        self._require_client()
        key = RESULTS.format(uuid=uuid)
        encoded = encode_payload(dict(info), codec=self._codec)
        self.r.hset(key, mapping={"payload": encoded})

    def _require_client(self) -> None:
        if self.r is None:
            raise RuntimeError("Redis client is not connected; call connect() or inject client")

    def _calc_free_key(self, name: str) -> str:
        return f"{name}:free"

    def _calc_busy_key(self, name: str) -> str:
        return f"{name}:busy"

    def _set_calc_counts(self, pipe: Any, name: str, *, free: int, busy: int) -> None:
        pipe.hset(
            CALC_STATUS,
            mapping={
                self._calc_free_key(name): max(0, int(free)),
                self._calc_busy_key(name): max(0, int(busy)),
            },
        )

    def _adjust_calc_counts(
        self,
        pipe: Any,
        name: str,
        *,
        free_delta: int,
        busy_delta: int,
    ) -> None:
        current = self.r.hgetall(CALC_STATUS) or {}
        free = int(current.get(self._calc_free_key(name), 0)) + free_delta
        busy = int(current.get(self._calc_busy_key(name), 0)) + busy_delta
        self._set_calc_counts(pipe, name, free=free, busy=busy)


def make_fakeredis_queue(**config: Any) -> RedisQueue:
    """Construct a RedisQueue backed by fakeredis (for unit tests)."""
    import fakeredis

    codec = str(config.get("codec", "json")).strip().lower()
    client = fakeredis.FakeStrictRedis(decode_responses=codec == "json")
    return RedisQueue(config, client=client)