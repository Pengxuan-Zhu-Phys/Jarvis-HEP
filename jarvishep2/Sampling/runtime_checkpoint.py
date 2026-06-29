#!/usr/bin/env python3
"""Distributed runtime checkpoint for Jarvis-HEP V2 (WP-D6.2)."""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import threading
from collections.abc import Mapping
from copy import deepcopy
from datetime import datetime, timezone
from typing import Any

import numpy as np

from jarvishep2.calculator_pools import register_calculator_pools
from jarvishep2.redis_queue import RedisQueue

CHECKPOINT_FORMAT = "jarvis-hep.v2-distributed"
CHECKPOINT_VERSION = 1

V1_CHECKPOINT_FORMAT = "jarvis-hep.statesaver"
V1_MCMc_FORMAT = "jarvis-hep.mcmc-runtime"
THROUGHPUT_CORE_FORMAT = "jarvis-hep.v2"

CHECKPOINT_HEARTBEAT_SEC = 30.0
RESUME_PROMPT = (
    "Detected checkpoint file. Re-run from scratch? [y/N] (default: resume in 30s): "
)

REFUSAL_MESSAGES = {
    V1_CHECKPOINT_FORMAT: (
        "Checkpoint format 'jarvis-hep.statesaver' is from Jarvis-HEP V1 and cannot be "
        "resumed under the distributed V2 runtime (jarvis-hep2). Start a fresh run."
    ),
    V1_MCMc_FORMAT: (
        "Checkpoint format 'jarvis-hep.mcmc-runtime' is from Jarvis-HEP V1 and cannot be "
        "resumed under the distributed V2 runtime (jarvis-hep2). Start a fresh run."
    ),
    THROUGHPUT_CORE_FORMAT: (
        "Checkpoint format 'jarvis-hep.v2' is from the retired throughput-core runtime and "
        "cannot be resumed under jarvis-hep2. Start a fresh distributed run."
    ),
}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_pickle_dump(path: str, payload: Mapping[str, Any]) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "wb") as handle:
        pickle.dump(dict(payload), handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, path)


def pickle_load(path: str) -> dict[str, Any]:
    with open(path, "rb") as handle:
        loaded = pickle.load(handle)
    if not isinstance(loaded, dict):
        raise ValueError("checkpoint payload must be a mapping")
    return loaded


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def stable_json_hash(payload: Any) -> str:
    blob = json.dumps(_json_safe(payload), sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


def serialize_seed_sequence(seq: np.random.SeedSequence) -> dict[str, Any]:
    state = seq.state
    return {
        "entropy": int(state["entropy"]),
        "spawn_key": [int(item) for item in state.get("spawn_key", ())],
        "pool_size": int(state.get("pool_size", 4)),
        "pool": [int(item) for item in state.get("pool", ())],
    }


def deserialize_seed_sequence(payload: Mapping[str, Any]) -> np.random.SeedSequence:
    spawn_key = tuple(int(item) for item in payload.get("spawn_key") or ())
    pool = tuple(int(item) for item in payload.get("pool") or ())
    kwargs: dict[str, Any] = {
        "entropy": int(payload.get("entropy", 0)),
        "spawn_key": spawn_key,
    }
    if pool:
        kwargs["pool"] = pool
    if payload.get("pool_size") is not None:
        kwargs["pool_size"] = int(payload["pool_size"])
    return np.random.SeedSequence(**kwargs)


def derive_sample_seed(master: np.random.SeedSequence, sample_index: int) -> np.random.SeedSequence:
    """Child stream for one Sample — independent of Worker count (WP-D6.2)."""
    index = max(0, int(sample_index))
    return master.spawn(1 + index)[index]


def checkpoint_path(
    *,
    task_root: str,
    scan_name: str,
    sampler_name: str,
) -> str:
    root = os.path.join(
        os.path.abspath(str(task_root)),
        "checkpoints",
        str(scan_name),
        str(sampler_name),
    )
    return os.path.join(root, "state.pkl")


def build_run_spec(
    *,
    config: Mapping[str, Any] | None,
    scan_name: str,
    task_root: str,
    task_result_dir: str,
    sampler_name: str,
    worker_parallel: int = 0,
) -> dict[str, Any]:
    return {
        "normalized_config": deepcopy(dict(config or {})),
        "scan_name": str(scan_name),
        "task_root": str(task_root),
        "task_result_dir": str(task_result_dir),
        "sampler_method": str(sampler_name),
        "worker_parallel": int(worker_parallel),
    }


def build_payload(
    *,
    run_spec: Mapping[str, Any],
    sampler_state: Mapping[str, Any],
    persistence: Mapping[str, Any] | None = None,
    reason: str = "",
    safe_barrier_confirmed: bool = True,
) -> dict[str, Any]:
    created_at = utc_now_iso()
    integrity = {
        "config_hash": stable_json_hash(run_spec.get("normalized_config", {})),
        "variable_signature": list(
            (run_spec.get("normalized_config") or {}).get("Sampling", {}).get("Variables", [])
            or []
        ),
        "safe_barrier_confirmed": bool(safe_barrier_confirmed),
        "checkpoint_reason": str(reason or ""),
    }
    return {
        "format": CHECKPOINT_FORMAT,
        "version": CHECKPOINT_VERSION,
        "created_at_utc": created_at,
        "timestamp_utc": created_at,
        "run_spec": dict(run_spec),
        "sampler_state": dict(sampler_state),
        "persistence": dict(persistence or {}),
        "integrity": integrity,
    }


def validate_checkpoint_payload(payload: Any) -> tuple[bool, str]:
    if not isinstance(payload, dict):
        return False, "checkpoint payload must be a mapping"
    payload_format = str(payload.get("format") or "")
    if payload_format in REFUSAL_MESSAGES:
        return False, REFUSAL_MESSAGES[payload_format]
    if payload_format and payload_format != CHECKPOINT_FORMAT:
        return False, f"unsupported checkpoint format: {payload_format!r}"
    version = payload.get("version")
    if version is not None and int(version) != CHECKPOINT_VERSION:
        return False, f"checkpoint version mismatch: expect {CHECKPOINT_VERSION}, got {version!r}"
    if not isinstance(payload.get("run_spec"), dict):
        return False, "checkpoint payload missing run_spec"
    if not isinstance(payload.get("sampler_state"), dict):
        return False, "checkpoint payload missing sampler_state"
    if "integrity" in payload and not isinstance(payload.get("integrity"), dict):
        return False, "checkpoint payload missing integrity"
    return True, "ok"


def load_checkpoint(path: str) -> dict[str, Any]:
    payload = pickle_load(path)
    ok, reason = validate_checkpoint_payload(payload)
    if not ok:
        raise ValueError(reason)
    return payload


def save_checkpoint(path: str, payload: Mapping[str, Any]) -> str:
    ok, reason = validate_checkpoint_payload(payload)
    if not ok:
        raise ValueError(reason)
    atomic_pickle_dump(path, payload)
    return path


def safe_barrier_ready(
    *,
    sampler_at_barrier: bool,
    submitted_uuids: set[str] | frozenset[str],
    archiver_persistence: Mapping[str, Any] | None,
) -> bool:
    if not sampler_at_barrier:
        return False
    acked = set(str(item) for item in (archiver_persistence or {}).get("acked_uuids") or [])
    if not submitted_uuids:
        return True
    return submitted_uuids <= acked


def prepare_resume(
    redis: RedisQueue,
    *,
    worker_config: Mapping[str, Any] | None,
) -> int:
    """Rebuild Redis pools and drain stale task queue (WP-D6.2)."""
    drained = redis.drain_task_queue()
    register_calculator_pools(redis, worker_config or {})
    return drained


class CheckpointHeartbeat:
    """30 s checkpoint heartbeat (frozen UX, WP-D6.2)."""

    def __init__(
        self,
        *,
        interval_sec: float = CHECKPOINT_HEARTBEAT_SEC,
        save_callback,
    ) -> None:
        self._interval_sec = max(1.0, float(interval_sec))
        self._save_callback = save_callback
        self._stop = threading.Event()
        self._thread: threading.Thread | None = None

    def start(self) -> None:
        if self._thread is not None and self._thread.is_alive():
            return
        self._stop.clear()
        self._thread = threading.Thread(
            target=self._loop,
            name="Jarvis2-CheckpointHeartbeat",
            daemon=True,
        )
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=2.0)
            self._thread = None

    def _loop(self) -> None:
        while not self._stop.wait(self._interval_sec):
            try:
                self._save_callback(reason="checkpoint_heartbeat")
            except Exception:
                pass


__all__ = [
    "CHECKPOINT_FORMAT",
    "CHECKPOINT_HEARTBEAT_SEC",
    "CHECKPOINT_VERSION",
    "CheckpointHeartbeat",
    "REFUSAL_MESSAGES",
    "RESUME_PROMPT",
    "build_payload",
    "build_run_spec",
    "checkpoint_path",
    "derive_sample_seed",
    "deserialize_seed_sequence",
    "load_checkpoint",
    "prepare_resume",
    "safe_barrier_ready",
    "save_checkpoint",
    "serialize_seed_sequence",
    "stable_json_hash",
    "validate_checkpoint_payload",
]