"""Small sampler-only pickle snapshots; no queues, factory or output database."""

from __future__ import annotations

import hashlib
import json
import os
import pickle
import tempfile
from collections.abc import Mapping
from typing import Any

# Terminal status only. AdaptiveBridson, Dynesty and MultiNest use the same
# pickle files but still need observable feedback (f, or the nested logL).
COMPLETION_METHODS = frozenset({"grid", "random", "bridson", "csv"})
METHODS = COMPLETION_METHODS | frozenset({"adaptivebridson", "dynesty", "multinest"})
FORMAT = "jarvis-hep.sampler-pickle"
VERSION = 1


def _method_name(config: Mapping[str, Any]) -> str:
    return str((config.get("Sampling") or {}).get("Method", "")).lower()


def uses_sample_checkpoint(config: Mapping[str, Any]) -> bool:
    return _method_name(config) in METHODS


def uses_completion_feedback(config: Mapping[str, Any]) -> bool:
    """True when workers should publish a terminal ack and nothing else."""
    return _method_name(config) in COMPLETION_METHODS


def config_identity(config: Mapping[str, Any]) -> str:
    """Compute once at setup, never on the sampling/checkpoint hot path."""
    sampling = dict(config.get("Sampling") or {})
    bounds = dict(sampling.get("Bounds") or {})
    bounds.pop("max_worker", None)  # Concurrency does not change the point source.
    sampling["Bounds"] = bounds
    data = {"sampling": sampling, "units": config.get("_unit_execution")}
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _atomic_replace(path: str, write: Any, *, prefix: str) -> None:
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=prefix, suffix=".tmp", dir=directory)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            write(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        if os.name == "posix":
            directory_fd = os.open(directory, os.O_RDONLY)
            try:
                os.fsync(directory_fd)
            finally:
                os.close(directory_fd)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def atomic_pickle(path: str, payload: Mapping[str, Any]) -> None:
    _atomic_replace(
        path,
        lambda handle: pickle.dump(dict(payload), handle, protocol=pickle.HIGHEST_PROTOCOL),
        prefix=".pickle-",
    )


def atomic_write_bytes(path: str, data: bytes) -> None:
    """Replace *path* with *data*; fsync the file and its directory."""
    _atomic_replace(path, lambda handle: handle.write(data), prefix=f".{os.path.basename(path)}-")


def read_pickle(path: str) -> dict[str, Any]:
    with open(path, "rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"{path}: expected a pickle mapping")
    return payload


def file_fingerprint(path: str) -> dict[str, Any]:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(4 * 1024 * 1024):
            digest.update(chunk)
    return {"size": os.path.getsize(path), "sha256": digest.hexdigest()}


def load_runtime(path: str, config: Mapping[str, Any]) -> dict[str, Any]:
    """Only the new format is accepted; rejection never starts a fresh run."""
    payload = read_pickle(path)
    if payload.get("format") != FORMAT or payload.get("version") != VERSION:
        raise ValueError("unsupported sampler checkpoint format; use a fresh scan")
    if payload.get("config_id") != config_identity(config):
        raise ValueError("sampler checkpoint configuration differs from the current scan")
    if payload.get("method") != str((config.get("Sampling") or {}).get("Method", "")).lower():
        raise ValueError("sampler checkpoint method differs")
    if not isinstance(payload.get("sampler_state"), dict):
        raise ValueError("sampler checkpoint has no state")
    return payload


