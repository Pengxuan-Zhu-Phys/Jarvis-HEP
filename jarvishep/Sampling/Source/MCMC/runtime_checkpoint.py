#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import os
import pickle
from copy import deepcopy
from datetime import datetime, timezone
from typing import Any, Dict, Sequence, Tuple

import numpy as np

from jarvishep.versioning import get_runtime_version


STATE_SAVER_FORMAT = "jarvis-hep.statesaver"
STATE_SAVER_VERSION = 1
CHECKPOINT_FORMAT = "jarvis-hep.mcmc-runtime"
CHECKPOINT_VERSION = 1
VALID_STATE_FORMATS = {STATE_SAVER_FORMAT, CHECKPOINT_FORMAT}

_ENGINE_STATE_EXCLUDE = {
    "_population_getter",
    "_gradient_provider",
    "gradient_provider",
    "logger",
}


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_pickle_dump(path: str, payload: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = path + ".tmp"
    with open(tmp_path, "wb") as handle:
        pickle.dump(payload, handle, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp_path, path)


def pickle_load(path: str) -> Dict[str, Any]:
    with open(path, "rb") as handle:
        return pickle.load(handle)


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
    return value


def atomic_json_dump(path: str, payload: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(_json_safe(payload), handle, indent=2, sort_keys=True)
    os.replace(tmp_path, path)


def stable_json_hash(payload: Any) -> str:
    blob = json.dumps(_json_safe(payload), sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()


class StateSaver:
    """Single-file checkpoint owner for resume payloads."""

    def __init__(self, checkpoint_file: str, *, logger=None, legacy_files: Sequence[str] | None = None) -> None:
        self.checkpoint_file = os.path.abspath(checkpoint_file)
        self.logger = logger
        self.legacy_files = [os.path.abspath(path) for path in (legacy_files or []) if path]
        os.makedirs(os.path.dirname(self.checkpoint_file), exist_ok=True)

    def _candidate_paths(self) -> list[str]:
        paths = [self.checkpoint_file]
        for path in self.legacy_files:
            if path not in paths:
                paths.append(path)
        return paths

    def exists(self) -> bool:
        return any(os.path.exists(path) for path in self._candidate_paths())

    def load(self, default: Any = None) -> Any:
        for path in self._candidate_paths():
            if not os.path.exists(path):
                continue
            try:
                with open(path, "rb") as handle:
                    payload = pickle.load(handle)
            except Exception as exc:
                raise ValueError(f"StateSaver checkpoint load failed -> {path}: {exc}") from exc
            if self.logger is not None:
                self.logger.info(f"StateSaver checkpoint loaded -> {path}")
            return payload
        return default

    def save(self, payload: Dict[str, Any]) -> str:
        try:
            atomic_pickle_dump(self.checkpoint_file, payload)
        except Exception as exc:
            raise RuntimeError(f"StateSaver checkpoint save failed -> {self.checkpoint_file}: {exc}") from exc
        if self.logger is not None:
            self.logger.info(f"StateSaver checkpoint saved -> {self.checkpoint_file}")
        return self.checkpoint_file


def build_sampler_signature(sampler) -> Dict[str, Any]:
    variable_names = []
    for var in getattr(sampler, "vars", ()) or ():
        name = getattr(var, "name", getattr(var, "_name", None))
        if name is not None:
            variable_names.append(str(name))

    sample_cfg = getattr(sampler, "info", {}) or {}
    sample_cfg = sample_cfg.get("sample", {}) if isinstance(sample_cfg, dict) else {}

    return {
        "method": str(getattr(sampler, "method", "MCMC")),
        "dimensions": int(getattr(sampler, "_dimensions", len(variable_names) or 0) or 0),
        "nchains": int(getattr(sampler, "_nchains", 1) or 1),
        "niters": int(getattr(sampler, "_niters", 0) or 0),
        "variable_names": list(variable_names),
        "task_result_dir": str(sample_cfg.get("task_result_dir", "")),
    }


def build_run_spec(core) -> Dict[str, Any]:
    config = deepcopy(getattr(getattr(core, "yaml", None), "config", None))
    config_file = getattr(core, "info", {}).get("config_file") if isinstance(getattr(core, "info", {}), dict) else None
    raw_yaml_text = None
    if isinstance(config_file, str) and os.path.exists(config_file):
        try:
            with open(config_file, "r", encoding="utf-8") as handle:
                raw_yaml_text = handle.read()
        except Exception:
            raw_yaml_text = None

    sample_cfg = dict(getattr(core, "info", {}).get("sample", {}) or {})
    task_result_dir = str(sample_cfg.get("task_result_dir", ""))
    task_root = str(getattr(core, "path", {}).get("task_root", ""))
    scan_name = str(getattr(core, "info", {}).get("scan_name", ""))
    logs_dir = str(getattr(core, "info", {}).get("logs_dir", ""))
    images_dir = str(getattr(core, "info", {}).get("images_dir", ""))
    try:
        worker_parallel = int(getattr(core.yaml, "get_worker_parallel")())
    except Exception:
        worker_parallel = int(getattr(getattr(core, "factory", None), "_max_workers", 0) or 0)

    workflow_layers = deepcopy(getattr(getattr(core, "workflow", None), "workflow", {}))
    workflow_calc_layer = deepcopy(getattr(getattr(core, "workflow", None), "calc_layer", {}))

    return {
        "raw_yaml_text": raw_yaml_text,
        "normalized_config": config,
        "scan_name": scan_name,
        "task_root": task_root,
        "task_result_dir": task_result_dir,
        "logs_dir": logs_dir,
        "images_dir": images_dir,
        "worker_parallel": worker_parallel,
        "sampler_method": str(getattr(getattr(core, "sampler", None), "method", "")),
        "workflow": workflow_layers,
        "workflow_layers": workflow_calc_layer,
    }


def build_factory_blueprint(core) -> Dict[str, Any]:
    factory = getattr(core, "factory", None)
    module_manager = getattr(core, "module_manager", None)
    workflow = getattr(core, "workflow", None)

    module_pools = {}
    if module_manager is not None and getattr(module_manager, "module_pools", None):
        for module_name, pool in module_manager.module_pools.items():
            if hasattr(pool, "export_blueprint"):
                module_pools[str(module_name)] = pool.export_blueprint()

    return {
        "max_workers": int(getattr(factory, "_max_workers", 0) or 0),
        "workflow": deepcopy(getattr(workflow, "workflow", {})),
        "calc_layer": deepcopy(getattr(workflow, "calc_layer", {})),
        "module_pools": module_pools,
    }


def build_state_payload(*, run_spec: Dict[str, Any], factory_blueprint: Dict[str, Any], sampler_state: Dict[str, Any], reason: str = "") -> Dict[str, Any]:
    sampler_signature = dict(sampler_state.get("sampler_signature", {}))
    created_at = utc_now_iso()
    integrity = {
        "config_hash": stable_json_hash(run_spec.get("normalized_config", {})),
        "sampler_signature": sampler_signature,
        "variable_signature": list(run_spec.get("normalized_config", {}).get("Sampling", {}).get("Variables", []) or []),
        "checkpoint_reason": str(reason or ""),
        "safe_barrier_confirmed": True,
    }
    return {
        "format": STATE_SAVER_FORMAT,
        "version": STATE_SAVER_VERSION,
        "created_at_utc": created_at,
        "timestamp_utc": created_at,
        "jarvishep_version": get_runtime_version(),
        "resume_policy": {
            "restore_inflight": False,
            "rebuild_worker_factory": True,
            "config_source": "checkpoint",
        },
        "run_spec": run_spec,
        "factory_blueprint": factory_blueprint,
        "sampler_state": sampler_state,
        "integrity": integrity,
    }


def validate_state_payload(
    payload: Any,
    *,
    require_run_spec: bool = True,
    require_factory_blueprint: bool = True,
) -> Tuple[bool, str]:
    if not isinstance(payload, dict):
        return False, "StateSaver checkpoint payload must be a mapping"

    payload_format = payload.get("format")
    if payload_format is not None and payload_format not in VALID_STATE_FORMATS:
        return False, f"StateSaver checkpoint unsupported format: {payload_format!r}"

    version = payload.get("version")
    if version is not None and not isinstance(version, int):
        return False, "StateSaver checkpoint version must be an integer"
    if payload_format == STATE_SAVER_FORMAT and version != STATE_SAVER_VERSION:
        return False, f"StateSaver checkpoint version mismatch: expect {STATE_SAVER_VERSION}, got {version!r}"
    if payload_format == CHECKPOINT_FORMAT and version != CHECKPOINT_VERSION:
        return False, f"StateSaver checkpoint version mismatch: expect {CHECKPOINT_VERSION}, got {version!r}"

    if require_run_spec and not isinstance(payload.get("run_spec"), dict):
        return False, "StateSaver checkpoint payload missing run_spec"

    if require_factory_blueprint and not isinstance(payload.get("factory_blueprint"), dict):
        return False, "StateSaver checkpoint payload missing factory_blueprint"

    if "sampler_state" in payload and not isinstance(payload.get("sampler_state"), dict):
        return False, "StateSaver checkpoint payload missing sampler_state"

    if "integrity" in payload and not isinstance(payload.get("integrity"), dict):
        return False, "StateSaver checkpoint payload missing integrity"

    return True, "ok"


def validate_sampler_signature(
    current: Dict[str, Any],
    saved: Dict[str, Any],
) -> Tuple[bool, str]:
    keys = ("method", "dimensions", "nchains", "niters", "variable_names")
    for key in keys:
        if current.get(key) != saved.get(key):
            return False, f"StateSaver checkpoint signature mismatch on {key}"
    return True, "ok"


def export_engine_state(engine: Any) -> Dict[str, Any]:
    raw = vars(engine)
    state = {}
    for key, value in raw.items():
        if key in _ENGINE_STATE_EXCLUDE:
            continue
        if callable(value):
            continue
        state[key] = value
    return {
        "module": engine.__class__.__module__,
        "class": engine.__class__.__name__,
        "state": state,
    }


def restore_engine_state(engine: Any, payload: Dict[str, Any]) -> None:
    expect_module = str(payload.get("module", ""))
    expect_class = str(payload.get("class", ""))
    if engine.__class__.__module__ != expect_module or engine.__class__.__name__ != expect_class:
        raise ValueError(
            "StateSaver checkpoint engine mismatch: "
            f"expect {expect_module}.{expect_class}, "
            f"got {engine.__class__.__module__}.{engine.__class__.__name__}"
        )
    for key, value in dict(payload.get("state", {})).items():
        if key in _ENGINE_STATE_EXCLUDE:
            continue
        setattr(engine, key, value)


def build_checkpoint_manifest(payload: Dict[str, Any], checkpoint_file: str) -> Dict[str, Any]:
    machine = payload.get("state_machine", {}) if isinstance(payload, dict) else {}
    chains = machine.get("chains", []) if isinstance(machine, dict) else []
    completed_steps = sum(int(chain.get("iter", 0)) for chain in chains if isinstance(chain, dict))
    return {
        "format": payload.get("format"),
        "version": payload.get("version"),
        "timestamp_utc": payload.get("timestamp_utc"),
        "checkpoint_file": os.path.abspath(checkpoint_file),
        "sampler_signature": payload.get("sampler_signature", {}),
        "state": machine.get("state"),
        "chains": len(chains),
        "completed_steps": int(completed_steps),
    }
