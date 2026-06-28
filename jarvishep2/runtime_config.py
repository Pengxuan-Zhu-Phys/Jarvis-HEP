#!/usr/bin/env python3
"""Runtime block normalization for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping


RUNTIME_DEFAULTS: dict[str, Any] = {
    "mode": "auto",
    "workers": 0,
    "batch_size": 256,
    "sample_artifacts": "auto",
}

_VALID_SAMPLE_ARTIFACTS = frozenset({"auto", "always", "never"})
_VALID_RUNTIME_MODES = frozenset({"auto", "redis"})


def _coerce_positive_int(value: Any, *, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed >= 0 else default


def normalize_runtime_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return a normalized Runtime block with V2 defaults."""
    runtime = dict(RUNTIME_DEFAULTS)
    if not isinstance(raw, Mapping):
        return runtime

    mode = str(raw.get("mode", runtime["mode"])).strip().lower()
    runtime["mode"] = mode if mode in _VALID_RUNTIME_MODES else RUNTIME_DEFAULTS["mode"]

    runtime["workers"] = _coerce_positive_int(raw.get("workers", runtime["workers"]), default=0)

    batch_size = _coerce_positive_int(raw.get("batch_size", runtime["batch_size"]), default=256)
    runtime["batch_size"] = batch_size if batch_size > 0 else RUNTIME_DEFAULTS["batch_size"]

    sample_artifacts = str(raw.get("sample_artifacts", runtime["sample_artifacts"])).strip().lower()
    if sample_artifacts not in _VALID_SAMPLE_ARTIFACTS:
        sample_artifacts = RUNTIME_DEFAULTS["sample_artifacts"]
    runtime["sample_artifacts"] = sample_artifacts

    redis_block = raw.get("redis")
    if isinstance(redis_block, Mapping):
        runtime["redis"] = dict(redis_block)

    subprocess = raw.get("Subprocess")
    if isinstance(subprocess, Mapping):
        runtime["Subprocess"] = dict(subprocess)

    return runtime


def get_runtime_block(config: Mapping[str, Any] | None) -> dict[str, Any]:
    if not isinstance(config, Mapping):
        return dict(RUNTIME_DEFAULTS)
    return normalize_runtime_block(config.get("Runtime"))


def workflow_has_calculator(config: Mapping[str, Any] | None) -> bool:
    if not isinstance(config, Mapping):
        return False
    calculators = (config.get("Calculators", {}) or {}).get("Modules", []) or []
    return bool(calculators)


def _mapping_contains_token(value: Any, token: str) -> bool:
    if isinstance(value, str):
        return token in value
    if isinstance(value, Mapping):
        return any(_mapping_contains_token(item, token) for item in value.values())
    if isinstance(value, (list, tuple, set)):
        return any(_mapping_contains_token(item, token) for item in value)
    return False


def workflow_references_sdir(config: Mapping[str, Any] | None) -> bool:
    if not isinstance(config, Mapping):
        return False
    operas = (config.get("Operas", {}) or {}).get("Modules", []) or []
    calculators = (config.get("Calculators", {}) or {}).get("Modules", []) or []
    return any(_mapping_contains_token(module, "@Sdir") for module in (*operas, *calculators))


def should_eager_materialize(sample_cfg: Mapping[str, Any] | None) -> bool:
    """Return True when per-sample dirs/logs must exist before workflow execution."""
    if not isinstance(sample_cfg, Mapping):
        return True

    mode = str(sample_cfg.get("sample_artifacts", RUNTIME_DEFAULTS["sample_artifacts"])).strip().lower()
    if mode == "always":
        return True
    if mode == "never":
        return False

    if bool(sample_cfg.get("workflow_has_calculator")):
        return True
    if bool(sample_cfg.get("workflow_references_sdir")):
        return True
    return False


def should_materialize_on_failure(sample_cfg: Mapping[str, Any] | None) -> bool:
    if not isinstance(sample_cfg, Mapping):
        return True
    mode = str(sample_cfg.get("sample_artifacts", RUNTIME_DEFAULTS["sample_artifacts"])).strip().lower()
    return mode != "never"