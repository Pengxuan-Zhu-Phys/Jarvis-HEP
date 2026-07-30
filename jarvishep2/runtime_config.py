#!/usr/bin/env python3
"""Runtime block normalization for Jarvis-HEP V2."""

from __future__ import annotations

from typing import Any, Mapping

from jarvishep2.archive_handoff import normalize_move_strategy, resolve_staging_dir
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, normalize_delete_method
from jarvishep2.sample_bucket import (
    SAMPLE_DIRECTORY_DEFAULTS,
    normalize_sample_directory,
)


RUNTIME_DEFAULTS: dict[str, Any] = {
    "mode": "auto",
    "workers": 0,
    "batch_size": 256,
    "sample_artifacts": "auto",
}

VALID_SAMPLE_ARTIFACTS = frozenset({"auto", "always", "never"})
VALID_RUNTIME_MODES = frozenset({"auto", "redis"})
VALID_CLEANUP_STRATEGIES = frozenset({"mv_to_staging", "direct"})
VALID_ARCHIVER_MODES = frozenset({"thread", "process"})
VALID_HANDOFF_MODES = frozenset({"direct", "staging"})
# Backward-compatible private aliases (tests / older imports).
_VALID_SAMPLE_ARTIFACTS = VALID_SAMPLE_ARTIFACTS
_VALID_RUNTIME_MODES = VALID_RUNTIME_MODES
_VALID_CLEANUP_STRATEGIES = VALID_CLEANUP_STRATEGIES
_VALID_ARCHIVER_MODES = VALID_ARCHIVER_MODES
_VALID_HANDOFF_MODES = VALID_HANDOFF_MODES

ARCHIVER_DEFAULTS: dict[str, Any] = {
    "mode": "process",
    "batch_size": 200,
    "flush_interval_sec": 1.0,
    # Seal samples.NNNN.hdf5 and generate its sibling CSV after this size.
    "max_hdf5_bytes": 1024 * 1024 * 1024,
    "strategy": "move",
    "delete_after_archive": True,
    # Default: no staging hop — Worker lands products under SAMPLE/<bucket>/<uuid>.
    "handoff": "direct",
    # Default: pack sealed SAMPLE buckets into <bucket>.tar.gz (V1 parity).
    "pack_buckets": True,
}
CLEANUP_DEFAULTS: dict[str, Any] = {
    # Default: skip staging; optional mv_to_staging retained for heavy-archive paths.
    "strategy": "direct",
    "staging_dir": None,
}
WATCHDOG_DEFAULTS: dict[str, Any] = {
    "enabled": True,
    "stale_sec": 30.0,
    "poll_interval_sec": 1.0,
    "max_sample_retries": 3,
}
FACTORY_DEFAULTS: dict[str, Any] = {
    "monitor_hz": 120.0,
}
# EnvReqs.V2 top-level keys accepted by the task loader (D12.4).
SUPPORTED_ENVREQS_V2_KEYS = frozenset(
    {
        "workers",
        "batch_size",
        "sample_directory",
        "cleanup",
        "archiver",
        "redis",
        "factory",
        "worker",
        "check_modules",
    }
)

# Defaults for ``Jarvis2 check`` / ``--check-modules`` (overridable via EnvReqs.V2).
CHECK_MODULES_DEFAULTS: dict[str, Any] = {
    # Prefer project data CSV; if missing at runtime, sampler draws n_samples points.
    "data": "&J/data/check_modules_points.csv",
    "n_samples": 10,
}


def _coerce_positive_int(value: Any, *, default: int) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError):
        return default
    return parsed if parsed >= 0 else default


def normalize_redis_config(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Merge optional ``EnvReqs.V2.redis`` over the internal broker defaults.

    Defaults remain local zero-config (``127.0.0.1:6379`` / db ``0``). Only
    ``host``, ``port``, and ``db`` are consumed; other keys are ignored.
    """
    from jarvishep2.redis_queue import INTERNAL_REDIS_CONFIG

    config = dict(INTERNAL_REDIS_CONFIG)
    if not isinstance(raw, Mapping):
        return config

    host = raw.get("host")
    if host is not None and str(host).strip():
        config["host"] = str(host).strip()

    if "port" in raw and raw.get("port") is not None:
        try:
            port = int(raw["port"])
        except (TypeError, ValueError):
            port = int(INTERNAL_REDIS_CONFIG["port"])
        if 1 <= port <= 65535:
            config["port"] = port

    if "db" in raw and raw.get("db") is not None:
        try:
            db = int(raw["db"])
        except (TypeError, ValueError):
            db = int(INTERNAL_REDIS_CONFIG["db"])
        if db >= 0:
            config["db"] = db
    return config


def normalize_factory_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize ``EnvReqs.V2.factory`` (monitor + watchdog knobs)."""
    factory = dict(FACTORY_DEFAULTS)
    if not isinstance(raw, Mapping):
        return factory

    monitor = raw.get("monitor")
    hz_raw = raw.get("monitor_hz")
    if hz_raw is None and isinstance(monitor, Mapping):
        hz_raw = monitor.get("hz", monitor.get("update_hz"))
    if hz_raw is not None:
        try:
            factory["monitor_hz"] = max(1.0, float(hz_raw))
        except (TypeError, ValueError):
            factory["monitor_hz"] = FACTORY_DEFAULTS["monitor_hz"]

    watchdog = raw.get("watchdog")
    if watchdog is None:
        watchdog = raw.get("Watchdog")
    if isinstance(watchdog, Mapping):
        factory["watchdog"] = normalize_watchdog_block(watchdog)
    return factory


def normalize_worker_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize ``EnvReqs.V2.worker`` (per-Worker execution policy)."""
    worker: dict[str, Any] = {
        "force_serial_layers": False,
        "sample_artifacts": RUNTIME_DEFAULTS["sample_artifacts"],
    }
    if not isinstance(raw, Mapping):
        return worker

    if "force_serial_layers" in raw:
        worker["force_serial_layers"] = bool(raw.get("force_serial_layers"))

    sample_artifacts = str(
        raw.get("sample_artifacts", worker["sample_artifacts"])
    ).strip().lower()
    if sample_artifacts in _VALID_SAMPLE_ARTIFACTS:
        worker["sample_artifacts"] = sample_artifacts
    return worker


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
        runtime["redis"] = normalize_redis_config(redis_block)

    if "force_serial_layers" in raw:
        runtime["force_serial_layers"] = bool(raw.get("force_serial_layers"))

    factory_block = raw.get("Factory")
    if factory_block is None:
        factory_block = raw.get("factory")
    if isinstance(factory_block, Mapping):
        runtime["Factory"] = normalize_factory_block(factory_block)

    subprocess = raw.get("Subprocess")
    if isinstance(subprocess, Mapping):
        runtime["Subprocess"] = dict(subprocess)

    file_operation = raw.get("FileOperation")
    if isinstance(file_operation, Mapping):
        runtime["FileOperation"] = normalize_file_operation(file_operation)

    watchdog = raw.get("Watchdog")
    if isinstance(watchdog, Mapping):
        runtime["Watchdog"] = normalize_watchdog_block(watchdog)

    return runtime


def normalize_file_operation(raw: Mapping[str, Any] | None) -> dict[str, str]:
    """Normalize ``Runtime.FileOperation`` settings."""
    if not isinstance(raw, Mapping):
        return {"delete_method": DEFAULT_DELETE_METHOD}
    return {
        "delete_method": normalize_delete_method(raw.get("delete_method")),
    }


def normalize_cleanup_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize ``Calculators.Cleanup`` settings."""
    cleanup = dict(CLEANUP_DEFAULTS)
    if not isinstance(raw, Mapping):
        return cleanup
    strategy = str(raw.get("strategy", cleanup["strategy"])).strip().lower()
    # Accept legacy handoff alias from Archiver.handoff when cleanup uses it.
    if strategy in {"staging", "mv_to_staging"}:
        strategy = "mv_to_staging"
    elif strategy in {"direct", "none", "off"}:
        strategy = "direct"
    cleanup["strategy"] = (
        strategy if strategy in _VALID_CLEANUP_STRATEGIES else CLEANUP_DEFAULTS["strategy"]
    )
    staging_dir = raw.get("staging_dir")
    cleanup["staging_dir"] = str(staging_dir).strip() if staging_dir else None
    return cleanup


def normalize_archiver_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize ``Calculators.Archiver`` settings."""
    archiver = dict(ARCHIVER_DEFAULTS)
    if not isinstance(raw, Mapping):
        return archiver
    mode = str(raw.get("mode", archiver["mode"])).strip().lower()
    archiver["mode"] = mode if mode in _VALID_ARCHIVER_MODES else ARCHIVER_DEFAULTS["mode"]
    batch_size = _coerce_positive_int(raw.get("batch_size", archiver["batch_size"]), default=200)
    archiver["batch_size"] = batch_size if batch_size > 0 else ARCHIVER_DEFAULTS["batch_size"]
    flush_interval = raw.get("flush_interval_sec", archiver["flush_interval_sec"])
    try:
        archiver["flush_interval_sec"] = max(0.05, float(flush_interval))
    except (TypeError, ValueError):
        archiver["flush_interval_sec"] = ARCHIVER_DEFAULTS["flush_interval_sec"]
    max_hdf5_bytes = _coerce_positive_int(
        raw.get("max_hdf5_bytes", archiver["max_hdf5_bytes"]),
        default=ARCHIVER_DEFAULTS["max_hdf5_bytes"],
    )
    archiver["max_hdf5_bytes"] = (
        max_hdf5_bytes if max_hdf5_bytes > 0 else ARCHIVER_DEFAULTS["max_hdf5_bytes"]
    )
    archiver["strategy"] = normalize_move_strategy(raw.get("strategy", archiver["strategy"]))
    archiver["delete_after_archive"] = bool(
        raw.get("delete_after_archive", archiver["delete_after_archive"])
    )
    handoff = str(raw.get("handoff", archiver["handoff"])).strip().lower()
    if handoff in {"staging", "mv_to_staging"}:
        handoff = "staging"
    elif handoff in {"direct", "none", "off"}:
        handoff = "direct"
    archiver["handoff"] = handoff if handoff in _VALID_HANDOFF_MODES else ARCHIVER_DEFAULTS["handoff"]
    if "pack_buckets" in raw:
        archiver["pack_buckets"] = bool(raw.get("pack_buckets"))
    return archiver


def get_calculators_block(config: Mapping[str, Any] | None) -> dict[str, Any]:
    if not isinstance(config, Mapping):
        return {}
    calculators = config.get("Calculators") or {}
    return dict(calculators) if isinstance(calculators, Mapping) else {}


def get_cleanup_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    calculators = get_calculators_block(config)
    return normalize_cleanup_block(calculators.get("Cleanup"))


def get_archiver_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    calculators = get_calculators_block(config)
    return normalize_archiver_block(calculators.get("Archiver"))


def get_staging_dir(
    config: Mapping[str, Any] | None,
    *,
    task_result_dir: str,
    create: bool | None = None,
) -> str:
    """Resolve staging root path.

    When ``create`` is None, create only if staging handoff is enabled for this
    config (never mkdir for default direct handoff).
    """
    cleanup = get_cleanup_config(config)
    if create is None:
        create = handoff_to_staging_enabled(config)
    return resolve_staging_dir(
        task_result_dir,
        configured=cleanup.get("staging_dir"),
        create=bool(create),
    )


def handoff_to_staging_enabled(config: Mapping[str, Any] | None) -> bool:
    """True only when Cleanup/Archiver explicitly request the staging hop."""
    cleanup = get_cleanup_config(config)
    if cleanup.get("strategy") == "mv_to_staging":
        return True
    archiver = get_archiver_config(config)
    return str(archiver.get("handoff", "direct")).strip().lower() == "staging"


def get_sample_directory_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return normalized SAMPLE bucket settings (Scan or EnvReqs.V2)."""
    if not isinstance(config, Mapping):
        return normalize_sample_directory(None)
    scan = config.get("Scan") if isinstance(config.get("Scan"), Mapping) else {}
    if isinstance(scan, Mapping) and isinstance(scan.get("sample_directory"), Mapping):
        return normalize_sample_directory(scan.get("sample_directory"))
    runtime = config.get("Runtime") if isinstance(config.get("Runtime"), Mapping) else {}
    if isinstance(runtime, Mapping) and isinstance(runtime.get("sample_directory"), Mapping):
        return normalize_sample_directory(runtime.get("sample_directory"))
    envreqs = config.get("EnvReqs") if isinstance(config.get("EnvReqs"), Mapping) else {}
    v2 = envreqs.get("V2") if isinstance(envreqs, Mapping) else None
    if isinstance(v2, Mapping) and isinstance(v2.get("sample_directory"), Mapping):
        return normalize_sample_directory(v2.get("sample_directory"))
    return normalize_sample_directory(dict(SAMPLE_DIRECTORY_DEFAULTS))


def pack_buckets_enabled(config: Mapping[str, Any] | None) -> bool:
    sample_dir = get_sample_directory_config(config)
    if not bool(sample_dir.get("enabled", True)):
        return False
    archiver = get_archiver_config(config)
    if "pack_buckets" in archiver:
        return bool(archiver.get("pack_buckets"))
    return bool(sample_dir.get("pack", True))


def get_delete_method(config: Mapping[str, Any] | None) -> str:
    """Return the configured delete backend (default ``shutil``)."""
    runtime = get_runtime_block(config)
    file_operation = runtime.get("FileOperation") or {}
    if isinstance(file_operation, Mapping):
        return normalize_delete_method(file_operation.get("delete_method"))
    return DEFAULT_DELETE_METHOD


def normalize_watchdog_block(raw: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize ``Runtime.Watchdog`` settings (WP-D6.1)."""
    watchdog = dict(WATCHDOG_DEFAULTS)
    if not isinstance(raw, Mapping):
        return watchdog
    watchdog["enabled"] = bool(raw.get("enabled", watchdog["enabled"]))
    stale_sec = raw.get("stale_sec", watchdog["stale_sec"])
    try:
        watchdog["stale_sec"] = max(1.0, float(stale_sec))
    except (TypeError, ValueError):
        watchdog["stale_sec"] = WATCHDOG_DEFAULTS["stale_sec"]
    poll_interval = raw.get("poll_interval_sec", watchdog["poll_interval_sec"])
    try:
        watchdog["poll_interval_sec"] = max(0.1, float(poll_interval))
    except (TypeError, ValueError):
        watchdog["poll_interval_sec"] = WATCHDOG_DEFAULTS["poll_interval_sec"]
    max_retries = raw.get("max_sample_retries", watchdog["max_sample_retries"])
    try:
        watchdog["max_sample_retries"] = max(0, int(max_retries))
    except (TypeError, ValueError):
        watchdog["max_sample_retries"] = WATCHDOG_DEFAULTS["max_sample_retries"]
    return watchdog


def get_watchdog_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    runtime = get_runtime_block(config)
    watchdog = runtime.get("Watchdog")
    if isinstance(watchdog, Mapping):
        return normalize_watchdog_block(watchdog)
    factory = runtime.get("Factory")
    if isinstance(factory, Mapping) and isinstance(factory.get("watchdog"), Mapping):
        return normalize_watchdog_block(factory.get("watchdog"))
    return dict(WATCHDOG_DEFAULTS)


def get_factory_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return monitor/watchdog factory knobs (``Runtime.Factory`` / EnvReqs.V2.factory)."""
    runtime = get_runtime_block(config)
    factory = runtime.get("Factory")
    if isinstance(factory, Mapping):
        return normalize_factory_block(factory)
    return dict(FACTORY_DEFAULTS)


def get_redis_config(config: Mapping[str, Any] | None) -> dict[str, Any]:
    """Return the broker connection (``Runtime.redis`` / EnvReqs.V2.redis / internal)."""
    runtime = get_runtime_block(config)
    redis_block = runtime.get("redis")
    if isinstance(redis_block, Mapping):
        return normalize_redis_config(redis_block)
    return normalize_redis_config(None)


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


def parse_registered_executables(config: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """Return raw ``LibDeps.registered_executables`` entries from a task config.

    Implementation lives in ``command_parser`` (sole consumer) to avoid import cycles;
    this re-export keeps the historical public path stable.
    """
    from jarvishep2.command_parser import parse_registered_executables as _parse

    return _parse(config)


def should_materialize_on_failure(sample_cfg: Mapping[str, Any] | None) -> bool:
    if not isinstance(sample_cfg, Mapping):
        return True
    mode = str(sample_cfg.get("sample_artifacts", RUNTIME_DEFAULTS["sample_artifacts"])).strip().lower()
    return mode != "never"
