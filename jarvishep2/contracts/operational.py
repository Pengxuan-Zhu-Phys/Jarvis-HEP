#!/usr/bin/env python3
"""EnvReqs.V2 / Archiver / Cleanup / sample_directory contracts (D14 L3)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING, Any

from jarvishep2.contracts.common import try_float, try_int, unknown_keys
from jarvishep2.runtime_config import (
    ARCHIVER_DEFAULTS,
    CLEANUP_DEFAULTS,
    SUPPORTED_ENVREQS_V2_KEYS,
    VALID_ARCHIVER_MODES,
    VALID_CLEANUP_STRATEGIES,
    VALID_HANDOFF_MODES,
    VALID_SAMPLE_ARTIFACTS,
)
from jarvishep2.sample_bucket import SAMPLE_DIRECTORY_DEFAULTS

if TYPE_CHECKING:
    from jarvishep2.task_validation import ValidationIssue

_WORKER_POLICY_KEYS = frozenset({"force_serial_layers", "sample_artifacts"})
_FACTORY_KEYS = frozenset({"monitor_hz", "monitor", "watchdog", "Watchdog"})
_REDIS_KEYS = frozenset({"host", "port", "db"})
_WATCHDOG_KEYS = frozenset(
    {"enabled", "stale_sec", "poll_interval_sec", "max_sample_retries"}
)
_ARCHIVER_KEYS = frozenset(ARCHIVER_DEFAULTS.keys())
_CLEANUP_KEYS = frozenset(CLEANUP_DEFAULTS.keys())
_SAMPLE_DIR_KEYS = frozenset(SAMPLE_DIRECTORY_DEFAULTS.keys())

_VALID_ARCHIVER_STRATEGY = frozenset({"move", "copy"})


def validate_operational_blocks(config: Mapping[str, Any]) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []

    envreqs = config.get("EnvReqs")
    if envreqs is not None and not isinstance(envreqs, Mapping):
        issues.append(
            issue(
                "error",
                "JV2-ENV-001",
                "EnvReqs",
                f"expected a mapping, got {type(envreqs).__name__}",
            )
        )
        return issues

    v2 = (envreqs or {}).get("V2") if isinstance(envreqs, Mapping) else None
    if v2 is not None:
        if not isinstance(v2, Mapping):
            issues.append(
                issue(
                    "error",
                    "JV2-ENV-002",
                    "EnvReqs.V2",
                    f"expected a mapping, got {type(v2).__name__}",
                )
            )
        else:
            # Unknown top-level V2 keys already rejected at load; re-check for safety.
            extra = unknown_keys(v2, SUPPORTED_ENVREQS_V2_KEYS)
            if extra:
                issues.append(
                    issue(
                        "error",
                        "JV2-ENV-003",
                        "EnvReqs.V2",
                        f"unsupported setting(s): {', '.join(extra)}; "
                        f"supported: {', '.join(sorted(SUPPORTED_ENVREQS_V2_KEYS))}",
                    )
                )

            if "workers" in v2:
                workers = try_int(v2.get("workers"))
                if workers is None or workers < 0:
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-010",
                            "EnvReqs.V2.workers",
                            f"expected non-negative integer, got {v2.get('workers')!r}",
                        )
                    )

            if "batch_size" in v2:
                batch = try_int(v2.get("batch_size"))
                if batch is None or batch < 1:
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-011",
                            "EnvReqs.V2.batch_size",
                            f"expected integer ≥ 1, got {v2.get('batch_size')!r}",
                        )
                    )

            if "checkpoint_heartbeat_sec" in v2:
                heartbeat = try_float(v2.get("checkpoint_heartbeat_sec"))
                if heartbeat is None or heartbeat < 1.0:
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-012",
                            "EnvReqs.V2.checkpoint_heartbeat_sec",
                            "expected number ≥ 1 second, got "
                            f"{v2.get('checkpoint_heartbeat_sec')!r}",
                        )
                    )

            worker = v2.get("worker")
            if worker is not None and not isinstance(worker, Mapping):
                issues.append(
                    issue(
                        "error",
                        "JV2-ENV-019",
                        "EnvReqs.V2.worker",
                        "worker is a policy mapping; use workers for the worker count",
                    )
                )
            elif isinstance(worker, Mapping):
                w_extra = unknown_keys(worker, _WORKER_POLICY_KEYS)
                if w_extra:
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-020",
                            "EnvReqs.V2.worker",
                            f"unknown key(s): {', '.join(w_extra)}; "
                            f"allowed: {', '.join(sorted(_WORKER_POLICY_KEYS))}",
                        )
                    )
                if "sample_artifacts" in worker:
                    sa = str(worker.get("sample_artifacts")).strip().lower()
                    if sa not in VALID_SAMPLE_ARTIFACTS:
                        issues.append(
                            issue(
                                "error",
                                "JV2-ENV-021",
                                "EnvReqs.V2.worker.sample_artifacts",
                                f"{sa!r} invalid; "
                                f"expected one of: {', '.join(sorted(VALID_SAMPLE_ARTIFACTS))}",
                            )
                        )

            factory = v2.get("factory")
            if factory is not None:
                if not isinstance(factory, Mapping):
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-030",
                            "EnvReqs.V2.factory",
                            f"expected a mapping, got {type(factory).__name__}",
                        )
                    )
                else:
                    f_extra = unknown_keys(factory, _FACTORY_KEYS)
                    if f_extra:
                        issues.append(
                            issue(
                                "error",
                                "JV2-ENV-031",
                                "EnvReqs.V2.factory",
                                f"unknown key(s): {', '.join(f_extra)}",
                            )
                        )
                    if "monitor_hz" in factory:
                        hz = try_float(factory.get("monitor_hz"))
                        if hz is None or hz < 1.0:
                            issues.append(
                                issue(
                                    "error",
                                    "JV2-ENV-032",
                                    "EnvReqs.V2.factory.monitor_hz",
                                    f"expected number ≥ 1, got {factory.get('monitor_hz')!r}",
                                )
                            )
                    for wd_key in ("watchdog", "Watchdog"):
                        wd = factory.get(wd_key)
                        if wd is None:
                            continue
                        if not isinstance(wd, Mapping):
                            issues.append(
                                issue(
                                    "error",
                                    "JV2-ENV-033",
                                    f"EnvReqs.V2.factory.{wd_key}",
                                    f"expected a mapping, got {type(wd).__name__}",
                                )
                            )
                        else:
                            wd_extra = unknown_keys(wd, _WATCHDOG_KEYS)
                            if wd_extra:
                                issues.append(
                                    issue(
                                        "error",
                                        "JV2-ENV-034",
                                        f"EnvReqs.V2.factory.{wd_key}",
                                        f"unknown key(s): {', '.join(wd_extra)}",
                                    )
                                )

            check_modules = v2.get("check_modules")
            if isinstance(check_modules, Mapping) and "timeout" in check_modules:
                issues.append(
                    issue(
                        "error",
                        "JV2-ENV-060",
                        "EnvReqs.V2.check_modules.timeout",
                        "timeout is a removed alias; use timeout_sec",
                        suggestion="Replace timeout with timeout_sec and remove the old key.",
                    )
                )
            redis = v2.get("redis")
            if redis is not None:
                if not isinstance(redis, Mapping):
                    issues.append(
                        issue(
                            "error",
                            "JV2-ENV-040",
                            "EnvReqs.V2.redis",
                            f"expected a mapping, got {type(redis).__name__}",
                        )
                    )
                else:
                    r_extra = unknown_keys(redis, _REDIS_KEYS)
                    if r_extra:
                        issues.append(
                            issue(
                                "error",
                                "JV2-ENV-041",
                                "EnvReqs.V2.redis",
                                f"unknown key(s): {', '.join(r_extra)}; "
                                f"allowed: {', '.join(sorted(_REDIS_KEYS))}",
                            )
                        )
                    if "port" in redis:
                        port = try_int(redis.get("port"))
                        if port is None or not (1 <= port <= 65535):
                            issues.append(
                                issue(
                                    "error",
                                    "JV2-ENV-042",
                                    "EnvReqs.V2.redis.port",
                                    f"expected integer 1..65535, got {redis.get('port')!r}",
                                )
                            )
                    if "db" in redis:
                        db = try_int(redis.get("db"))
                        if db is None or db < 0:
                            issues.append(
                                issue(
                                    "error",
                                    "JV2-ENV-043",
                                    "EnvReqs.V2.redis.db",
                                    f"expected non-negative integer, got {redis.get('db')!r}",
                                )
                            )

            sample_dir = v2.get("sample_directory")
            if sample_dir is not None:
                issues.extend(_validate_sample_directory(sample_dir, "EnvReqs.V2.sample_directory"))

    calculators = config.get("Calculators")
    if isinstance(calculators, Mapping):
        archiver = calculators.get("Archiver")
        if archiver is not None:
            issues.extend(_validate_archiver(archiver))
        cleanup = calculators.get("Cleanup")
        if cleanup is not None:
            issues.extend(_validate_cleanup(cleanup))

    return issues


def _validate_sample_directory(raw: Any, path: str) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    if not isinstance(raw, Mapping):
        return [
            issue(
                "error",
                "JV2-ENV-050",
                path,
                f"expected a mapping, got {type(raw).__name__}",
            )
        ]
    extra = unknown_keys(raw, _SAMPLE_DIR_KEYS)
    if extra:
        issues.append(
            issue(
                "error",
                "JV2-ENV-051",
                path,
                f"unknown key(s): {', '.join(extra)}; "
                f"allowed: {', '.join(sorted(_SAMPLE_DIR_KEYS))}",
            )
        )
    for key in ("limit", "width", "start_bucket"):
        if key not in raw:
            continue
        val = try_int(raw.get(key))
        if val is None or val < 1:
            issues.append(
                issue(
                    "error",
                    "JV2-ENV-052",
                    f"{path}.{key}",
                    f"expected integer ≥ 1, got {raw.get(key)!r}",
                )
            )
    return issues


def _validate_archiver(raw: Any) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    path = "Calculators.Archiver"
    if not isinstance(raw, Mapping):
        return [
            issue(
                "error",
                "JV2-ARC-001",
                path,
                f"expected a mapping, got {type(raw).__name__}",
            )
        ]
    extra = unknown_keys(raw, _ARCHIVER_KEYS)
    if extra:
        issues.append(
            issue(
                "error",
                "JV2-ARC-002",
                path,
                f"unknown key(s): {', '.join(extra)}; "
                f"allowed: {', '.join(sorted(_ARCHIVER_KEYS))}",
            )
        )
    if "mode" in raw:
        mode = str(raw.get("mode")).strip().lower()
        if mode not in VALID_ARCHIVER_MODES:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-010",
                    f"{path}.mode",
                    f"{mode!r} invalid; expected one of: "
                    f"{', '.join(sorted(VALID_ARCHIVER_MODES))}",
                )
            )
    if "batch_size" in raw:
        batch = try_int(raw.get("batch_size"))
        if batch is None or batch < 1:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-011",
                    f"{path}.batch_size",
                    f"expected integer ≥ 1, got {raw.get('batch_size')!r}",
                )
            )
    if "flush_interval_sec" in raw:
        flush = try_float(raw.get("flush_interval_sec"))
        if flush is None or flush < 0.05:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-012",
                    f"{path}.flush_interval_sec",
                    f"expected number ≥ 0.05, got {raw.get('flush_interval_sec')!r}",
                )
            )
    if "max_hdf5_bytes" in raw:
        max_bytes = try_int(raw.get("max_hdf5_bytes"))
        if max_bytes is None or max_bytes < 1:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-014",
                    f"{path}.max_hdf5_bytes",
                    f"expected integer ≥ 1, got {raw.get('max_hdf5_bytes')!r}",
                )
            )
    if "strategy" in raw:
        strategy = str(raw.get("strategy")).strip().lower()
        if strategy not in _VALID_ARCHIVER_STRATEGY:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-013",
                    f"{path}.strategy",
                    f"{strategy!r} invalid; expected one of: "
                    f"{', '.join(sorted(_VALID_ARCHIVER_STRATEGY))}",
                )
            )
    if "handoff" in raw:
        handoff = str(raw.get("handoff")).strip().lower()
        # Accept documented aliases, then map to canonical.
        if handoff in {"staging", "mv_to_staging"}:
            handoff = "staging"
        elif handoff in {"direct", "none", "off"}:
            handoff = "direct"
        if handoff not in VALID_HANDOFF_MODES:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-014",
                    f"{path}.handoff",
                    f"{raw.get('handoff')!r} invalid; expected one of: "
                    f"{', '.join(sorted(VALID_HANDOFF_MODES))} "
                    f"(aliases: staging/mv_to_staging, direct/none/off)",
                )
            )
    return issues


def _validate_cleanup(raw: Any) -> list[ValidationIssue]:
    from jarvishep2.task_validation import issue

    issues: list[ValidationIssue] = []
    path = "Calculators.Cleanup"
    if not isinstance(raw, Mapping):
        return [
            issue(
                "error",
                "JV2-ARC-020",
                path,
                f"expected a mapping, got {type(raw).__name__}",
            )
        ]
    extra = unknown_keys(raw, _CLEANUP_KEYS)
    if extra:
        issues.append(
            issue(
                "error",
                "JV2-ARC-021",
                path,
                f"unknown key(s): {', '.join(extra)}; "
                f"allowed: {', '.join(sorted(_CLEANUP_KEYS))}",
            )
        )
    if "strategy" in raw:
        strategy = str(raw.get("strategy")).strip().lower()
        if strategy in {"staging", "mv_to_staging"}:
            strategy = "mv_to_staging"
        elif strategy in {"direct", "none", "off"}:
            strategy = "direct"
        if strategy not in VALID_CLEANUP_STRATEGIES:
            issues.append(
                issue(
                    "error",
                    "JV2-ARC-022",
                    f"{path}.strategy",
                    f"{raw.get('strategy')!r} invalid; expected one of: "
                    f"{', '.join(sorted(VALID_CLEANUP_STRATEGIES))}",
                )
            )
    return issues


__all__ = ["validate_operational_blocks"]
