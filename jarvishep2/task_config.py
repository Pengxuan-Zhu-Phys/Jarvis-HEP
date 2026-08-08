#!/usr/bin/env python3
"""Task YAML loading and path normalization for Jarvis CLI runs."""

from __future__ import annotations

import os
import re
import tempfile
from collections.abc import Mapping
from copy import deepcopy
from typing import Any

import yaml

from jarvishep2.base import decode_path, infer_project_root, scan_output_root
from jarvishep2.runtime_config import (
    CHECK_MODULES_DEFAULTS,
    SUPPORTED_ENVREQS_V2_KEYS,
    normalize_factory_block,
    normalize_redis_config,
    normalize_runtime_block,
    normalize_worker_block,
)


class TaskCardLoadError(ValueError):
    """A task-card loading error with a stable diagnostic and correction."""

    def __init__(
        self, code: str, message: str, suggestion: str, example: str | None = None
    ) -> None:
        self.code = code
        self.suggestion = suggestion
        self.example = example
        super().__init__(message)


class LoadedTaskConfig(dict[str, Any]):
    """Normalized task configuration retaining its pre-injection YAML document."""

    raw_task_card: dict[str, Any]


def _load_yaml_document(path: str, *, label: str) -> Any:
    """Load one user-selected YAML file with a consistent syntax diagnostic."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return yaml.safe_load(handle)
    except yaml.YAMLError as exc:
        mark = getattr(exc, "problem_mark", None)
        location = ""
        if mark is not None:
            location = f" at line {mark.line + 1}, column {mark.column + 1}"
        problem = str(getattr(exc, "problem", "invalid YAML syntax"))
        raise TaskCardLoadError(
            "JV2-YAML-001",
            f"invalid YAML syntax in {label} {path}{location}: {problem}",
            "Check indentation, list markers, quoting, and ':' placement near the reported location.",
            "Sampling:\n  Method: Random\n  Point number: 100",
        ) from exc

def _deep_merge(defaults: Mapping[str, Any], overrides: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge task values over externally supplied defaults."""
    merged = deepcopy(dict(defaults))
    for key, value in overrides.items():
        current = merged.get(key)
        if isinstance(current, Mapping) and isinstance(value, Mapping):
            merged[key] = _deep_merge(current, value)
        else:
            merged[key] = deepcopy(value)
    return merged


def _is_enabled(value: Any) -> bool:
    """Interpret the conventional YAML boolean spellings used by V1 projects."""
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _canonicalize_v2_settings(settings: Mapping[str, Any]) -> dict[str, Any]:
    """Normalize the singular ``worker`` count spelling before merging.

    ``EnvReqs.V2.worker`` is dual-purpose (D12.4):

    * **scalar** (``worker: 4``) — compatibility alias for ``workers``
    * **mapping** (``worker: {force_serial_layers: …}``) — Worker policy group

    A scalar ``worker`` cannot coexist with ``workers``. A mapping ``worker``
    group may sit alongside ``workers``.
    """
    normalized = deepcopy(dict(settings))
    if "worker" in normalized:
        worker_value = normalized["worker"]
        if isinstance(worker_value, Mapping):
            # Keep the Worker policy group under EnvReqs.V2.worker.
            pass
        else:
            if "workers" in normalized:
                raise ValueError(
                    "EnvReqs.V2 accepts either workers or worker (count), not both; "
                    "use workers: N plus worker: {…} for the policy group"
                )
            normalized["workers"] = normalized.pop("worker")
    return normalized

def _v2_defaults_from_envreqs(
    config: Mapping[str, Any], *, project_root: str, yaml_dir: str
) -> dict[str, Any]:
    """Read optional V2 defaults through V1's ``EnvReqs`` entry point.

    The task-level V1 shape is deliberately retained::

        EnvReqs.Check_default_dependencies.default_yaml_path: "&J/deps/..."

    Only ``EnvReqs.V2`` from that external YAML is interpreted by V2.  It is
    merged with the task's ``EnvReqs.V2`` block; task values take precedence.
    """
    envreqs = config.get("EnvReqs")
    if not isinstance(envreqs, Mapping):
        return {}
    dependency_check = envreqs.get("Check_default_dependencies")
    if not isinstance(dependency_check, Mapping) or not _is_enabled(
        dependency_check.get("required")
    ):
        return {}

    raw_path = dependency_check.get("default_yaml_path")
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError(
            "EnvReqs.Check_default_dependencies.required is true, but "
            "default_yaml_path is missing"
        )
    defaults_path = decode_path(raw_path, project_root=project_root, base_dir=yaml_dir)
    if not os.path.isfile(defaults_path):
        raise FileNotFoundError(
            "EnvReqs.Check_default_dependencies default YAML not found: "
            f"{defaults_path}"
        )
    defaults_document = _load_yaml_document(defaults_path, label="default environment YAML")
    if defaults_document is None:
        return {}
    if not isinstance(defaults_document, Mapping):
        raise ValueError(f"default environment YAML must contain a mapping: {defaults_path}")

    defaults_envreqs = defaults_document.get("EnvReqs")
    if not isinstance(defaults_envreqs, Mapping):
        return {}
    v2_defaults = defaults_envreqs.get("V2")
    if v2_defaults is None:
        return {}
    if not isinstance(v2_defaults, Mapping):
        raise ValueError(f"EnvReqs.V2 must be a mapping: {defaults_path}")
    return _canonicalize_v2_settings(v2_defaults)


def _legacy_environment_requirements_from_defaults(
    config: Mapping[str, Any], *, project_root: str, yaml_dir: str
) -> dict[str, Any]:
    """Load V1 Python/ROOT requirements from the selected default environment YAML."""
    defaults_path = _default_environment_yaml_path(
        config, project_root=project_root, yaml_dir=yaml_dir
    )
    if defaults_path is None:
        return {}
    document = _load_yaml_document(defaults_path, label="default environment YAML")
    if not isinstance(document, Mapping):
        return {}
    envreqs = document.get("EnvReqs")
    if not isinstance(envreqs, Mapping):
        return {}
    return {
        key: deepcopy(envreqs[key])
        for key in ("Python", "CERN_ROOT")
        if isinstance(envreqs.get(key), Mapping)
    }


def _default_environment_yaml_path(
    config: Mapping[str, Any], *, project_root: str, yaml_dir: str
) -> str | None:
    """Return the project default-environment YAML selected by a task, if any."""
    envreqs = config.get("EnvReqs")
    if not isinstance(envreqs, Mapping):
        return None
    dependency_check = envreqs.get("Check_default_dependencies")
    if not isinstance(dependency_check, Mapping) or not _is_enabled(
        dependency_check.get("required")
    ):
        return None
    raw_path = dependency_check.get("default_yaml_path")
    if not isinstance(raw_path, str) or not raw_path.strip():
        return None
    return decode_path(raw_path, project_root=project_root, base_dir=yaml_dir)


def update_default_redis_port(defaults_path: str, port: int) -> None:
    """Atomically set ``EnvReqs.V2.redis.port`` without discarding YAML comments.

    Project environment files are user-maintained documents, so a full
    ``yaml.safe_dump`` rewrite would be unnecessarily destructive.  This small
    structural edit retains their layout and comments while adding the Redis
    block when a legacy defaults file has not declared one yet.
    """
    path = os.path.abspath(str(defaults_path))
    if not os.path.isfile(path):
        raise FileNotFoundError(f"default environment YAML not found: {path}")
    try:
        selected_port = int(port)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid Redis port: {port!r}") from exc
    if not 1 <= selected_port <= 65535:
        raise ValueError(f"invalid Redis port: {selected_port}")

    with open(path, "r", encoding="utf-8") as handle:
        text = handle.read()
    document = yaml.safe_load(text)
    if not isinstance(document, Mapping):
        raise ValueError(f"default environment YAML must contain a mapping: {path}")
    envreqs = document.get("EnvReqs")
    if not isinstance(envreqs, Mapping) or not isinstance(envreqs.get("V2"), Mapping):
        raise ValueError(f"default environment YAML must contain EnvReqs.V2: {path}")

    lines = text.splitlines(keepends=True)
    v2_line = next(
        (
            index
            for index, line in enumerate(lines)
            if re.match(r"^\s{2}V2:\s*(?:#.*)?(?:\r?\n)?$", line)
        ),
        None,
    )
    if v2_line is None:
        raise ValueError(f"could not locate EnvReqs.V2 in default YAML: {path}")
    end_v2 = len(lines)
    for index in range(v2_line + 1, len(lines)):
        if lines[index].strip() and not lines[index].startswith((" ", "\t", "#")):
            end_v2 = index
            break

    redis_line = next(
        (
            index
            for index in range(v2_line + 1, end_v2)
            if re.match(r"^\s{4}redis:\s*(?:#.*)?(?:\r?\n)?$", lines[index])
        ),
        None,
    )
    if redis_line is None:
        newline = "\r\n" if "\r\n" in text else "\n"
        lines[v2_line + 1:v2_line + 1] = [
            f"    redis:{newline}",
            f"      host: 127.0.0.1{newline}",
            f"      port: {selected_port}{newline}",
            f"      db: 0{newline}",
        ]
    else:
        end_redis = end_v2
        for index in range(redis_line + 1, end_v2):
            if lines[index].strip() and not lines[index].startswith(("      ", "\t", "#")):
                end_redis = index
                break
        port_line = next(
            (
                index
                for index in range(redis_line + 1, end_redis)
                if re.match(r"^\s{6}port:\s*", lines[index])
            ),
            None,
        )
        if port_line is None:
            newline = "\r\n" if "\r\n" in text else "\n"
            lines.insert(redis_line + 1, f"      port: {selected_port}{newline}")
        else:
            lines[port_line] = re.sub(
                r"^(\s{6}port:\s*)\S*(.*)$",
                rf"\g<1>{selected_port}\2",
                lines[port_line],
            )

    directory = os.path.dirname(path) or "."
    fd, temporary_path = tempfile.mkstemp(prefix=".environment_default.", suffix=".yaml", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.writelines(lines)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_path, path)
    except Exception:
        try:
            os.unlink(temporary_path)
        except OSError:
            pass
        raise


def _runtime_defaults_from_envreqs(
    config: Mapping[str, Any], *, project_root: str, yaml_dir: str
) -> dict[str, Any]:
    """Load V1-shaped ``EnvReqs.Runtime.default_runtime_settings`` when present."""
    envreqs = config.get("EnvReqs")
    if not isinstance(envreqs, Mapping):
        return {}
    runtime_reference = envreqs.get("Runtime")
    if runtime_reference is None:
        return {}
    if not isinstance(runtime_reference, Mapping):
        raise ValueError("EnvReqs.Runtime must be a mapping")
    raw_path = runtime_reference.get("default_runtime_settings")
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError("EnvReqs.Runtime.default_runtime_settings is required")

    defaults_path = decode_path(raw_path, project_root=project_root, base_dir=yaml_dir)
    if not os.path.isfile(defaults_path):
        raise FileNotFoundError(f"runtime default YAML not found: {defaults_path}")
    defaults_document = _load_yaml_document(defaults_path, label="runtime default YAML")
    if not isinstance(defaults_document, Mapping):
        raise ValueError(f"runtime default YAML must contain a mapping: {defaults_path}")
    runtime_defaults = defaults_document.get("Runtime", defaults_document)
    if not isinstance(runtime_defaults, Mapping):
        raise ValueError(f"Runtime in default YAML must be a mapping: {defaults_path}")
    return _canonicalize_v2_settings(runtime_defaults)


def load_task_yaml(path: str) -> dict[str, Any]:
    """Load a task YAML file and attach normalized layout metadata."""
    task_path = os.path.abspath(str(path))
    if not os.path.isfile(task_path):
        raise FileNotFoundError(f"task YAML not found: {task_path}")

    loaded = _load_yaml_document(task_path, label="task YAML")
    if not isinstance(loaded, dict):
        raise TaskCardLoadError(
            "JV2-YAML-002",
            f"task YAML must contain a mapping at top level: {task_path}",
            "Start the file with top-level section names such as Scan and Sampling, not a list or scalar.",
            "Scan:\n  name: my_scan\nSampling:\n  Method: Random",
        )

    yaml_dir = os.path.dirname(task_path)
    project_root = infer_project_root(yaml_dir)
    if "Runtime" in loaded:
        raise ValueError(
            "top-level Runtime is no longer a V2 YAML interface; move workers and "
            "batch_size to EnvReqs.V2"
        )

    task_envreqs = loaded.get("EnvReqs")
    if task_envreqs is not None and not isinstance(task_envreqs, Mapping):
        raise ValueError("EnvReqs must be a mapping")
    task_envreqs = task_envreqs if isinstance(task_envreqs, Mapping) else {}
    task_v2 = task_envreqs.get("V2")
    if task_v2 is not None and not isinstance(task_v2, Mapping):
        raise ValueError("EnvReqs.V2 must be a mapping")

    v2_defaults = _v2_defaults_from_envreqs(
        loaded, project_root=project_root, yaml_dir=yaml_dir
    )
    defaults_path = _default_environment_yaml_path(
        loaded, project_root=project_root, yaml_dir=yaml_dir
    )
    legacy_environment_defaults = _legacy_environment_requirements_from_defaults(
        loaded, project_root=project_root, yaml_dir=yaml_dir
    )
    # Legacy EnvReqs.Runtime defaults may include V1 Runtime keys (mode, …);
    # only V2 knobs are retained so old defaults files do not hard-fail.
    runtime_defaults_raw = _runtime_defaults_from_envreqs(
        loaded, project_root=project_root, yaml_dir=yaml_dir
    )
    runtime_defaults = {
        key: value
        for key, value in runtime_defaults_raw.items()
        if key in SUPPORTED_ENVREQS_V2_KEYS
    }
    task_v2_settings = _canonicalize_v2_settings(task_v2) if isinstance(task_v2, Mapping) else {}

    for source_name, source in (
        ("default environment YAML EnvReqs.V2", v2_defaults),
        ("task EnvReqs.V2", task_v2_settings),
    ):
        unsupported = set(source) - SUPPORTED_ENVREQS_V2_KEYS
        if unsupported:
            keys = ", ".join(sorted(str(key) for key in unsupported))
            supported = ", ".join(sorted(SUPPORTED_ENVREQS_V2_KEYS))
            raise ValueError(
                f"unsupported {source_name} setting(s): {keys}; "
                f"supported settings are {supported}"
            )

    v2_settings = _deep_merge(_deep_merge(v2_defaults, runtime_defaults), task_v2_settings)

    config = LoadedTaskConfig(deepcopy(loaded))
    # Validation must identify text the user actually wrote, not values derived
    # below such as ``scan_name`` and absolute result-directory paths.
    config.raw_task_card = deepcopy(loaded)
    resolved_envreqs = deepcopy(dict(task_envreqs))
    for key, default in legacy_environment_defaults.items():
        existing = resolved_envreqs.get(key)
        resolved_envreqs[key] = (
            _deep_merge(default, existing) if isinstance(existing, Mapping) else default
        )
    if v2_settings:
        resolved_envreqs["V2"] = v2_settings
    if resolved_envreqs:
        config["EnvReqs"] = resolved_envreqs
    scan_block = config.get("Scan") if isinstance(config.get("Scan"), Mapping) else {}
    scan_block = dict(scan_block)
    # Apply EnvReqs.V2 sample_directory defaults unless Scan already defines them.
    sample_directory = v2_settings.get("sample_directory")
    if isinstance(sample_directory, Mapping):
        existing = scan_block.get("sample_directory")
        if isinstance(existing, Mapping):
            scan_block["sample_directory"] = _deep_merge(sample_directory, existing)
        else:
            scan_block["sample_directory"] = deepcopy(sample_directory)
        config["Scan"] = scan_block
    # Map Cleanup / Archiver defaults into Calculators.* (task YAML wins).
    calculators = config.get("Calculators")
    calculators = dict(calculators) if isinstance(calculators, Mapping) else {}
    cleanup_defaults = v2_settings.get("cleanup")
    if isinstance(cleanup_defaults, Mapping):
        existing_cleanup = calculators.get("Cleanup")
        if isinstance(existing_cleanup, Mapping):
            calculators["Cleanup"] = _deep_merge(cleanup_defaults, existing_cleanup)
        else:
            calculators["Cleanup"] = deepcopy(cleanup_defaults)
    archiver_defaults = v2_settings.get("archiver")
    if isinstance(archiver_defaults, Mapping):
        existing_archiver = calculators.get("Archiver")
        if isinstance(existing_archiver, Mapping):
            calculators["Archiver"] = _deep_merge(archiver_defaults, existing_archiver)
        else:
            calculators["Archiver"] = deepcopy(archiver_defaults)
    if calculators:
        config["Calculators"] = calculators

    scan_name = str(scan_block.get("name") or config.get("scan_name") or "default").strip() or "default"
    task_result_dir = str(
        config.get("task_result_dir")
        or scan_output_root(project_root=project_root, scan_name=scan_name)
    )
    # Runtime adapter: scheduling scalars + optional redis/factory/worker groups.
    runtime_payload: dict[str, Any] = {
        "mode": "redis",
        **{
            key: value
            for key, value in v2_settings.items()
            if key in {"workers", "batch_size", "checkpoint_heartbeat_sec"}
        },
    }
    if isinstance(v2_settings.get("redis"), Mapping):
        runtime_payload["redis"] = normalize_redis_config(v2_settings["redis"])

    factory_settings = v2_settings.get("factory")
    if isinstance(factory_settings, Mapping):
        factory_norm = normalize_factory_block(factory_settings)
        runtime_payload["Factory"] = factory_norm
        watchdog = factory_norm.get("watchdog")
        if isinstance(watchdog, Mapping):
            runtime_payload["Watchdog"] = deepcopy(watchdog)

    worker_settings = v2_settings.get("worker")
    if isinstance(worker_settings, Mapping):
        worker_norm = normalize_worker_block(worker_settings)
        runtime_payload["sample_artifacts"] = worker_norm["sample_artifacts"]
        runtime_payload["force_serial_layers"] = worker_norm["force_serial_layers"]

    runtime = normalize_runtime_block(runtime_payload)
    if isinstance(sample_directory, Mapping):
        runtime["sample_directory"] = deepcopy(
            scan_block.get("sample_directory") or sample_directory
        )
    config["Runtime"] = runtime
    config["task_yaml"] = task_path
    if defaults_path is not None:
        config["environment_defaults_path"] = os.path.abspath(defaults_path)
    task_redis = task_v2.get("redis") if isinstance(task_v2, Mapping) else None
    config["redis_port_task_override"] = bool(
        isinstance(task_redis, Mapping) and task_redis.get("port") is not None
    )
    config["task_root"] = project_root
    config["project_root"] = project_root
    config["task_result_dir"] = os.path.abspath(task_result_dir)
    config["scan_name"] = scan_name
    return config


def sampling_method(config: Mapping[str, Any]) -> str:
    sampling = dict(config.get("Sampling") or {})
    return str(sampling.get("Method") or "").strip()


def is_check_modules_task(config: Mapping[str, Any]) -> bool:
    sampling = dict(config.get("Sampling") or {})
    mode = str(sampling.get("mode") or "").strip().lower()
    return mode in {"check_modules", "check-modules"}


def resolve_sampling_path(config: Mapping[str, Any], raw: str) -> str:
    project_root = str(config.get("project_root") or config.get("task_root") or os.getcwd())
    if raw.startswith("&J/") or raw.startswith("&J"):
        from jarvishep2.base import expand_j

        return os.path.abspath(expand_j(raw, project_root=project_root))
    if os.path.isabs(raw):
        return raw
    task_yaml = str(config.get("task_yaml") or "")
    anchor = os.path.dirname(task_yaml) if task_yaml else project_root
    return os.path.abspath(os.path.join(anchor, raw))


def _coerce_check_timeout_sec(value: Any, *, default: float) -> float:
    """Return a positive timeout in seconds (minimum 0.1s)."""
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float(default)
    if parsed <= 0:
        return float(default)
    return max(0.1, parsed)


def get_check_modules_settings(config: Mapping[str, Any]) -> dict[str, Any]:
    """Return normalized check-modules knobs (task Sampling + EnvReqs.V2 + defaults)."""
    settings = dict(CHECK_MODULES_DEFAULTS)
    envreqs = config.get("EnvReqs") if isinstance(config.get("EnvReqs"), Mapping) else {}
    v2 = envreqs.get("V2") if isinstance(envreqs, Mapping) else None
    block = v2.get("check_modules") if isinstance(v2, Mapping) else None
    if isinstance(block, Mapping):
        if block.get("data") is not None and str(block.get("data")).strip():
            settings["data"] = str(block.get("data")).strip()
        if "n_samples" in block and block.get("n_samples") is not None:
            try:
                settings["n_samples"] = max(1, int(block.get("n_samples")))
            except (TypeError, ValueError):
                pass
        # Accept timeout_sec (preferred) or timeout as alias.
        raw_timeout = block.get("timeout_sec", block.get("timeout"))
        if raw_timeout is not None:
            settings["timeout_sec"] = _coerce_check_timeout_sec(
                raw_timeout,
                default=float(CHECK_MODULES_DEFAULTS["timeout_sec"]),
            )
    # Task Sampling.data always wins for the CSV path when present.
    sampling = config.get("Sampling") if isinstance(config.get("Sampling"), Mapping) else {}
    if isinstance(sampling, Mapping):
        raw = sampling.get("data") or sampling.get("points_csv")
        if raw is not None and str(raw).strip():
            settings["data"] = str(raw).strip()
    return settings


def resolve_check_modules_csv(config: Mapping[str, Any]) -> tuple[str | None, str]:
    """Resolve the check-modules CSV path.

    Returns
    -------
    (path_or_none, raw_spec)
        *path_or_none* is an absolute path when the file exists; ``None`` when
        the configured path is missing (caller should fall back to sampler).
        *raw_spec* is the unresolved YAML string (for logging).
    """
    settings = get_check_modules_settings(config)
    raw = str(settings.get("data") or CHECK_MODULES_DEFAULTS["data"]).strip()
    if not raw:
        return None, ""
    path = resolve_sampling_path(config, raw)
    if os.path.isfile(path):
        return path, raw
    return None, raw


def check_modules_points_path(config: Mapping[str, Any]) -> str:
    """Return an existing check-modules CSV path or raise.

    Prefer :func:`resolve_check_modules_csv` when a missing file should fall back
    to sampler-drawn smoke points.
    """
    path, raw = resolve_check_modules_csv(config)
    if path is None:
        raise ValueError(
            "check-modules CSV not found"
            + (f" (configured: {raw!r})" if raw else "")
            + "; set Sampling.data / EnvReqs.V2.check_modules.data or rely on "
            "sampler fallback (n_samples)"
        )
    return path


def check_modules_n_samples(config: Mapping[str, Any]) -> int:
    """How many sampler-drawn smoke points to use when CSV is missing."""
    settings = get_check_modules_settings(config)
    try:
        return max(1, int(settings.get("n_samples", 10)))
    except (TypeError, ValueError):
        return int(CHECK_MODULES_DEFAULTS["n_samples"])


def check_modules_timeout_sec(config: Mapping[str, Any]) -> float:
    """Max seconds to wait for check-modules samples to archive (deadline)."""
    settings = get_check_modules_settings(config)
    return _coerce_check_timeout_sec(
        settings.get("timeout_sec"),
        default=float(CHECK_MODULES_DEFAULTS["timeout_sec"]),
    )


__all__ = [
    "LoadedTaskConfig",
    "TaskCardLoadError",
    "check_modules_n_samples",
    "check_modules_points_path",
    "check_modules_timeout_sec",
    "get_check_modules_settings",
    "is_check_modules_task",
    "load_task_yaml",
    "resolve_check_modules_csv",
    "resolve_sampling_path",
    "sampling_method",
]
