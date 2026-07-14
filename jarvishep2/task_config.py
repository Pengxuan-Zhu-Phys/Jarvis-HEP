#!/usr/bin/env python3
"""Task YAML loading and path normalization for Jarvis2 CLI runs."""

from __future__ import annotations

import os
from collections.abc import Mapping
from copy import deepcopy
from typing import Any

import yaml

from jarvishep2.base import decode_path, infer_project_root, scan_output_root
from jarvishep2.runtime_config import normalize_runtime_block


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
    """Normalize the singular ``worker`` compatibility spelling before merging."""
    normalized = deepcopy(dict(settings))
    if "worker" in normalized:
        if "workers" in normalized:
            raise ValueError("EnvReqs.V2 accepts either workers or worker, not both")
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
    with open(defaults_path, "r", encoding="utf-8") as handle:
        defaults_document = yaml.safe_load(handle)
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
    with open(defaults_path, "r", encoding="utf-8") as handle:
        defaults_document = yaml.safe_load(handle)
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

    with open(task_path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"task YAML must contain a mapping at top level: {task_path}")

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
    runtime_defaults = _runtime_defaults_from_envreqs(
        loaded, project_root=project_root, yaml_dir=yaml_dir
    )
    task_v2_settings = _canonicalize_v2_settings(task_v2) if isinstance(task_v2, Mapping) else {}
    v2_settings = _deep_merge(_deep_merge(v2_defaults, runtime_defaults), task_v2_settings)
    # Scheduling knobs + SAMPLE archive policy (EnvReqs.V2 / default_yaml).
    _SUPPORTED_V2 = {
        "workers",
        "batch_size",
        "sample_directory",
        "cleanup",
        "archiver",
    }
    unsupported_v2_settings = set(v2_settings) - _SUPPORTED_V2
    if unsupported_v2_settings:
        unsupported = ", ".join(sorted(str(key) for key in unsupported_v2_settings))
        raise ValueError(
            "unsupported EnvReqs.V2 setting(s): "
            f"{unsupported}; supported settings are "
            + ", ".join(sorted(_SUPPORTED_V2))
        )

    config = deepcopy(loaded)
    resolved_envreqs = deepcopy(dict(task_envreqs))
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
    # Runtime adapter keeps only scheduling scalars; nested archive policy lives
    # on Scan / Calculators (and EnvReqs.V2 for inspection).
    runtime_scalars = {
        key: value
        for key, value in v2_settings.items()
        if key in {"workers", "batch_size"}
    }
    runtime = normalize_runtime_block({"mode": "redis", **runtime_scalars})
    if isinstance(sample_directory, Mapping):
        runtime["sample_directory"] = deepcopy(
            scan_block.get("sample_directory") or sample_directory
        )
    config["Runtime"] = runtime
    config["task_yaml"] = task_path
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


def check_modules_points_path(config: Mapping[str, Any]) -> str:
    sampling = dict(config.get("Sampling") or {})
    raw = str(sampling.get("data") or sampling.get("points_csv") or "").strip()
    if not raw:
        raise ValueError("check-modules task requires Sampling.data pointing to a CSV file")
    return resolve_sampling_path(config, raw)


__all__ = [
    "check_modules_points_path",
    "is_check_modules_task",
    "load_task_yaml",
    "resolve_sampling_path",
    "sampling_method",
]
