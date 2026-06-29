#!/usr/bin/env python3
"""Task YAML loading and path normalization for Jarvis2 CLI runs."""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

import yaml

from jarvishep2.base import infer_project_root, scan_output_root
from jarvishep2.runtime_config import normalize_runtime_block


def load_task_yaml(path: str) -> dict[str, Any]:
    """Load a task YAML file and attach normalized layout metadata."""
    task_path = os.path.abspath(str(path))
    if not os.path.isfile(task_path):
        raise FileNotFoundError(f"task YAML not found: {task_path}")

    with open(task_path, "r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"task YAML must contain a mapping at top level: {task_path}")

    config = dict(loaded)
    yaml_dir = os.path.dirname(task_path)
    project_root = infer_project_root(yaml_dir)
    scan_block = config.get("Scan") if isinstance(config.get("Scan"), Mapping) else {}
    scan_name = str(scan_block.get("name") or config.get("scan_name") or "default").strip() or "default"
    task_result_dir = str(
        config.get("task_result_dir")
        or scan_output_root(project_root=project_root, scan_name=scan_name)
    )
    runtime = normalize_runtime_block(config.get("Runtime"))
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


def check_modules_points_path(config: Mapping[str, Any]) -> str:
    sampling = dict(config.get("Sampling") or {})
    raw = str(sampling.get("data") or sampling.get("points_csv") or "").strip()
    if not raw:
        raise ValueError("check-modules task requires Sampling.data pointing to a CSV file")
    project_root = str(config.get("project_root") or config.get("task_root") or os.getcwd())
    if raw.startswith("&J/") or raw.startswith("&J"):
        from jarvishep2.base import expand_j

        return os.path.abspath(expand_j(raw, project_root=project_root))
    if os.path.isabs(raw):
        return raw
    task_yaml = str(config.get("task_yaml") or "")
    anchor = os.path.dirname(task_yaml) if task_yaml else project_root
    return os.path.abspath(os.path.join(anchor, raw))


__all__ = [
    "check_modules_points_path",
    "is_check_modules_task",
    "load_task_yaml",
    "sampling_method",
]