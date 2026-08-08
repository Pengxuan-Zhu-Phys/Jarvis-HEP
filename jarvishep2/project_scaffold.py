#!/usr/bin/env python3
"""Standalone project scaffold for ``Jarvis project create`` (D12.5 / V1 layout)."""

from __future__ import annotations

import json
import os
from importlib import resources as importlib_resources
from typing import Optional

PROJECT_SUBDIRS = (
    "bin",
    "data",
    "deps",
)
PROJECT_NESTED_SUBDIRS: tuple[str, ...] = ()
LEGACY_COMPAT_SUBDIRS: tuple[str, ...] = ()
PROJECT_MARKER_NAME = ".jarvis-project.json"
PROJECT_DESCRIPTOR_NAME = "jarvis.project.yaml"
_TEMPLATE_PACKAGE = "jarvishep2"
_TEMPLATE_SUBDIR = "project_template"


def _contains_path_separator(name: str) -> bool:
    separators = [os.sep]
    if os.altsep:
        separators.append(os.altsep)
    return any(sep in name for sep in separators)


def _copy_template_tree(src, dst: str) -> bool:
    copied_any = False
    for entry in src.iterdir():
        if entry.name in {"__pycache__", ".DS_Store"}:
            continue
        target = os.path.join(dst, entry.name)
        if entry.is_dir():
            if _copy_template_tree(entry, target):
                copied_any = True
            continue
        os.makedirs(dst, exist_ok=True)
        with entry.open("rb") as f_src, open(target, "wb") as f_dst:
            f_dst.write(f_src.read())
        copied_any = True
    return copied_any


def _write_project_marker(project_root: str) -> None:
    marker_payload = {
        "format": "jarvis-hep-standalone-project",
        "version": 1,
        "layout": list(PROJECT_SUBDIRS),
        "runtime": "jarvis2",
    }
    marker_path = os.path.join(project_root, PROJECT_MARKER_NAME)
    with open(marker_path, "w", encoding="utf-8") as handle:
        json.dump(marker_payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def _write_project_descriptor(project_root: str) -> None:
    descriptor_path = os.path.join(project_root, PROJECT_DESCRIPTOR_NAME)
    descriptor = [
        "project:",
        "  format: jarvis-hep-standalone-project",
        "  version: 1",
        "  runtime: jarvis2",
        "  output_root: outputs",
        "  artifact_roots:",
        "    outputs: outputs",
        "    logs: logs",
        "    images: images",
        "    checkpoints: checkpoints",
        "  path_markers:",
        "    task_root: '&J'",
        "  defaults:",
        "    main_config: bin/quickstart_bridson_operas.yaml",
    ]
    with open(descriptor_path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(descriptor) + "\n")


def _bootstrap_minimal_files(project_root: str) -> None:
    readme_path = os.path.join(project_root, "README.md")
    if not os.path.exists(readme_path):
        with open(readme_path, "w", encoding="utf-8") as handle:
            handle.write(
                "# Jarvis-HEP V2 Project\n\n"
                "Create this project via `Jarvis project create <name>`.\n"
                "Put runnable YAML cards under `bin/` and run with:\n\n"
                "```bash\n"
                "Jarvis run bin/quickstart_bridson_operas.yaml\n"
                "```\n\n"
                "Runtime artifact directories such as `outputs/`, `logs/`, and "
                "`images/` are created automatically on first use.\n"
            )


def _apply_project_template(project_root: str) -> None:
    try:
        template_root = importlib_resources.files(_TEMPLATE_PACKAGE).joinpath(_TEMPLATE_SUBDIR)
    except Exception:
        _bootstrap_minimal_files(project_root)
        return
    if not template_root.is_dir():
        _bootstrap_minimal_files(project_root)
        return
    _copy_template_tree(template_root, project_root)


def create_project_scaffold(project_name: str, cwd: Optional[str] = None) -> str:
    """Create a standalone project directory with the V1 layout markers."""
    name = project_name.strip()
    if not name:
        raise ValueError("project create requires a non-empty project name.")
    if _contains_path_separator(name):
        raise ValueError("Project name should be a folder name, not a path.")

    base_dir = os.getcwd() if cwd is None else cwd
    project_root = os.path.abspath(os.path.join(base_dir, name))
    if os.path.exists(project_root):
        raise FileExistsError(project_root)

    os.makedirs(project_root, exist_ok=False)
    for subdir in PROJECT_SUBDIRS:
        os.makedirs(os.path.join(project_root, subdir), exist_ok=False)
    for subdir in PROJECT_NESTED_SUBDIRS:
        os.makedirs(os.path.join(project_root, subdir), exist_ok=True)
    for subdir in LEGACY_COMPAT_SUBDIRS:
        os.makedirs(os.path.join(project_root, subdir), exist_ok=False)
    _apply_project_template(project_root)
    _write_project_marker(project_root)
    _write_project_descriptor(project_root)

    return project_root


__all__ = [
    "LEGACY_COMPAT_SUBDIRS",
    "PROJECT_DESCRIPTOR_NAME",
    "PROJECT_MARKER_NAME",
    "PROJECT_NESTED_SUBDIRS",
    "PROJECT_SUBDIRS",
    "create_project_scaffold",
]
