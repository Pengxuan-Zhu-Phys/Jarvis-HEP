#!/usr/bin/env python3
from __future__ import annotations

import json
import os
from typing import Optional

from importlib import resources as importlib_resources


PROJECT_SUBDIRS = (
    "bin",
    "data",
    "outputs",
    "calculators",
    "deps",
    "images",
    "logs",
    "checkpoints",
)
LEGACY_COMPAT_SUBDIRS = ()
PROJECT_MARKER_NAME = ".jarvis-project.json"
PROJECT_DESCRIPTOR_NAME = "jarvis.project.yaml"
_TEMPLATE_PACKAGE = "jarvishep"
_TEMPLATE_SUBDIR = "project_template"


def _contains_path_separator(name: str) -> bool:
    separators = [os.sep]
    if os.altsep:
        separators.append(os.altsep)
    return any(sep in name for sep in separators)


def _copy_template_tree(src, dst: str) -> None:
    for entry in src.iterdir():
        if entry.name in {"__pycache__", ".DS_Store"}:
            continue
        target = os.path.join(dst, entry.name)
        if entry.is_dir():
            os.makedirs(target, exist_ok=True)
            _copy_template_tree(entry, target)
            continue
        with entry.open("rb") as f_src, open(target, "wb") as f_dst:
            f_dst.write(f_src.read())


def _write_project_marker(project_root: str) -> None:
    marker_payload = {
        "format": "jarvis-hep-standalone-project",
        "version": 1,
        "layout": list(PROJECT_SUBDIRS),
    }
    marker_path = os.path.join(project_root, PROJECT_MARKER_NAME)
    with open(marker_path, "w", encoding="utf-8") as f1:
        json.dump(marker_payload, f1, indent=2, ensure_ascii=False)


def _write_project_descriptor(project_root: str) -> None:
    descriptor_path = os.path.join(project_root, PROJECT_DESCRIPTOR_NAME)
    descriptor = [
        "project:",
        "  format: jarvis-hep-standalone-project",
        "  version: 1",
        "  output_root: outputs",
        "  path_markers:",
        "    task_root: '&J'",
        "    package_root: '&SRC'",
        "  defaults:",
        "    main_config: bin/quickstart_mcmc_operas.yaml",
    ]
    with open(descriptor_path, "w", encoding="utf-8") as f1:
        f1.write("\n".join(descriptor) + "\n")


def _bootstrap_minimal_files(project_root: str) -> None:
    readme_path = os.path.join(project_root, "README.md")
    if not os.path.exists(readme_path):
        with open(readme_path, "w", encoding="utf-8") as f1:
            f1.write(
                "# Jarvis-HEP Project\n\n"
                "Create this project via `Jarvis --mkproject <name>`.\n"
                "Put runnable YAML cards under `bin/` and run with:\n\n"
                "```bash\n"
                "Jarvis bin/quickstart_mcmc_operas.yaml\n"
                "```\n"
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
    name = project_name.strip()
    if not name:
        raise ValueError("--mkproject requires a non-empty project name.")
    if _contains_path_separator(name):
        raise ValueError("Project name should be a folder name, not a path.")

    base_dir = os.getcwd() if cwd is None else cwd
    project_root = os.path.abspath(os.path.join(base_dir, name))
    if os.path.exists(project_root):
        raise FileExistsError(project_root)

    os.makedirs(project_root, exist_ok=False)
    for subdir in PROJECT_SUBDIRS:
        os.makedirs(os.path.join(project_root, subdir), exist_ok=False)
    for subdir in LEGACY_COMPAT_SUBDIRS:
        os.makedirs(os.path.join(project_root, subdir), exist_ok=False)
    _apply_project_template(project_root)
    _write_project_marker(project_root)
    _write_project_descriptor(project_root)

    return project_root
