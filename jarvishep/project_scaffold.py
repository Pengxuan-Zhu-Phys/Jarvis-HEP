#!/usr/bin/env python3
from __future__ import annotations

import os
from typing import Optional

PROJECT_SUBDIRS = ("bin", "Library", "Workshop", "Result")


def _contains_path_separator(name: str) -> bool:
    separators = [os.sep]
    if os.altsep:
        separators.append(os.altsep)
    return any(sep in name for sep in separators)


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

    return project_root
