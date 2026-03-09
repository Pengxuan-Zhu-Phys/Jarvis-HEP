#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
import sys
import tarfile
import tempfile
from importlib import metadata as importlib_metadata


PROJECT_MARKERS = (".jarvis-project.json", "jarvis.project.yaml")
ALLOWED_PROFILES = ("share", "repro", "full")

ALWAYS_EXCLUDE_DIRS = {
    "__pycache__",
    ".git",
    ".mypy_cache",
    ".pytest_cache",
    ".ruff_cache",
}
ALWAYS_EXCLUDE_BASENAMES = {
    ".DS_Store",
}
ALWAYS_EXCLUDE_SUFFIXES = (
    ".pyc",
    ".pyo",
    ".tmp",
    ".swp",
    ".lock",
)

SHARE_TOPLEVEL_DIRS = {
    "bin",
    "data",
    "outputs",
    "images",
}
SHARE_TOPLEVEL_FILES = {
    ".jarvis-project.json",
    "jarvis.project.yaml",
    "README.md",
}


class ProjectPackError(RuntimeError):
    pass


class ProjectNotFoundError(ProjectPackError):
    pass


class ProjectProfileError(ProjectPackError):
    pass


@dataclass
class ProjectPackReport:
    project_root: str
    archive_path: str
    profile: str
    included_files: int
    excluded_files: int
    total_bytes: int


def _jarvis_version() -> str:
    for candidate in ("Jarvis-HEP", "jarvis-hep"):
        try:
            return str(importlib_metadata.version(candidate))
        except importlib_metadata.PackageNotFoundError:
            continue
        except Exception:
            continue
    return "unknown"


def _is_project_root(path: str) -> bool:
    for marker in PROJECT_MARKERS:
        if os.path.isfile(os.path.join(path, marker)):
            return True
    return False


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f1:
        while True:
            chunk = f1.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _normalize_profile(profile: str) -> str:
    p = str(profile or "").strip().lower()
    if not p:
        p = "repro"
    if p not in ALLOWED_PROFILES:
        raise ProjectProfileError(
            f"Unsupported profile '{profile}'. Expected one of: "
            + ", ".join(ALLOWED_PROFILES)
        )
    return p


def _common_exclude(parts: list[str], basename: str) -> str | None:
    for part in parts[:-1]:
        if part in ALWAYS_EXCLUDE_DIRS:
            return f"exclude-dir:{part}"
    if basename in ALWAYS_EXCLUDE_BASENAMES:
        return f"exclude-name:{basename}"
    for suffix in ALWAYS_EXCLUDE_SUFFIXES:
        if basename.endswith(suffix):
            return f"exclude-suffix:{suffix}"
    return None


def _select_profile(profile: str, relpath: str) -> tuple[bool, str]:
    parts = relpath.split("/")
    top = parts[0]
    basename = parts[-1]

    common = _common_exclude(parts, basename)
    if common is not None:
        return False, common

    if profile == "full":
        return True, "include:full"

    if profile == "share":
        if len(parts) == 1:
            if top in SHARE_TOPLEVEL_FILES:
                return True, "include:share-root"
            return False, "share-root-filter"
        if top in SHARE_TOPLEVEL_DIRS:
            return True, "include:share-dir"
        return False, "share-dir-filter"

    # repro (default): keep all non-hidden top-level entries and explicit markers.
    if len(parts) == 1:
        if top in PROJECT_MARKERS:
            return True, "include:repro-marker"
        if top.startswith("."):
            return False, "repro-hidden-root-file"
        return True, "include:repro-root"

    if top.startswith("."):
        return False, "repro-hidden-root-dir"
    return True, "include:repro"


def _collect_project_files(project_root: str, profile: str) -> tuple[list[str], list[tuple[str, str]], int]:
    included: list[str] = []
    excluded: list[tuple[str, str]] = []
    total_bytes = 0

    for root, dirs, files in os.walk(project_root):
        dirs.sort()
        files.sort()

        rel_root = os.path.relpath(root, project_root)
        if rel_root == ".":
            prefix_parts: list[str] = []
        else:
            prefix_parts = rel_root.replace("\\", "/").split("/")

        # Prune always-excluded directories early.
        kept_dirs = []
        for dname in dirs:
            if dname in ALWAYS_EXCLUDE_DIRS:
                rel = "/".join(prefix_parts + [dname])
                excluded.append((rel, f"exclude-dir:{dname}"))
                continue
            kept_dirs.append(dname)
        dirs[:] = kept_dirs

        for fname in files:
            rel = "/".join(prefix_parts + [fname]) if prefix_parts else fname
            keep, reason = _select_profile(profile, rel)
            if not keep:
                excluded.append((rel, reason))
                continue
            fpath = os.path.join(project_root, rel)
            try:
                size = int(os.path.getsize(fpath))
            except OSError:
                size = 0
            total_bytes += max(0, size)
            included.append(rel)

    included.sort()
    excluded.sort(key=lambda x: x[0])
    return included, excluded, total_bytes


def _write_metadata_files(
    meta_root: str,
    *,
    project_name: str,
    profile: str,
    archive_name: str,
    included: list[str],
    excluded: list[tuple[str, str]],
    checksums: list[str],
    total_bytes: int,
) -> dict[str, str]:
    os.makedirs(meta_root, exist_ok=True)
    manifest_path = os.path.join(meta_root, "manifest.json")
    checksums_path = os.path.join(meta_root, "checksums.sha256")
    included_path = os.path.join(meta_root, "included_files.txt")
    excluded_path = os.path.join(meta_root, "excluded_files.txt")

    manifest = {
        "format": "jarvis-hep-project-package",
        "version": 1,
        "project_name": project_name,
        "profile": profile,
        "archive_name": archive_name,
        "created_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "jarvis_hep_version": _jarvis_version(),
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "included_files": len(included),
        "excluded_files": len(excluded),
        "payload_bytes": int(total_bytes),
    }
    with open(manifest_path, "w", encoding="utf-8") as f1:
        json.dump(manifest, f1, indent=2, ensure_ascii=False)
        f1.write("\n")

    with open(checksums_path, "w", encoding="utf-8") as f1:
        for line in checksums:
            f1.write(line + "\n")

    with open(included_path, "w", encoding="utf-8") as f1:
        for rel in included:
            f1.write(rel + "\n")

    with open(excluded_path, "w", encoding="utf-8") as f1:
        for rel, reason in excluded:
            f1.write(f"{reason}\t{rel}\n")

    return {
        "manifest": manifest_path,
        "checksums": checksums_path,
        "included": included_path,
        "excluded": excluded_path,
    }


def create_project_package(project_root: str | None, profile: str = "repro") -> ProjectPackReport:
    target = os.getcwd() if project_root is None else str(project_root)
    target = os.path.abspath(os.path.expanduser(target))
    profile = _normalize_profile(profile)

    if not os.path.isdir(target):
        raise ProjectNotFoundError(f"Project directory not found: {target}")
    if not _is_project_root(target):
        raise ProjectNotFoundError("Jarvis project not found in target path.")

    included, excluded, total_bytes = _collect_project_files(target, profile=profile)
    if not included:
        raise ProjectPackError("No files selected for packaging.")

    project_name = os.path.basename(os.path.normpath(target))
    output_root = os.path.dirname(target)
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    archive_name = f"{project_name}_{profile}_{timestamp}.tar.gz"
    archive_path = os.path.join(output_root, archive_name)

    checksums: list[str] = []
    with tarfile.open(archive_path, "w:gz") as tar:
        for rel in included:
            src = os.path.join(target, rel)
            checksums.append(f"{_sha256(src)}  {rel}")
            arcname = os.path.join(project_name, rel)
            tar.add(src, arcname=arcname, recursive=False)

        with tempfile.TemporaryDirectory() as tmpdir:
            meta_paths = _write_metadata_files(
                os.path.join(tmpdir, ".jarvis-pack"),
                project_name=project_name,
                profile=profile,
                archive_name=archive_name,
                included=included,
                excluded=excluded,
                checksums=checksums,
                total_bytes=total_bytes,
            )
            for _, mpath in sorted(meta_paths.items()):
                rel = os.path.relpath(mpath, tmpdir).replace("\\", "/")
                arcname = os.path.join(project_name, rel)
                tar.add(mpath, arcname=arcname, recursive=False)

    return ProjectPackReport(
        project_root=target,
        archive_path=archive_path,
        profile=profile,
        included_files=len(included),
        excluded_files=len(excluded),
        total_bytes=total_bytes,
    )

