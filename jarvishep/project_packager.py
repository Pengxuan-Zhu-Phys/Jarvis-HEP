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

import yaml


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


@dataclass
class ProjectPackManifestReport:
    project_root: str
    manifest_path: str
    pack_id: str
    profile: str
    included_files: int
    excluded_files: int
    output: str


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


def _next_pack_id(output_dir: str) -> tuple[str, str]:
    date_part = datetime.now().strftime("%Y%m%d")
    for idx in range(1, 1000):
        pack_id = f"pack_{date_part}_{idx:03d}"
        manifest_path = os.path.join(output_dir, f"{pack_id}.yaml")
        if not os.path.exists(manifest_path):
            return pack_id, manifest_path
    raise ProjectPackError(f"No available pack manifest name in: {output_dir}")


def _manifest_exclude_entries(excluded: list[tuple[str, str]]) -> list[str]:
    entries: list[str] = []
    seen: set[str] = set()
    for rel, reason in excluded:
        if reason.startswith("exclude-dir:"):
            entry = rel.rstrip("/") + "/"
        else:
            entry = rel
        if entry not in seen:
            seen.add(entry)
            entries.append(entry)
    return entries


def create_project_pack_manifest(
    project_root: str | None,
    profile: str = "repro",
    manifest_dir: str | None = None,
) -> ProjectPackManifestReport:
    target = os.getcwd() if project_root is None else str(project_root)
    target = os.path.abspath(os.path.expanduser(target))
    profile = _normalize_profile(profile)

    if not os.path.isdir(target):
        raise ProjectNotFoundError(f"Project directory not found: {target}")
    if not _is_project_root(target):
        raise ProjectNotFoundError("Jarvis project not found in target path.")

    included, excluded, _total_bytes = _collect_project_files(target, profile=profile)
    if not included:
        raise ProjectPackError("No files selected for manifest.")

    output_dir = os.getcwd() if manifest_dir is None else str(manifest_dir)
    output_dir = os.path.abspath(os.path.expanduser(output_dir))
    os.makedirs(output_dir, exist_ok=True)
    pack_id, manifest_path = _next_pack_id(output_dir)

    project_name = os.path.basename(os.path.normpath(target))
    output_name = f"{project_name}_{profile}_{pack_id}.tar.gz"
    exclude_entries = _manifest_exclude_entries(excluded)
    payload = {
        "pack_id": pack_id,
        "mode": profile,
        "project_root": target,
        "output": output_name,
        "include": included,
        "exclude": exclude_entries,
    }

    with open(manifest_path, "w", encoding="utf-8") as f1:
        yaml.safe_dump(payload, f1, sort_keys=False)

    return ProjectPackManifestReport(
        project_root=target,
        manifest_path=manifest_path,
        pack_id=pack_id,
        profile=profile,
        included_files=len(included),
        excluded_files=len(exclude_entries),
        output=output_name,
    )


def _load_manifest_yaml(manifest_path: str) -> dict:
    try:
        with open(manifest_path, "r", encoding="utf-8") as f1:
            payload = yaml.safe_load(f1)
    except OSError as exc:
        raise ProjectPackError(f"Failed to read pack manifest: {exc}") from exc
    except yaml.YAMLError as exc:
        raise ProjectPackError(f"Invalid pack manifest YAML: {exc}") from exc

    if not isinstance(payload, dict):
        raise ProjectPackError("Pack manifest must contain a YAML mapping.")
    return payload


def _manifest_list(payload: dict, key: str) -> list[str]:
    value = payload.get(key, [])
    if value is None:
        return []
    if not isinstance(value, list):
        raise ProjectPackError(f"Pack manifest '{key}' must be a list.")
    entries: list[str] = []
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ProjectPackError(f"Pack manifest '{key}' contains an invalid path.")
        entries.append(item.strip().replace("\\", "/"))
    return entries


def _validate_manifest_relpath(relpath: str, key: str) -> str:
    if os.path.isabs(relpath) or (len(relpath) >= 2 and relpath[1] == ":"):
        raise ProjectPackError(f"Pack manifest '{key}' path must be relative: {relpath}")
    is_dir = relpath.endswith("/")
    parts = relpath.replace("\\", "/").split("/")
    if any(part == ".." for part in parts):
        raise ProjectPackError(f"Pack manifest '{key}' path cannot traverse upward: {relpath}")
    cleaned = "/".join(part for part in parts if part not in ("", "."))
    if not cleaned:
        raise ProjectPackError(f"Pack manifest '{key}' contains an empty path.")
    if is_dir:
        return cleaned.rstrip("/") + "/"
    return cleaned


def _is_manifest_excluded(relpath: str, excludes: list[str]) -> bool:
    rel = relpath.rstrip("/")
    for raw in excludes:
        is_dir = raw.endswith("/")
        ex = raw.rstrip("/")
        if rel == ex:
            return True
        if is_dir and rel.startswith(ex + "/"):
            return True
    return False


def create_project_package_from_manifest(manifest_path: str) -> ProjectPackReport:
    manifest = os.path.abspath(os.path.expanduser(str(manifest_path)))
    payload = _load_manifest_yaml(manifest)

    project_root = payload.get("project_root")
    if not isinstance(project_root, str) or not project_root.strip():
        raise ProjectPackError("Pack manifest requires 'project_root'.")
    project_root = os.path.abspath(os.path.expanduser(project_root.strip()))
    if not os.path.isdir(project_root):
        raise ProjectNotFoundError(f"Project directory not found: {project_root}")

    mode = _normalize_profile(str(payload.get("mode") or "repro"))
    output = payload.get("output")
    if not isinstance(output, str) or not output.strip():
        raise ProjectPackError("Pack manifest requires 'output'.")
    output = os.path.expanduser(output.strip())
    archive_path = output if os.path.isabs(output) else os.path.abspath(output)

    includes = [
        _validate_manifest_relpath(rel, "include")
        for rel in _manifest_list(payload, "include")
    ]
    excludes = [
        _validate_manifest_relpath(rel, "exclude")
        for rel in _manifest_list(payload, "exclude")
    ]

    for rel in includes:
        src = os.path.join(project_root, rel)
        if not os.path.isfile(src):
            raise ProjectPackError(f"Included file does not exist: {rel}")

    final_files = [
        rel for rel in includes
        if not _is_manifest_excluded(rel, excludes)
    ]
    if not final_files:
        raise ProjectPackError("No files selected for packaging.")

    output_parent = os.path.dirname(archive_path)
    if output_parent:
        os.makedirs(output_parent, exist_ok=True)

    total_bytes = 0
    with tarfile.open(archive_path, "w:gz") as tar:
        for rel in final_files:
            src = os.path.join(project_root, rel)
            total_bytes += max(0, int(os.path.getsize(src)))
            tar.add(src, arcname=rel, recursive=False)

    return ProjectPackReport(
        project_root=project_root,
        archive_path=archive_path,
        profile=mode,
        included_files=len(final_files),
        excluded_files=len(includes) - len(final_files),
        total_bytes=total_bytes,
    )


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
