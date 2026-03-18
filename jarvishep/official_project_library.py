#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from importlib import resources as importlib_resources
import json
import os
import shutil
import tarfile
import tempfile
from urllib.parse import urlparse
from urllib.request import urlopen
import zipfile


DEFAULT_OFFICIAL_LIBRARY_NAME = "official Jarvis library"
DEFAULT_OFFICIAL_LIBRARY_INDEX_URL = (
    "https://raw.githubusercontent.com/Pengxuan-Zhu-Phys/Jarvis-HEP/main/"
    "jarvishep/card/official_project_library.json"
)
OFFICIAL_LIBRARY_INDEX_ENV = "JARVIS_OFFICIAL_LIBRARY_INDEX_URL"
OFFICIAL_LIBRARY_TIMEOUT_ENV = "JARVIS_OFFICIAL_LIBRARY_TIMEOUT_SEC"
DEFAULT_TIMEOUT_SEC = 20.0
OFFICIAL_LIBRARY_CARD_PATH = ("card", "official_project_library.json")


class OfficialLibraryError(RuntimeError):
    pass


class OfficialCatalogError(OfficialLibraryError):
    pass


class OfficialProjectNotFoundError(OfficialLibraryError):
    pass


class OfficialProjectFetchError(OfficialLibraryError):
    pass


@dataclass(frozen=True)
class OfficialProjectFetchReport:
    project_name: str
    target_dir: str
    entrypoint: str


def _coerce_timeout(timeout_sec: float | None) -> float:
    if timeout_sec is not None:
        try:
            value = float(timeout_sec)
            if value > 0:
                return value
        except Exception:
            pass
    env_value = os.environ.get(OFFICIAL_LIBRARY_TIMEOUT_ENV, "")
    if env_value:
        try:
            value = float(env_value)
            if value > 0:
                return value
        except Exception:
            pass
    return DEFAULT_TIMEOUT_SEC


def _read_packaged_catalog() -> dict:
    try:
        card = importlib_resources.files("jarvishep")
        for part in OFFICIAL_LIBRARY_CARD_PATH:
            card = card.joinpath(part)
        with card.open("r", encoding="utf-8") as f1:
            payload = json.load(f1)
    except Exception as exc:
        raise OfficialCatalogError(
            f"Cannot load packaged official Jarvis library metadata: {exc}"
        ) from exc
    return payload


def _read_catalog_from_url(index_url: str, timeout_sec: float) -> dict:
    parsed = urlparse(index_url)
    if parsed.scheme not in {"https", "http", "file"}:
        raise OfficialCatalogError(
            f"Unsupported official Jarvis library source URL scheme: {parsed.scheme or '<none>'}"
        )

    try:
        with urlopen(index_url, timeout=timeout_sec) as response:
            payload = json.loads(response.read().decode("utf-8"))
    except Exception as exc:
        raise OfficialCatalogError(
            f"Cannot read official Jarvis library source: {index_url} ({exc})"
        ) from exc

    return payload


def _normalize_project_entry(entry: object) -> dict | None:
    if not isinstance(entry, dict):
        return None

    name = str(entry.get("name") or "").strip()
    if not name:
        return None

    return {
        "name": name,
        "category": str(entry.get("category") or "").strip(),
        "summary": str(entry.get("summary") or "").strip(),
        "entrypoint": str(entry.get("entrypoint") or "").strip(),
        "archive_url": str(entry.get("archive_url") or "").strip(),
        "archive_root": str(entry.get("archive_root") or "").strip(),
        "compatibility_notes": str(entry.get("compatibility_notes") or "").strip(),
    }


def _normalize_catalog(payload: object) -> dict:
    if not isinstance(payload, dict):
        raise OfficialCatalogError("Official Jarvis library metadata must be a JSON object.")

    projects = payload.get("projects")
    if not isinstance(projects, list):
        raise OfficialCatalogError("Official Jarvis library metadata is missing the 'projects' list.")

    normalized = []
    for item in projects:
        project = _normalize_project_entry(item)
        if project is None:
            continue
        normalized.append(project)

    normalized.sort(key=lambda entry: entry["name"].lower())

    return {
        "library_name": str(payload.get("library_name") or DEFAULT_OFFICIAL_LIBRARY_NAME),
        "projects": normalized,
    }


def _load_catalog_payload(index_url: str | None, timeout_sec: float) -> dict:
    effective_url = (
        str(index_url).strip()
        if index_url is not None
        else str(os.environ.get(OFFICIAL_LIBRARY_INDEX_ENV) or "").strip()
    )
    if not effective_url:
        effective_url = DEFAULT_OFFICIAL_LIBRARY_INDEX_URL

    try:
        return _read_catalog_from_url(effective_url, timeout_sec=timeout_sec)
    except OfficialCatalogError:
        return _read_packaged_catalog()


def load_official_project_catalog(
    index_url: str | None = None,
    *,
    timeout_sec: float | None = None,
) -> dict:
    timeout = _coerce_timeout(timeout_sec)
    payload = _load_catalog_payload(index_url=index_url, timeout_sec=timeout)
    return _normalize_catalog(payload)


def list_official_projects(index_url: str | None = None) -> list[dict]:
    catalog = load_official_project_catalog(index_url=index_url)
    return list(catalog["projects"])


def get_official_project(project_name: str, index_url: str | None = None) -> dict:
    name = str(project_name or "").strip()
    if not name:
        raise OfficialProjectNotFoundError("Missing official project name.")

    for project in list_official_projects(index_url=index_url):
        if project["name"].lower() == name.lower():
            return project

    raise OfficialProjectNotFoundError(
        f"Official project not found in the official Jarvis library: {name}"
    )


def _download_archive(archive_url: str, target_path: str, timeout_sec: float) -> None:
    parsed = urlparse(archive_url)
    if parsed.scheme not in {"https", "http", "file"}:
        raise OfficialProjectFetchError(
            f"Unsupported archive URL scheme: {parsed.scheme or '<none>'}"
        )

    try:
        with urlopen(archive_url, timeout=timeout_sec) as response, open(
            target_path, "wb"
        ) as f1:
            shutil.copyfileobj(response, f1)
    except Exception as exc:
        raise OfficialProjectFetchError(
            f"Cannot download project archive: {archive_url} ({exc})"
        ) from exc


def _is_within_directory(base: str, target: str) -> bool:
    base_real = os.path.realpath(base)
    target_real = os.path.realpath(target)
    if base_real == target_real:
        return True
    return target_real.startswith(base_real + os.sep)


def _safe_extract_tar(archive_path: str, destination: str) -> None:
    with tarfile.open(archive_path, "r:*") as tf:
        for member in tf.getmembers():
            member_path = os.path.join(destination, member.name)
            if not _is_within_directory(destination, member_path):
                raise OfficialProjectFetchError(
                    f"Unsafe archive path encountered: {member.name}"
                )
            if member.islnk() or member.issym():
                raise OfficialProjectFetchError(
                    f"Symbolic links are not allowed in official project archives: {member.name}"
                )
        try:
            tf.extractall(destination, filter="data")
        except TypeError:
            tf.extractall(destination)


def _safe_extract_zip(archive_path: str, destination: str) -> None:
    with zipfile.ZipFile(archive_path, "r") as zf:
        for name in zf.namelist():
            member_path = os.path.join(destination, name)
            if not _is_within_directory(destination, member_path):
                raise OfficialProjectFetchError(
                    f"Unsafe archive path encountered: {name}"
                )
        zf.extractall(destination)


def _extract_archive(archive_path: str, destination: str) -> None:
    if tarfile.is_tarfile(archive_path):
        _safe_extract_tar(archive_path, destination)
        return
    if zipfile.is_zipfile(archive_path):
        _safe_extract_zip(archive_path, destination)
        return
    raise OfficialProjectFetchError("Unsupported archive format. Expected tar or zip archive.")


def _resolve_archive_project_root(extract_root: str, project: dict) -> str:
    archive_root = str(project.get("archive_root") or "").strip()
    if archive_root:
        candidate = os.path.join(extract_root, archive_root)
        if os.path.isdir(candidate):
            return candidate
        raise OfficialProjectFetchError(
            f"Configured archive_root not found for project {project['name']}: {archive_root}"
        )

    direct_name_root = os.path.join(extract_root, project["name"])
    if os.path.isdir(direct_name_root):
        return direct_name_root

    top_entries = sorted(
        entry
        for entry in os.listdir(extract_root)
        if entry not in {".", "..", "__MACOSX"}
    )
    top_dirs = [entry for entry in top_entries if os.path.isdir(os.path.join(extract_root, entry))]
    if len(top_dirs) == 1:
        return os.path.join(extract_root, top_dirs[0])

    if os.path.exists(os.path.join(extract_root, ".jarvis-project.json")):
        return extract_root

    raise OfficialProjectFetchError(
        "Cannot locate project root in archive. The official Jarvis library metadata must set archive_root."
    )


def fetch_official_project(
    project_name: str,
    target_dir: str | None = None,
    *,
    index_url: str | None = None,
    timeout_sec: float | None = None,
) -> OfficialProjectFetchReport:
    project = get_official_project(project_name, index_url=index_url)
    archive_url = str(project.get("archive_url") or "").strip()
    if not archive_url:
        raise OfficialProjectFetchError(
            f"No archive_url configured for official project: {project['name']}"
        )

    resolved_target = (
        os.path.abspath(os.path.join(os.getcwd(), project["name"]))
        if target_dir is None
        else os.path.abspath(os.path.expanduser(str(target_dir)))
    )
    if os.path.exists(resolved_target):
        raise OfficialProjectFetchError(
            f"Target directory already exists: {resolved_target}"
        )

    target_parent = os.path.dirname(resolved_target)
    if target_parent:
        os.makedirs(target_parent, exist_ok=True)

    timeout = _coerce_timeout(timeout_sec)

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            archive_path = os.path.join(tmpdir, "official_project.archive")
            extract_root = os.path.join(tmpdir, "extract")
            staging_root = os.path.join(tmpdir, "staged")
            os.makedirs(extract_root, exist_ok=True)

            _download_archive(archive_url, archive_path, timeout_sec=timeout)
            _extract_archive(archive_path, extract_root)
            source_root = _resolve_archive_project_root(extract_root, project)
            shutil.copytree(source_root, staging_root)
            shutil.move(staging_root, resolved_target)
    except OfficialProjectFetchError:
        if os.path.isdir(resolved_target):
            shutil.rmtree(resolved_target, ignore_errors=True)
        raise
    except Exception as exc:
        if os.path.isdir(resolved_target):
            shutil.rmtree(resolved_target, ignore_errors=True)
        raise OfficialProjectFetchError(
            f"Failed to fetch official project '{project['name']}': {exc}"
        ) from exc

    return OfficialProjectFetchReport(
        project_name=project["name"],
        target_dir=resolved_target,
        entrypoint=project.get("entrypoint") or "",
    )
