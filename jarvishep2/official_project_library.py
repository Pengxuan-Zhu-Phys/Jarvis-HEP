#!/usr/bin/env python3
"""Official project catalog client (D12.5 / D12.6).

True source of the catalog is a **GitHub JSON** in Jarvis-Examples (no PyPI package)::

    https://raw.githubusercontent.com/Pengxuan-Zhu-Phys/Jarvis-Examples/main/catalog/official_project_library.json

Override with ``JARVIS_OFFICIAL_LIBRARY_INDEX_URL`` (https or file://). Offline fallback:
packaged snapshot under ``jarvishep2/card/``, then ``~/.jarvis/cache/official_catalog.json``.
"""

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

from jarvishep2.project_crypto import (
    ENCRYPTION_SCHEME_NONE,
    ENCRYPTION_SCHEME_OPENSSL,
    PROJECT_FETCH_KEY_ENV,
    ProjectCryptoError,
    maybe_decrypt_archive,
)

DEFAULT_OFFICIAL_LIBRARY_NAME = "official Jarvis library"
# D12.6: catalog lives in Jarvis-Examples (GitHub JSON only — no PyPI package).
DEFAULT_OFFICIAL_LIBRARY_INDEX_URL = (
    "https://raw.githubusercontent.com/Pengxuan-Zhu-Phys/Jarvis-Examples/main/"
    "catalog/official_project_library.json"
)
OFFICIAL_LIBRARY_INDEX_ENV = "JARVIS_OFFICIAL_LIBRARY_INDEX_URL"
OFFICIAL_LIBRARY_TIMEOUT_ENV = "JARVIS_OFFICIAL_LIBRARY_TIMEOUT_SEC"
DEFAULT_TIMEOUT_SEC = 20.0
OFFICIAL_LIBRARY_CARD_PATH = ("card", "official_project_library.json")
# Last successful remote pull (offline soft cache).
USER_CATALOG_CACHE_PATH = os.path.join(
    os.path.expanduser("~"), ".jarvis", "cache", "official_catalog.json"
)
SUPPORTED_SCHEMA_MAJOR = 1


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
    access: str = "public"
    required_key: bool = False


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
        card = importlib_resources.files("jarvishep2")
        for part in OFFICIAL_LIBRARY_CARD_PATH:
            card = card.joinpath(part)
        with card.open("r", encoding="utf-8") as f1:
            payload = json.load(f1)
    except Exception as exc:
        raise OfficialCatalogError(
            f"Cannot load packaged official Jarvis library metadata: {exc}"
        ) from exc
    return payload


def _read_user_cache_catalog() -> dict | None:
    path = USER_CATALOG_CACHE_PATH
    if not os.path.isfile(path):
        return None
    try:
        with open(path, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        if isinstance(payload, dict):
            return payload
    except Exception:
        return None
    return None


def _write_user_cache_catalog(payload: dict) -> None:
    try:
        parent = os.path.dirname(USER_CATALOG_CACHE_PATH)
        os.makedirs(parent, exist_ok=True)
        tmp = USER_CATALOG_CACHE_PATH + ".tmp"
        with open(tmp, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(tmp, USER_CATALOG_CACHE_PATH)
    except Exception:
        pass


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


def _normalize_access(entry: dict) -> tuple[str, bool, str, str]:
    """Return (access, requires_key, encryption_scheme, encryption_hint)."""
    encryption = entry.get("encryption") if isinstance(entry.get("encryption"), dict) else {}
    scheme = str(
        encryption.get("scheme")
        or entry.get("encryption_scheme")
        or ENCRYPTION_SCHEME_NONE
    ).strip().lower()
    if scheme in {"", "null"}:
        scheme = ENCRYPTION_SCHEME_NONE

    raw_access = str(entry.get("access") or "").strip().lower()
    if raw_access in {"public", "open"}:
        access = "public"
    elif raw_access in {"restricted", "private", "encrypted", "closed"}:
        access = "restricted"
    elif scheme not in {ENCRYPTION_SCHEME_NONE, "none"}:
        access = "restricted"
    else:
        access = "public"

    if "requires_key" in entry:
        requires_key = bool(entry.get("requires_key"))
    else:
        requires_key = access == "restricted" or scheme not in {
            ENCRYPTION_SCHEME_NONE,
            "none",
        }

    if requires_key and scheme in {ENCRYPTION_SCHEME_NONE, "none"}:
        scheme = ENCRYPTION_SCHEME_OPENSSL

    if access == "public":
        requires_key = False
        if scheme not in {ENCRYPTION_SCHEME_NONE, "none"}:
            # Explicit encryption on a "public" row still needs a key.
            requires_key = True
            access = "restricted"

    hint = str(encryption.get("hint") or entry.get("key_hint") or "").strip()
    return access, requires_key, scheme, hint


def _normalize_project_entry(entry: object) -> dict | None:
    if not isinstance(entry, dict):
        return None

    name = str(entry.get("name") or "").strip()
    if not name:
        return None

    access, requires_key, scheme, hint = _normalize_access(entry)

    return {
        "name": name,
        "category": str(entry.get("category") or "").strip(),
        "summary": str(entry.get("summary") or "").strip(),
        "entrypoint": str(entry.get("entrypoint") or "").strip(),
        "archive_url": str(entry.get("archive_url") or "").strip(),
        "archive_root": str(entry.get("archive_root") or "").strip(),
        "compatibility_notes": str(entry.get("compatibility_notes") or "").strip(),
        "min_hep_version": str(entry.get("min_hep_version") or "").strip(),
        "access": access,
        "requires_key": requires_key,
        "encryption_scheme": scheme,
        "encryption_hint": hint,
    }


def _normalize_catalog(payload: object) -> dict:
    if not isinstance(payload, dict):
        raise OfficialCatalogError("Official Jarvis library metadata must be a JSON object.")

    schema_version = payload.get("schema_version")
    if schema_version is not None:
        try:
            major = int(str(schema_version).split(".", 1)[0])
        except (TypeError, ValueError) as exc:
            raise OfficialCatalogError(
                f"Invalid catalog schema_version: {schema_version!r}"
            ) from exc
        if major > SUPPORTED_SCHEMA_MAJOR:
            raise OfficialCatalogError(
                f"Official library catalog schema_version {schema_version} is newer "
                "than this Jarvis supports; please upgrade Jarvis-HEP."
            )

    projects = payload.get("projects")
    if not isinstance(projects, list):
        raise OfficialCatalogError(
            "Official Jarvis library metadata is missing the 'projects' list."
        )

    normalized = []
    for item in projects:
        project = _normalize_project_entry(item)
        if project is None:
            continue
        normalized.append(project)

    normalized.sort(key=lambda entry: entry["name"].lower())

    return {
        "library_name": str(payload.get("library_name") or DEFAULT_OFFICIAL_LIBRARY_NAME),
        "schema_version": schema_version if schema_version is not None else 1,
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
        payload = _read_catalog_from_url(effective_url, timeout_sec=timeout_sec)
        if isinstance(payload, dict):
            _write_user_cache_catalog(payload)
        return payload
    except OfficialCatalogError:
        cached = _read_user_cache_catalog()
        if cached is not None:
            return cached
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


def format_project_list_table(projects: list[dict]) -> str:
    """Plain-text project table for non-Rich consumers (shows key requirements)."""
    headers = ("Name", "Access", "Key", "Category", "Summary")
    rows: list[tuple[str, str, str, str, str]] = []
    for project in projects:
        access = str(project.get("access") or "public")
        needs_key = bool(project.get("requires_key"))
        key_label = "required" if needs_key else "no"
        rows.append(
            (
                str(project.get("name") or ""),
                access,
                key_label,
                str(project.get("category") or "-"),
                str(project.get("summary") or ""),
            )
        )
    if not rows:
        return "No verified projects are listed in the official library."

    widths = [len(h) for h in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], min(len(cell), 48 if idx == 4 else 24))

    def fmt_row(cells: tuple[str, ...] | list[str]) -> str:
        parts = []
        for idx, cell in enumerate(cells):
            text = cell if len(cell) <= widths[idx] else cell[: widths[idx] - 1] + "…"
            parts.append(f"{text:<{widths[idx]}}")
        return "  ".join(parts)

    lines = [
        fmt_row(headers),
        "  ".join("-" * w for w in widths),
        *[fmt_row(row) for row in rows],
        "",
        "Key: required → fetch needs --key or "
        f"{PROJECT_FETCH_KEY_ENV}; no → open download.",
    ]
    return "\n".join(lines)


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
    raise OfficialProjectFetchError(
        "Unsupported archive format. Expected tar/zip (or encrypted openssl blob)."
    )


def _resolve_archive_project_root(extract_root: str, project: dict) -> str:
    archive_root = str(project.get("archive_root") or "").strip()
    # "." means "auto-detect" (not literal extract_root), matching public catalog entries.
    if archive_root and archive_root not in {".", "./"}:
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
    top_dirs = [
        entry for entry in top_entries if os.path.isdir(os.path.join(extract_root, entry))
    ]
    if len(top_dirs) == 1:
        return os.path.join(extract_root, top_dirs[0])

    if os.path.exists(os.path.join(extract_root, ".jarvis-project.json")):
        return extract_root

    raise OfficialProjectFetchError(
        "Cannot locate project root in archive. "
        "The official Jarvis library metadata must set archive_root."
    )


def fetch_official_project(
    project_name: str,
    target_dir: str | None = None,
    *,
    index_url: str | None = None,
    timeout_sec: float | None = None,
    key: str | None = None,
) -> OfficialProjectFetchReport:
    project = get_official_project(project_name, index_url=index_url)
    archive_url = str(project.get("archive_url") or "").strip()
    if not archive_url:
        raise OfficialProjectFetchError(
            f"No archive_url configured for official project: {project['name']}"
        )

    requires_key = bool(project.get("requires_key"))
    if requires_key:
        from jarvishep2.project_crypto import resolve_fetch_key

        if not resolve_fetch_key(key):
            hint = project.get("encryption_hint") or ""
            msg = (
                f"Project '{project['name']}' is restricted. Fetch with:\n"
                f"  Jarvis project fetch {project['name']} --key YOUR_KEY\n"
                f"or set {PROJECT_FETCH_KEY_ENV}."
            )
            if hint:
                msg += f"\nKey hint: {hint}"
            raise OfficialProjectFetchError(msg)

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
    plain_tmp: str | None = None

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            archive_path = os.path.join(tmpdir, "official_project.archive")
            extract_root = os.path.join(tmpdir, "extract")
            staging_root = os.path.join(tmpdir, "staged")
            os.makedirs(extract_root, exist_ok=True)

            _download_archive(archive_url, archive_path, timeout_sec=timeout)
            try:
                plain_path = maybe_decrypt_archive(
                    archive_path,
                    requires_key=requires_key,
                    encryption_scheme=str(project.get("encryption_scheme") or "none"),
                    key=key,
                )
            except ProjectCryptoError as exc:
                raise OfficialProjectFetchError(str(exc)) from exc

            if plain_path != archive_path:
                plain_tmp = plain_path
            _extract_archive(plain_path, extract_root)
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
    finally:
        if plain_tmp and os.path.isfile(plain_tmp):
            try:
                os.unlink(plain_tmp)
            except OSError:
                pass

    return OfficialProjectFetchReport(
        project_name=project["name"],
        target_dir=resolved_target,
        entrypoint=project.get("entrypoint") or "",
        access=str(project.get("access") or "public"),
        required_key=requires_key,
    )


__all__ = [
    "DEFAULT_OFFICIAL_LIBRARY_INDEX_URL",
    "OFFICIAL_LIBRARY_INDEX_ENV",
    "PROJECT_FETCH_KEY_ENV",
    "OfficialCatalogError",
    "OfficialLibraryError",
    "OfficialProjectFetchError",
    "OfficialProjectFetchReport",
    "OfficialProjectNotFoundError",
    "fetch_official_project",
    "format_project_list_table",
    "get_official_project",
    "list_official_projects",
    "load_official_project_catalog",
]
