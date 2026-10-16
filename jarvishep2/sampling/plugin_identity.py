"""Metadata-only scan source binding for external samplers (S2a)."""

from __future__ import annotations

import hashlib
import json
import os
import platform
import tempfile
from collections.abc import Mapping
from importlib.metadata import version
from pathlib import Path
from typing import Any

from jarvis_hep_sampler.errors import PluginError, PluginResumeMismatch
from jarvis_hep_sampler.metadata import parse_metadata

from jarvishep2.sampler_catalog import PluginRecord, SOURCE_SCAN, builtin_plugin_path, builtin_specs
from jarvishep2.sampler_sdk import require_sampler_sdk

BINDING_VERSION = 1


def trusted_builtin_binding_method(
    method: str, original_source: object, source_sha256: object,
) -> str | None:
    """Recognize a shipped plug-in binding by its package-relative source path.

    The saved scan copy is checked against ``source_sha256`` in preflight. The
    installed package path and bytes may have changed since the scan was first
    bound, so neither can decide whether the old scan copy is usable.
    """
    package_source = builtin_plugin_path(method)
    if (
        package_source is None
        or not isinstance(original_source, str)
        or not isinstance(source_sha256, str)
    ):
        return None
    try:
        original_path = Path(original_source)
        package_path = package_source.resolve()
        package_suffix = ("jarvishep2", "sampling", "plugins", package_path.name)
        if (
            not original_path.is_absolute()
            or package_path.parts[-4:] != package_suffix
            or original_path.parts[-4:] != package_suffix
        ):
            return None
    except OSError:
        return None
    return method


def identity_path(config: Mapping[str, Any]) -> Path:
    root = Path(str(config.get("task_root") or os.getcwd())).resolve()
    scan = str(
        config.get("scan_name") or (config.get("Scan") or {}).get("name") or "scan"
    )
    result = Path(str(config.get("task_result_dir") or root / "outputs" / scan))
    database = result / "DATABASE"
    if config.get("_check_modules_sample_layout"):
        database /= "test"
    return database / "sampler_identity.json"


def canonical_hash(fields: Mapping[str, Any]) -> str:
    data = json.dumps(
        dict(fields),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(data).hexdigest()


def read_binding(config: Mapping[str, Any]) -> dict[str, Any] | None:
    path = identity_path(config)
    if not path.exists():
        return None
    method = str((config.get("Sampling") or {}).get("Method") or "")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        if data["binding_version"] != BINDING_VERSION:
            raise ValueError("unsupported binding version")
        fields = data["identity_fields"]
        if canonical_hash(fields) != data["hard_identity"]:
            raise ValueError("identity fields do not match their hash")
        digest = data["source_sha256"]
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(c not in "0123456789abcdef" for c in digest)
        ):
            raise ValueError("invalid source SHA-256")
        source_path(config, data)
        trusted = trusted_builtin_binding_method(method, data.get("original_source"), digest)
        meta = parse_metadata(data["metadata"], trusted_builtin_method=trusted)
        if meta.method != method or fields["method"] != method:
            raise ValueError("saved sampler method differs from Sampling.Method")
        if method in {s.name for s in builtin_specs()}:
            raise ValueError("a plug-in binding cannot replace a built-in sampler")
        if not isinstance(data["hooks"], list) or not isinstance(
            data["environment"], dict
        ):
            raise ValueError("invalid saved metadata")
        schema = data["options_schema"]
        if schema is not None and not isinstance(schema, dict):
            raise ValueError("invalid saved options schema")
        return data
    except (OSError, ValueError, KeyError, TypeError, PluginError) as exc:
        raise PluginResumeMismatch(
            f"invalid sampler identity record: {exc}; use a new scan name",
            method=method,
            source=str(path),
        ) from exc


def source_path(config: Mapping[str, Any], binding: Mapping[str, Any]) -> Path:
    root = Path(str(config.get("task_root") or os.getcwd())).resolve()
    rel = Path(binding["source_file"])
    path = (root / rel).resolve()
    if (
        rel.is_absolute()
        or ".." in rel.parts
        or not path.is_relative_to(root / "checkpoints")
    ):
        raise ValueError("source binding must stay inside this project's checkpoints")
    return path


def binding_record(config: Mapping[str, Any]) -> PluginRecord | None:
    """Read only metadata; never import source or depend on the current registry."""
    data = read_binding(config)
    if data is None:
        return None
    trusted = trusted_builtin_binding_method(
        str((config.get("Sampling") or {}).get("Method") or ""),
        data.get("original_source"), data.get("source_sha256"),
    )
    meta = parse_metadata(data["metadata"], trusted_builtin_method=trusted)
    return PluginRecord(
        method=meta.method,
        source=SOURCE_SCAN,
        version=meta.version,
        profile=meta.profile,
        description=meta.description,
        statistical=meta.statistical,
        options_schema=data["options_schema"],
        hooks=tuple(data["hooks"]),
        source_path=data["original_source"],
        snapshot_path=str(source_path(config, data)),
        source_sha256=data["source_sha256"],
    )


def environment_versions() -> dict[str, str]:
    return {
        "sdk": require_sampler_sdk(),
        "python": platform.python_version(),
        "numpy": version("numpy"),
    }


def write_binding(path: Path, data: Mapping[str, Any]) -> None:
    """Atomically publish only after successful initialization."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".sampler_identity.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(dict(data), handle, sort_keys=True, indent=2, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def check_identity(
    saved: Mapping[str, Any], current: Mapping[str, Any], *, method: str
) -> None:
    differences = [
        key
        for key in sorted(set(saved) | set(current))
        if saved.get(key) != current.get(key)
    ]
    if differences:
        raise PluginResumeMismatch(
            "science settings differ: "
            + ", ".join(differences)
            + "; use a new scan name",
            method=method,
        )


def validate_checkpoint(
    state: Mapping[str, Any],
    *,
    identity: str,
    source: str,
    fields: Mapping[str, Any],
    method: str,
) -> None:
    if not isinstance(state, Mapping):
        raise PluginResumeMismatch(
            "checkpoint sampler state must be a mapping", method=method
        )
    saved = state.get("plugin") or {}
    if not isinstance(saved, Mapping):
        raise PluginResumeMismatch(
            "checkpoint plug-in identity must be a mapping", method=method
        )
    if saved.get("identity") != identity or saved.get("source_sha256") != source:
        raise PluginResumeMismatch(
            "checkpoint identity/source differs from the scan binding", method=method
        )
    if (
        state.get("maxp") != fields["point_number"]
        or state.get("seed") != fields["seed"]
    ):
        raise PluginResumeMismatch(
            "checkpoint budget/seed differs from the scan binding", method=method
        )
    if state.get("selectionexp") != fields["selection"]:
        raise PluginResumeMismatch(
            "checkpoint selection differs from the scan binding", method=method
        )
    try:
        index, accepted = state["index"], state["accepted_index"]
        if (
            isinstance(index, bool)
            or isinstance(accepted, bool)
            or not isinstance(index, int)
            or not isinstance(accepted, int)
            or not 0 <= accepted <= index <= fields["point_number"]
        ):
            raise ValueError("invalid candidate/accepted cursor")
    except (ValueError, KeyError, TypeError) as exc:
        raise PluginResumeMismatch(
            f"invalid checkpoint cursor: {exc}", method=method
        ) from exc
