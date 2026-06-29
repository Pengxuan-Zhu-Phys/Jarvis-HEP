#!/usr/bin/env python3
"""Two-phase command/path resolution for Jarvis-HEP V2 (WP-D3.1)."""

from __future__ import annotations

import copy
import os
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

from jarvishep2.base import decode_path, expand_j, infer_project_root, scan_output_root
from jarvishep2.sample import ensure_sample_materialized

SAMPLE_TOKENS: tuple[str, ...] = ("@SampleID", "@Sdir", "@PackID")
_LIBDEPS_PATTERN = re.compile(r"\$\{LibDeps:([^}]+)\}")
_SCAN_PATTERN = re.compile(r"\$\{Scan:([^}]+)\}")
def _registered_name_pattern(name: str) -> re.Pattern[str]:
    return re.compile(rf"(?<![\w./-]){re.escape(name)}(?![\w./-])")


@dataclass(frozen=True)
class ResolvedExecutable:
    name: str
    path: str
    resolution: str


def _contains_sample_tokens(text: str) -> bool:
    return any(token in text for token in SAMPLE_TOKENS)


def _walk_strings(value: Any, callback) -> Any:
    if isinstance(value, str):
        return callback(value)
    if isinstance(value, Mapping):
        return {key: _walk_strings(item, callback) for key, item in value.items()}
    if isinstance(value, list):
        return [_walk_strings(item, callback) for item in value]
    if isinstance(value, tuple):
        return tuple(_walk_strings(item, callback) for item in value)
    return value


class CommandParser:
    """Phase 1 (control) + Phase 2 (Worker) command resolution."""

    def __init__(
        self,
        *,
        project_root: str,
        scan_name: str = "default",
        libdeps_paths: Mapping[str, str] | None = None,
        registered: Mapping[str, ResolvedExecutable] | None = None,
        registered_symlink_root: str | None = None,
    ) -> None:
        self.project_root = os.path.abspath(str(project_root))
        self.scan_name = str(scan_name or "default").strip() or "default"
        self.libdeps_paths = {
            str(name): os.path.abspath(str(path))
            for name, path in (libdeps_paths or {}).items()
            if str(name).strip()
        }
        self.registered: dict[str, ResolvedExecutable] = dict(registered or {})
        self.registered_symlink_root = os.path.abspath(
            str(registered_symlink_root or os.path.join(self.project_root, "deps", "registered"))
        )

    @classmethod
    def from_picklable(cls, payload: Mapping[str, Any] | None) -> CommandParser | None:
        if not isinstance(payload, Mapping) or not payload:
            return None
        registered: dict[str, ResolvedExecutable] = {}
        raw_registered = payload.get("registered") or {}
        if isinstance(raw_registered, Mapping):
            for name, item in raw_registered.items():
                if isinstance(item, Mapping):
                    registered[str(name)] = ResolvedExecutable(
                        name=str(item.get("name", name)),
                        path=str(item.get("path", "")),
                        resolution=str(item.get("resolution", "direct_path")),
                    )
        return cls(
            project_root=str(payload.get("project_root") or infer_project_root()),
            scan_name=str(payload.get("scan_name") or "default"),
            libdeps_paths=dict(payload.get("libdeps_paths") or {}),
            registered=registered,
            registered_symlink_root=str(
                payload.get("registered_symlink_root")
                or os.path.join(str(payload.get("project_root") or infer_project_root()), "deps", "registered")
            ),
        )

    @classmethod
    def from_config(
        cls,
        config: Mapping[str, Any] | None,
        *,
        project_root: str | None = None,
        task_root: str | None = None,
    ) -> CommandParser:
        cfg = dict(config or {})
        root = str(
            project_root
            or task_root
            or cfg.get("task_result_dir")
            or cfg.get("project_root")
            or infer_project_root()
        )
        scan_name = str((cfg.get("Scan") or {}).get("name") or cfg.get("scan_name") or "default")
        libdeps_paths = _build_libdeps_paths(cfg, project_root=root)
        parser = cls(project_root=root, scan_name=scan_name, libdeps_paths=libdeps_paths)
        for spec in _registered_specs(cfg):
            resolved = parser.register(spec)
            parser.registered[resolved.name] = resolved
        return parser

    def register(self, spec: Mapping[str, Any]) -> ResolvedExecutable:
        name = str(spec.get("name", "")).strip()
        if not name:
            raise ValueError("registered_executables entry requires a non-empty name")
        resolution = str(spec.get("resolution", "direct_path")).strip().lower() or "direct_path"
        if resolution not in {"direct_path", "symlink"}:
            raise ValueError(
                f"registered executable '{name}' has invalid resolution '{resolution}'"
            )
        source = self.resolve_static(str(spec.get("source", "")).strip())
        if not source:
            raise ValueError(f"registered executable '{name}' requires a source path")
        source_path = decode_path(source, project_root=self.project_root)
        if not os.path.exists(source_path):
            raise FileNotFoundError(
                f"registered executable '{name}' source does not exist: {source_path}"
            )
        if resolution == "symlink":
            os.makedirs(self.registered_symlink_root, exist_ok=True)
            link_path = os.path.join(self.registered_symlink_root, name)
            if not os.path.lexists(link_path):
                os.symlink(source_path, link_path)
            path = os.path.abspath(link_path)
        else:
            path = os.path.abspath(source_path)
        return ResolvedExecutable(name=name, path=path, resolution=resolution)

    def resolve_static(self, text: str) -> str:
        if text is None:
            return ""
        raw = str(text)
        if not raw:
            return raw
        resolved = expand_j(raw, project_root=self.project_root)
        resolved = _SCAN_PATTERN.sub(self._replace_scan_token, resolved)
        resolved = _LIBDEPS_PATTERN.sub(self._replace_libdeps_token, resolved)
        resolved = self._replace_registered_names(resolved)
        if "~" in resolved:
            resolved = os.path.expanduser(resolved)
        if resolved.startswith("&"):
            raise ValueError(f"Unresolved static token remains after Phase 1: {resolved}")
        if _LIBDEPS_PATTERN.search(resolved) or _SCAN_PATTERN.search(resolved) or "&J" in resolved:
            raise ValueError(f"Unresolved static token remains after Phase 1: {resolved}")
        return resolved

    def resolve_static_config(self, config: Mapping[str, Any]) -> dict[str, Any]:
        payload = copy.deepcopy(dict(config))

        def _phase1(text: str) -> str:
            if _contains_sample_tokens(text):
                return _walk_sample_aware_static(text)
            return self.resolve_static(text)

        def _walk_sample_aware_static(text: str) -> str:
            # Resolve static prefixes while preserving per-Sample tokens for Phase 2.
            parts: list[str] = []
            cursor = 0
            while cursor < len(text):
                next_index = len(text)
                next_token = ""
                for token in SAMPLE_TOKENS:
                    index = text.find(token, cursor)
                    if index >= 0 and index < next_index:
                        next_index = index
                        next_token = token
                if next_index > cursor:
                    parts.append(self.resolve_static(text[cursor:next_index]))
                if next_token:
                    parts.append(next_token)
                    cursor = next_index + len(next_token)
                else:
                    break
            return "".join(parts) if parts else self.resolve_static(text)

        return _walk_strings(payload, _phase1)

    def resolve_sample(
        self,
        text: str,
        *,
        sample_info: Mapping[str, Any],
        pack_id: str | None = None,
        stage: str = "execution",
        field: str = "text",
    ) -> str:
        if text is None:
            return ""
        raw = str(text)
        if stage == "install" or not _contains_sample_tokens(raw):
            if self.has_static_tokens(raw):
                raise ValueError(
                    f"Phase 2 cannot resolve static tokens in stage '{stage}' field '{field}': {raw}"
                )
            return raw
        if self.has_static_tokens(raw):
            raise ValueError(
                f"Phase 2 cannot resolve static tokens in stage '{stage}' field '{field}': {raw}"
            )
        if not isinstance(sample_info, Mapping):
            raise RuntimeError(
                f"Runtime token requires sample_info during stage '{stage}' for field '{field}'"
            )
        resolved = raw
        if "@PackID" in resolved:
            effective_pack = pack_id or sample_info.get("pack_id")
            if not effective_pack:
                raise RuntimeError(
                    f"@PackID requires pack_id during stage '{stage}' for field '{field}'"
                )
            resolved = resolved.replace("@PackID", str(effective_pack))
        if "@SampleID" in resolved:
            sample_uuid = sample_info.get("uuid")
            if not sample_uuid:
                raise RuntimeError(
                    f"@SampleID requires uuid during stage '{stage}' for field '{field}'"
                )
            resolved = resolved.replace("@SampleID", str(sample_uuid))
        if "@Sdir" in resolved:
            save_dir = ensure_sample_materialized(dict(sample_info))
            if save_dir is None:
                raise RuntimeError(
                    f"@Sdir requires save_dir during stage '{stage}' for field '{field}'"
                )
            resolved = resolved.replace("@Sdir", str(save_dir))
        return resolved

    def has_static_tokens(self, text: str) -> bool:
        if text is None:
            return False
        raw = str(text)
        if "&J" in raw:
            return True
        if _LIBDEPS_PATTERN.search(raw):
            return True
        if _SCAN_PATTERN.search(raw):
            return True
        for name, entry in self.registered.items():
            if entry.path in raw:
                continue
            if _registered_name_pattern(name).search(raw):
                return True
        return False

    def _replace_scan_token(self, match: re.Match[str]) -> str:
        token_name = str(match.group(1) or "").strip().lower()
        if token_name == "name":
            return self.scan_name
        if token_name in {"root", "output", "outputs"}:
            return scan_output_root(project_root=self.project_root, scan_name=self.scan_name)
        return scan_output_root(project_root=self.project_root, scan_name=self.scan_name)

    def _replace_libdeps_token(self, match: re.Match[str]) -> str:
        name = str(match.group(1) or "").strip()
        if name not in self.libdeps_paths:
            raise KeyError(f"Unknown LibDeps reference '${{LibDeps:{name}}}'")
        return self.libdeps_paths[name]

    def _replace_registered_names(self, text: str) -> str:
        resolved = text
        for name in sorted(self.registered, key=len, reverse=True):
            path = self.registered[name].path
            resolved = _registered_name_pattern(name).sub(path, resolved)
        return resolved


def _registered_specs(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    from jarvishep2.runtime_config import parse_registered_executables

    return parse_registered_executables(config)


def _build_libdeps_paths(config: Mapping[str, Any], *, project_root: str) -> dict[str, str]:
    libdeps = config.get("LibDeps") or {}
    if not isinstance(libdeps, Mapping):
        return {}
    paths: dict[str, str] = {}
    base = decode_path(str(libdeps.get("path", project_root)), project_root=project_root)
    modules = libdeps.get("Modules") or []
    if not isinstance(modules, Sequence):
        return paths
    for module in modules:
        if not isinstance(module, Mapping):
            continue
        name = str(module.get("name", "")).strip()
        if not name:
            continue
        install = module.get("installation")
        candidate = None
        if isinstance(install, Mapping):
            candidate = install.get("path") or install.get("source")
        candidate = candidate or module.get("path") or module.get("source")
        if candidate:
            paths[name] = decode_path(str(candidate), project_root=project_root, base_dir=base)
        else:
            paths[name] = os.path.join(base, name)
    return paths


def prepare_calculator_modules(
    modules: Sequence[Mapping[str, Any]] | None,
    parser: CommandParser,
) -> list[dict[str, Any]]:
    """Apply Phase 1 resolution to calculator module configs."""
    prepared: list[dict[str, Any]] = []
    for module in modules or []:
        if isinstance(module, Mapping):
            prepared.append(parser.resolve_static_config(module))
    return prepared


__all__ = [
    "CommandParser",
    "ResolvedExecutable",
    "SAMPLE_TOKENS",
    "prepare_calculator_modules",
]