#!/usr/bin/env python3
"""LibDeps installation, reuse, and calculator-isolation helpers."""

from __future__ import annotations

import hashlib
import json
import os
import shlex
import subprocess
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any

from jarvishep2.base import decode_path
from jarvishep2.command_parser import CommandParser, ResolvedExecutable
from jarvishep2.workflow import resolve_module_layers


LIBRARY_INSTALL_STAMP_BASENAME = ".jarvis_install_stamp.json"
LIBRARY_INSTALL_CONTROL_BASENAME = "jarvis_install.json"
# Deliberately identical to the calculator install records (D13.11): operators
# can inspect one stable file format for both calculator packs and LibDeps.
LIBRARY_INSTALL_STAMP_SCHEMA = "jarvishep2.calc_install/v1"
LIBRARY_INSTALL_CONTROL_SCHEMA = "jarvishep2.calc_install/v2"


class LibraryInstallError(RuntimeError):
    """A control-process library installation failed before Workers started."""


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _source_stat_payload(source: str) -> dict[str, Any] | None:
    text = str(source or "").strip()
    if not text:
        return None
    path = os.path.abspath(text)
    if not os.path.exists(path):
        return {"path": path, "exists": False}
    try:
        stat = os.stat(path)
    except OSError:
        return {"path": path, "exists": False}
    return {
        "path": path,
        "exists": True,
        "mtime_ns": int(getattr(stat, "st_mtime_ns", int(stat.st_mtime * 1e9))),
        "size": int(stat.st_size),
        "is_dir": os.path.isdir(path),
    }


def _read_json(path: str) -> dict[str, Any] | None:
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, json.JSONDecodeError):
        return None
    return payload if isinstance(payload, dict) else None


def _write_json(path: str, payload: Mapping[str, Any]) -> str:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    temporary = path + ".tmp"
    with open(temporary, "w", encoding="utf-8") as handle:
        json.dump(dict(payload), handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(temporary, path)
    return path


def library_install_control_path(libdeps_path: str) -> str:
    return os.path.join(os.path.abspath(str(libdeps_path)), LIBRARY_INSTALL_CONTROL_BASENAME)


def library_install_stamp_path(installation_path: str) -> str:
    return os.path.join(os.path.abspath(str(installation_path)), LIBRARY_INSTALL_STAMP_BASENAME)


@dataclass(frozen=True)
class LibraryModuleSpec:
    """Normalized V1-shaped LibDeps module for the control process."""

    name: str
    required_modules: tuple[str, ...]
    installation_path: str
    source: str
    commands: tuple[dict[str, str], ...]

    @classmethod
    def from_config(
        cls,
        config: Mapping[str, Any],
        *,
        libdeps_path: str,
        parser: CommandParser,
    ) -> "LibraryModuleSpec":
        name = str(config.get("name") or "").strip()
        if not name:
            raise LibraryInstallError("LibDeps.Modules entry requires a non-empty name")
        installation = config.get("installation")
        installation = installation if isinstance(installation, Mapping) else {}
        raw_path = (
            installation.get("path")
            or config.get("path")
            or os.path.join(libdeps_path, name)
        )
        raw_source = installation.get("source") or config.get("source") or ""
        installation_path = decode_path(
            parser.resolve_static(str(raw_path)),
            project_root=parser.project_root,
            base_dir=libdeps_path,
        )
        source = ""
        if raw_source:
            source = decode_path(
                parser.resolve_static(str(raw_source)),
                project_root=parser.project_root,
                base_dir=libdeps_path,
            )
        raw_commands = installation.get("commands")
        if raw_commands is None:
            command_items: list[Any] = []
        elif isinstance(raw_commands, Sequence) and not isinstance(raw_commands, (str, bytes)):
            command_items = list(raw_commands)
        else:
            command_items = [raw_commands]
        commands = _normalize_library_commands(
            command_items,
            initial_cwd=libdeps_path,
            installation_path=installation_path,
            source=source,
            parser=parser,
        )
        dependencies = config.get("required_modules") or []
        if isinstance(dependencies, str):
            dependencies = [dependencies]
        if not isinstance(dependencies, Sequence):
            dependencies = []
        return cls(
            name=name,
            required_modules=tuple(str(item).strip() for item in dependencies if str(item).strip()),
            installation_path=os.path.abspath(installation_path),
            source=os.path.abspath(source) if source else "",
            commands=tuple(commands),
        )

    def fingerprint(self) -> str:
        payload = {
            "schema": LIBRARY_INSTALL_STAMP_SCHEMA,
            "module": self.name,
            "installation_path": self.installation_path,
            "source": self.source,
            "source_stat": _source_stat_payload(self.source),
            "commands": list(self.commands),
        }
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def _expand_library_placeholders(text: str, *, installation_path: str, source: str) -> str:
    return str(text).replace("${path}", installation_path).replace("${source}", source)


def _next_cwd_from_command(command: str, current_cwd: str, *, project_root: str) -> str:
    """Mirror V1's sequential standalone ``cd`` handling for string commands."""
    try:
        parts = shlex.split(str(command).strip(), posix=True)
    except ValueError:
        return current_cwd
    if len(parts) != 2 or parts[0] != "cd":
        return current_cwd
    return decode_path(parts[1], project_root=project_root, base_dir=current_cwd)


def _normalize_library_commands(
    command_items: Sequence[Any],
    *,
    initial_cwd: str,
    installation_path: str,
    source: str,
    parser: CommandParser,
) -> list[dict[str, str]]:
    """Turn V1 strings into explicit command/cwd records without losing ``cd``."""
    current_cwd = os.path.abspath(initial_cwd)
    normalized: list[dict[str, str]] = []
    for item in command_items:
        if isinstance(item, Mapping):
            command = item.get("cmd", item.get("command", ""))
            raw_cwd = item.get("cwd", current_cwd)
        else:
            command = item
            raw_cwd = current_cwd
        command_text = _expand_library_placeholders(
            str(command or ""), installation_path=installation_path, source=source
        )
        command_text = parser.resolve_static(command_text)
        cwd_text = _expand_library_placeholders(
            str(raw_cwd or current_cwd), installation_path=installation_path, source=source
        )
        cwd = decode_path(
            parser.resolve_static(cwd_text),
            project_root=parser.project_root,
            base_dir=current_cwd,
        )
        normalized.append({"cmd": command_text, "cwd": os.path.abspath(cwd)})
        if not isinstance(item, Mapping):
            current_cwd = _next_cwd_from_command(
                command_text, os.path.abspath(cwd), project_root=parser.project_root
            )
    return normalized


class LibraryInstaller:
    """Install one shared LibDeps graph before Redis and Workers exist."""

    def __init__(
        self,
        config: Mapping[str, Any] | None,
        *,
        parser: CommandParser,
        logs_dir: str,
        logger: Any | None = None,
        skip_installation: bool = False,
    ) -> None:
        self.config = dict(config or {})
        self.parser = parser
        self.logs_dir = os.path.abspath(str(logs_dir))
        self.logger = logger
        self.skip_installation = bool(skip_installation)
        libdeps = self.config.get("LibDeps")
        self.libdeps = dict(libdeps) if isinstance(libdeps, Mapping) else {}
        raw_base = self.libdeps.get("path") or parser.project_root
        self.libdeps_path = decode_path(
            parser.resolve_static(str(raw_base)), project_root=parser.project_root
        )
        raw_modules = self.libdeps.get("Modules") or []
        self.modules = [
            LibraryModuleSpec.from_config(item, libdeps_path=self.libdeps_path, parser=parser)
            for item in raw_modules
            if isinstance(item, Mapping)
        ]

    @property
    def control_path(self) -> str:
        return library_install_control_path(self.libdeps_path)

    def _log(self, message: str, *args: Any) -> None:
        if self.logger is not None:
            self.logger.info(message, *args)

    def _warn(self, message: str, *args: Any) -> None:
        if self.logger is not None:
            self.logger.warning(message, *args)

    def _control_state(self) -> tuple[dict[str, Any], int, bool]:
        state = dict(_read_json(self.control_path) or {})
        try:
            epoch = max(0, int(state.get("reinstall_epoch", 0) or 0))
        except (TypeError, ValueError):
            epoch = 0
        requested = bool(state.get("reinstall"))
        if requested:
            epoch += 1
        return state, epoch, requested

    def _write_control(self, *, epoch: int, modules: Mapping[str, Any]) -> None:
        _write_json(
            self.control_path,
            {
                "schema": LIBRARY_INSTALL_CONTROL_SCHEMA,
                "reinstall": False,
                "reinstall_epoch": epoch,
                "modules": dict(modules),
                "_help": "Set \"reinstall\": true and run again to rebuild all libraries.",
            },
        )

    def _layers(self) -> dict[int, list[LibraryModuleSpec]]:
        names = [spec.name for spec in self.modules]
        if len(names) != len(set(names)):
            raise LibraryInstallError("LibDeps.Modules contains duplicate module names")
        known = set(names)
        for spec in self.modules:
            unknown = sorted(set(spec.required_modules) - known)
            if unknown:
                raise LibraryInstallError(
                    f"Library '{spec.name}' requires unknown module(s): {', '.join(unknown)}"
                )
        layer_map = resolve_module_layers(
            [{"name": spec.name, "required_modules": list(spec.required_modules)} for spec in self.modules]
        )
        for spec in self.modules:
            if any(layer_map[dependency] >= layer_map[spec.name] for dependency in spec.required_modules):
                raise LibraryInstallError(
                    f"LibDeps.Modules has a dependency cycle involving '{spec.name}'"
                )
        grouped: dict[int, list[LibraryModuleSpec]] = {}
        for spec in self.modules:
            grouped.setdefault(layer_map[spec.name], []).append(spec)
        return grouped

    def _can_reuse(self, spec: LibraryModuleSpec, *, epoch: int, force: bool) -> bool:
        if force or not os.path.isdir(spec.installation_path):
            return False
        stamp = _read_json(library_install_stamp_path(spec.installation_path))
        if not stamp or str(stamp.get("fingerprint") or "") != spec.fingerprint():
            return False
        try:
            stamp_epoch = max(0, int(stamp.get("epoch", 0) or 0))
        except (TypeError, ValueError):
            return False
        return stamp_epoch >= epoch

    def _run_module(self, spec: LibraryModuleSpec, *, epoch: int, force: bool) -> dict[str, Any]:
        log_path = os.path.join(self.logs_dir, f"library-{spec.name}.log")
        if self.skip_installation:
            if not os.path.isdir(spec.installation_path):
                raise LibraryInstallError(
                    f"{spec.name}: --skip-library-installation requires an existing path "
                    f"{spec.installation_path}"
                )
            self._warn(
                "LibDeps %s: skipping installation; verified %s",
                spec.name,
                spec.installation_path,
            )
            return {"status": "skipped", "path": spec.installation_path, "log_path": log_path}
        if self._can_reuse(spec, epoch=epoch, force=force):
            self._log("LibDeps %s: reusing installation from %s", spec.name, spec.installation_path)
            return {"status": "reused", "path": spec.installation_path, "log_path": log_path}

        os.makedirs(self.logs_dir, exist_ok=True)
        # V1 cards often create ``installation.path`` themselves (for example
        # ``mkdir Pythia8``).  Creating it here would turn those valid commands
        # into an immediate “File exists” failure, so only ensure its parent.
        os.makedirs(os.path.dirname(spec.installation_path) or ".", exist_ok=True)
        self._log("LibDeps %s: building (log: %s)", spec.name, log_path)
        try:
            with open(log_path, "w", encoding="utf-8") as log_handle:
                for command in spec.commands:
                    cmd = command["cmd"]
                    cwd = command["cwd"]
                    log_handle.write(f"$ {cmd}\n[cwd] {cwd}\n")
                    log_handle.flush()
                    result = subprocess.run(
                        cmd,
                        shell=True,
                        cwd=cwd,
                        stdout=log_handle,
                        stderr=subprocess.STDOUT,
                        text=True,
                        check=False,
                    )
                    if result.returncode != 0:
                        raise LibraryInstallError(
                            f"Library '{spec.name}' command failed (exit {result.returncode}); "
                            f"see {log_path}"
                        )
        except OSError as exc:
            raise LibraryInstallError(
                f"Library '{spec.name}' could not start; see {log_path}: {exc}"
            ) from exc
        stamp_path = library_install_stamp_path(spec.installation_path)
        _write_json(
            stamp_path,
            {
                "schema": LIBRARY_INSTALL_STAMP_SCHEMA,
                "module": spec.name,
                "fingerprint": spec.fingerprint(),
                "epoch": epoch,
                "installed_at_utc": _utc_now_iso(),
            },
        )
        self._log("LibDeps %s: build complete", spec.name)
        return {"status": "installed", "path": spec.installation_path, "log_path": log_path}

    def prepare(self) -> dict[str, dict[str, Any]]:
        """Install/reuse all libraries, returning the state recorded in control JSON."""
        if not self.modules:
            return {}
        os.makedirs(self.libdeps_path, exist_ok=True)
        state, epoch, reinstall_requested = self._control_state()
        del state
        max_workers = max(1, int(self.libdeps.get("make_paraller", 1) or 1))
        results: dict[str, dict[str, Any]] = {}
        for _layer, specs in sorted(self._layers().items()):
            with ThreadPoolExecutor(max_workers=min(max_workers, len(specs))) as executor:
                futures = {
                    executor.submit(
                        self._run_module, spec, epoch=epoch, force=reinstall_requested
                    ): spec
                    for spec in specs
                }
                for future in as_completed(futures):
                    spec = futures[future]
                    results[spec.name] = future.result()
        if self.skip_installation:
            return results
        self._write_control(epoch=epoch, modules=results)
        return results


class LibraryManager:
    """Minimal LibDeps helper: symlink safe tools into Sample dirs."""

    def __init__(self, config: Mapping[str, Any] | None = None, *, task_root: str | None = None) -> None:
        self.config = dict(config or {})
        self.task_root = str(task_root or os.getcwd())

    @staticmethod
    def requires_shadow(clone_shadow: bool) -> bool:
        return bool(clone_shadow)

    def link_into_sample(self, source_path: str, sample_dir: str, link_name: str) -> str:
        """Symlink a concurrency-safe tool into a Sample directory (zero-copy)."""
        source = os.path.abspath(str(source_path))
        if not os.path.exists(source):
            raise FileNotFoundError(f"LibDeps source does not exist: {source}")
        sample_root = os.path.abspath(str(sample_dir))
        os.makedirs(sample_root, exist_ok=True)
        link_path = os.path.join(sample_root, str(link_name))
        if os.path.lexists(link_path):
            return link_path
        os.symlink(source, link_path)
        return link_path

    def resolve_registered(self, parser: CommandParser) -> dict[str, ResolvedExecutable]:
        """Return the Phase-1 registered executable map from a CommandParser."""
        return dict(parser.registered)


__all__ = [
    "LIBRARY_INSTALL_CONTROL_BASENAME",
    "LIBRARY_INSTALL_STAMP_BASENAME",
    "LibraryInstallError",
    "LibraryInstaller",
    "LibraryManager",
    "LibraryModuleSpec",
    "library_install_control_path",
    "library_install_stamp_path",
]
