#!/usr/bin/env python3
"""Picklable calculator configuration (D9.3 CalculatorSpec + D12.1 V1 YAML parity)."""

from __future__ import annotations

import os
import shlex
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any


def _normalize_timeout(value: Any) -> float | None:
    if value is None:
        return None
    timeout = float(value)
    return timeout if timeout > 0 else None


def next_cwd_from_command(command: str, current_cwd: str) -> str:
    """V1 command-chain semantics: standalone ``cd <path>`` updates inherited cwd.

    Only exact two-token forms are treated as directory switches
    (``cd path``, ``cd "path with spaces"``). Combined shell lines such as
    ``cd foo && make`` do **not** change the inherited cwd (same as V1).

    Relative targets are joined against ``current_cwd``. Paths may still contain
    runtime tokens (``@PackID``); do not force ``abspath`` so those survive until
    Worker/CommandParser resolution.
    """
    line = str(command or "").strip()
    if not line:
        return current_cwd
    try:
        parts = shlex.split(line, posix=True)
    except ValueError:
        return current_cwd
    if not parts or parts[0] != "cd" or len(parts) != 2:
        return current_cwd
    target = parts[1]
    if not target:
        return current_cwd
    if os.path.isabs(target):
        return os.path.normpath(target)
    anchor = str(current_cwd or ".")
    return os.path.normpath(os.path.join(anchor, target))


def is_cwd_switch_command(command: str) -> bool:
    """True when ``command`` is a pure standalone ``cd <path>`` (no other argv)."""
    line = str(command or "").strip()
    if not line:
        return False
    try:
        parts = shlex.split(line, posix=True)
    except ValueError:
        return False
    return bool(parts) and parts[0] == "cd" and len(parts) == 2


def _normalize_command_entry(item: Any, *, default_cwd: str) -> dict[str, Any] | None:
    """Normalize a V1 string or V2 mapping command into ``{cmd, cwd, ...}``.

    V1 cards use plain string lists under ``installation`` / ``initialization`` /
    ``execution.commands``. V2 historically dropped non-mapping entries silently.
    """
    if isinstance(item, str):
        cmd = item.strip()
        if not cmd:
            return None
        return {"cmd": cmd, "cwd": default_cwd}
    if not isinstance(item, Mapping):
        return None
    entry = dict(item)
    if "cmd" not in entry and "command" in entry:
        entry["cmd"] = entry.pop("command")
    if "cmd" not in entry:
        # Keep unexpected mapping shapes only if they already look executable.
        return None
    entry["cmd"] = str(entry["cmd"])
    if "cwd" not in entry or entry.get("cwd") in (None, ""):
        entry["cwd"] = default_cwd
    else:
        entry["cwd"] = str(entry["cwd"])
    return entry


def normalize_command_list(
    items: Any,
    *,
    default_cwd: str = ".",
) -> tuple[dict[str, Any], ...]:
    """Normalize a command list with V1 ``cd`` → inherited-cwd chaining.

    Empty/invalid entries are skipped. Each entry receives the running cwd;
    a pure ``cd <path>`` updates that cwd for subsequent entries (V1
    ``ConfigLoader._next_cwd_from_command`` contract).
    """
    if items is None:
        return ()
    if isinstance(items, (str, bytes)) or not isinstance(items, Sequence):
        # A single bare string is a one-element list.
        if isinstance(items, str):
            items = [items]
        else:
            return ()
    current_cwd = str(default_cwd or ".")
    normalized: list[dict[str, Any]] = []
    for item in items:
        entry = _normalize_command_entry(item, default_cwd=current_cwd)
        if entry is None:
            continue
        # Explicit mapping ``cwd`` wins for this entry only; otherwise inherit.
        if isinstance(item, Mapping) and item.get("cwd") not in (None, ""):
            entry["cwd"] = str(item["cwd"])
        else:
            entry["cwd"] = current_cwd
        current_cwd = next_cwd_from_command(str(entry["cmd"]), str(entry["cwd"]))
        normalized.append(entry)
    return tuple(normalized)


def _mapping_list(items: Any) -> tuple[dict[str, Any], ...]:
    if not isinstance(items, Sequence) or isinstance(items, (str, bytes)):
        return ()
    return tuple(dict(item) for item in items if isinstance(item, Mapping))


@dataclass(frozen=True)
class CalculatorSpec:
    """Immutable view of a Calculators.Modules[] entry after parse."""

    name: str
    clone_shadow: bool = False
    installation: tuple[Mapping[str, Any], ...] = ()
    initialization: tuple[Mapping[str, Any], ...] = ()
    commands: tuple[Mapping[str, Any], ...] = ()
    input_specs: tuple[Mapping[str, Any], ...] = ()
    output_specs: tuple[Mapping[str, Any], ...] = ()
    timeout: float | None = None
    basepath: str = "."
    source: str = ""
    symlink_name: str = ""
    env_setup: tuple[Mapping[str, Any], ...] = ()
    selection: str | None = None
    raw: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_config(cls, name: str, config: Mapping[str, Any]) -> CalculatorSpec:
        cfg = dict(config or {})
        execution = dict(cfg.get("execution") or {})
        module_name = str(name or cfg.get("name", "")).strip() or "Calculator"
        symlink = str(cfg.get("symlink_name", module_name)).strip() or module_name
        basepath = str(cfg.get("path", execution.get("path", ".")) or ".")
        execution_path = str(execution.get("path") or basepath)
        # V1 tolerates unused keys: modes, make_paraller (module-level), required_modules, …
        # They remain in ``raw``; pools consume make_paraller via worker_config.
        selection_raw = cfg.get("selection")
        selection: str | None
        if selection_raw is None or selection_raw is False:
            selection = None
        else:
            text = str(selection_raw).strip()
            selection = text or None
        return cls(
            name=module_name,
            clone_shadow=bool(cfg.get("clone_shadow", False)),
            installation=normalize_command_list(
                cfg.get("installation"), default_cwd=basepath
            ),
            initialization=normalize_command_list(
                cfg.get("initialization"), default_cwd=basepath
            ),
            commands=normalize_command_list(
                execution.get("commands"), default_cwd=execution_path
            ),
            input_specs=_mapping_list(execution.get("input")),
            output_specs=_mapping_list(execution.get("output")),
            timeout=_normalize_timeout(cfg.get("timeout")),
            basepath=basepath,
            source=str(cfg.get("source", "")).strip(),
            symlink_name=symlink,
            env_setup=_mapping_list(cfg.get("env_setup")),
            selection=selection,
            raw=cfg,
        )


__all__ = [
    "CalculatorSpec",
    "is_cwd_switch_command",
    "next_cwd_from_command",
    "normalize_command_list",
]
