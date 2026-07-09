#!/usr/bin/env python3
"""Picklable calculator configuration (D9.3 CalculatorSpec)."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any


def _normalize_timeout(value: Any) -> float | None:
    if value is None:
        return None
    timeout = float(value)
    return timeout if timeout > 0 else None


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
    raw: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def from_config(cls, name: str, config: Mapping[str, Any]) -> CalculatorSpec:
        cfg = dict(config or {})
        execution = dict(cfg.get("execution") or {})
        module_name = str(name or cfg.get("name", "")).strip() or "Calculator"
        symlink = str(cfg.get("symlink_name", module_name)).strip() or module_name
        return cls(
            name=module_name,
            clone_shadow=bool(cfg.get("clone_shadow", False)),
            installation=tuple(dict(item) for item in (cfg.get("installation") or []) if isinstance(item, Mapping)),
            initialization=tuple(
                dict(item) for item in (cfg.get("initialization") or []) if isinstance(item, Mapping)
            ),
            commands=tuple(dict(item) for item in (execution.get("commands") or []) if isinstance(item, Mapping)),
            input_specs=tuple(dict(item) for item in (execution.get("input") or []) if isinstance(item, Mapping)),
            output_specs=tuple(dict(item) for item in (execution.get("output") or []) if isinstance(item, Mapping)),
            timeout=_normalize_timeout(cfg.get("timeout")),
            basepath=str(cfg.get("path", execution.get("path", "."))),
            source=str(cfg.get("source", "")).strip(),
            symlink_name=symlink,
            env_setup=tuple(dict(item) for item in (cfg.get("env_setup") or []) if isinstance(item, Mapping)),
            raw=cfg,
        )


__all__ = ["CalculatorSpec"]
