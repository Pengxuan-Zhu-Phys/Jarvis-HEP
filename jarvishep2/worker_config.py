#!/usr/bin/env python3
"""Worker blueprint builder with Phase-1 command resolution (WP-D3.1)."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_cleanup_config,
    get_delete_method,
    get_staging_dir,
    handoff_to_staging_enabled,
    workflow_has_calculator,
    workflow_references_sdir,
)


def _default_mapper(cfg: Mapping[str, Any]) -> dict[str, Any]:
    mapper = cfg.get("Mapper")
    if isinstance(mapper, Mapping):
        return dict(mapper)
    variables = (cfg.get("Sampling") or {}).get("Variables") if isinstance(cfg.get("Sampling"), Mapping) else None
    if variables:
        return {"type": "distribution", "variables": list(variables)}
    return {"type": "identity", "keys": ["x", "y"]}


def _config_references_sdir(modules: list[dict[str, Any]]) -> bool:
    import json

    return "@Sdir" in json.dumps(modules, sort_keys=True)


def build_command_parser(config: Mapping[str, Any] | None) -> CommandParser:
    """Build a CommandParser from a task config mapping."""
    cfg = dict(config or {})
    project_root = str(cfg.get("project_root") or cfg.get("task_root") or "")
    return CommandParser.from_config(cfg, project_root=project_root or None)


def build_worker_config(
    config: Mapping[str, Any] | None,
    *,
    task_result_dir: str,
    sample_dirs: str | None = None,
    opera_modules: list[dict[str, Any]] | None = None,
    calculator_modules: list[dict[str, Any]] | None = None,
    likelihood_expressions: list[dict[str, Any]] | None = None,
    parser: CommandParser | None = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Normalize a picklable Worker config with Phase-1 static resolution applied."""
    cfg = dict(config or {})
    cfg["task_result_dir"] = task_result_dir
    command_parser = parser or build_command_parser(cfg)

    extra_payload = dict(extra or {})
    calc_modules = list(calculator_modules or extra_payload.pop("calculator_modules", []) or [])
    if not calc_modules:
        calculators = (cfg.get("Calculators") or {}).get("Modules") or []
        if isinstance(calculators, list):
            calc_modules = [dict(item) for item in calculators if isinstance(item, Mapping)]

    opera = list(opera_modules or extra_payload.pop("opera_modules", []) or [])
    if not opera:
        operas = (cfg.get("Operas") or {}).get("Modules") or []
        if isinstance(operas, list):
            opera = [dict(item) for item in operas if isinstance(item, Mapping)]

    resolved_calculators = prepare_calculator_modules(calc_modules, command_parser)
    sample_root = sample_dirs or __import__("os").path.join(task_result_dir, "SAMPLE")
    sample_config = dict(extra_payload.pop("sample_config", {}) or {})
    sample_config.setdefault("task_result_dir", task_result_dir)
    sample_config.setdefault("sample_dirs", sample_root)
    sample_config.setdefault(
        "sample_artifacts",
        str((cfg.get("Runtime") or {}).get("sample_artifacts", "auto")),
    )
    sample_config.setdefault(
        "workflow_has_calculator",
        workflow_has_calculator(cfg) or bool(resolved_calculators),
    )
    sample_config.setdefault(
        "workflow_references_sdir",
        workflow_references_sdir(cfg) or _config_references_sdir(resolved_calculators),
    )
    worker_config: dict[str, Any] = {
        "sample_config": sample_config,
        "mapper": extra_payload.pop("mapper", _default_mapper(cfg)),
        "opera_modules": opera,
        "calculator_modules": resolved_calculators,
        "likelihood_expressions": list(
            likelihood_expressions
            or (cfg.get("Likelihood") or {}).get("expressions")
            or (cfg.get("Sampling") or {}).get("LogLikelihood")
            or []
        ),
        "pull_timeout": 1,
        "delete_method": get_delete_method(cfg),
        "staging_dir": get_staging_dir(cfg, task_result_dir=task_result_dir),
        "handoff_to_staging": handoff_to_staging_enabled(cfg),
        "cleanup_config": get_cleanup_config(cfg),
        "archiver_config": get_archiver_config(cfg),
        "command_parser": {
            "project_root": command_parser.project_root,
            "scan_name": command_parser.scan_name,
            "libdeps_paths": dict(command_parser.libdeps_paths),
            "registered": {
                name: {
                    "name": item.name,
                    "path": item.path,
                    "resolution": item.resolution,
                }
                for name, item in command_parser.registered.items()
            },
            "registered_symlink_root": command_parser.registered_symlink_root,
        },
    }
    calc_block = cfg.get("Calculators") if isinstance(cfg.get("Calculators"), Mapping) else {}
    pools = None
    if isinstance(calc_block, Mapping):
        pools = calc_block.get("Pools") or calc_block.get("pools")
    if isinstance(pools, Mapping) and pools and "calculator_pools" not in extra_payload:
        worker_config["calculator_pools"] = {
            str(name): max(1, int(count or 1)) for name, count in pools.items()
        }

    if extra_payload:
        worker_config.update(extra_payload)
    return worker_config


__all__ = ["build_command_parser", "build_worker_config"]