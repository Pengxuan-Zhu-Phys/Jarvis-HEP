#!/usr/bin/env python3
"""Worker blueprint builder with Phase-1 command resolution (WP-D3.1)."""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.calculator_modes import expand_calculator_modes
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_cleanup_config,
    get_delete_method,
    get_sample_directory_config,
    get_staging_dir,
    handoff_to_staging_enabled,
    workflow_has_calculator,
    workflow_references_sdir,
)


def _default_mapper(cfg: Mapping[str, Any]) -> dict[str, Any]:
    """Build picklable Worker mapper config (D22: MapperPipeline / MapperSpec).

    Top-level ``Mapper`` is not a V2 interface (schema-rejected).  Optional
    Optional ``Sampling.Mapper`` (name → expression) is folded into a pipeline
    spec so Workers share the same u→params implementation as the control process.
    """
    sampling = cfg.get("Sampling") if isinstance(cfg.get("Sampling"), Mapping) else {}
    method = str(sampling.get("Method") or "").strip()
    if method == "CSV":
        return {"type": "none"}
    variables = sampling.get("Variables")
    if not variables:
        return {"type": "identity", "keys": ["x", "y"]}
    try:
        from jarvishep2.mapper import MapperSpec, build_mapper_spec_from_config

        spec: MapperSpec = build_mapper_spec_from_config(cfg)
        return {"type": "pipeline", "spec": spec.to_dict()}
    except Exception:
        # Fall back to distribution-only so a Worker can still start if the
        # card was not re-validated (tests inject partial configs).
        return {"type": "distribution", "variables": list(variables)}


def _config_references_sdir(modules: list[dict[str, Any]]) -> bool:
    import json

    return "@Sdir" in json.dumps(modules, sort_keys=True)


def build_command_parser(
    config: Mapping[str, Any] | None,
    *,
    root_path: str | None = None,
    register_executables: bool = True,
) -> CommandParser:
    """Build a CommandParser from a task config mapping."""
    cfg = dict(config or {})
    project_root = str(cfg.get("project_root") or cfg.get("task_root") or "")
    return CommandParser.from_config(
        cfg,
        project_root=project_root or None,
        root_path=root_path,
        register_executables=register_executables,
    )


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
    calc_modules = expand_calculator_modes(calc_modules)

    opera = list(opera_modules or extra_payload.pop("opera_modules", []) or [])
    if not opera:
        operas = (cfg.get("Operas") or {}).get("Modules") or []
        if isinstance(operas, list):
            opera = [dict(item) for item in operas if isinstance(item, Mapping)]

    resolved_calculators = prepare_calculator_modules(calc_modules, command_parser)
    sample_root = sample_dirs or os.path.join(task_result_dir, "SAMPLE")
    sample_config = dict(extra_payload.pop("sample_config", {}) or {})
    sample_config.setdefault("task_result_dir", task_result_dir)
    sample_config.setdefault("sample_dirs", sample_root)
    # Used by Worker OS process title (Jarvis2-Worker-N:scan).
    scan_name = str(cfg.get("scan_name") or sample_config.get("scan_name") or "").strip()
    if scan_name:
        sample_config.setdefault("scan_name", scan_name)
    runtime_block = cfg.get("Runtime") if isinstance(cfg.get("Runtime"), Mapping) else {}
    sample_config.setdefault(
        "sample_artifacts",
        str((runtime_block or {}).get("sample_artifacts", "auto")),
    )
    sample_config.setdefault(
        "workflow_has_calculator",
        workflow_has_calculator(cfg) or bool(resolved_calculators),
    )
    sample_config.setdefault(
        "workflow_references_sdir",
        workflow_references_sdir(cfg) or _config_references_sdir(resolved_calculators),
    )
    sampling_method = str((cfg.get("Sampling") or {}).get("Method", "")).strip()
    publish_feedback = bool(extra_payload.pop("publish_feedback", False))
    # Feedback channel for stateful / feedback-driven samplers (ALS, MCMC family).
    if not publish_feedback and sampling_method:
        try:
            from jarvishep2.distributor import STATELESS_METHODS

            if sampling_method not in STATELESS_METHODS:
                publish_feedback = True
        except Exception:
            if sampling_method in (
                "AdaptiveBridson",
                "MCMC",
                "AMMCMC",
                "AM",
                "DRAM",
            ):
                publish_feedback = True

    from jarvishep2.Module.nuisance import extract_nuisance_config
    from jarvishep2.feedback_return import resolve_feedback_return

    nuisance_cfg = extract_nuisance_config(cfg)
    # D13.8: flat feedback projection policy (default {uuid, logL}).
    feedback_return = extra_payload.pop("feedback_return", None)
    if not isinstance(feedback_return, dict):
        feedback_return = resolve_feedback_return(
            cfg,
            sampler=extra_payload.pop("sampler", None),
            method=sampling_method,
        )
    else:
        extra_payload.pop("sampler", None)
    worker_config: dict[str, Any] = {
        "sample_config": sample_config,
        "mapper": extra_payload.pop("mapper", _default_mapper(cfg)),
        "opera_modules": opera,
        "calculator_modules": resolved_calculators,
        "likelihood_expressions": list(
            likelihood_expressions
            or (cfg.get("Sampling") or {}).get("LogLikelihood")
            or []
        ),
        "nuisance_config": nuisance_cfg,
        "pull_timeout": 1,
        "delete_method": get_delete_method(cfg),
        "cleanup_config": get_cleanup_config(cfg),
        "archiver_config": get_archiver_config(cfg),
        "sample_directory": get_sample_directory_config(cfg),
        "command_parser": command_parser.to_picklable(),
        "publish_feedback": publish_feedback,
        "feedback_return": feedback_return,
    }
    # Only resolve/create staging when the staging hop is actually enabled.
    handoff = handoff_to_staging_enabled(cfg)
    worker_config["handoff_to_staging"] = handoff
    if handoff:
        worker_config["staging_dir"] = get_staging_dir(
            cfg, task_result_dir=task_result_dir, create=True
        )
    else:
        worker_config["staging_dir"] = ""
    if scan_name:
        worker_config.setdefault("scan_name", scan_name)
    # Component logs: always <task_root>/logs/<scan>/worker-NN.log (never cwd jarvis_worker_PID).
    from jarvishep2.logging import scan_logs_dir

    scan_for_logs = scan_name or str(
        worker_config.get("scan_name") or sample_config.get("scan_name") or "scan"
    ).strip() or "scan"
    task_root = str(
        cfg.get("task_root")
        or cfg.get("project_root")
        or worker_config.get("task_root")
        or ""
    ).strip()
    if not task_root:
        trd = os.path.abspath(str(task_result_dir or "."))
        parent = os.path.dirname(trd)
        # Usual layout: <task_root>/outputs/<scan>
        if os.path.basename(parent) == "outputs":
            task_root = os.path.dirname(parent)
        else:
            task_root = parent or os.getcwd()
    worker_config.setdefault("logs_dir", scan_logs_dir(task_root, scan_for_logs))
    worker_config.setdefault("task_root", task_root)
    # Console policy (CLI --console-level / --silence) propagated from control.
    if "log_silence" in extra_payload:
        worker_config["log_silence"] = bool(extra_payload.pop("log_silence"))
    if "console_level" in extra_payload:
        worker_config["console_level"] = str(extra_payload.pop("console_level"))
    # EnvReqs.V2.worker.force_serial_layers → Runtime → Worker blueprint.
    if isinstance(runtime_block, Mapping) and "force_serial_layers" in runtime_block:
        worker_config["force_serial_layers"] = bool(runtime_block.get("force_serial_layers"))
    calc_block = cfg.get("Calculators") if isinstance(cfg.get("Calculators"), Mapping) else {}
    pools = None
    if isinstance(calc_block, Mapping):
        pools = calc_block.get("Pools") or calc_block.get("pools")
        # V1 invariant 12: keep the make_paraller spelling; consume for pool sizing.
        if "calculator_make_paraller" not in worker_config and "make_paraller" in calc_block:
            try:
                worker_config["calculator_make_paraller"] = max(
                    1, int(calc_block.get("make_paraller") or 1)
                )
            except (TypeError, ValueError):
                worker_config["calculator_make_paraller"] = 1
        # Top-level Calculators.path is V1 layout metadata; modules carry their own path.
        # Tolerate it (do not error); stamp resolved path for diagnostics only.
        calc_root = calc_block.get("path")
        if isinstance(calc_root, str) and calc_root.strip():
            try:
                resolved_root = command_parser.resolve_static(calc_root.strip())
            except (ValueError, KeyError):
                resolved_root = calc_root.strip()
            sample_config.setdefault("calculators_path", resolved_root)
    if isinstance(pools, Mapping) and pools and "calculator_pools" not in extra_payload:
        worker_config["calculator_pools"] = {
            str(name): max(1, int(count or 1)) for name, count in pools.items()
        }

    if extra_payload:
        worker_config.update(extra_payload)
    return worker_config


__all__ = ["build_command_parser", "build_worker_config"]
