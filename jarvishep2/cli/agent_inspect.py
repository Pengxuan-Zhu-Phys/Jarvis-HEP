"""Explain loaded task settings without bringing up a scan or user operators."""

from __future__ import annotations

import argparse
import json
import math
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from jarvishep2.cli.run import _manual_navigation, _validation_metadata
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE
from jarvishep2.versioning import get_runtime_version


CONTRACT_VERSION = 1
_MISSING = object()
_SECTION_FIELDS = {
    "paths": ("paths", "path_notes"),
    "runtime": ("runtime", "runtime_decisions"),
    "sources": ("sources",),
    "validation": (),
    "card": ("raw_task_card", "merged_envreqs"),
    "policies": ("policies", "policy_notes"),
}


class _UsageError(ValueError):
    pass


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise _UsageError(message)


class _Preview:
    """Bound large user documents explicitly; never claim omitted data is absent."""

    def __init__(self) -> None:
        self.nodes = 1500
        self.characters = 20000
        self.omissions: list[dict[str, Any]] = []
        self.omission_count = 0

    def omitted(self, path: str, reason: str) -> None:
        self.omission_count += 1
        if len(self.omissions) < 20:
            self.omissions.append({"path": path[:200], "reason": reason})

    def take(self, value: Any, path: str, depth: int = 0) -> Any:
        if self.nodes <= 0 or depth > 16:
            self.omitted(path, "node/depth preview limit")
            return {"preview_omitted": True}
        self.nodes -= 1
        if isinstance(value, Mapping):
            result = {}
            for index, (key, child) in enumerate(value.items()):
                if index >= 50 or self.nodes <= 0:
                    self.omitted(path, "remaining mapping entries omitted")
                    break
                name = str(key)
                if not isinstance(key, str):
                    name = f"[{type(key).__name__}] {key}"
                    self.omitted(path, "non-string YAML key represented with its type")
                if len(name) > 200:
                    self.omitted(path, "mapping key exceeds 200 characters")
                    name = name[:200] + f" [entry {index}]"
                result[name] = self.take(child, f"{path}.{name}", depth + 1)
            return result
        if isinstance(value, (list, tuple)):
            result = []
            for index, child in enumerate(value):
                if index >= 50 or self.nodes <= 0:
                    self.omitted(path, "remaining list items omitted")
                    break
                result.append(self.take(child, f"{path}[{index}]", depth + 1))
            return result
        if isinstance(value, str):
            limit = min(3000, max(0, self.characters))
            result = value[:limit]
            self.characters -= len(result)
            if len(result) < len(value):
                self.omitted(path, "string/text preview limit")
                result += " [preview truncated]"
            return result
        if isinstance(value, float) and not math.isfinite(value):
            return {"yaml_scalar_type": "float", "value": str(value)}
        if value is None or isinstance(value, (bool, int, float)):
            return value
        # YAML can construct dates and other non-JSON scalar types.
        return {"yaml_scalar_type": type(value).__name__, "value": self.take(str(value), path, depth + 1)}


def _base_payload(task_yaml: str) -> dict[str, Any]:
    metadata = _validation_metadata(completed=False)
    metadata.update({
        "inspection_scope": "Loaded defaults, configuration normalization, planned paths, and static validation for run mode",
        "user_operator_imports_performed": False,
        "resource_preflight_performed": False,
        "limitations": [
            *metadata["limitations"],
            "No user operator availability/signature checks, shell execution, installations, port probes, or scan resource creation are performed.",
            "Paths describe loaded or planned locations; inspect does not create them or verify runtime outputs.",
            "When validation fails, normalized values explain loader behavior and are not a proposed repair.",
            "Jarvis check has separate temporary worker and output overrides; this response describes run configuration.",
        ],
    })
    return {
        "ok": False, "kind": "agent.inspect", "contract_version": CONTRACT_VERSION,
        "installed_version": get_runtime_version(), "task_yaml": task_yaml,
        "scan_name": "", "metadata": metadata, "issues": [],
    }


def _emit(payload: dict[str, Any], *, section: str | None = None) -> None:
    if section is not None:
        selected_fields = set(_SECTION_FIELDS[section])
        data_fields = {field for fields in _SECTION_FIELDS.values() for field in fields}
        for key in data_fields - selected_fields:
            payload.pop(key, None)
        payload["metadata"].update({
            "requested_section": section,
            "omitted_sections": [name for name in _SECTION_FIELDS if name != section],
            "section_scope": "The complete card is loaded and validated. Other sections are omitted by request, not absent or unvalidated.",
        })
    preview = _Preview()
    # Preserve response/navigation metadata; bound user-controlled document data.
    # Effective settings take precedence over potentially huge user kwargs.
    for key in ("paths", "scan_name", "runtime", "runtime_decisions", "policies", "policy_notes", "sources", "issues", "explanation_unavailable", "merged_envreqs", "raw_task_card"):
        if key in payload:
            payload[key] = preview.take(payload[key], f"$.{key}")
    payload["preview"] = {
        "truncated": bool(preview.omission_count),
        "omission_count": preview.omission_count, "omissions": preview.omissions,
        "limits": {"nodes": 1500, "text_characters": 20000, "items_per_container": 50, "depth": 16},
        "note": "Truncation affects display only. Validation uses the complete loaded card; inspect the original files for omitted data.",
    }
    print(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False))


def _failure(payload: dict[str, Any], *, code: str, message: str, suggestion: str,
             exception: Exception | None = None) -> None:
    payload["ok"] = False
    payload["metadata"]["validation_completed"] = False
    diagnostic = {
        "level": "error", "code": code, "path": getattr(exception, "path", "$"),
        "message": message, "hint": None, "suggestion": suggestion,
        "example": getattr(exception, "example", None), **_manual_navigation(code),
    }
    for name in ("location", "related_locations"):
        value = getattr(exception, name, None)
        if value is not None:
            diagnostic[name] = value
    payload["issues"] = [diagnostic]
    payload["issue_count"] = 1


def _at(document: Any, parts: list[str]) -> Any:
    for part in parts:
        if not isinstance(document, Mapping) or part not in document:
            return _MISSING
        document = document[part]
    return document


def _leaves(document: Any, parts: tuple[str, ...] = ()):
    if isinstance(document, Mapping) and document:
        for key, value in document.items():
            yield from _leaves(value, (*parts, str(key)))
    else:
        yield parts


def _origin(parts: list[str], task: Mapping[str, Any], defaults: Mapping[str, Any],
            task_path: str, defaults_path: str | None) -> dict[str, Any]:
    if _at(task, parts) is not _MISSING:
        return {"layer": "task", "file": task_path, "input_path": "$.EnvReqs." + ".".join(parts)}
    if _at(defaults, parts) is not _MISSING:
        return {"layer": "default_environment", "file": defaults_path, "input_path": "$.EnvReqs." + ".".join(parts)}
    return {"layer": "internal_normalization", "file": None, "input_path": None}


def _runtime_input(parts: tuple[str, ...]) -> list[str]:
    key, *rest = parts
    if key in {"workers", "batch_size", "store_samples", "checkpoint", "redis"}:
        return ["V2", key, *rest]
    if key in {"sample_artifacts", "force_serial_layers"}:
        return ["V2", "worker", key, *rest]
    if key == "Factory":
        return ["V2", "factory", *rest]
    if key == "Watchdog":
        return ["V2", "factory", "watchdog", *rest]
    return []


def _explain(config: Mapping[str, Any]) -> dict[str, Any]:
    from jarvishep2.base import project_images_dir
    from jarvishep2.logging.toplevel import scan_logs_dir
    from jarvishep2.runtime_config import (
        get_archiver_config, get_checkpoint_config, get_factory_config,
        get_redis_config, get_runtime_block, get_sample_directory_config,
    )
    from jarvishep2.sampling.runtime_checkpoint import checkpoint_path
    from jarvishep2.task_config import (
        _load_yaml_document, get_check_modules_settings, resolve_sampling_path,
        sampling_method,
    )

    task_path = str(config["task_yaml"])
    root, scan = str(config["project_root"]), str(config["scan_name"])
    output = str(config["task_result_dir"])
    defaults_path = config.get("environment_defaults_path")
    raw = config.raw_task_card
    task_envreqs = raw.get("EnvReqs") or {}
    defaults = {}
    if defaults_path:
        document = _load_yaml_document(str(defaults_path), label="default environment YAML")
        if isinstance(document, Mapping) and isinstance(document.get("EnvReqs"), Mapping):
            defaults = document["EnvReqs"]
    envreqs = config.get("EnvReqs") or {}
    runtime = get_runtime_block(config)
    method = sampling_method(config) if isinstance(config.get("Sampling"), Mapping) else ""
    check_settings = get_check_modules_settings(config)
    check_data = str(check_settings.get("data") or "")
    sources = {
        "task_yaml": task_path, "default_environment_yaml": defaults_path,
        "merge_order": ["default environment EnvReqs", "task EnvReqs overrides", "internal runtime normalization"],
        "envreqs_fields": {
            "$.EnvReqs." + ".".join(parts): _origin(list(parts), task_envreqs, defaults, task_path, defaults_path)
            for parts in _leaves(envreqs)
        },
        "runtime_fields": {
            "$.Runtime." + ".".join(parts): {
                **(_origin(_runtime_input(parts), task_envreqs, defaults, task_path, defaults_path)
                   if _runtime_input(parts) else {"layer": "internal_normalization", "file": None, "input_path": None}),
                "normalizer": "get_runtime_block",
            }
            for parts in _leaves(runtime)
        },
        "policy_functions": {
            "checkpoint": "get_checkpoint_config", "redis": "get_redis_config",
            "archiver": "get_archiver_config", "sample_directory": "get_sample_directory_config",
            "factory": "get_factory_config", "check_modules": "get_check_modules_settings",
        },
    }
    return {
        "scan_name": scan,
        "paths": {
            "project_root": root, "environment_defaults_yaml": defaults_path,
            "scan_output_root": output, "database_dir": os.path.join(output, "DATABASE"),
            "sample_dir": os.path.join(output, "SAMPLE"), "logs_dir": scan_logs_dir(root, scan),
            "images_dir": project_images_dir(project_root=root, scan_name=scan),
            "declared_checkpoint_file": checkpoint_path(task_root=root, scan_name=scan, sampler_name=method) if method else None,
            "check_modules_data": resolve_sampling_path(config, check_data) if check_data else None,
        },
        "path_notes": {
            "project_root": "Resolved by the task loader: configured task-root environment variable, nearest project marker/bin parent, or task directory.",
            "output_paths": "Planned run locations derived from the loaded scan output root; files/directories are not created or verified.",
            "sample_storage_enabled": runtime["store_samples"],
            "sample_dir": "Conventional SAMPLE target. When store_samples is false, this target is inactive and per-sample work uses disposable scratch space.",
            "declared_checkpoint_file": "Pre-run path for the declared Sampling.Method using checkpoint_path; the running sampler's canonical method can change this path. No checkpoint is read.",
            "check_modules_data": "Resolved configured/default CSV location, without checking its contents or selecting the runtime fallback.",
        },
        "raw_task_card": raw, "merged_envreqs": envreqs, "runtime": runtime,
        "policies": {
            "checkpoint": get_checkpoint_config(config), "redis": get_redis_config(config),
            "archiver": get_archiver_config(config), "sample_directory": get_sample_directory_config(config),
            "factory": get_factory_config(config), "check_modules": check_settings,
        },
        "policy_notes": {
            "scope": "Normalized run configuration only; policies are not preflight results or evidence that execution succeeded.",
            "checkpoint": "heartbeat is in seconds. Invalid input may normalize to a default/minimum while still failing validation; this is not an automatic repair.",
            "redis": "The configured port has not been probed or reserved; startup can choose a different port when permitted.",
            "check_modules": "These settings apply only when Jarvis check executes; inspect does not activate check mode or select its CSV fallback.",
        },
        "sources": sources,
        "runtime_decisions": {
            "configured_workers": runtime["workers"], "selected_workers": None,
            "workers_note": "workers=0 requests runtime automatic sizing. inspect does not select a Worker count; check mode separately forces one worker.",
            "configured_redis_port": get_redis_config(config)["port"], "selected_redis_port": None,
            "redis_port_task_override": bool(config.get("redis_port_task_override")),
            "redis_note": "The port is configured, not probed or reserved. Startup preflight decides whether it is usable and may update project defaults when allowed.",
        },
    }


def dispatch_inspect(argv: list[str] | None = None) -> int:
    """Emit one JSON response for an explicit card; never initialize execution."""
    parser = _Parser(prog="Jarvis agent inspect", add_help=False, allow_abbrev=False,
                     description="Explain loaded task defaults, planned paths and validation without starting a scan. JSON is the default.")
    parser.add_argument("task_yaml", nargs="?", help="Explicit path to the user's task YAML")
    parser.add_argument("--json", action="store_true", help="JSON output (the default)")
    parser.add_argument("--strict", action="store_true", help="Treat validation warnings as errors")
    parser.add_argument("--section", choices=tuple(_SECTION_FIELDS),
                        help="Return only this section plus validation/navigation metadata; the full card is still validated")
    parser.add_argument("-h", "--help", action="store_true", help="Show this help")
    raw_argv = list(argv or [])
    payload = _base_payload("")
    section = None
    try:
        args = parser.parse_args(raw_argv)
        section = args.section
        if args.help:
            if args.json:
                print(json.dumps({"contract_version": CONTRACT_VERSION, "help": parser.format_help()}, indent=2))
            else:
                parser.print_help()
            return EXIT_OK
        if not args.task_yaml:
            raise _UsageError("An explicit task YAML path is required.")
        payload["task_yaml"] = str(Path(args.task_yaml).absolute())
    except _UsageError as exc:
        _failure(payload, code="JV2-CLI-001", message=str(exc), suggestion=parser.format_usage().strip())
        _emit(payload, section=section)
        return EXIT_USAGE

    try:
        from jarvishep2.task_config import load_task_yaml
        config = load_task_yaml(args.task_yaml)
    except FileNotFoundError as exc:
        _failure(payload, code="JV2-LOAD-001", message=str(exc), suggestion="Provide an existing task YAML and readable referenced default files.", exception=exc)
        _emit(payload, section=section)
        return EXIT_RUN_FAILED
    except ValueError as exc:
        _failure(payload, code=getattr(exc, "code", "JV2-LOAD-002"), message=str(exc),
                 suggestion=getattr(exc, "suggestion", "Correct the reported task/default configuration and inspect again."), exception=exc)
        _emit(payload, section=section)
        return EXIT_USAGE
    except Exception as exc:
        _failure(payload, code="JV2-LOAD-003" if isinstance(exc, OSError) else "JV2-LOAD-004",
                 message=f"Inspection loading failed ({type(exc).__name__}): {exc}",
                 suggestion="Resolve the file/dependency failure and inspect again; validation did not complete.", exception=exc)
        _emit(payload, section=section)
        return EXIT_RUN_FAILED

    try:
        from jarvishep2.task_validation import validate_task_config
        report = validate_task_config(config, strict=args.strict)
        payload["ok"] = report.ok
        payload["metadata"]["validation_completed"] = True
        payload["issue_count"] = len(report.issues)
        payload["issues"] = [{
            "level": item.level, "code": item.code, "path": item.path,
            "message": item.message, "hint": item.hint, "suggestion": item.suggestion,
            "example": item.example, **_manual_navigation(item.code),
        } for item in report.issues]
        try:
            payload.update(_explain(config))
        except (ValueError, TypeError, AttributeError) as exc:
            if report.ok:
                raise
            # Invalid types may prevent normalized policy readers from running.
            # Keep the complete validator diagnostics instead of masking them.
            payload["raw_task_card"] = config.raw_task_card
            payload["merged_envreqs"] = config.get("EnvReqs")
            payload["scan_name"] = str(config.get("scan_name") or "")
            payload["paths"] = {
                "project_root": config.get("project_root"),
                "environment_defaults_yaml": config.get("environment_defaults_path"),
                "scan_output_root": config.get("task_result_dir"),
            }
            payload["explanation_unavailable"] = f"Invalid configuration prevents policy normalization ({type(exc).__name__}): {exc}"
        _emit(payload, section=section)
        return EXIT_OK if report.ok else EXIT_USAGE
    except Exception as exc:
        _failure(payload, code="JV2-VAL-001", message=f"Inspection validation failed ({type(exc).__name__}): {exc}",
                 suggestion="Resolve the validator/dependency failure and inspect again; the card has not been validated.", exception=exc)
        _emit(payload, section=section)
        return EXIT_RUN_FAILED


__all__ = ["dispatch_inspect"]
