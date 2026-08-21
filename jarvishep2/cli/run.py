#!/usr/bin/env python3
"""CLI adapters for validate / run / check / convert (D25.9)."""

from __future__ import annotations

import json
import sys

from jarvishep2.logging import get_jarvis_logger
from jarvishep2.run_outcome import (
    EXIT_OK,
    EXIT_RUN_FAILED,
    EXIT_USAGE,
    RunOutcome,
)

def _log_outcome(outcome: RunOutcome) -> None:
    """Write the final run summary through the configured process logger."""
    logger = get_jarvis_logger("core")
    ncall = None
    if isinstance(outcome.extras, dict):
        ncall = outcome.extras.get("ncall")
    ncall_bit = f" ncall={ncall}" if ncall is not None else ""
    logger.warning(
        f"RunOutcome status={outcome.status} "
        f"submitted={outcome.submitted} completed={outcome.completed} "
        f"failed={outcome.failed} archived={outcome.archived}"
        f"{ncall_bit} "
        f"exit={outcome.exit_code}",
    )
    if outcome.error:
        logger.error(
            "RunOutcome error_type=%s error=%s",
            outcome.error_type,
            outcome.error,
        )


def _emit_load_diagnostic(
    *,
    task_yaml: str,
    code: str,
    message: str,
    suggestion: str,
    example: str | None,
    as_json: bool,
) -> None:
    """Render pre-validation YAML/loading failures like normal diagnostics."""
    if as_json:
        print(json.dumps({
            "ok": False,
            "task_yaml": task_yaml,
            "scan_name": "",
            "issues": [{
                "level": "error", "code": code, "path": "$",
                "message": message, "hint": None,
                "suggestion": suggestion, "example": example,
            }],
        }, indent=2, sort_keys=True))
        return
    print(f"Config validation failed (1 error, 0 warnings):\n"
          f"  [error] {code}  $\n"
          f"          {message}\n"
          f"          suggestion: {suggestion}", file=sys.stderr)
    if example:
        print("          example:\n" + "\n".join(
            f"            {line}" for line in example.splitlines()
        ), file=sys.stderr)


def dispatch_validate(
    task_yaml: str,
    *,
    strict: bool = False,
    as_json: bool = False,
    check_modules: bool | None = None,
) -> int:
    """Load + validate a task card; never start Redis or workers."""
    from jarvishep2.task_config import TaskCardLoadError
    from jarvishep2.task_validation import (
        ConfigValidationError,
        format_report,
        format_report_success,
        validate_task_config,
    )

    if not task_yaml:
        print("Task YAML is required.", file=sys.stderr)
        return EXIT_USAGE

    try:
        from jarvishep2.task_config import load_task_yaml

        config = load_task_yaml(task_yaml)
    except FileNotFoundError as exc:
        _emit_load_diagnostic(
            task_yaml=task_yaml, code="JV2-LOAD-001", message=str(exc),
            suggestion="Provide the path to an existing task YAML file.",
            example="Jarvis validate path/to/task.yaml", as_json=as_json,
        )
        return EXIT_RUN_FAILED
    except ValueError as exc:
        _emit_load_diagnostic(
            task_yaml=task_yaml,
            code=getattr(exc, "code", "JV2-LOAD-002"),
            message=str(exc),
            suggestion=getattr(exc, "suggestion", "Correct the task-card structure described in the error and validate again."),
            example=getattr(exc, "example", None),
            as_json=as_json,
        )
        return EXIT_USAGE

    report = validate_task_config(
        config, strict=strict, check_modules=check_modules
    )
    if as_json:
        payload = {
            "ok": report.ok,
            "task_yaml": str(config.get("task_yaml") or task_yaml),
            "scan_name": str(config.get("scan_name") or ""),
            "issues": [
                {
                    "level": i.level,
                    "code": i.code,
                    "path": i.path,
                    "message": i.message,
                    "hint": i.hint,
                    "suggestion": i.suggestion,
                    "example": i.example,
                }
                for i in report.issues
            ],
        }
        print(json.dumps(payload, indent=2, sort_keys=True))
        return EXIT_OK if report.ok else EXIT_USAGE

    if report.ok:
        # Mirror run-path logging: explicit success line on stderr for humans.
        print(
            "Config validation successful. The task YAML meets the V2 contract.",
            file=sys.stderr,
        )
        if report.warnings():
            print(format_report(report), file=sys.stderr)
        else:
            # Keep format_report_success available for tests / quieter tooling.
            _ = format_report_success(report)
        return EXIT_OK

    print(format_report(report), file=sys.stderr)
    return EXIT_USAGE


def dispatch_convert(task_yaml: str) -> int:
    """Refresh DATABASE HDF5 → CSV snapshots when their source content changes."""
    if not task_yaml:
        print("Task YAML is required for convert.", file=sys.stderr)
        return EXIT_USAGE

    from jarvishep2.core import Jarvis2Core
    from jarvishep2.task_config import TaskCardLoadError

    try:
        core = Jarvis2Core()
        # Post-run utility: path resolution only; skip full config gate.
        core.load_task_yaml(task_yaml, validate=False)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except TaskCardLoadError as exc:
        _emit_load_diagnostic(
            task_yaml=task_yaml, code=exc.code, message=str(exc),
            suggestion=exc.suggestion, example=exc.example, as_json=False,
        )
        return EXIT_USAGE
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    results = core.convert()
    if not results:
        return EXIT_RUN_FAILED

    if any(item.get("status") in {"converted", "skipped_unchanged"} for item in results):
        return EXIT_OK
    return EXIT_RUN_FAILED


def dispatch_run(
    task_yaml: str,
    *,
    resume: bool = False,
    check_modules: bool = False,
    skip_draw_flowchart: bool = False,
    console_level: str = "WARNING",
    silence: bool = False,
    debug: bool = False,
    strict: bool = False,
    skip_library_installation: bool = False,
    check_timeout: float | None = None,
) -> int:
    if not task_yaml:
        print("Task YAML is required.", file=sys.stderr)
        return EXIT_USAGE

    from jarvishep2.core import Jarvis2Core
    from jarvishep2.task_config import TaskCardLoadError

    try:
        core = Jarvis2Core()
        core.load_task_yaml(
            task_yaml,
            validate=True,
            strict=strict,
            check_modules=check_modules if check_modules else None,
        )
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_RUN_FAILED
    except TaskCardLoadError as exc:
        _emit_load_diagnostic(
            task_yaml=task_yaml, code=exc.code, message=str(exc),
            suggestion=exc.suggestion, example=exc.example, as_json=False,
        )
        return EXIT_USAGE
    except ValueError as exc:
        # Includes ConfigValidationError (subclass of ValueError).
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE

    if debug:
        console_level = "DEBUG"

    core.set_logging_options(
        console_level=console_level,
        silence=silence,
    )
    core.set_skip_library_installation(skip_library_installation)

    if skip_draw_flowchart:
        core.config["skip_draw_flowchart"] = True
        core.info["skip_draw_flowchart"] = True

    try:
        if check_modules:
            outcome = core.check_modules(timeout=check_timeout)
        else:
            outcome = core.run(resume=resume)
    except NotImplementedError as exc:
        print(str(exc), file=sys.stderr)
        return EXIT_USAGE
    except (RuntimeError, TimeoutError, OSError) as exc:
        print(f"Run failed: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    if not isinstance(outcome, RunOutcome):
        # Backward safety if a test double returns a bare int.
        try:
            count = int(outcome)  # type: ignore[arg-type]
        except Exception:
            count = 0
        return EXIT_OK if count > 0 else EXIT_RUN_FAILED

    _log_outcome(outcome)
    return int(outcome.exit_code)
