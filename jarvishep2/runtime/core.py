#!/usr/bin/env python3
"""Jarvis control-process orchestrator for the distributed runtime.

D25.3: ``Jarvis2Core`` is the public façade. Redis/Factory/Archiver bring-up
lives on ``_RuntimeSupervisor``, scan loops on ``_ScanDriver``, and checkpoint
resume on ``_ResumeService``. Public ``run()`` / ``shutdown()`` signatures are
unchanged; tests and ``client.dispatch_run`` still construct this class.
"""

from __future__ import annotations

import atexit
import csv
import os
import signal
import sys
import threading
import time
import uuid
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.archiver import ArchiverProcess, SimpleArchiver
from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.calculator_modes import (
    expand_calculator_modes,
    expand_calculator_modes_in_config,
    mode_info,
)
from jarvishep2.library import LibraryInstaller
from jarvishep2.dashboard import SnapshotReader, format_monitor_view
from jarvishep2.database import (
    convert_database_dir,
    read_persisted_outcome_counts,
    read_persisted_sample_index_state,
    prepare_hdf5_database_for_writer,
)
from jarvishep2.environment_requirements import (
    EnvironmentRequirementError,
    check_environment_requirements,
)
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.runtime.factory import TaskFactory
from jarvishep2.log_kv import PermilleProgress, format_duration
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.monitoring.run_summary import (
    RunSummaryRenderer,
    build_run_summary,
    format_scan_performance_log,
)
from jarvishep2.versioning import render_logo_with_version
from jarvishep2.redis_queue import (
    CONTROL_LOCK_TTL_SEC,
    RedisQueue,
)
from jarvishep2.redis_server import (
    ManagedRedisServer,
    find_available_redis_port,
    redis_port_open,
)
from jarvishep2.task_config import update_default_redis_port
from jarvishep2.runtime_metadata import write_scan_metadata
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_delete_method,
    get_redis_config,
    get_runtime_block,
    get_sample_directory_config,
    get_watchdog_config,
    pack_buckets_enabled,
)
from jarvishep2.sample import Sample
from jarvishep2.Sampling.runtime_checkpoint import (
    RESUME_PROMPT,
    build_payload,
    build_run_spec,
    check_mapper_fingerprint,
    check_operas_constants_fingerprint,
    checkpoint_path,
    load_checkpoint,
    prepare_resume,
    save_checkpoint,
)
from jarvishep2.Sampling.sampler import SamplingVirtial
from jarvishep2.check_modules import (
    build_check_module_samples,
    verify_check_modules_golden,
)
from jarvishep2.run_outcome import RunOutcome
from jarvishep2.runtime.worker_config import build_command_parser, build_worker_config
from jarvishep2.Module.runtime_preparer import (
    prepare_install_controls,
    refresh_install_control_summaries,
)
from jarvishep2.task_config import (
    check_modules_n_samples,
    check_modules_timeout_sec,
    load_task_yaml,
    resolve_check_modules_csv,
    sampling_method,
)
from jarvishep2.task_validation import (
    ValidationReport,
    raise_if_errors,
    validate_task_config,
)


from jarvishep2.runtime._resume_service import _ResumeService
from jarvishep2.runtime._runtime_supervisor import _RuntimeSupervisor
from jarvishep2.runtime._scan_driver import _ScanDriver


class Jarvis2Core:
    """Distributed run façade: load/validate/run, holding private collaborators."""

    def __init__(self, config: Mapping[str, Any] | None = None) -> None:
        self.config = dict(config or {})
        self.runtime = get_runtime_block(self.config)
        self.info: dict[str, Any] = {}
        if self.config:
            self._populate_info_from_config()
        self.redis: RedisQueue | None = None
        self.factory: TaskFactory | None = None
        self.archiver: SimpleArchiver | ArchiverProcess | None = None
        self.sampler: Any = None
        self.command_parser: CommandParser | None = None
        self._environment_summary: dict[str, Any] = {}
        self._skip_library_installation = False
        self._resume_checkpoint_payload: dict[str, Any] | None = None
        self._resume_policy: str = "auto"
        self._control_lock_owner: str | None = None
        self._control_lease_stop: threading.Event | None = None
        self._control_lease_thread: threading.Thread | None = None
        self._persisted_uuids: set[str] = set()
        self._persisted_index_prefix = 0
        self._persisted_records_count = 0
        # D21.13: durable Failed rows must seed SAMPLE_STATS.failed on resume.
        self._persisted_failed_count = 0
        self._managed_redis: ManagedRedisServer | None = None
        self._shutdown_done = False
        # Captured after Archiver.stop() so the final CLI RunOutcome can use
        # the post-flush count rather than the pre-shutdown snapshot.
        self._final_archived_records: int | None = None
        self._interrupt_requested = False
        self._signal_handlers_installed = False
        self._previous_signal_handlers: dict[int, Any] = {}
        self._logger = get_jarvis_logger("core")
        # Console / file logging options (CLI sets before bootstrap).
        self._log_console: bool = True
        self._log_silence: bool = False
        self._console_level: str = "WARNING"
        # D13.9: last gate report (re-logged after init_logger so console/file see it).
        self._validation_report: ValidationReport | None = None
        self._runtime = _RuntimeSupervisor(self)
        self._scan = _ScanDriver(self)
        self._resume = _ResumeService(self)
        # Last-resort cleanup if the process exits without an orderly finally.
        atexit.register(self._atexit_cleanup)

    def _ensure_collaborators(self) -> None:
        if getattr(self, "_runtime", None) is None:
            self._runtime = _RuntimeSupervisor(self)
        if getattr(self, "_scan", None) is None:
            self._scan = _ScanDriver(self)
        if getattr(self, "_resume", None) is None:
            self._resume = _ResumeService(self)

    def _get_runtime(self) -> _RuntimeSupervisor:
        self._ensure_collaborators()
        return self._runtime

    def _get_scan(self) -> _ScanDriver:
        self._ensure_collaborators()
        return self._scan

    def _get_resume(self) -> _ResumeService:
        self._ensure_collaborators()
        return self._resume

    def set_logging_options(
        self,
        *,
        console_level: str | None = None,
        silence: bool | None = None,
    ) -> None:
        """CLI-facing console/file logging policy (applied in ``init_logger``)."""
        if console_level is not None:
            self._console_level = str(console_level).strip().upper() or "WARNING"
        if silence is not None:
            self._log_silence = bool(silence)
            self._log_console = not self._log_silence

    def set_skip_library_installation(self, skip: bool) -> None:
        """Control-process policy for the non-interactive LibDeps installer."""
        self._skip_library_installation = bool(skip)

    def init_logger(self) -> None:
        """Configure top-level logging with V1 visual format and scan-scoped file path."""
        from jarvishep2.logging import component_log_path, scan_logs_dir
        from jarvishep2.proc_title import control_title, set_process_title

        scan_name = str(self.info.get("scan_name") or "scan")
        # Refresh process title once the scan name is known.
        set_process_title(control_title(scan_name=scan_name))
        task_root = str(self.info.get("task_root") or os.getcwd())
        logs_dir = scan_logs_dir(task_root, scan_name)
        # Control multi-sink: core.log + factory.log + sampler.log + archiver.log
        # (omit log_path so setup opens filtered handlers per component).
        jarvis_log = component_log_path(logs_dir, "core")
        self.info["logs_dir"] = logs_dir
        self.info["jarvis_log"] = jarvis_log
        setup_jarvis_logging(
            role="core",
            component="core",
            scan_logs_dir=logs_dir,
            console=self._log_console,
            console_level=self._console_level,
            silence=self._log_silence,
            use_queue=True,
            multi_sink=True,
        )
        self._logger = get_jarvis_logger("core")
        try:
            self._logger.warning("\n" + render_logo_with_version())
            self._logger.warning("Jarvis-HEP V2 logging system initialized successful!")
            self._logger.info(
                "component logs under %s "
                "(core.log, factory.log, sampler.log, archiver.log, "
                "datarecorder.log; workers: worker-NN.log)",
                logs_dir,
            )
            self._logger.info("Jarvis-HEP write into main log file -> %s", jarvis_log)
            if self._log_silence:
                self._logger.info("console output silenced (--silence)")
            # Validation runs before logging is wired; replay so users see the gate.
            self._log_validation_report(self._validation_report)
        except Exception:
            # Logging must never block bootstrap; banner is best-effort.
            pass

    def check_environment_requirements(self) -> None:
        """Log and enforce V1-compatible OS, Python, and CERN ROOT requirements."""
        report = check_environment_requirements(self.config)
        self._environment_summary = dict(report.summary)
        if not report.summary and not report.warnings and not report.errors:
            return
        envreqs = self.config.get("EnvReqs")
        envreqs = envreqs if isinstance(envreqs, Mapping) else {}
        lines = [
            "Environment requirements preflight",
            "  Check                      | Required        | Detected        | Status",
            "  -------------------------- | --------------- | --------------- | --------",
        ]

        def format_check(
            label: str,
            required: str,
            detected: str,
            status: str,
            *,
            child: bool = False,
        ) -> str:
            """Render every environment check against the same column grid."""
            # Children have two extra spaces only *inside* the Check column;
            # shrink their label field by the same amount so every separator
            # remains on the parent row's column boundary.
            indent = "    " if child else "  "
            label_width = 24 if child else 26
            return (
                f"{indent}{label:<{label_width}} | {required:<15} | {detected:<15} | {status}"
            )

        python_requirement = envreqs.get("Python")
        if isinstance(python_requirement, Mapping):
            required = str(python_requirement.get("version") or "not specified")
            found = str(report.summary.get("Python version") or "not detected")
            passed = not any(error.startswith("Python ") for error in report.errors)
            lines.append(format_check("Python", required, found, "PASS" if passed else "FAIL"))
            dependencies = python_requirement.get("Dependencies")
            if isinstance(dependencies, Sequence) and not isinstance(dependencies, (str, bytes)):
                for dependency in dependencies:
                    if not isinstance(dependency, Mapping):
                        continue
                    name = str(dependency.get("name") or "unnamed")
                    required_version = str(dependency.get("version") or "any")
                    found_version = report.summary.get(name)
                    is_required = bool(dependency.get("required", False))
                    failed = any(
                        error.startswith(f"Python package {name}") for error in report.errors
                    )
                    warned = any(
                        warning.startswith(f"Python package {name}") for warning in report.warnings
                    )
                    status = "FAIL" if failed else "WARN" if warned else "PASS"
                    found_text = str(found_version) if found_version is not None else "not installed"
                    requirement_label = f">={required_version.removeprefix('>=').strip()}" if required_version != "any" else "any"
                    label = f"PyPI·{name}" + (" (optional)" if not is_required else "")
                    lines.append(
                        format_check(
                            label,
                            requirement_label,
                            found_text,
                            status,
                            child=True,
                        )
                    )

        os_requirements = envreqs.get("OS")
        if isinstance(os_requirements, Sequence) and not isinstance(os_requirements, (str, bytes)):
            detected = str(report.summary.get("OS") or "not detected")
            passed = not any(error.startswith(("Operating system ", "OS ")) for error in report.errors)
            lines.append(
                format_check(
                    "OS",
                    "configured" if os_requirements else "none",
                    detected,
                    "PASS" if passed else "FAIL",
                )
            )

        root_requirement = envreqs.get("CERN_ROOT")
        if isinstance(root_requirement, Mapping):
            required = bool(root_requirement.get("required", False))
            version_required = str(root_requirement.get("version") or "not specified")
            if not required:
                lines.append(format_check("CERN ROOT", "false", "not checked", "SKIPPED"))
            else:
                found = str(report.summary.get("ROOT version") or "not detected")
                passed = not any(error.startswith("CERN ROOT") for error in report.errors)
                lines.append(
                    format_check("CERN ROOT", version_required, found, "PASS" if passed else "FAIL")
                )
                dependencies = root_requirement.get("Dependencies")
                if isinstance(dependencies, Sequence) and not isinstance(dependencies, (str, bytes)):
                    for dependency in dependencies:
                        if not isinstance(dependency, Mapping):
                            continue
                        name = str(dependency.get("name") or "unnamed")
                        feature_found = bool(report.summary.get(f"ROOT-{name}", False))
                        feature_required = bool(dependency.get("required", False))
                        status = "PASS" if feature_found else "FAIL" if feature_required else "WARN"
                        lines.append(
                            format_check(
                                f"feature {name}",
                                "yes" if feature_required else "optional",
                                "yes" if feature_found else "no",
                                status,
                                child=True,
                            )
                        )
        if report.errors:
            lines.append(format_check("Result", "-", "-", "FAIL"))
            self._logger.error("\n%s", "\n".join(lines))
            raise EnvironmentRequirementError(report)
        lines.append(format_check("Result", "-", "-", "PASS"))
        # Use warning so the complete one-record summary is visible at V2's
        # default terminal threshold; it is one logging call, not one per check.
        self._logger.warning("\n%s", "\n".join(lines))

    def load_task_yaml(
        self,
        path: str,
        *,
        validate: bool = True,
        strict: bool = False,
        check_modules: bool | None = None,
    ) -> dict[str, Any]:
        """Load task YAML and merge normalized layout into ``self.config``.

        By default runs the D13.9 config gate (:func:`validate_task_config`) before
        any Redis/Worker bootstrap.  Pass ``validate=False`` only for internal
        tests that intentionally exercise malformed cards past the loader.
        """
        self.config = load_task_yaml(path)
        self.runtime = get_runtime_block(self.config)
        self._populate_info_from_config()
        if validate:
            self.validate_loaded_config(strict=strict, check_modules=check_modules)
        return self.config

    def validate_loaded_config(
        self,
        *,
        strict: bool = False,
        check_modules: bool | None = None,
    ) -> ValidationReport:
        """Run the pure D13.9 gate on ``self.config``; raise on errors.

        The report is stored and re-emitted after :meth:`init_logger` so the
        success line appears in console/file logs (validation itself runs
        before the logging system is fully configured).
        """
        report = validate_task_config(
            self.config,
            strict=strict,
            check_modules=check_modules,
        )
        self._validation_report = report
        # Success is logged after init_logger (no sinks yet here). Failures raise below.
        raise_if_errors(report)
        # D20's user-facing parent/mode syntax is a configuration-time macro.
        # Everything below this boundary (workflow, Redis, Workers) receives
        # ordinary dotted calculator module names only.
        self.config = expand_calculator_modes_in_config(self.config)
        return report

    def _log_validation_report(self, report: ValidationReport | None) -> None:
        """Emit human-visible validation status (V1-style success line).

        Called from :meth:`init_logger` after the logging system is wired so the
        message appears on the console and in ``core.log``.
        """
        if report is None:
            return
        n_err = len(report.errors())
        n_warn = len(report.warnings())
        try:
            logger = get_jarvis_logger("config")
        except Exception:
            return

        if n_err:
            logger.error(
                "Config validation failed (%d error(s), %d warning(s))",
                n_err,
                n_warn,
            )
            for item in report.errors():
                logger.error("  [%s] %s: %s", item.code, item.path, item.message)
            for item in report.warnings():
                logger.warning("  [%s] %s: %s", item.code, item.path, item.message)
            return

        if n_warn:
            # warning level: same visibility as the logo banner under default levels
            logger.warning(
                "Config validation successful with %d warning(s). "
                "The task YAML meets the V2 contract (non-fatal issues below).",
                n_warn,
            )
            for item in report.warnings():
                logger.warning("  [%s] %s: %s", item.code, item.path, item.message)
            return

        logger.warning(
            "Config validation successful. The task YAML meets the V2 contract."
        )

    def _populate_info_from_config(self) -> None:
        from jarvishep2.logging import component_log_path, scan_logs_dir

        task_root = str(self.config.get("task_root") or os.getcwd())
        scan = self.config.get("Scan")
        scan_mapping = scan if isinstance(scan, Mapping) else {}
        # ``load_task_yaml`` projects Scan.name to the legacy top-level alias,
        # but direct Python callers may pass a raw task-card mapping.  Preserve
        # the task-card value in both paths so scan isolation never silently
        # collapses to the default name ``scan``.
        scan_name = str(
            self.config.get("scan_name") or scan_mapping.get("name") or "scan"
        )
        logs_dir = scan_logs_dir(task_root, scan_name)
        self.info = {
            "task_root": task_root,
            "task_result_dir": str(self.config.get("task_result_dir") or os.getcwd()),
            "scan_name": scan_name,
            "sampler_name": "SamplingVirtial",
            "run_id": str(self.config.get("run_id") or uuid.uuid4()),
            "task_yaml": self.config.get("task_yaml"),
            "logs_dir": logs_dir,
            "jarvis_log": component_log_path(logs_dir, "core"),
        }

    def prepare_resume(self, *, resume: bool = False, fresh: bool = False) -> None:
        """Preload checkpoint state before sampler bring-up (CLI ``--resume``)."""
        return self._get_resume().prepare_resume(resume=resume, fresh=fresh)

    def init_sampler_from_config(self) -> SamplingVirtial:
        method = sampling_method(self.config)
        if method:
            sampler = Distributor.set_method(method)
        else:
            sampler = SamplingVirtial()
        sampler.set_config(self.config)
        calc_modules = list(
            (self.config.get("Calculators") or {}).get("Modules") or []
        )
        opera_modules = list((self.config.get("Operas") or {}).get("Modules") or [])
        sampling_block = self.config.get("Sampling") or {}
        likelihood_exprs = list(sampling_block.get("LogLikelihood") or [])
        from jarvishep2.Module.nuisance import extract_nuisance_config

        include_nuisance = extract_nuisance_config(self.config) is not None
        sampler.set_execution_plan_template(
            calculator_modules=calc_modules,
            opera_modules=opera_modules,
            include_likelihood=bool(likelihood_exprs),
            include_nuisance=include_nuisance,
        )
        self.set_sampler(sampler)
        self.info["sampler_name"] = str(getattr(sampler, "method", type(sampler).__name__))
        # Safety net: prepare_resume may have looked up the wrong path before
        # Method was known.  Re-resolve under the real sampler name and import.
        if (
            self._resume_policy == "resume"
            and self._resume_checkpoint_payload is None
            and getattr(self, "_resume_checkpoint_preloaded", False)
        ):
            path = self.checkpoint_file()
            if os.path.isfile(path):
                try:
                    self._resume_checkpoint_payload = load_checkpoint(path)
                except ValueError as exc:
                    self._logger.warning(
                        "Late checkpoint load rejected (starting fresh) → %s",
                        exc,
                    )
                    self._discard_resume_checkpoint(
                        reason=f"late checkpoint load rejected: {exc}",
                        level="warning",
                    )
                else:
                    if self._import_sampler_checkpoint_state(sampler):
                        self._logger.info(
                            "Loaded resume checkpoint after sampler init → %s", path
                        )
        return sampler

    def bootstrap_distributed_runtime(self) -> None:
        """Bring up Redis, Workers, and Archiver for a distributed run."""
        return self._get_runtime().bootstrap_distributed_runtime()

    @staticmethod
    def _load_check_module_rows(csv_path: str) -> tuple[list[dict[str, str]], list[str]]:
        return _ScanDriver._load_check_module_rows(csv_path)

    def _build_check_module_samples(self) -> list[Sample]:
        """Build smoke samples: prefer fixed CSV, else draw from the sampler.

        Resolution order for the CSV path:

        1. ``EnvReqs.V2.check_modules.data`` (task or default environment YAML)
        2. Built-in default ``&J/data/check_modules_points.csv``

        When the resolved file is missing, fall back to
        ``EnvReqs.V2.check_modules.n_samples`` (default 10) points drawn from
        ``Sampling.Method`` (V1-like assembly-line smoke).
        """
        return self._get_scan()._build_check_module_samples()

    def _sample_check_module_points_from_sampler(self, n_samples: int) -> list[Sample]:
        """Draw up to *n_samples* smoke points from the active sampler (V1-like)."""
        return self._get_scan()._sample_check_module_points_from_sampler(n_samples)

    def run_distributed_scan(self, *, timeout: float = 3600.0) -> int:
        """Drive a stateless sampler through propose → Redis → Archiver."""
        return self._get_scan().run_distributed_scan(timeout=timeout)

    def run_adaptive_scan(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Drive feedback samplers; timeout is a per-generation wait, not a run budget."""
        return self._get_scan().run_adaptive_scan(generation_timeout=generation_timeout, timeout=timeout)

    def _wait_for_archive_caught_up(self, *, timeout: float = 3600.0) -> None:
        """Block until Archiver has persisted every Redis-completed sample.

        Process-mode :class:`ArchiverProcess.drain` is a parent-side no-op; the
        child keeps consuming ``hep:archive``. Plot/CSV export must wait until
        ``records_written >= ok + failed`` (and workers are idle), otherwise
        ``DATABASE/samples.csv`` misses the last archive batches.
        """
        return self._get_scan()._wait_for_archive_caught_up(timeout=timeout)

    def run_check_modules(
        self,
        *,
        timeout: float | None = None,
        verify_golden: Mapping[str, Any] | None = None,
    ) -> int:
        """Run the calculator/opera smoke path (CSV fixed points or N sampler draws).

        ``timeout`` is a wait deadline (seconds) for samples to archive, not a
        fixed sleep. Default comes from ``EnvReqs.V2.check_modules.timeout_sec``
        (120s); CLI ``Jarvis check --timeout SEC`` can override.
        """
        return self._get_scan().run_check_modules(timeout=timeout, verify_golden=verify_golden)

    def _apply_check_modules_runtime_policy(self) -> None:
        """Force smoke-friendly layout: 1 worker, SAMPLE/test, no tar pack.

        Applied before bootstrap so Factory/Archiver/Workers all see the policy.
        """
        return self._get_scan()._apply_check_modules_runtime_policy()

    def _resolve_sample_root(self) -> str:
        """Return SAMPLE root; ``Jarvis check`` uses ``SAMPLE/test`` (no tar pack)."""
        return self._get_scan()._resolve_sample_root()

    def _export_workflow_flowchart(self) -> None:
        """Write flowchart.json (+ optional PNG) under images/ (D12.3 / V1 paths).

        Semantic graph is V1-compatible ``jarvisplot.scene/v1`` so stock JarvisPLOT
        can render with no V2 adapter. ``--skip-draw-flowchart`` skips both JSON and PNG.
        """
        if bool(self.config.get("skip_draw_flowchart")) or bool(
            self.info.get("skip_draw_flowchart")
        ):
            return
        try:
            from jarvishep2.flowchart import (
                build_flowchart_scene_from_config,
                export_flowchart_semantics,
                render_flowchart_png,
            )
        except Exception as exc:
            self._logger.warning("flowchart module unavailable -> %s", exc)
            return
        task_result_dir = str(
            self.info.get("task_result_dir") or self.config.get("task_result_dir") or os.getcwd()
        )
        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        project_root = str(
            self.config.get("project_root")
            or self.config.get("task_root")
            or self.info.get("project_root")
            or ""
        ).strip()
        if not project_root:
            from jarvishep2.base import infer_project_root_from_task_result_dir

            project_root = infer_project_root_from_task_result_dir(task_result_dir)
        # Project layout: <project>/images/<scan>/flowchart.{json,png}
        from jarvishep2.base import project_images_dir

        images_dir = project_images_dir(project_root=project_root, scan_name=scan_name)
        os.makedirs(images_dir, exist_ok=True)
        json_path = os.path.join(images_dir, "flowchart.json")
        png_path = os.path.join(images_dir, "flowchart.png")
        try:
            # Stamp project/scan names into config for scene metadata.
            cfg = dict(self.config)
            if not cfg.get("project_name") and self.info.get("project_name"):
                cfg["project_name"] = self.info.get("project_name")
            if not cfg.get("scan_name") and self.info.get("scan_name"):
                cfg["scan_name"] = self.info.get("scan_name")
            scene = build_flowchart_scene_from_config(cfg)
            export_flowchart_semantics(scene, json_path)
            self.info["flowchart_semantic_path"] = json_path
            self._logger.info(
                "Export workflow semantic graph into %s", json_path
            )
            rendered = render_flowchart_png(scene, png_path)
            if rendered:
                self.info["flowchart_path"] = rendered
                self._logger.info("workflow flowchart rendered → %s", rendered)
            else:
                self._logger.info(
                    "JarvisPLOT not installed; skipped flowchart.png "
                    "(scene JSON still written)"
                )
        except Exception as exc:
            self._logger.warning("workflow flowchart export failed -> %s", exc)

    def _finalize_nested_result_csv(self) -> None:
        """Write clean nested-sampler CSV (dead points only, not all logL calls).

        Dynesty/MultiNest keep a much smaller set of nested samples than the
        full Worker likelihood call stream. That table lives at
        ``DATABASE/dynesty_result.csv`` or ``DATABASE/multinest_result.csv`` and
        is the scientific artifact for evidence/weights — distinct from
        ``DATABASE/samples.csv`` (every archived evaluation).
        """
        sampler = self.sampler
        if sampler is None:
            return
        save = getattr(sampler, "save_dynesty_results_to_csv", None)
        if not callable(save):
            return
        try:
            path = save()
        except Exception as exc:
            self._logger.warning("nested result CSV finalize failed -> %s", exc)
            return
        if not path:
            return
        path_text = str(path)
        method = str(getattr(sampler, "method", "") or "")
        if "multinest" in path_text.lower() or method.lower() == "multinest":
            self.info["multinest_result_csv"] = path_text
            self.info.pop("dynesty_result_csv", None)
        else:
            self.info["dynesty_result_csv"] = path_text
        self._logger.info(
            "nested clean result CSV ready → %s (%s)",
            path_text,
            method or "nested",
        )

    def _emit_plot_scenes(self) -> None:
        """Emit scan/levelset jplot YAML under ``<project>/images/<scan>/``."""
        if bool(self.config.get("skip_emit_plot_scenes")):
            return
        try:
            from jarvishep2.plot_scene import emit_plot_scenes_from_run
        except Exception:
            return
        task_result_dir = str(
            self.info.get("task_result_dir") or self.config.get("task_result_dir") or os.getcwd()
        )
        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        project_root = str(
            self.config.get("project_root")
            or self.config.get("task_root")
            or self.info.get("project_root")
            or ""
        ).strip()
        if not project_root:
            from jarvishep2.base import infer_project_root_from_task_result_dir

            project_root = infer_project_root_from_task_result_dir(task_result_dir)
        # Auto-render when JarvisPLOT is available unless explicitly disabled.
        auto_render = not bool(
            self.config.get("skip_render_plots") or self.info.get("skip_render_plots")
        )
        try:
            written = emit_plot_scenes_from_run(
                task_result_dir,
                scan_name=scan_name,
                project_root=project_root,
                auto_render=auto_render,
                config=self.config,
            )
            if written:
                self.info["plot_scenes"] = written
                self._logger.info(
                    "plot scenes emitted → %s",
                    ", ".join(f"{key}={path}" for key, path in written.items()),
                )
                for key, label in (
                    ("jplot_levelset", "jplot_levelset_png"),
                    ("jplot_dynesty_runplot", "jplot_dynesty_runplot_png"),
                    ("jplot_multinest_runplot", "jplot_multinest_runplot_png"),
                ):
                    jplot = written.get(key)
                    if jplot and label not in written:
                        self._logger.info(
                            "jplot YAML ready (render with: Jarvis plot %s)", jplot
                        )
        except Exception as exc:
            self._logger.warning("plot scene emit failed -> %s", exc)

    def _capture_run_outcome(
        self,
        *,
        submitted: int,
        interrupted: bool = False,
        error: str | None = None,
        error_type: str | None = None,
    ) -> RunOutcome:
        counters = self._live_sample_counters()
        archived = 0
        try:
            archived = int(self._archiver_records_written())
        except Exception:
            archived = 0
        run_id = str(self.info.get("run_id") or "")
        ok = int(counters.get("ok") or 0)
        failed = int(counters.get("failed") or 0)
        queued = int(counters.get("queued") or 0)
        running = int(counters.get("running") or 0)
        redis_total = ok + failed + queued + running
        extras: dict[str, Any] = {
            "queued": queued,
            "running": running,
            "archive_q": int(counters.get("archive_q") or 0),
        }
        # Nested samplers return *ncall* from run_adaptive; Redis SAMPLE_STATS
        # count real Worker tasks. Prefer Redis for submitted/completed and
        # expose ncall separately so CLI is not misleading (D13.8 hygiene).
        ncall: int | None = None
        reported = int(submitted or 0)
        if reported > 0 and redis_total > 0 and reported > max(1, redis_total) * 2:
            ncall = reported
            submitted_out = max(redis_total, ok + failed)
        else:
            submitted_out = reported if reported > 0 else redis_total
            # Still surface dynesty ncall when available on the sampler summary.
            try:
                summary = (
                    self.sampler.summary()
                    if self.sampler is not None and hasattr(self.sampler, "summary")
                    else None
                )
                if isinstance(summary, Mapping) and summary.get("ncall") is not None:
                    ncall = int(summary["ncall"])
            except Exception:
                pass
        if ncall is not None:
            extras["ncall"] = int(ncall)
        return RunOutcome.from_counters(
            submitted=int(submitted_out or 0),
            completed=ok,
            failed=failed,
            archived=archived,
            run_id=run_id,
            interrupted=interrupted,
            error=error,
            error_type=error_type,
            extras=extras,
        )

    def run(
        self,
        *,
        resume: bool = False,
        check_modules: bool = False,
        verify_golden: Mapping[str, Any] | None = None,
        write_run_summary: bool = True,
        check_timeout: float | None = None,
    ) -> RunOutcome:
        """Execute a distributed scan; return a truthful :class:`RunOutcome` (D11.1)."""
        self.prepare_resume(resume=resume, fresh=False)
        self._install_control_signal_handlers()
        submitted = 0
        outcome: RunOutcome | None = None
        is_check = bool(check_modules)
        if is_check:
            self._apply_check_modules_runtime_policy()
        try:
            self.bootstrap_distributed_runtime()
            method = sampling_method(self.config)
            if is_check:
                submitted = self.run_check_modules(
                    timeout=check_timeout,
                    verify_golden=verify_golden,
                )
            elif method in STATELESS_METHODS:
                submitted = self.run_distributed_scan()
            elif method:
                # Stateful / feedback-driven methods (e.g. AdaptiveBridson, MCMC, nested).
                submitted = self.run_adaptive_scan()
            else:
                raise NotImplementedError(
                    "Unsupported task: configure Sampling.Method "
                    "(Bridson|AdaptiveBridson|Dynesty|…), or run Jarvis check TASK.yaml."
                )
            outcome = self._capture_run_outcome(submitted=submitted)
            if not is_check and not self._interrupt_requested:
                # Every normal scan, including stateless runs, must finish only
                # after its Archiver has persisted the final partial batch.  In
                # process mode the child owns that flush; waiting here keeps the
                # outcome, CSV export, and plots consistent with DATABASE.
                self._logger.info("final archive verification before normal exit")
                self._wait_for_archive_caught_up(timeout=120.0)
                outcome = self._capture_run_outcome(submitted=submitted)

            # D19.2: partial_failure is common in physics (failed points are
            # useful science). Emit plots/CSV from surviving samples; only a
            # clean success, or a mixed outcome, should produce scenes. Full
            # failure / interrupt still skip with an explicit reason.
            if not self._interrupt_requested:
                if is_check:
                    self._logger.warning(
                        "Jarvis check smoke complete (submitted=%d); "
                        "skipping nested result CSV export and plot scenes",
                        submitted,
                    )
                elif outcome.ok or outcome.status == "partial_failure":
                    if outcome.status == "partial_failure":
                        total = max(1, outcome.completed + outcome.failed)
                        fail_pct = 100.0 * outcome.failed / total
                        self._logger.warning(
                            "Scan finished with partial_failure "
                            "(completed=%d failed=%d, %.1f%% failed); "
                            "emitting nested CSV / plot scenes from surviving samples",
                            outcome.completed,
                            outcome.failed,
                            fail_pct,
                        )
                    self._finalize_nested_result_csv()
                    self._emit_plot_scenes()
                else:
                    self._logger.info(
                        "Skipping nested CSV export and plot scenes: "
                        "run status=%s (completed=%d failed=%d)",
                        outcome.status,
                        outcome.completed,
                        outcome.failed,
                    )
            return outcome
        except KeyboardInterrupt:
            self._interrupt_requested = True
            try:
                self._logger.warning(
                    "KeyboardInterrupt / stop signal — shutting down Workers, "
                    "Archiver, and managed Redis (prefer Ctrl+C over Ctrl+Z)"
                )
            except Exception:
                pass
            # D8.3 human path: force a resume checkpoint *before* teardown
            # (agent run_state / stop-ack remain parked with D8).
            self._save_interrupt_checkpoint()
            outcome = self._capture_run_outcome(
                submitted=submitted,
                interrupted=True,
                error="interrupted",
                error_type="KeyboardInterrupt",
            )
            return outcome
        finally:
            try:
                self.shutdown(
                    wait=True,
                    write_run_summary=write_run_summary and not self._interrupt_requested,
                )
            finally:
                if outcome is not None and self._final_archived_records is not None:
                    outcome.archived = max(
                        int(outcome.archived),
                        int(self._final_archived_records),
                    )
                self._restore_control_signal_handlers()

    def check_modules(
        self,
        *,
        verify_golden: Mapping[str, Any] | None = None,
        timeout: float | None = None,
    ) -> RunOutcome:
        """CLI entry for ``Jarvis check <task>.yaml``."""
        return self.run(
            check_modules=True,
            verify_golden=verify_golden,
            check_timeout=timeout,
        )

    def is_redis_runtime(self) -> bool:
        """Return True when the distributed Redis path should be used."""
        return str(self.runtime.get("mode", "redis")).strip().lower() == "redis"

    def _resolve_checkpoint_sampler_name(self) -> str:
        """Directory name under ``checkpoints/<scan>/``.

        Must match what a running scan writes (``sampler.method``, e.g. ``Dynesty``).
        ``prepare_resume`` runs **before** ``init_sampler_from_config``, so
        ``info["sampler_name"]`` is still the placeholder ``SamplingVirtial`` —
        fall through to YAML ``Sampling.Method`` instead of looking for a
        non-existent ``…/SamplingVirtial/state.pkl`` (which silently skipped
        nested engine restore and made Dynesty restart from scratch).
        """
        return self._get_resume()._resolve_checkpoint_sampler_name()

    def checkpoint_file(self) -> str:
        return self._get_resume().checkpoint_file()

    @staticmethod
    def prompt_resume_from_checkpoint(
        checkpoint_file: str,
        *,
        timeout_seconds: float = 30.0,
    ) -> bool:
        return _ResumeService.prompt_resume_from_checkpoint(checkpoint_file, timeout_seconds=timeout_seconds)

    def _preload_resume_checkpoint(
        self,
        *,
        resume: bool = False,
        fresh: bool = False,
    ) -> None:
        return self._get_resume()._preload_resume_checkpoint(resume=resume, fresh=fresh)

    def _load_persisted_database_state(self, database_dir: str) -> None:
        """Refresh index prefix, UUID set, and completed/failed counts from DATABASE.

        Index prefix drives resume fast-forward; outcome counts seed Redis
        SAMPLE_STATS so a resume does not re-label durable Failed rows as
        successful (D21.13).
        """
        return self._get_resume()._load_persisted_database_state(database_dir)

    def apply_resume_checkpoint(self, worker_config: Mapping[str, Any] | None = None) -> int:
        return self._get_resume().apply_resume_checkpoint(worker_config)

    def save_runtime_checkpoint(self, *, reason: str = "") -> str | None:
        return self._get_resume().save_runtime_checkpoint(reason=reason)

    def _save_interrupt_checkpoint(self) -> str | None:
        """Force a resume checkpoint on SIGINT/SIGTERM (D8.3 human remainder).

        Agent-facing pieces (run_state interrupted, stop-ack) stay parked with D8.
        Failures are logged and never re-raised — teardown must still run.
        """
        return self._get_resume()._save_interrupt_checkpoint()

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        # EnvReqs.V2.redis (optional) overlays INTERNAL_REDIS_CONFIG defaults.
        return self._get_runtime().init_redis(client=client)

    def _ensure_managed_redis(self, redis_config: Mapping[str, Any]) -> None:
        """Start redis-server as Jarvis-Redis:<scan> when the port is free."""
        return self._get_runtime()._ensure_managed_redis(redis_config)

    def _set_runtime_redis_config(self, redis_config: Mapping[str, Any]) -> None:
        """Keep the loaded task in sync with an automatically reassigned port."""
        return self._get_runtime()._set_runtime_redis_config(redis_config)

    def _control_lock_owner_id(self) -> str:
        return self._get_runtime()._control_lock_owner_id()

    def _calculator_module_names(self) -> list[str]:
        return self._get_runtime()._calculator_module_names()

    def _claim_redis_control_lock(self) -> None:
        """Refuse to start if another Jarvis still holds the Redis control lease."""
        return self._get_runtime()._claim_redis_control_lock()

    def _control_lease_loop(self, stop: threading.Event, owner: str) -> None:
        return self._get_runtime()._control_lease_loop(stop, owner)

    def _start_control_lease_refresh(self) -> None:
        return self._get_runtime()._start_control_lease_refresh()

    def _stop_control_lease_refresh(self) -> None:
        return self._get_runtime()._stop_control_lease_refresh()

    def _reset_redis_for_fresh_run(self) -> None:
        """Drop ephemeral queues/stats/calc pools before Workers are spawned."""
        return self._get_runtime()._reset_redis_for_fresh_run()

    def _release_redis_control_lock(self) -> None:
        return self._get_runtime()._release_redis_control_lock()

    def _publish_runtime_metadata(self) -> None:
        return self._get_runtime()._publish_runtime_metadata()

    def init_command_parser(self) -> CommandParser:
        """Run Phase-1 static command resolution for the loaded task config."""
        return self._get_runtime().init_command_parser()

    def prepare_libraries(self, *, parser: CommandParser | None = None) -> dict[str, dict[str, Any]]:
        """Install or reuse LibDeps before Redis, Workers, or an Archiver start."""
        return self._get_runtime().prepare_libraries(parser=parser)

    def _apply_command_parser_to_worker_config(self, worker_config: dict[str, Any]) -> dict[str, Any]:
        return self._get_runtime()._apply_command_parser_to_worker_config(worker_config)

    def build_worker_config(self, **overrides: Any) -> dict[str, Any]:
        """Build a picklable Worker blueprint with Phase-1 command resolution applied."""
        return self._get_runtime().build_worker_config(**overrides)

    def init_factory(self, worker_config: Mapping[str, Any] | None = None) -> TaskFactory | None:
        return self._get_runtime().init_factory(worker_config)

    @staticmethod
    def _format_convert_report(
        *,
        title: str,
        rows: Sequence[tuple[str, Any]],
    ) -> str:
        """Multi-line convert status with space-aligned key/value fields."""
        text = str(title or "").strip()
        label_w = max((len(str(k)) for k, _ in rows), default=0)
        lines = [f"{text}"]
        for key, value in rows:
            lines.append(f"  {str(key):<{label_w}}  →  {value}")
        return "\n".join(lines)

    def convert(self) -> list[dict[str, Any]]:
        """Refresh CSV snapshots when their DATABASE HDF5 source has changed."""
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        results = convert_database_dir(database_dir)
        logger = getattr(self, "_logger", None)
        if logger is None:
            try:
                logger = get_jarvis_logger()
            except Exception:
                logger = None

        def _emit(msg: str, *, error: bool = False) -> None:
            if logger is not None:
                if error:
                    logger.error(msg)
                else:
                    logger.warning(msg)
            else:
                print(msg, file=sys.stderr if error else sys.stdout)

        if not results:
            _emit(
                self._format_convert_report(
                    title="HDF5 → CSV  (nothing to convert)",
                    rows=[
                        ("database", database_dir),
                        ("reason", "no *.hdf5 files found"),
                    ],
                ),
                error=True,
            )
            return results

        for item in results:
            status = str(item.get("status") or "")
            csv_path = str(item.get("csv") or "")
            hdf5_path = str(item.get("hdf5") or "")
            rows_n = int(item.get("rows") or 0)
            if status == "converted":
                msg = self._format_convert_report(
                    title="HDF5 → CSV  (converted)",
                    rows=[
                        ("rows", f"{rows_n:,}"),
                        ("from", hdf5_path),
                        ("to", csv_path),
                    ],
                )
            elif status == "skipped_unchanged":
                msg = self._format_convert_report(
                    title="HDF5 → CSV  (skipped)",
                    rows=[
                        ("reason", "HDF5 is unchanged since the last export"),
                        ("csv", csv_path),
                        ("hdf5", hdf5_path),
                    ],
                )
            elif status == "empty":
                msg = self._format_convert_report(
                    title="HDF5 → CSV  (empty)",
                    rows=[
                        ("reason", "HDF5 has no records"),
                        ("hdf5", hdf5_path),
                    ],
                )
            else:
                msg = self._format_convert_report(
                    title="HDF5 → CSV  (missing)",
                    rows=[
                        ("reason", "HDF5 missing or unreadable"),
                        ("hdf5", hdf5_path),
                    ],
                )
            _emit(msg, error=status in {"empty", "missing"})
        return results

    def init_archiver(self, db_path: str | None = None) -> SimpleArchiver | ArchiverProcess:
        return self._get_runtime().init_archiver(db_path)

    def _restore_archiver_persistence(self) -> None:
        """Compatibility hook; DATABASE seeding now happens in SimpleArchiver.

        Never promote checkpoint/Redis acknowledgements into completion truth:
        only UUIDs actually read from HDF5 may seed the deduplication set.
        """
        return self._get_runtime()._restore_archiver_persistence()

    def _archiver_records_written(self) -> int:
        return self._get_runtime()._archiver_records_written()

    def _discard_resume_checkpoint(
        self,
        *,
        reason: str,
        level: str = "error",
    ) -> None:
        """Drop a rejected resume snapshot and force a fresh sampling start."""
        return self._get_resume()._discard_resume_checkpoint(reason=reason, level=level)

    def _import_sampler_checkpoint_state(self, sampler: Any | None = None) -> bool:
        """Import sampler_state from the preloaded checkpoint, or start fresh.

        Shared by ``set_sampler``, late Method resolution, bootstrap, and
        ``apply_resume_checkpoint``. On Bounds/shape mismatch (``ValueError``)
        the stale checkpoint is discarded so a finished or drifted snapshot
        cannot silently yield ``submitted=0`` under a new card.

        Returns True when import succeeded.
        """
        return self._get_resume()._import_sampler_checkpoint_state(sampler)

    def set_sampler(self, sampler: Any) -> None:
        self.sampler = sampler
        if self._resume_policy == "auto" and not getattr(self, "_resume_checkpoint_preloaded", False):
            self._preload_resume_checkpoint()
            self._resume_checkpoint_preloaded = True
        if self.redis is not None and hasattr(sampler, "set_redis"):
            sampler.set_redis(self.redis)
        self._import_sampler_checkpoint_state(sampler)

    def start_runtime_checkpoint(self) -> None:
        """Enable the configured checkpoint heartbeat once sampler and archiver exist."""
        return self._get_resume().start_runtime_checkpoint()

    def submit_samples(self, samples: Sequence[Sample]) -> None:
        return self._get_scan().submit_samples(samples)

    def _live_sample_counters(self) -> dict[str, int]:
        """Snapshot Redis sample/queue counters for completion progress lines."""
        return self._get_scan()._live_sample_counters()

    def wait_for_results(
        self,
        expected: int,
        *,
        timeout: float = 30.0,
        poll_interval: float = 0.1,
        progress_total: int | None = None,
        progress_base: int | None = None,
        require_worker_completion: bool = False,
    ) -> None:
        """Block until this batch is complete and Archiver has written its rows.

        When ``progress_total`` is set, emit V1-style ‰ completion heartbeats
        using archived count relative to ``progress_base`` (default 0).
        """
        return self._get_scan().wait_for_results(expected, timeout=timeout, poll_interval=poll_interval, progress_total=progress_total, progress_base=progress_base, require_worker_completion=require_worker_completion)

    def get_monitor_snapshot(self) -> dict[str, Any]:
        if self.factory is None:
            return {}
        return self.factory.get_monitor_snapshot()

    def monitor_once(self) -> str:
        if self.factory is None:
            return ""
        view = SnapshotReader(self.factory).read()
        return format_monitor_view(view)

    def write_run_summary(self, output_dir: str | None = None) -> dict[str, str]:
        if self.factory is None:
            raise RuntimeError("factory is not configured")
        task_result_dir = str(
            output_dir
            or self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        metrics = dict(self.factory.get_run_metrics())
        started = metrics.pop("run_started_at", None)
        scan = self.config.get("Scan") or {}
        scan_name = scan.get("name") if isinstance(scan, Mapping) else None
        sampler_name = type(self.sampler).__name__ if self.sampler is not None else None
        summary = build_run_summary(
            factory_metrics=metrics,
            project_name=str(self.config.get("project_name") or scan_name or ""),
            sampler_name=sampler_name,
            run_id=str(self.info.get("run_id") or metrics.get("run_id") or "jarvis2-run"),
            start_epoch=float(started) if started is not None else None,
            configured_workers=int(self.runtime.get("workers", 0) or len(self.factory.workers)),
        )
        # D13.6: keep frozen run_summary schema pure; mirror sampler.summary()
        # into DATABASE/sampler_summary.json when present (additive contract).
        if self.sampler is not None and hasattr(self.sampler, "summary"):
            try:
                sampler_summary = self.sampler.summary()
                if isinstance(sampler_summary, Mapping) and sampler_summary:
                    from jarvishep2.Sampling.diagnostics_export import (
                        write_sampler_summary_json,
                    )

                    diag = write_sampler_summary_json(task_result_dir, sampler_summary)
                    self.info["sampler_summary"] = diag
                    # Nested samplers: surface ncall next to SAMPLE counts in the
                    # performance block (not the same as Redis sample submits).
                    if sampler_summary.get("ncall") is not None:
                        summary["ncall"] = int(sampler_summary["ncall"])
                    if sampler_summary.get("niter") is not None:
                        summary["niter"] = int(sampler_summary["niter"])
            except Exception as exc:
                self._logger.warning("sampler_summary export failed -> %s", exc)
        paths = RunSummaryRenderer().write_outputs(summary, task_result_dir)
        try:
            # End-of-scan performance block in the main Jarvis-HEP log / console.
            self._logger.warning("\n" + format_scan_performance_log(summary).rstrip())
            self._logger.info(
                "run_summary written -> json=%s csv=%s txt=%s",
                paths.get("json"),
                paths.get("csv"),
                paths.get("txt"),
            )
        except Exception:
            pass
        return paths

    def _init_sample_buckets(self) -> None:
        """Register Redis SAMPLE bucket meta (numbering + active/completed state)."""
        return self._get_runtime()._init_sample_buckets()

    def _finalize_sample_buckets(self) -> None:
        """Seal the open bucket after all samples finish; Archiver packs when idle."""
        return self._get_scan()._finalize_sample_buckets()

    def _install_control_signal_handlers(self) -> None:
        """SIGINT/SIGTERM → clean shutdown; refuse Ctrl+Z suspend."""
        return self._get_runtime()._install_control_signal_handlers()

    def _restore_control_signal_handlers(self) -> None:
        return self._get_runtime()._restore_control_signal_handlers()

    def _stop_managed_redis(self) -> None:
        return self._get_runtime()._stop_managed_redis()

    def _atexit_cleanup(self) -> None:
        """Best-effort stop of managed Redis if the process dies unexpectedly."""
        return self._get_runtime()._atexit_cleanup()

    def shutdown(self, *, wait: bool = True, write_run_summary: bool = False) -> None:
        """Stop Archiver, Workers, release the control lock, and managed Redis.

        Idempotent: safe to call from ``finally``, signal paths, and atexit.
        Managed Redis is always stopped last in an inner ``finally`` so a failure
        earlier in the chain cannot leave ``Jarvis-Redis:*`` orphaned.
        """
        return self._get_runtime().shutdown(wait=wait, write_run_summary=write_run_summary)


__all__ = ["Jarvis2Core"]
