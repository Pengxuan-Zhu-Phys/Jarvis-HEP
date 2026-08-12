#!/usr/bin/env python3
"""Jarvis control-process orchestrator for the distributed runtime."""

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
from jarvishep2.factory import TaskFactory
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
from jarvishep2.testing.check_modules import (
    build_check_module_samples,
    verify_check_modules_golden,
)
from jarvishep2.run_outcome import RunOutcome
from jarvishep2.worker_config import build_command_parser, build_worker_config
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


class Jarvis2Core:
    """Minimal distributed run orchestrator for the D1.1 opera-only MVP."""

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
        # Last-resort cleanup if the process exits without an orderly finally.
        atexit.register(self._atexit_cleanup)

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
        self._preload_resume_checkpoint(resume=resume, fresh=fresh)
        self._resume_checkpoint_preloaded = True

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
        if not self.is_redis_runtime():
            raise RuntimeError("distributed runtime requires the internal Redis runtime")
        from jarvishep2.process_cleanup import ensure_scan_name_available

        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        ensure_scan_name_available(scan_name, cleanup_stale=self._resume_policy == "resume")
        db_path = os.path.join(
            str(self.info.get("task_result_dir") or os.getcwd()),
            "DATABASE",
            "samples.hdf5",
        )
        if self._resume_policy == "resume":
            try:
                # Control must not open samples.hdf5 for append (write lock would
                # block Archiver). Only clear SWMR flags / reclaim stale Archivers.
                # Logs a single summary line internally when recovery runs.
                prepare_hdf5_database_for_writer(
                    db_path,
                    logger=self._logger,
                    scan_name=scan_name,
                    probe_append=False,
                )
            except Exception as exc:
                # Non-fatal at this early stage — Archiver will retry recover on open.
                self._logger.warning(
                    "Early HDF5 prepare failed (will retry at Archiver start): %s",
                    exc,
                )
        if self._resume_policy == "resume":
            self._load_persisted_database_state(os.path.dirname(db_path))
        self.init_logger()
        # Preflight remains before command setup, Redis, Workers, and Archiver,
        # while init_logger above ensures all outcomes reach core.log.
        self.check_environment_requirements()
        # LibDeps is globally shared.  Resolve its V1 tokens with no registered
        # executable side effects, install/reuse it once, then register tools.
        root_path = self._environment_summary.get("ROOT path")
        library_parser = build_command_parser(
            self.config,
            root_path=str(root_path) if root_path else None,
            register_executables=False,
        )
        self.prepare_libraries(parser=library_parser)
        self.init_command_parser()
        self.init_redis()
        self._claim_redis_control_lock()
        self._publish_runtime_metadata()
        if self._resume_policy != "resume":
            self._reset_redis_for_fresh_run()
        self.init_sampler_from_config()
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_uuids"):
            self.sampler.set_persisted_uuids(self._persisted_uuids)
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_prefix"):
            self.sampler.set_persisted_prefix(self._persisted_index_prefix)
        if (
            self._resume_policy == "resume"
            and self.sampler is not None
            and hasattr(self.sampler, "advance_to_persisted_prefix")
        ):
            self.sampler.advance_to_persisted_prefix(self._persisted_index_prefix)
        self._init_sample_buckets()
        self.init_factory()
        # Immediately before Archiver spawn: reclaim orphan Archivers only —
        # still no append-open in Control (probe_append=False).
        if self._resume_policy == "resume" and os.path.isfile(db_path):
            try:
                # Only re-runs work if the file is still dirty; one summary log max.
                prepare_hdf5_database_for_writer(
                    db_path,
                    logger=self._logger,
                    scan_name=scan_name,
                    probe_append=False,
                )
            except Exception as exc:
                self._logger.warning(
                    "Pre-Archiver HDF5 prepare failed (Archiver will retry): %s",
                    exc,
                )
        self.init_archiver(db_path)
        self._load_persisted_database_state(os.path.dirname(db_path))
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_uuids"):
            self.sampler.set_persisted_uuids(self._persisted_uuids)
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_prefix"):
            self.sampler.set_persisted_prefix(self._persisted_index_prefix)
        if (
            self._resume_policy == "resume"
            and self.sampler is not None
            and hasattr(self.sampler, "advance_to_persisted_prefix")
        ):
            self.sampler.advance_to_persisted_prefix(self._persisted_index_prefix)
        if self._resume_policy == "resume":
            self._logger.warning(
                "Resume DATABASE reconciliation found prefix=%d, "
                "completed=%d failed=%d, %d legacy UUID(s)",
                self._persisted_index_prefix,
                max(0, self._persisted_records_count - self._persisted_failed_count),
                self._persisted_failed_count,
                len(self._persisted_uuids),
            )
            # Nested / feedback engines: ensure runtime state was applied once Redis
            # and DATABASE facts exist (idempotent if set_sampler already imported).
            self._import_sampler_checkpoint_state(self.sampler)
            # Always arm nested engine resume under --resume, even if state.pkl
            # is missing: nested_engine.pkl alone is enough to continue.
            if self.sampler is not None and hasattr(self.sampler, "arm_engine_resume"):
                self.sampler.arm_engine_resume()
            ckpt = self.checkpoint_file()
            self._logger.warning(
                "Resume checkpoint path → %s (payload=%s)",
                ckpt,
                "yes" if self._resume_checkpoint_payload is not None else "MISSING",
            )
        self.start_runtime_checkpoint()
        self._export_workflow_flowchart()

    @staticmethod
    def _load_check_module_rows(csv_path: str) -> tuple[list[dict[str, str]], list[str]]:
        with open(csv_path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            return list(reader), fieldnames

    def _build_check_module_samples(self) -> list[Sample]:
        """Build smoke samples: prefer fixed CSV, else draw from the sampler.

        Resolution order for the CSV path:

        1. ``EnvReqs.V2.check_modules.data`` (task or default environment YAML)
        2. Built-in default ``&J/data/check_modules_points.csv``

        When the resolved file is missing, fall back to
        ``EnvReqs.V2.check_modules.n_samples`` (default 10) points drawn from
        ``Sampling.Method`` (V1-like assembly-line smoke).
        """
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        csv_path, raw_spec = resolve_check_modules_csv(self.config)
        if csv_path is not None:
            self._logger.warning(
                "Jarvis check: using fixed points from CSV -> %s (configured as %r)",
                csv_path,
                raw_spec,
            )
            rows, fieldnames = self._load_check_module_rows(csv_path)
            samples = build_check_module_samples(
                sampler=self.sampler,
                config=self.config,
                rows=rows,
                csv_fieldnames=fieldnames,
            )
            self._logger.warning(
                "Jarvis check: loaded %d fixed point(s) from CSV", len(samples)
            )
            return samples

        n_samples = check_modules_n_samples(self.config)
        method = sampling_method(self.config) or type(self.sampler).__name__
        self._logger.warning(
            "Jarvis check: CSV not found (configured as %r); "
            "drawing %d smoke point(s) from Sampling.Method=%r",
            raw_spec or "(none)",
            n_samples,
            method,
        )
        samples = self._sample_check_module_points_from_sampler(n_samples)
        self._logger.warning(
            "Jarvis check: prepared %d sampler-drawn smoke point(s)", len(samples)
        )
        return samples

    def _sample_check_module_points_from_sampler(self, n_samples: int) -> list[Sample]:
        """Draw up to *n_samples* smoke points from the active sampler (V1-like)."""
        n_samples = max(1, int(n_samples))
        samples: list[Sample] = []

        # Prefer propose_next for Bridson/Random/Grid/CSV when available.
        propose = getattr(self.sampler, "propose_next", None)
        if callable(propose):
            for _ in range(n_samples):
                try:
                    sample = propose()
                except Exception as exc:
                    self._logger.warning(
                        "Jarvis check: propose_next failed after %d sample(s) -> %s; "
                        "falling back to unit-cube draws",
                        len(samples),
                        exc,
                    )
                    break
                if sample is None:
                    break
                samples.append(sample)

        if len(samples) >= n_samples:
            return samples[:n_samples]

        # Nested / feedback methods (Dynesty, MultiNest, MCMC, …) have no
        # propose_next: draw unit-cube u and let the Worker mapper + plan run.
        try:
            from jarvishep2.Sampling.variables import load_variables

            variables = load_variables(self.config)
            ndim = max(1, len(variables))
        except Exception:
            ndim = 2
        remaining = n_samples - len(samples)
        seed = 0
        try:
            sampling = dict(self.config.get("Sampling") or {})
            bounds = (
                sampling.get("Bounds")
                if isinstance(sampling.get("Bounds"), dict)
                else {}
            )
            # seed lives under Sampling.Bounds only (V2).
            seed = int(bounds.get("seed", 0) or 0)
        except (TypeError, ValueError):
            seed = 0
        rng = np.random.default_rng(seed if seed else None)
        for _ in range(remaining):
            u_coords = rng.random(ndim).astype(np.float64)
            samples.append(self.sampler._build_sample(u_coords))

        if not samples:
            raise RuntimeError(
                "Jarvis check could not build any smoke samples: "
                "CSV missing and sampler produced no points"
            )
        return samples

    def run_distributed_scan(self, *, timeout: float = 3600.0) -> int:
        """Drive a stateless sampler through propose → Redis → Archiver."""
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        start_records = self._archiver_records_written()
        requeued = 0
        if self._resume_policy == "resume" and hasattr(self.sampler, "repropose_unfinished"):
            requeued = len(self.sampler.repropose_unfinished())
        if hasattr(self.sampler, "run_distributed"):
            pushed = int(self.sampler.run_distributed())
        elif hasattr(self.sampler, "submit_all_remaining"):
            pushed = len(self.sampler.submit_all_remaining())
        else:
            raise RuntimeError(
                f"Sampler {type(self.sampler).__name__} does not implement a distributed run API"
            )
        expected_records = start_records + requeued + pushed
        total_work = requeued + pushed
        if expected_records > start_records:
            self.wait_for_results(
                expected_records,
                timeout=timeout,
                progress_total=total_work,
                progress_base=start_records,
            )
        self._finalize_sample_buckets()
        return requeued + pushed

    def run_adaptive_scan(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Drive feedback samplers; timeout is a per-generation wait, not a run budget."""
        if timeout is not None:
            generation_timeout = timeout
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        if self.redis is None:
            raise RuntimeError("adaptive scan requires redis")
        if hasattr(self.sampler, "set_redis"):
            self.sampler.set_redis(self.redis)
        # Stale feedback from a prior run must not poison barriers.
        drained_fb = int(self.redis.drain_feedback_queue())
        if drained_fb:
            self._logger.info("drained %d stale feedback record(s)", drained_fb)
        requeued = 0
        if self._resume_policy == "resume" and hasattr(self.sampler, "repropose_unfinished"):
            requeued = len(self.sampler.repropose_unfinished())
        if hasattr(self.sampler, "run_adaptive"):
            pushed = int(
                self.sampler.run_adaptive(generation_timeout=generation_timeout)
            )
        elif hasattr(self.sampler, "run_distributed"):
            pushed = int(self.sampler.run_distributed())
        else:
            raise RuntimeError(
                f"Sampler {type(self.sampler).__name__} does not implement run_adaptive"
            )
        # Nested/adaptive samplers barrier on *feedback* first; Archiver lags
        # behind on the archive queue. Wait until DATABASE rows catch Redis
        # completed counts so samples.csv export is not a partial HDF5 snapshot.
        self._wait_for_archive_caught_up(timeout=generation_timeout)
        self._finalize_sample_buckets()
        return requeued + pushed

    def _wait_for_archive_caught_up(self, *, timeout: float = 3600.0) -> None:
        """Block until Archiver has persisted every Redis-completed sample.

        Process-mode :class:`ArchiverProcess.drain` is a parent-side no-op; the
        child keeps consuming ``hep:archive``. Plot/CSV export must wait until
        ``records_written >= ok + failed`` (and workers are idle), otherwise
        ``DATABASE/samples.csv`` misses the last archive batches.
        """
        if self.archiver is None:
            return
        deadline = time.monotonic() + max(1.0, float(timeout))
        # Stabilize Redis counters once workers finish the feedback barrier.
        # Resume may deliberately replay a partially persisted adaptive
        # generation to reconstruct its feedback.  Those Worker completions
        # are not new DATABASE rows because the Archiver deduplicates UUIDs.
        # Prefer the sampler's identity set; counters remain the fallback for
        # nested/external samplers that do not expose submitted UUIDs.
        expected = 0
        while time.monotonic() < deadline:
            counters = self._live_sample_counters()
            expected = int(counters["ok"]) + int(counters["failed"])
            inflight = int(counters["running"]) + int(counters["queued"])
            if inflight == 0:
                break
            time.sleep(0.05)
        submitted_uuids = {
            str(item)
            for item in getattr(self.sampler, "submitted_uuids", frozenset())
            if str(item)
        }
        if bool(getattr(self.sampler, "uses_indexed_resume", False)):
            # Redis sample counters are run-local and can include work from a
            # killed control process.  Indexed samplers instead have one
            # durable logical target: the highest proposed sample index.
            expected = max(0, int(getattr(self.sampler, "_accepted_index", expected) or expected))
        elif submitted_uuids or self._persisted_uuids:
            expected = len(submitted_uuids | self._persisted_uuids)
        # Parent drain is meaningful only for in-process SimpleArchiver.
        try:
            self.archiver.drain(idle_timeout=2.0)
        except Exception:
            pass
        if expected <= 0:
            return
        remaining = max(1.0, deadline - time.monotonic())
        self.wait_for_results(
            expected,
            timeout=remaining,
            progress_total=expected,
            progress_base=0,
        )

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
        wait_timeout = (
            float(timeout)
            if timeout is not None
            else check_modules_timeout_sec(self.config)
        )
        if wait_timeout <= 0:
            wait_timeout = check_modules_timeout_sec(self.config)
        sample_root = self._resolve_sample_root()
        self._logger.warning(
            "Start Jarvis check smoke (assembly-line test) "
            "workers=1 sample_root=%s layout=flat-uuid pack=off timeout_sec=%.1f",
            sample_root,
            wait_timeout,
        )
        samples = self._build_check_module_samples()
        if not samples:
            raise RuntimeError("Jarvis check produced an empty sample list")
        self.submit_samples(samples)
        self.wait_for_results(len(samples), timeout=wait_timeout)
        self._finalize_sample_buckets()
        if verify_golden is not None:
            task_result_dir = str(self.info.get("task_result_dir") or os.getcwd())
            verify_check_modules_golden(task_result_dir=task_result_dir, golden=verify_golden)
        self._logger.warning(
            "Jarvis check finished: submitted %d sample(s) "
            "(pipeline smoke only — not a full MultiNest/Dynesty evidence run; "
            "inspect artifacts under %s)",
            len(samples),
            self._resolve_sample_root(),
        )
        return len(samples)

    def _apply_check_modules_runtime_policy(self) -> None:
        """Force smoke-friendly layout: 1 worker, SAMPLE/test, no tar pack.

        Applied before bootstrap so Factory/Archiver/Workers all see the policy.
        """
        # --- single worker ---
        runtime = dict(get_runtime_block(self.config))
        runtime["workers"] = 1
        self.config["Runtime"] = runtime
        self.runtime = runtime

        envreqs = self.config.get("EnvReqs")
        envreqs = dict(envreqs) if isinstance(envreqs, Mapping) else {}
        v2 = envreqs.get("V2")
        v2 = dict(v2) if isinstance(v2, Mapping) else {}
        v2["workers"] = 1
        # SAMPLE/test layout: no numbered buckets (flat uuid dirs), no tar pack.
        sample_dir = dict(get_sample_directory_config(self.config))
        sample_dir["pack"] = False
        sample_dir["enabled"] = False  # Worker skips Redis 00000N allocator
        v2_sd = v2.get("sample_directory")
        if isinstance(v2_sd, Mapping):
            merged_sd = dict(v2_sd)
            merged_sd["pack"] = False
            merged_sd["enabled"] = False
            v2["sample_directory"] = merged_sd
        else:
            v2["sample_directory"] = {"pack": False, "enabled": False}
        archiver = v2.get("archiver")
        archiver = dict(archiver) if isinstance(archiver, Mapping) else {}
        archiver["pack_buckets"] = False
        v2["archiver"] = archiver
        envreqs["V2"] = v2
        self.config["EnvReqs"] = envreqs

        # --- One calculator PackID only ---
        calculators = self.config.get("Calculators")
        calculators = dict(calculators) if isinstance(calculators, Mapping) else {}
        # Workers=1 alone is not enough: Redis still registers N slots from
        # ``make_parallel`` / Pools. Free-list rotation then hands sample 1 → 001,
        # sample 2 → 002, … and each new PackID runs a full clone_shadow install
        # (e.g. full micrOMEGAs copy+make). Check is a smoke: pin pool size to 1
        # so every point reuses ``001`` after the first install.
        calculators["make_parallel"] = 1
        if "Pools" in calculators:
            pools = calculators.get("Pools")
            if isinstance(pools, Mapping):
                calculators["Pools"] = {str(k): 1 for k in pools}
            else:
                calculators.pop("Pools", None)
        if "pools" in calculators:
            pools_l = calculators.get("pools")
            if isinstance(pools_l, Mapping):
                calculators["pools"] = {str(k): 1 for k in pools_l}
            else:
                calculators.pop("pools", None)
        modules = calculators.get("Modules")
        if isinstance(modules, list):
            pinned_modules: list[Any] = []
            for item in modules:
                if isinstance(item, Mapping):
                    mod = dict(item)
                    mod["make_parallel"] = 1
                    pinned_modules.append(mod)
                else:
                    pinned_modules.append(item)
            calculators["Modules"] = pinned_modules
        self.config["Calculators"] = calculators

        # Layout flag: SAMPLE/test instead of SAMPLE/
        self.config["_check_modules_sample_layout"] = True
        # Eager path so logging / helpers work even before _init_sample_buckets.
        root = self._resolve_sample_root()
        os.makedirs(root, exist_ok=True)
        self.info["sample_root"] = root
        self.info["check_modules"] = True

    def _resolve_sample_root(self) -> str:
        """Return SAMPLE root; ``Jarvis check`` uses ``SAMPLE/test`` (no tar pack)."""
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        base = os.path.join(task_result_dir, "SAMPLE")
        if bool(self.config.get("_check_modules_sample_layout")) or bool(
            self.info.get("check_modules")
        ):
            return os.path.join(base, "test")
        cached = self.info.get("sample_root")
        if isinstance(cached, str) and cached.strip():
            return os.path.abspath(cached)
        return base

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
        return str(self.runtime.get("mode", "auto")).strip().lower() == "redis"

    def _resolve_checkpoint_sampler_name(self) -> str:
        """Directory name under ``checkpoints/<scan>/``.

        Must match what a running scan writes (``sampler.method``, e.g. ``Dynesty``).
        ``prepare_resume`` runs **before** ``init_sampler_from_config``, so
        ``info["sampler_name"]`` is still the placeholder ``SamplingVirtial`` —
        fall through to YAML ``Sampling.Method`` instead of looking for a
        non-existent ``…/SamplingVirtial/state.pkl`` (which silently skipped
        nested engine restore and made Dynesty restart from scratch).
        """
        placeholders = {"", "sampler", "samplingvirtial", "samplingvirtual"}
        candidates: list[str] = []
        if self.sampler is not None:
            method = getattr(self.sampler, "method", None)
            if method:
                candidates.append(str(method))
            candidates.append(type(self.sampler).__name__)
        info_name = str(self.info.get("sampler_name") or "")
        if info_name:
            candidates.append(info_name)
        try:
            method_from_yaml = sampling_method(self.config)
        except Exception:
            method_from_yaml = ""
        if method_from_yaml:
            candidates.append(str(method_from_yaml))
        for name in candidates:
            text = str(name or "").strip()
            if text and text.lower() not in placeholders:
                return text
        return "sampler"

    def checkpoint_file(self) -> str:
        task_root = str(
            self.info.get("task_root")
            or self.config.get("task_root")
            or os.getcwd()
        )
        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        sampler_name = self._resolve_checkpoint_sampler_name()
        return checkpoint_path(
            task_root=task_root,
            scan_name=scan_name,
            sampler_name=sampler_name,
        )

    @staticmethod
    def prompt_resume_from_checkpoint(
        checkpoint_file: str,
        *,
        timeout_seconds: float = 30.0,
    ) -> bool:
        if not getattr(sys.stdin, "isatty", lambda: False)():
            return True
        response: dict[str, str | None] = {"text": None}

        def _read() -> None:
            try:
                response["text"] = sys.stdin.readline()
            except Exception:
                response["text"] = ""

        print(RESUME_PROMPT, end="", flush=True)
        reader = threading.Thread(target=_read, daemon=True)
        reader.start()
        reader.join(timeout_seconds)
        if reader.is_alive():
            return True
        answer = str(response["text"] or "").strip().lower()
        return answer not in {"y", "yes"}

    def _preload_resume_checkpoint(
        self,
        *,
        resume: bool = False,
        fresh: bool = False,
    ) -> None:
        if fresh:
            self._resume_policy = "fresh"
            self._resume_checkpoint_payload = None
            return

        path = self.checkpoint_file()
        if not os.path.exists(path):
            # Explicit --resume also supports DATABASE-only recovery when the
            # checkpoint is missing or deliberately removed.
            self._resume_policy = "resume" if resume else "fresh"
            self._resume_checkpoint_payload = None
            return

        if not resume:
            if not self.prompt_resume_from_checkpoint(path):
                self._resume_policy = "fresh"
                self._resume_checkpoint_payload = None
                try:
                    os.remove(path)
                except OSError:
                    pass
                self._logger.warning(
                    "Starting a fresh run from user confirmation; existing checkpoint was discarded."
                )
                return
        self._resume_policy = "resume"
        try:
            payload = load_checkpoint(path)
        except ValueError as exc:
            self._logger.warning("Checkpoint rejected -> %s", exc)
            self._resume_checkpoint_payload = None
            self._resume_policy = "fresh"
            return
        # D22.5: refuse resume when Sampling.Mapper / Variables fingerprint drifts.
        ok_mapper, mapper_reason = check_mapper_fingerprint(payload, self.config)
        if not ok_mapper:
            self._logger.error("Checkpoint rejected (starting fresh) -> %s", mapper_reason)
            self._resume_checkpoint_payload = None
            self._resume_policy = "fresh"
            # Drop the stale checkpoint so a subsequent auto-resume cannot
            # silently rebind uuid → coordinates under a changed Mapper.
            try:
                os.remove(path)
            except OSError:
                pass
            return
        # D23.8: refuse when Operas constant values drifted (outside the card).
        # D23.13: if Operas is merely unavailable, skip with WARNING — never treat
        # that as drift (would delete a good checkpoint and re-run from index 0).
        ok_consts, consts_reason = check_operas_constants_fingerprint(payload)
        if not ok_consts:
            self._logger.error("Checkpoint rejected (starting fresh) -> %s", consts_reason)
            self._resume_checkpoint_payload = None
            self._resume_policy = "fresh"
            try:
                os.remove(path)
            except OSError:
                pass
            return
        if str(consts_reason).startswith("skip:"):
            self._logger.warning("%s", consts_reason)
        self._resume_checkpoint_payload = payload

    def _load_persisted_database_state(self, database_dir: str) -> None:
        """Refresh index prefix, UUID set, and completed/failed counts from DATABASE.

        Index prefix drives resume fast-forward; outcome counts seed Redis
        SAMPLE_STATS so a resume does not re-label durable Failed rows as
        successful (D21.13).
        """
        (
            self._persisted_index_prefix,
            pending_indices,
            self._persisted_uuids,
        ) = read_persisted_sample_index_state(database_dir)
        completed, failed = read_persisted_outcome_counts(database_dir)
        total_from_status = completed + failed
        # Prefer status counts when any durable row exists; fall back to the
        # index/legacy UUID reconstruction for pre-status shards.
        if total_from_status > 0:
            self._persisted_records_count = total_from_status
            self._persisted_failed_count = failed
        else:
            self._persisted_records_count = (
                self._persisted_index_prefix
                + len(pending_indices)
                + len(self._persisted_uuids)
            )
            self._persisted_failed_count = 0

    def apply_resume_checkpoint(self, worker_config: Mapping[str, Any] | None = None) -> int:
        if self._resume_policy != "resume":
            return 0
        if self.redis is None:
            raise RuntimeError("init_redis() must run before apply_resume_checkpoint()")
        drained = prepare_resume(
            self.redis,
            worker_config=worker_config,
            persisted_count=self._persisted_records_count,
            persisted_failed=self._persisted_failed_count,
        )
        imported = self._import_sampler_checkpoint_state(self.sampler)
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_uuids"):
            self.sampler.set_persisted_uuids(self._persisted_uuids)
        if self.sampler is not None and hasattr(self.sampler, "set_persisted_prefix"):
            self.sampler.set_persisted_prefix(self._persisted_index_prefix)
        if (
            bool(getattr(self.sampler, "uses_indexed_resume", False))
            and hasattr(self.sampler, "advance_to_persisted_prefix")
        ):
            self.sampler.advance_to_persisted_prefix(self._persisted_index_prefix)
        if (
            imported
            and self.sampler is not None
            and not bool(getattr(self.sampler, "uses_indexed_resume", False))
            and hasattr(self.sampler, "set_resume_repropose_hint")
        ):
            self.sampler.set_resume_repropose_hint(True)
        self._logger.info(
            "Resumed from checkpoint; drained %d stale task(s) from hep:task_queue",
            drained,
        )
        return drained

    def save_runtime_checkpoint(self, *, reason: str = "") -> str | None:
        if self.sampler is None or not hasattr(self.sampler, "export_runtime_state"):
            return None
        persistence: dict[str, Any] = {}
        if self.archiver is not None and hasattr(self.archiver, "persistence_state"):
            persistence = dict(self.archiver.persistence_state())
        prefix = max(
            self._persisted_index_prefix,
            int(persistence.get("persisted_prefix", 0) or 0),
        )
        self._persisted_index_prefix = prefix
        if hasattr(self.sampler, "set_persisted_prefix"):
            self.sampler.set_persisted_prefix(prefix)
        task_root = str(self.info.get("task_root") or self.config.get("task_root") or os.getcwd())
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or task_root
        )
        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        sampler_name = str(self.info.get("sampler_name") or type(self.sampler).__name__)
        run_spec = build_run_spec(
            config=self.config,
            scan_name=scan_name,
            task_root=task_root,
            task_result_dir=task_result_dir,
            sampler_name=sampler_name,
            worker_parallel=int(self.runtime.get("workers", 0) or 0),
        )
        indexed = bool(getattr(self.sampler, "uses_indexed_resume", False))
        if indexed:
            # ``prefix`` is the sole hot-path completion fact.  It arrives
            # from the Archiver after fsync and is independent of scan size.
            # Do not let ``checkpoint_runtime_state`` capture the current
            # generator state here: only a batch-boundary snapshot may become
            # a resume point.
            safe = False
        else:
            submitted = frozenset(getattr(self.sampler, "submitted_uuids", frozenset()))
            at_barrier = bool(getattr(self.sampler, "at_safe_barrier", lambda: False)())
            safe = at_barrier and submitted <= frozenset(
                getattr(self.sampler, "_persisted_uuids", set())
            )
        # Nested samplers override checkpoint_runtime_state to always re-export
        # the live dynesty engine (ignore stale FeedbackSampler barrier copies).
        if hasattr(self.sampler, "checkpoint_runtime_state"):
            # Interrupt / dynesty_engine_checkpoint: always force a live export
            # for methods that do not use generation barriers.
            force_live = str(reason or "") in {
                "interrupt",
                "dynesty_engine_checkpoint",
                "dynesty_finished",
            } or bool(getattr(self.sampler, "method", "") in {"Dynesty", "MultiNest"})
            sampler_state = self.sampler.checkpoint_runtime_state(
                safe=True if force_live else safe
            )
        else:
            sampler_state = self.sampler.export_runtime_state()
        if indexed:
            checkpoint_index = int(
                sampler_state.get("checkpoint_index", sampler_state.get("accepted_index", 0)) or 0
            )
            safe = checkpoint_index <= prefix
        payload = build_payload(
            run_spec=run_spec,
            sampler_state=sampler_state,
            persistence=persistence,
            reason=reason,
            safe_barrier_confirmed=safe,
        )
        path = self.checkpoint_file()
        save_checkpoint(path, payload)
        return path

    def _save_interrupt_checkpoint(self) -> str | None:
        """Force a resume checkpoint on SIGINT/SIGTERM (D8.3 human remainder).

        Agent-facing pieces (run_state interrupted, stop-ack) stay parked with D8.
        Failures are logged and never re-raised — teardown must still run.
        """
        try:
            path = self.save_runtime_checkpoint(reason="interrupt")
            if path:
                try:
                    self._logger.warning(
                        "interrupt checkpoint written → %s "
                        "(resume: Jarvis run <task.yaml> --resume)",
                        path,
                    )
                except Exception:
                    pass
            return path
        except Exception as exc:
            try:
                self._logger.warning("interrupt checkpoint failed → %s", exc)
            except Exception:
                pass
            return None

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        # EnvReqs.V2.redis (optional) overlays INTERNAL_REDIS_CONFIG defaults.
        redis_config = get_redis_config(self.config)
        raw_runtime = self.config.get("Runtime")
        raw_redis = raw_runtime.get("redis") if isinstance(raw_runtime, Mapping) else None
        # Task cards always use a Jarvis-managed local broker.  Programmatic
        # callers (including the multiprocessing fakeredis integration suite)
        # may deliberately attach to an already-running broker instead.  This
        # internal Runtime marker is not a task-card YAML setting.
        use_managed_redis = not (
            isinstance(raw_redis, Mapping) and raw_redis.get("managed") is False
        )
        # Prefer a Jarvis-managed redis-server so `ps` shows Jarvis-Redis:<scan>.
        # A port collision is reassigned into the project's default YAML first.
        if client is None and use_managed_redis:
            self._ensure_managed_redis(redis_config)
            redis_config = get_redis_config(self.config)
        if client is not None:
            self.redis = RedisQueue(redis_config, client=client)
        else:
            self.redis = RedisQueue(redis_config)
        self.redis.connect()
        try:
            self.redis.ping()
        except Exception as exc:
            host = redis_config.get("host", "127.0.0.1")
            port = redis_config.get("port", 6379)
            raise RuntimeError(
                f"internal Redis is unavailable at {host}:{port}; "
                "install redis-server or start a local Redis service"
            ) from exc
        return self.redis

    def _ensure_managed_redis(self, redis_config: Mapping[str, Any]) -> None:
        """Start redis-server as Jarvis-Redis:<scan> when the port is free."""
        if self._managed_redis is not None and self._managed_redis.started_by_us:
            return
        redis_config = dict(redis_config)
        reassigned = False
        host = str(redis_config.get("host") or "127.0.0.1")
        requested_port = int(redis_config.get("port") or 6379)
        if redis_port_open(host, requested_port):
            defaults_path = self.config.get("environment_defaults_path")
            if not defaults_path or self.config.get("redis_port_task_override"):
                raise RuntimeError(
                    f"Redis port {host}:{requested_port} is already in use. "
                    "Set EnvReqs.V2.redis.port in the project default environment YAML "
                    "(not in the task YAML) to enable automatic port reassignment."
                )
            selected_port = find_available_redis_port(host, requested_port + 1)
            update_default_redis_port(str(defaults_path), selected_port)
            redis_config["port"] = selected_port
            self._set_runtime_redis_config(redis_config)
            reassigned = True
            self._logger.warning(
                "Redis port %s:%s is occupied; reassigned this project to %s and updated %s",
                host,
                requested_port,
                selected_port,
                defaults_path,
            )
        scan_name = str(
            self.info.get("scan_name") or self.config.get("scan_name") or "scan"
        ).strip()
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        work_dir = os.path.join(task_result_dir, ".redis")
        managed = ManagedRedisServer.from_redis_config(
            redis_config,
            scan_name=scan_name,
            work_dir=work_dir,
        )
        try:
            started = managed.ensure(
                scan_name=scan_name,
                work_dir=work_dir,
                force_start=reassigned,
            )
        except Exception as exc:
            if reassigned:
                raise RuntimeError(
                    "the replacement Redis port became occupied before Jarvis could start it; "
                    "run again to select another port"
                ) from exc
            # Fall through to connect attempt; init_redis will raise a clear error.
            self._logger.warning("managed redis-server ensure failed -> %s", exc)
            return
        self._managed_redis = managed
        if started:
            self._logger.info(
                "started managed redis-server as %s on %s:%s",
                managed.title,
                managed.host,
                managed.port,
            )
        else:
            self._logger.info(
                "using existing Redis at %s:%s "
                "(process title is only set when Jarvis starts redis-server)",
                managed.host,
                managed.port,
            )

    def _set_runtime_redis_config(self, redis_config: Mapping[str, Any]) -> None:
        """Keep the loaded task in sync with an automatically reassigned port."""
        normalized = dict(redis_config)
        envreqs = dict(self.config.get("EnvReqs") or {})
        v2 = dict(envreqs.get("V2") or {})
        v2["redis"] = normalized
        envreqs["V2"] = v2
        self.config["EnvReqs"] = envreqs
        runtime = dict(self.config.get("Runtime") or {})
        runtime["redis"] = normalized
        self.config["Runtime"] = runtime
        self.runtime = runtime

    def _control_lock_owner_id(self) -> str:
        run_id = str(self.info.get("run_id") or "jarvis2-run")
        return f"{os.uname().nodename}:{os.getpid()}:{run_id}"

    def _calculator_module_names(self) -> list[str]:
        modules = (self.config.get("Calculators") or {}).get("Modules") or []
        names: list[str] = []
        if isinstance(modules, list):
            for item in expand_calculator_modes(modules):
                if isinstance(item, Mapping):
                    info = mode_info(item)
                    name = info[0] if info is not None else str(item.get("name") or "").strip()
                    if name:
                        names.append(name)
        return list(dict.fromkeys(names))

    def _claim_redis_control_lock(self) -> None:
        """Refuse to start if another Jarvis still holds the Redis control lease."""
        if self.redis is None:
            raise RuntimeError("init_redis() must run before claiming the control lock")
        owner = self._control_lock_owner_id()
        if self.redis.claim_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
            self._control_lock_owner = owner
            self._start_control_lease_refresh()
            self._logger.info("claimed Redis control lock (%s)", owner)
            return
        current = self.redis.get_control_lock_owner()
        if self._resume_policy == "resume" and current:
            parts = str(current).split(":", 2)
            same_host = bool(parts and parts[0] == os.uname().nodename)
            stale_pid = False
            if same_host and len(parts) >= 2:
                try:
                    os.kill(int(parts[1]), 0)
                except ProcessLookupError:
                    stale_pid = True
                except (OSError, ValueError):
                    stale_pid = False
            if stale_pid:
                self.redis.release_control_lock(str(current))
                if self.redis.claim_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
                    self._control_lock_owner = owner
                    self._start_control_lease_refresh()
                    self._logger.warning(
                        "took over stale Redis control lease from dead pid %s", parts[1]
                    )
                    return
        raise RuntimeError(
            "another Jarvis instance already holds the Redis control lock "
            f"({current!r}). Stop residual Jarvis/HEP2-Worker processes "
            "(do not leave runs suspended with ^Z), then retry. "
            "Emergency: redis-cli DEL hep:control:lock"
        )

    def _control_lease_loop(self, stop: threading.Event, owner: str) -> None:
        interval = max(1.0, CONTROL_LOCK_TTL_SEC / 3.0)
        while not stop.wait(interval):
            redis = self.redis
            if redis is None:
                return
            try:
                if not redis.refresh_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
                    self._logger.error("lost Redis control lease; requesting shutdown")
                    self._interrupt_requested = True
                    return
            except Exception as exc:
                self._logger.warning("control lease refresh failed -> %s", exc)

    def _start_control_lease_refresh(self) -> None:
        owner = str(self._control_lock_owner or "")
        if not owner:
            return
        self._stop_control_lease_refresh()
        stop = threading.Event()
        thread = threading.Thread(
            target=self._control_lease_loop,
            args=(stop, owner),
            name="Jarvis-ControlLease",
            daemon=True,
        )
        self._control_lease_stop = stop
        self._control_lease_thread = thread
        thread.start()

    def _stop_control_lease_refresh(self) -> None:
        if self._control_lease_stop is not None:
            self._control_lease_stop.set()
        if self._control_lease_thread is not None:
            self._control_lease_thread.join(timeout=2.0)
        self._control_lease_stop = None
        self._control_lease_thread = None

    def _reset_redis_for_fresh_run(self) -> None:
        """Drop ephemeral queues/stats/calc pools before Workers are spawned."""
        if self.redis is None:
            return
        workers = int(self.runtime.get("workers", 0) or 0)
        result = self.redis.reset_run_ephemeral_keys(
            calculator_names=self._calculator_module_names(),
            worker_ids=max(workers, 32),
        )
        self.redis.clear_archived_uuids(str(self.info.get("scan_name") or "scan"))
        self.redis.clear_archived_index_prefix(str(self.info.get("scan_name") or "scan"))
        self._logger.info(
            "reset Redis run keys (deleted=%s sample_stats_reset=%s)",
            result.get("deleted_keys"),
            result.get("sample_stats_reset"),
        )

    def _release_redis_control_lock(self) -> None:
        self._stop_control_lease_refresh()
        owner = getattr(self, "_control_lock_owner", None)
        if not owner or self.redis is None:
            return
        try:
            self.redis.release_control_lock(str(owner))
            self._logger.info("released Redis control lock (%s)", owner)
        except Exception as exc:
            self._logger.warning("failed to release Redis control lock -> %s", exc)
        self._control_lock_owner = None

    def _publish_runtime_metadata(self) -> None:
        if self.redis is None:
            return
        redis_config = self.redis.connection_config()
        path = write_scan_metadata(config=self.config, info=self.info, redis=redis_config)
        self.redis.set_runtime_metadata_path(path)
        self.info["runtime_metadata_path"] = path

    def init_command_parser(self) -> CommandParser:
        """Run Phase-1 static command resolution for the loaded task config."""
        root_path = self._environment_summary.get("ROOT path")
        self.command_parser = build_command_parser(
            self.config,
            root_path=str(root_path) if root_path else None,
        )
        return self.command_parser

    def prepare_libraries(self, *, parser: CommandParser | None = None) -> dict[str, dict[str, Any]]:
        """Install or reuse LibDeps before Redis, Workers, or an Archiver start."""
        root_path = self._environment_summary.get("ROOT path")
        library_parser = parser or build_command_parser(
            self.config,
            root_path=str(root_path) if root_path else None,
            register_executables=False,
        )
        logs_dir = str(
            self.info.get("logs_dir")
            or os.path.join(
                str(self.info.get("task_root") or os.getcwd()),
                "logs",
                str(self.info.get("scan_name") or "scan"),
            )
        )
        result = LibraryInstaller(
            self.config,
            parser=library_parser,
            logs_dir=logs_dir,
            logger=self._logger,
            skip_installation=self._skip_library_installation,
        ).prepare()
        if result:
            self.info["library_installations"] = result
        return result

    def _apply_command_parser_to_worker_config(self, worker_config: dict[str, Any]) -> dict[str, Any]:
        if self.command_parser is None:
            self.init_command_parser()
        assert self.command_parser is not None
        merged = dict(worker_config)
        calculator_modules = merged.get("calculator_modules")
        if calculator_modules:
            merged["calculator_modules"] = prepare_calculator_modules(
                calculator_modules,
                self.command_parser,
            )
        merged["command_parser"] = self.command_parser.to_picklable()
        return merged

    def build_worker_config(self, **overrides: Any) -> dict[str, Any]:
        """Build a picklable Worker blueprint with Phase-1 command resolution applied."""
        task_result_dir = str(
            overrides.pop("task_result_dir", None)
            or self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        # Propagate CLI console policy to every Worker process.
        overrides.setdefault("log_silence", bool(getattr(self, "_log_silence", False)))
        console_level = str(getattr(self, "_console_level", "WARNING")).strip().upper()
        overrides.setdefault("console_level", console_level)
        sample_config = dict(overrides.get("sample_config") or {})
        sample_config.setdefault(
            "sample_console_level",
            "DEBUG" if console_level == "DEBUG" else "ERROR",
        )
        sample_config.setdefault(
            "sample_log_silence",
            bool(getattr(self, "_log_silence", False)),
        )
        overrides["sample_config"] = sample_config
        # Prefer control-resolved scan logs dir so Workers never land in cwd/logs/jarvis_worker_*.
        if self.info.get("logs_dir"):
            overrides.setdefault("logs_dir", str(self.info["logs_dir"]))
        if self.info.get("scan_name"):
            overrides.setdefault("scan_name", str(self.info["scan_name"]))
        if self.info.get("task_root"):
            overrides.setdefault("task_root", str(self.info["task_root"]))
        if self._control_lock_owner:
            overrides.setdefault("control_lock_owner", str(self._control_lock_owner))
        extra = dict(overrides or {})
        # Stamp active sampler into feedback_return resolution (D13.8).
        if "sampler" not in extra and self.sampler is not None:
            extra["sampler"] = self.sampler
        # Prefer Jarvis check SAMPLE/test (or any resolved sample_root) over SAMPLE/.
        sample_dirs = extra.pop("sample_dirs", None)
        if sample_dirs is None:
            sample_dirs = self._resolve_sample_root()
        return build_worker_config(
            self.config,
            task_result_dir=task_result_dir,
            parser=self.command_parser,
            calculator_modules=extra.pop("calculator_modules", None),
            likelihood_expressions=extra.pop("likelihood_expressions", None),
            opera_modules=extra.pop("opera_modules", None),
            sample_dirs=sample_dirs,
            extra=extra or None,
        )

    def init_factory(self, worker_config: Mapping[str, Any] | None = None) -> TaskFactory | None:
        if not self.is_redis_runtime():
            self._logger.info("Runtime.mode != redis; skipping TaskFactory bring-up")
            return None

        redis_config = get_redis_config(self.config)
        # Explicit instance ownership (D9.4) — not the deprecated singleton shell.
        if self.factory is None:
            self.factory = TaskFactory(redis_config)
        elif redis_config and not self.factory.redis_config:
            self.factory.redis_config.update(redis_config)
        if self.redis is not None:
            self.factory.redis = self.redis
        else:
            self.factory.init_redis()

        workers = int(self.runtime.get("workers", 1) or 1)
        if workers <= 0:
            workers = 1

        if worker_config is None:
            if self.command_parser is None:
                self.init_command_parser()
            merged_config = self.build_worker_config()
        else:
            merged_config = dict(worker_config)
            if "command_parser" not in merged_config:
                merged_config = self._apply_command_parser_to_worker_config(merged_config)
        if self._control_lock_owner:
            merged_config.setdefault("control_lock_owner", str(self._control_lock_owner))
        # Direct Python callers may supply the public parent/mode shape as a
        # prebuilt worker config.  Normalize here as well as in the standard
        # config builder so Redis and the Worker always agree on dotted names.
        merged_config["calculator_modules"] = expand_calculator_modes(
            merged_config.get("calculator_modules") or []
        )
        # Calculator install control is single-writer: bump epochs and write
        # jarvis_install.json in this control process before Workers spawn.
        merged_config["calculator_modules"] = prepare_install_controls(
            merged_config.get("calculator_modules") or [],
            logger=self._logger,
        )
        self._install_control_modules = merged_config.get("calculator_modules") or []
        if self._resume_policy == "resume":
            prepare_resume(
                self.redis,
                worker_config=merged_config,
                persisted_count=self._persisted_records_count,
                persisted_failed=self._persisted_failed_count,
            )
        self.factory.start_workers(workers, **merged_config)
        self.factory.start_monitor()
        watchdog = get_watchdog_config(self.config)
        self.factory.start_watchdog(**watchdog)
        self._logger.info("TaskFactory started with %d worker(s)", workers)
        return self.factory

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
        if self.redis is None:
            raise RuntimeError("init_redis() must run before init_archiver()")
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        sample_root = self._resolve_sample_root()
        os.makedirs(database_dir, exist_ok=True)
        os.makedirs(sample_root, exist_ok=True)
        resolved_db_path = db_path or os.path.join(database_dir, "samples.hdf5")
        archiver_config = dict(get_archiver_config(self.config))
        if bool(self.config.get("_check_modules_sample_layout")):
            # Belt-and-suspenders: never tar SAMPLE/test during check smoke.
            archiver_config["pack_buckets"] = False
        # Same CLI console policy as Workers.
        archiver_config.setdefault("log_silence", bool(getattr(self, "_log_silence", False)))
        archiver_config.setdefault("console_level", str(getattr(self, "_console_level", "WARNING")))
        if self._control_lock_owner:
            archiver_config["control_lock_owner"] = str(self._control_lock_owner)
        delete_method = get_delete_method(self.config)
        redis_config = get_redis_config(self.config)
        if self.redis is not None:
            # Prefer the live connection if init_redis already bound a client.
            try:
                redis_config = dict(self.redis.connection_config())
            except Exception:
                pass

        from jarvishep2.logging import component_log_path

        scan_name = str(self.info.get("scan_name") or "").strip() or None
        logs_dir = str(self.info.get("logs_dir") or "").strip() or None
        if str(archiver_config.get("mode", "process")).strip().lower() == "process":
            archiver_log = component_log_path(logs_dir, "archiver") if logs_dir else None
            self.archiver = ArchiverProcess(
                redis_config,
                db_path=resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
                scan_name=scan_name,
                log_dir=logs_dir,
                log_path=archiver_log,
            )
            self.archiver.start()
            self._logger.info(
                "Archiver process started (own logger → %s)",
                archiver_log or "logs/<scan>/archiver.log",
            )
        else:
            # Thread mode: same process as control; Archiver logger → archiver.log.
            self.archiver = SimpleArchiver(
                self.redis,
                resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
                scan_name=scan_name,
                logger=get_jarvis_logger("archiver"),
            )
            self.archiver.start()
            self._logger.info("Archiver thread started (Jarvis-HEP.Archiver)")
        self._restore_archiver_persistence()
        return self.archiver

    def _restore_archiver_persistence(self) -> None:
        """Compatibility hook; DATABASE seeding now happens in SimpleArchiver.

        Never promote checkpoint/Redis acknowledgements into completion truth:
        only UUIDs actually read from HDF5 may seed the deduplication set.
        """
        return

    def _archiver_records_written(self) -> int:
        archiver = self.archiver
        if archiver is None:
            return 0
        counter = getattr(archiver, "records_written", 0)
        if hasattr(counter, "value"):
            return int(counter.value)
        return int(counter)

    def _discard_resume_checkpoint(
        self,
        *,
        reason: str,
        level: str = "error",
    ) -> None:
        """Drop a rejected resume snapshot and force a fresh sampling start."""
        log = self._logger.error if level == "error" else self._logger.warning
        log("Checkpoint sampler state rejected (starting fresh) -> %s", reason)
        self._resume_policy = "fresh"
        self._resume_checkpoint_payload = None
        path = self.checkpoint_file()
        try:
            if os.path.isfile(path):
                os.remove(path)
        except OSError:
            pass

    def _import_sampler_checkpoint_state(self, sampler: Any | None = None) -> bool:
        """Import sampler_state from the preloaded checkpoint, or start fresh.

        Shared by ``set_sampler``, late Method resolution, bootstrap, and
        ``apply_resume_checkpoint``. On Bounds/shape mismatch (``ValueError``)
        the stale checkpoint is discarded so a finished or drifted snapshot
        cannot silently yield ``submitted=0`` under a new card.

        Returns True when import succeeded.
        """
        target = sampler if sampler is not None else self.sampler
        if self._resume_policy != "resume":
            return False
        if self._resume_checkpoint_payload is None or target is None:
            return False
        if not hasattr(target, "import_runtime_state"):
            return False
        sampler_state = self._resume_checkpoint_payload.get("sampler_state") or {}
        if not isinstance(sampler_state, Mapping):
            sampler_state = {}
        try:
            target.import_runtime_state(sampler_state)
        except ValueError as exc:
            # MCMC Bounds drift (num_chains / num_iters / scale / seed) and
            # similar sampler-shape mismatches: refuse the stale snapshot.
            self._discard_resume_checkpoint(reason=str(exc), level="error")
            return False
        if (
            not bool(getattr(target, "uses_indexed_resume", False))
            and hasattr(target, "set_resume_repropose_hint")
        ):
            target.set_resume_repropose_hint(True)
        return True

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
        if self.sampler is None or not hasattr(self.sampler, "configure_checkpoint"):
            return
        self.sampler.configure_checkpoint(
            checkpoint_file=self.checkpoint_file(),
            save_callback=lambda reason="": self.save_runtime_checkpoint(reason=reason),
        )

    def submit_samples(self, samples: Sequence[Sample]) -> None:
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        if hasattr(self.sampler, "_submit_group"):
            self.sampler._submit_group(list(samples))
        else:
            for sample in samples:
                self.sampler._submit(sample)

    def _live_sample_counters(self) -> dict[str, int]:
        """Snapshot Redis sample/queue counters for completion progress lines."""
        ok = failed = running = queued = archive_q = 0
        if self.redis is not None:
            try:
                stats = self.redis.fetch_sample_stats()
                lengths = self.redis.get_queue_lengths()
                ok = int(stats.get("completed", 0) or 0)
                failed = int(stats.get("failed", 0) or 0)
                running = max(0, int(stats.get("running", 0) or 0))
                queued = int(lengths.get("task_queue_length", 0) or 0)
                archive_q = int(lengths.get("archive_queue_length", 0) or 0)
            except Exception:
                pass
        return {
            "ok": ok,
            "failed": failed,
            "running": running,
            "queued": queued,
            "archive_q": archive_q,
        }

    def wait_for_results(
        self,
        expected: int,
        *,
        timeout: float = 30.0,
        poll_interval: float = 0.1,
        progress_total: int | None = None,
        progress_base: int | None = None,
    ) -> None:
        """Block until Archiver has written ``expected`` records.

        When ``progress_total`` is set, emit V1-style ‰ completion heartbeats
        using archived count relative to ``progress_base`` (default 0).
        """
        if self.archiver is None:
            raise RuntimeError("archiver is not configured")
        deadline = time.monotonic() + max(0.1, float(timeout))
        base = int(progress_base or 0)
        total = int(progress_total) if progress_total is not None else max(0, int(expected) - base)
        progress: PermilleProgress | None = None
        if total > 0:
            # Control-plane archive ‰ heartbeats are noisy vs DataRecorder/Archiver;
            # keep them at DEBUG. Final "sample drain complete" stays INFO below.
            progress = PermilleProgress(
                self._logger,
                total=total,
                label="samples archived",
                level="debug",
                milestone_level="debug",
            )
            progress.update(
                0,
                extra="ok=? failed=? queued=? running=? archive_q=? archived=0/?",
                force=True,
            )

        last_written = -1
        stall_since: float | None = None
        while time.monotonic() < deadline:
            written = self._archiver_records_written()
            counters = self._live_sample_counters()
            if progress is not None:
                done = max(0, written - base)
                extra = (
                    f"ok={counters['ok']} failed={counters['failed']} "
                    f"queued={counters['queued']} running={counters['running']} "
                    f"archive_q={counters['archive_q']} "
                    f"archived={written}/{expected}"
                )
                progress.update(done, extra=extra)
            if written >= expected:
                if progress is not None:
                    progress.update(
                        total,
                        extra=(
                            f"ok={counters['ok']} failed={counters['failed']} "
                            f"queued={counters['queued']} running={counters['running']} "
                            f"archive_q={counters['archive_q']}"
                        ),
                        force=True,
                    )
                    self._logger.info(
                        "sample drain complete: %d/%d archived in %s",
                        total,
                        total,
                        format_duration(time.time() - progress.t0),
                    )
                return
            # Detect permanent stall: workers done, archive queue empty, count frozen.
            workers_done = (
                counters["running"] == 0
                and counters["queued"] == 0
                and (counters["ok"] + counters["failed"]) >= total
            )
            if workers_done and counters["archive_q"] == 0 and written == last_written:
                if stall_since is None:
                    stall_since = time.monotonic()
                elif time.monotonic() - stall_since >= 5.0:
                    raise TimeoutError(
                        f"archive drain stalled: workers finished "
                        f"(ok={counters['ok']} failed={counters['failed']}) but only "
                        f"{written}/{expected} DATABASE rows written and archive_q=0. "
                        f"Archiver made no progress; inspect the Archiver log and "
                        f"DATABASE persistence configuration."
                    )
            else:
                stall_since = None
            last_written = written
            time.sleep(poll_interval)
        raise TimeoutError(
            f"timed out waiting for {expected} archived results; "
            f"got {self._archiver_records_written()}"
        )

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
        if self.redis is None:
            return
        sample_dir = get_sample_directory_config(self.config)
        sample_root = self._resolve_sample_root()
        os.makedirs(sample_root, exist_ok=True)
        meta = self.redis.init_sample_buckets(
            sample_root=sample_root,
            limit=int(sample_dir.get("limit", 200)),
            width=int(sample_dir.get("width", 6)),
            start_bucket=int(sample_dir.get("start_bucket", 1)),
            pack=pack_buckets_enabled(self.config),
            enabled=bool(sample_dir.get("enabled", True)),
        )
        self.info["sample_root"] = sample_root
        self.info["sample_directory"] = dict(sample_dir)
        enabled = bool(sample_dir.get("enabled", True))
        self._logger.info(
            "SAMPLE ready -> root=%s enabled=%s limit=%s width=%s pack=%s",
            sample_root,
            enabled,
            meta.get("limit"),
            meta.get("width"),
            bool(int(meta.get("pack", 0))),
        )

    def _finalize_sample_buckets(self) -> None:
        """Seal the open bucket after all samples finish; Archiver packs when idle."""
        if self.redis is None:
            return
        try:
            self.redis.seal_current_sample_bucket()
        except Exception as exc:
            self._logger.warning("seal current SAMPLE bucket failed -> %s", exc)
        if self.archiver is not None:
            try:
                self.archiver.drain(idle_timeout=2.0)
            except Exception as exc:
                self._logger.warning("archiver bucket drain failed -> %s", exc)

    def _install_control_signal_handlers(self) -> None:
        """SIGINT/SIGTERM → clean shutdown; refuse Ctrl+Z suspend."""
        if self._signal_handlers_installed:
            return

        def _stop_handler(signum: int, _frame: Any) -> None:
            try:
                sig_name = signal.Signals(signum).name
            except Exception:
                sig_name = str(signum)
            self._interrupt_requested = True
            try:
                self._logger.warning(
                    "received %s — initiating clean shutdown "
                    "(Workers / Archiver / managed Redis)",
                    sig_name,
                )
            except Exception:
                pass
            # Unwind the main thread so run()'s finally → shutdown() always runs.
            raise KeyboardInterrupt(f"interrupted by {sig_name}")

        def _refuse_suspend(signum: int, _frame: Any) -> None:
            try:
                self._logger.warning(
                    "Ctrl+Z / SIGTSTP ignored — process suspend leaves Workers and "
                    "Redis half-alive. Use Ctrl+C to stop the scan cleanly."
                )
            except Exception:
                pass

        for sig in (signal.SIGINT, signal.SIGTERM):
            try:
                self._previous_signal_handlers[sig] = signal.getsignal(sig)
                signal.signal(sig, _stop_handler)
            except (ValueError, OSError):
                # Not in main thread or signal unsupported.
                pass
        if hasattr(signal, "SIGTSTP"):
            try:
                self._previous_signal_handlers[signal.SIGTSTP] = signal.getsignal(
                    signal.SIGTSTP
                )
                signal.signal(signal.SIGTSTP, _refuse_suspend)
            except (ValueError, OSError):
                pass
        self._signal_handlers_installed = True

    def _restore_control_signal_handlers(self) -> None:
        if not self._signal_handlers_installed:
            return
        for sig, previous in list(self._previous_signal_handlers.items()):
            try:
                signal.signal(sig, previous)
            except (ValueError, OSError):
                pass
        self._previous_signal_handlers.clear()
        self._signal_handlers_installed = False

    def _stop_managed_redis(self) -> None:
        managed = self._managed_redis
        if managed is None:
            return
        self._managed_redis = None
        try:
            if managed.started_by_us:
                try:
                    self._logger.info(
                        "stopping managed redis-server (%s)",
                        managed.title or "Jarvis-Redis",
                    )
                except Exception:
                    pass
            managed.stop()
        except Exception as exc:
            try:
                self._logger.warning("managed redis stop failed -> %s", exc)
            except Exception:
                pass

    def _atexit_cleanup(self) -> None:
        """Best-effort stop of managed Redis if the process dies unexpectedly."""
        if self._shutdown_done:
            return
        try:
            self._stop_managed_redis()
        except Exception:
            pass

    def shutdown(self, *, wait: bool = True, write_run_summary: bool = False) -> None:
        """Stop Archiver, Workers, release the control lock, and managed Redis.

        Idempotent: safe to call from ``finally``, signal paths, and atexit.
        Managed Redis is always stopped last in an inner ``finally`` so a failure
        earlier in the chain cannot leave ``Jarvis-Redis:*`` orphaned.
        """
        if self._shutdown_done:
            # Still try Redis in case a previous partial shutdown left it.
            self._stop_managed_redis()
            return
        self._shutdown_done = True
        try:
            if self.sampler is not None and hasattr(self.sampler, "shutdown_checkpointing"):
                try:
                    self.sampler.shutdown_checkpointing()
                except Exception as exc:
                    self._logger.warning("sampler checkpoint shutdown failed -> %s", exc)
            if not self._interrupt_requested:
                try:
                    self._finalize_sample_buckets()
                except Exception:
                    pass
            if self.archiver is not None:
                try:
                    if isinstance(self.archiver, ArchiverProcess):
                        # Interrupt path: don't hang forever on tar packing.
                        timeout = 5.0 if self._interrupt_requested else 30.0
                        self.archiver.stop(wait=wait, timeout=timeout)
                    else:
                        self.archiver.stop(
                            wait=wait,
                            drain=not self._interrupt_requested,
                        )
                except Exception as exc:
                    self._logger.warning("archiver stop failed -> %s", exc)
                self.archiver = None
            if write_run_summary and self.factory is not None:
                try:
                    self.write_run_summary()
                except Exception as exc:
                    self._logger.warning("run_summary write failed -> %s", exc)
            try:
                if self.redis is not None:
                    self.redis.clear_runtime_metadata_path()
                self._release_redis_control_lock()
            except Exception:
                pass
            if self.factory is not None:
                try:
                    self.factory.shutdown(wait=wait)
                except Exception as exc:
                    self._logger.warning("factory shutdown failed -> %s", exc)
                self.factory = None
            elif self.redis is not None:
                try:
                    self.redis.close()
                except Exception:
                    pass
            self.redis = None
            try:
                refresh_install_control_summaries(getattr(self, "_install_control_modules", []))
            except Exception as exc:
                self._logger.debug("install-control summary refresh skipped -> %s", exc)
        finally:
            self._stop_managed_redis()


__all__ = ["Jarvis2Core"]
