#!/usr/bin/env python3
"""Jarvis2 control-process orchestrator for the distributed runtime."""

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
from jarvishep2.dashboard import SnapshotReader, format_monitor_view
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor, STATELESS_METHODS
from jarvishep2.factory import TaskFactory
from jarvishep2.log_kv import PermilleProgress, format_duration
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.monitoring.run_summary import RunSummaryRenderer, build_run_summary
from jarvishep2.versioning import render_logo_with_version
from jarvishep2.redis_queue import (
    CONTROL_LOCK_TTL_SEC,
    RedisQueue,
)
from jarvishep2.redis_server import ManagedRedisServer
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_delete_method,
    get_factory_config,
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
from jarvishep2.task_config import (
    check_modules_points_path,
    is_check_modules_task,
    load_task_yaml,
    sampling_method,
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
        self._resume_checkpoint_payload: dict[str, Any] | None = None
        self._resume_policy: str = "auto"
        self._control_lock_owner: str | None = None
        self._managed_redis: ManagedRedisServer | None = None
        self._shutdown_done = False
        self._interrupt_requested = False
        self._signal_handlers_installed = False
        self._previous_signal_handlers: dict[int, Any] = {}
        self._logger = get_jarvis_logger("core")
        # Last-resort cleanup if the process exits without an orderly finally.
        atexit.register(self._atexit_cleanup)

    def init_logger(self) -> None:
        """Configure top-level logging with V1 visual format and scan-scoped file path."""
        from jarvishep2.proc_title import control_title, set_process_title

        scan_name = str(self.info.get("scan_name") or "scan")
        # Refresh process title once the scan name is known.
        set_process_title(control_title(scan_name=scan_name))
        task_root = str(self.info.get("task_root") or os.getcwd())
        logs_dir = os.path.join(task_root, "logs", scan_name)
        jarvis_log = os.path.join(logs_dir, f"{scan_name}.log")
        self.info["logs_dir"] = logs_dir
        self.info["jarvis_log"] = jarvis_log
        setup_jarvis_logging(
            role="core",
            log_dir=logs_dir,
            log_path=jarvis_log,
            level="INFO",
            console=True,
            use_queue=True,
        )
        self._logger = get_jarvis_logger("core").bind(module="Jarvis-HEP")
        try:
            self._logger.warning("\n" + render_logo_with_version())
            self._logger.warning("Jarvis-HEP V2 logging system initialized successful!")
            self._logger.info(f"Jarvis-HEP write into main log file -> {jarvis_log}")
        except Exception:
            # Logging must never block bootstrap; banner is best-effort.
            pass

    def load_task_yaml(self, path: str) -> dict[str, Any]:
        """Load task YAML and merge normalized layout into ``self.config``."""
        self.config = load_task_yaml(path)
        self.runtime = get_runtime_block(self.config)
        self._populate_info_from_config()
        return self.config

    def _populate_info_from_config(self) -> None:
        task_root = str(self.config.get("task_root") or os.getcwd())
        scan_name = str(self.config.get("scan_name") or "scan")
        logs_dir = os.path.join(task_root, "logs", scan_name)
        self.info = {
            "task_root": task_root,
            "task_result_dir": str(self.config.get("task_result_dir") or os.getcwd()),
            "scan_name": scan_name,
            "sampler_name": "SamplingVirtial",
            "run_id": str(self.config.get("run_id") or uuid.uuid4()),
            "task_yaml": self.config.get("task_yaml"),
            "logs_dir": logs_dir,
            "jarvis_log": os.path.join(logs_dir, f"{scan_name}.log"),
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
        likelihood_block = self.config.get("Likelihood") or {}
        sampling_block = self.config.get("Sampling") or {}
        likelihood_exprs = list(
            likelihood_block.get("expressions")
            or sampling_block.get("LogLikelihood")
            or []
        )
        sampler.set_execution_plan_template(
            calculator_modules=calc_modules,
            opera_modules=opera_modules,
            include_likelihood=bool(likelihood_exprs),
        )
        self.set_sampler(sampler)
        self.info["sampler_name"] = str(getattr(sampler, "method", type(sampler).__name__))
        return sampler

    def bootstrap_distributed_runtime(self) -> None:
        """Bring up Redis, Workers, and Archiver for a distributed run."""
        if not self.is_redis_runtime():
            raise RuntimeError("distributed runtime requires the internal Redis runtime")
        self.init_logger()
        self.init_command_parser()
        self.init_redis()
        self._claim_redis_control_lock()
        if self._resume_policy != "resume":
            self._reset_redis_for_fresh_run()
        self.init_sampler_from_config()
        self._init_sample_buckets()
        self.init_factory()
        db_path = os.path.join(
            str(self.info.get("task_result_dir") or os.getcwd()),
            "DATABASE",
            "samples.hdf5",
        )
        self.init_archiver(db_path)
        self.start_runtime_checkpoint()
        self._export_workflow_flowchart()

    @staticmethod
    def _load_check_module_rows(csv_path: str) -> tuple[list[dict[str, str]], list[str]]:
        with open(csv_path, "r", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            return list(reader), fieldnames

    def _build_check_module_samples(self) -> list[Sample]:
        if self.sampler is None:
            raise RuntimeError("sampler is not configured")
        csv_path = check_modules_points_path(self.config)
        rows, fieldnames = self._load_check_module_rows(csv_path)
        return build_check_module_samples(
            sampler=self.sampler,
            config=self.config,
            rows=rows,
            csv_fieldnames=fieldnames,
        )

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

    def run_adaptive_scan(self, *, timeout: float = 3600.0) -> int:
        """Drive a feedback sampler (AdaptiveLevelSet) through generation barriers."""
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
            pushed = int(self.sampler.run_adaptive(timeout=timeout))
        elif hasattr(self.sampler, "run_distributed"):
            pushed = int(self.sampler.run_distributed())
        else:
            raise RuntimeError(
                f"Sampler {type(self.sampler).__name__} does not implement run_adaptive"
            )
        # Archiver still drains the archive queue independently; wait for persisted rows.
        expected = self._archiver_records_written()
        # Best-effort: wait until archiver has caught up to total submitted.
        # Adaptive samplers already barriered on feedback; a short drain is enough.
        if self.archiver is not None:
            try:
                self.archiver.drain(idle_timeout=2.0)
            except Exception:
                pass
        self._finalize_sample_buckets()
        return requeued + pushed

    def run_check_modules(
        self,
        *,
        timeout: float = 120.0,
        verify_golden: Mapping[str, Any] | None = None,
    ) -> int:
        """Run the fixed-point calculator smoke path and return submitted sample count."""
        samples = self._build_check_module_samples()
        self.submit_samples(samples)
        self.wait_for_results(len(samples), timeout=timeout)
        self._finalize_sample_buckets()
        if verify_golden is not None:
            task_result_dir = str(self.info.get("task_result_dir") or os.getcwd())
            verify_check_modules_golden(task_result_dir=task_result_dir, golden=verify_golden)
        return len(samples)

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
            )
            if written:
                self.info["plot_scenes"] = written
                self._logger.info(
                    "plot scenes emitted → %s",
                    ", ".join(f"{key}={path}" for key, path in written.items()),
                )
                jplot = written.get("jplot_levelset")
                if jplot and "jplot_levelset_png" not in written:
                    self._logger.info(
                        "jplot YAML ready (render with: Jarvis2 plot %s)", jplot
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
        return RunOutcome.from_counters(
            submitted=int(submitted or 0),
            completed=int(counters.get("ok") or 0),
            failed=int(counters.get("failed") or 0),
            archived=archived,
            run_id=run_id,
            interrupted=interrupted,
            error=error,
            error_type=error_type,
            extras={
                "queued": int(counters.get("queued") or 0),
                "running": int(counters.get("running") or 0),
                "archive_q": int(counters.get("archive_q") or 0),
            },
        )

    def run(
        self,
        *,
        resume: bool = False,
        check_modules: bool = False,
        verify_golden: Mapping[str, Any] | None = None,
        write_run_summary: bool = True,
    ) -> RunOutcome:
        """Execute a distributed scan; return a truthful :class:`RunOutcome` (D11.1)."""
        self.prepare_resume(resume=resume, fresh=False)
        self._install_control_signal_handlers()
        submitted = 0
        outcome: RunOutcome | None = None
        try:
            self.bootstrap_distributed_runtime()
            method = sampling_method(self.config)
            if check_modules or is_check_modules_task(self.config):
                submitted = self.run_check_modules(verify_golden=verify_golden)
            elif method in STATELESS_METHODS:
                submitted = self.run_distributed_scan()
            elif method:
                # Stateful / feedback-driven methods (e.g. AdaptiveLevelSet).
                submitted = self.run_adaptive_scan()
            else:
                raise NotImplementedError(
                    "Unsupported task: configure Sampling.mode: check_modules, "
                    "Sampling.Method: Bridson|AdaptiveLevelSet, or pass --check-modules."
                )
            outcome = self._capture_run_outcome(submitted=submitted)
            if outcome.ok and not self._interrupt_requested:
                self._emit_plot_scenes()
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

    def check_modules(self, *, verify_golden: Mapping[str, Any] | None = None) -> RunOutcome:
        """CLI entry for ``Jarvis2 check <task>.yaml`` / ``--check-modules``."""
        return self.run(check_modules=True, verify_golden=verify_golden)

    def is_redis_runtime(self) -> bool:
        """Return True when the distributed Redis path should be used."""
        return str(self.runtime.get("mode", "auto")).strip().lower() == "redis"

    def checkpoint_file(self) -> str:
        task_root = str(
            self.info.get("task_root")
            or self.config.get("task_root")
            or os.getcwd()
        )
        scan_name = str(self.info.get("scan_name") or self.config.get("scan_name") or "scan")
        sampler_name = str(
            self.info.get("sampler_name")
            or (type(self.sampler).__name__ if self.sampler is not None else "sampler")
        )
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
            self._resume_policy = "fresh"
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
        self._resume_checkpoint_payload = payload

    def apply_resume_checkpoint(self, worker_config: Mapping[str, Any] | None = None) -> int:
        if self._resume_policy != "resume" or self._resume_checkpoint_payload is None:
            return 0
        if self.redis is None:
            raise RuntimeError("init_redis() must run before apply_resume_checkpoint()")
        drained = prepare_resume(self.redis, worker_config=worker_config)
        sampler_state = dict(self._resume_checkpoint_payload.get("sampler_state") or {})
        if self.sampler is not None and hasattr(self.sampler, "import_runtime_state"):
            self.sampler.import_runtime_state(sampler_state)
        if self.sampler is not None and hasattr(self.sampler, "set_resume_repropose_hint"):
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
        payload = build_payload(
            run_spec=run_spec,
            sampler_state=self.sampler.export_runtime_state(),
            persistence=persistence,
            reason=reason,
            safe_barrier_confirmed=True,
        )
        path = self.checkpoint_file()
        save_checkpoint(path, payload)
        return path

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        # EnvReqs.V2.redis (optional) overlays INTERNAL_REDIS_CONFIG defaults.
        redis_config = get_redis_config(self.config)
        # Prefer a Jarvis-managed redis-server so `ps` shows Jarvis-Redis:<scan>.
        # If something is already listening, we connect and leave its title alone.
        if client is None:
            self._ensure_managed_redis(redis_config)
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
            started = managed.ensure(scan_name=scan_name, work_dir=work_dir)
        except Exception as exc:
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

    def _control_lock_owner_id(self) -> str:
        run_id = str(self.info.get("run_id") or "jarvis2-run")
        return f"{os.uname().nodename}:{os.getpid()}:{run_id}"

    def _calculator_module_names(self) -> list[str]:
        modules = (self.config.get("Calculators") or {}).get("Modules") or []
        names: list[str] = []
        if isinstance(modules, list):
            for item in modules:
                if isinstance(item, Mapping):
                    name = str(item.get("name") or "").strip()
                    if name:
                        names.append(name)
        return names

    def _claim_redis_control_lock(self) -> None:
        """Refuse to start if another Jarvis2 still holds the Redis control lease."""
        if self.redis is None:
            raise RuntimeError("init_redis() must run before claiming the control lock")
        owner = self._control_lock_owner_id()
        if self.redis.claim_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
            self._control_lock_owner = owner
            self._logger.info("claimed Redis control lock (%s)", owner)
            return
        current = self.redis.get_control_lock_owner()
        raise RuntimeError(
            "another Jarvis2 instance already holds the Redis control lock "
            f"({current!r}). Stop residual Jarvis2/HEP2-Worker processes "
            "(do not leave runs suspended with ^Z), then retry. "
            "Emergency: redis-cli DEL hep:control:lock"
        )

    def _reset_redis_for_fresh_run(self) -> None:
        """Drop ephemeral queues/stats/calc pools before Workers are spawned."""
        if self.redis is None:
            return
        workers = int(self.runtime.get("workers", 0) or 0)
        result = self.redis.reset_run_ephemeral_keys(
            calculator_names=self._calculator_module_names(),
            worker_ids=max(workers, 32),
        )
        self._logger.info(
            "reset Redis run keys (deleted=%s sample_stats_reset=%s)",
            result.get("deleted_keys"),
            result.get("sample_stats_reset"),
        )

    def _release_redis_control_lock(self) -> None:
        owner = getattr(self, "_control_lock_owner", None)
        if not owner or self.redis is None:
            return
        try:
            self.redis.release_control_lock(str(owner))
            self._logger.info("released Redis control lock (%s)", owner)
        except Exception as exc:
            self._logger.warning("failed to release Redis control lock -> %s", exc)
        self._control_lock_owner = None

    def init_command_parser(self) -> CommandParser:
        """Run Phase-1 static command resolution for the loaded task config."""
        self.command_parser = build_command_parser(self.config)
        return self.command_parser

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
        return build_worker_config(
            self.config,
            task_result_dir=task_result_dir,
            parser=self.command_parser,
            calculator_modules=overrides.pop("calculator_modules", None),
            likelihood_expressions=overrides.pop("likelihood_expressions", None),
            opera_modules=overrides.pop("opera_modules", None),
            sample_dirs=overrides.pop("sample_dirs", None),
            extra=overrides or None,
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
        if self._resume_policy == "resume":
            prepare_resume(self.redis, worker_config=merged_config)
        self.factory.start_workers(workers, **merged_config)
        factory_cfg = get_factory_config(self.config)
        self.factory.start_monitor(update_hz=float(factory_cfg.get("monitor_hz", 120.0)))
        watchdog = get_watchdog_config(self.config)
        self.factory.start_watchdog(**watchdog)
        self._logger.info("TaskFactory started with %d worker(s)", workers)
        return self.factory

    def init_archiver(self, db_path: str | None = None) -> SimpleArchiver | ArchiverProcess:
        if self.redis is None:
            raise RuntimeError("init_redis() must run before init_archiver()")
        task_result_dir = str(
            self.info.get("task_result_dir")
            or self.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        sample_root = os.path.join(task_result_dir, "SAMPLE")
        os.makedirs(database_dir, exist_ok=True)
        os.makedirs(sample_root, exist_ok=True)
        resolved_db_path = db_path or os.path.join(database_dir, "samples.hdf5")
        archiver_config = get_archiver_config(self.config)
        delete_method = get_delete_method(self.config)
        redis_config = get_redis_config(self.config)
        if self.redis is not None:
            # Prefer the live connection if init_redis already bound a client.
            try:
                redis_config = dict(self.redis.connection_config())
            except Exception:
                pass

        scan_name = str(self.info.get("scan_name") or "").strip() or None
        if str(archiver_config.get("mode", "process")).strip().lower() == "process":
            self.archiver = ArchiverProcess(
                redis_config,
                db_path=resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
                scan_name=scan_name,
            )
            self.archiver.start()
        else:
            self.archiver = SimpleArchiver(
                self.redis,
                resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
            )
            self.archiver.start()
        self._restore_archiver_persistence()
        return self.archiver

    def _restore_archiver_persistence(self) -> None:
        if self._resume_policy != "resume" or self._resume_checkpoint_payload is None:
            return
        if self.archiver is None:
            return
        persistence = dict(self._resume_checkpoint_payload.get("persistence") or {})
        acked = persistence.get("acked_uuids") or []
        processor = getattr(self.archiver, "processor", None)
        if processor is not None:
            processor.acked_uuids.update(str(item) for item in acked)

    def _archiver_records_written(self) -> int:
        archiver = self.archiver
        if archiver is None:
            return 0
        counter = getattr(archiver, "records_written", 0)
        if hasattr(counter, "value"):
            return int(counter.value)
        return int(counter)

    def set_sampler(self, sampler: Any) -> None:
        self.sampler = sampler
        if self._resume_policy == "auto" and not getattr(self, "_resume_checkpoint_preloaded", False):
            self._preload_resume_checkpoint()
            self._resume_checkpoint_preloaded = True
        if self.redis is not None and hasattr(sampler, "set_redis"):
            sampler.set_redis(self.redis)
        if (
            self._resume_policy == "resume"
            and self._resume_checkpoint_payload is not None
            and hasattr(sampler, "import_runtime_state")
        ):
            sampler.import_runtime_state(
                self._resume_checkpoint_payload.get("sampler_state") or {}
            )
            if hasattr(sampler, "set_resume_repropose_hint"):
                sampler.set_resume_repropose_hint(True)

    def start_runtime_checkpoint(self) -> None:
        """Enable the 30 s checkpoint heartbeat once sampler + archiver exist."""
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
            progress = PermilleProgress(
                self._logger,
                total=total,
                label="samples archived",
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
                        f"Likely early SAMPLE bucket pack pruned dirs before Archiver "
                        f"wrote rows (fixed by packing only after archived==assigned)."
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
        return RunSummaryRenderer().write_outputs(summary, task_result_dir)

    def _init_sample_buckets(self) -> None:
        """Register Redis SAMPLE bucket meta (numbering + active/completed state)."""
        if self.redis is None:
            return
        sample_dir = get_sample_directory_config(self.config)
        task_result_dir = str(self.info.get("task_result_dir") or os.getcwd())
        sample_root = os.path.join(task_result_dir, "SAMPLE")
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
        self._logger.info(
            "SAMPLE buckets ready -> root=%s limit=%s width=%s pack=%s",
            sample_root,
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
        finally:
            self._stop_managed_redis()


__all__ = ["Jarvis2Core"]
