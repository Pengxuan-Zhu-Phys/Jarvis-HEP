#!/usr/bin/env python3
"""Redis-backed Worker process for Jarvis-HEP V2."""

from __future__ import annotations

import concurrent.futures
import json
import os
import signal
import threading
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.async_subprocess import AsyncSubprocessScheduler, SubprocessRuntimeConfig
from jarvishep2.command_parser import CommandParser
from jarvishep2.calculator_modes import mode_info, shared_mode_groups
from jarvishep2.env_setup import EnvCapture, resolve_env_setup_sources
from jarvishep2.expression import ExpressionContext
from jarvishep2.archive_handoff import list_product_names, resolve_staging_dir, stage_sample_dir
from jarvishep2.file_ops import DEFAULT_DELETE_METHOD, normalize_delete_method
from jarvishep2.file_operation_service import FileOperationService
from jarvishep2.likelihood import LogLikelihoodEvaluator
from jarvishep2.Module.calculator import CalculatorModule
from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging
from jarvishep2.mapper import build_mapper
from jarvishep2.mp_context import get_spawn_context
from jarvishep2.operas import preload_operas
from jarvishep2.operas_functions import (
    build_operas_expression_context,
    operas_expression_functions_required,
)
from jarvishep2.redis_queue import RedisQueue, SHARED_HELD_PREFIX
from jarvishep2.sample import Sample, materialize_failure_artifacts
from jarvishep2.sample import ExecutionStep
from jarvishep2.workflow import group_by_layer, max_layer_width, resolve_module_layers

_SPAWN_CTX = get_spawn_context()
Process = _SPAWN_CTX.Process


class Worker(Process):
    """Long-lived process that pulls one Sample at a time from Redis.

    ``__init__`` accepts either a live :class:`RedisQueue` (from the Factory) or
    a picklable connection mapping. Only the connection **settings** are stored
    for spawn; the child process builds its own Redis client in ``run()``.
    """

    def __init__(
        self,
        worker_id: int,
        redis: RedisQueue | Mapping[str, Any],
        worker_config: Mapping[str, Any],
    ) -> None:
        super().__init__(name=f"Jarvis-Worker-{worker_id}", daemon=False)
        self.worker_id = int(worker_id)
        self.redis_config = RedisQueue.extract_connection_config(redis)
        self.worker_config = dict(worker_config)
        self._redis: RedisQueue | None = None
        self._mapper = None
        self._operas: dict[str, Any] = {}
        self._calculators: dict[str, CalculatorModule] = {}
        self._shared_mode_sets: dict[str, list[str]] = {}
        self._likelihood: LogLikelihoodEvaluator | None = None
        self._nuisance_profiler = None  # Profile1DProfiler | None
        self._expression_context: ExpressionContext | None = None
        self._scheduler: AsyncSubprocessScheduler | None = None
        self._command_parser: CommandParser | None = None
        self._delete_method = DEFAULT_DELETE_METHOD
        self._file_ops: FileOperationService | None = None
        self._staging_dir = ""
        self._handoff_to_staging = False
        self._sample_buckets_enabled = True
        self._observables_lock: threading.Lock | None = None
        self._is_running = True
        self._current_sample_uuid: str | None = None
        self._current_task: dict[str, Any] | None = None
        self._held_calc_packs: dict[str, str] = {}
        # Created child-side in _init_redis: locks/events are not spawn-picklable.
        self._heartbeat_lock: threading.Lock | None = None
        self._last_status = "starting"
        self._heartbeat_stop: threading.Event | None = None
        self._heartbeat_thread: threading.Thread | None = None

    def _handle_signal(self, signum: int, _frame: Any) -> None:
        self._is_running = False

    def _hb_lock(self) -> threading.Lock:
        """Heartbeat-state lock, created lazily child-side (not spawn-picklable).

        First touch happens single-threaded (``_init_redis`` in the child, or
        the test harness driving Worker methods in-process), so lazy creation
        is race-free.
        """
        lock = self._heartbeat_lock
        if lock is None:
            lock = threading.Lock()
            self._heartbeat_lock = lock
        return lock

    def _init_redis(self) -> None:
        """Connect to Redis in the child process (spawn-safe)."""
        self._hb_lock()
        self._redis = RedisQueue(self.redis_config)
        self._redis.connect()
        self._redis.heartbeat(
            str(self.worker_id),
            status="idle",
            pid=self.pid,
        )

    def _init_runtime(self) -> None:
        self._observables_lock = threading.Lock()
        self._delete_method = normalize_delete_method(
            self.worker_config.get("delete_method", DEFAULT_DELETE_METHOD)
        )
        # Dedicated FileOperation process (save/copy/delete). Inline mode for tests.
        file_ops_mode = str(
            self.worker_config.get("file_operation_mode")
            or self.worker_config.get("file_ops_mode")
            or "process"
        ).strip().lower()
        if file_ops_mode not in {"process", "inline"}:
            file_ops_mode = "process"
        scan_name = str(self.worker_config.get("scan_name") or "").strip()
        if not scan_name:
            sample_config = self.worker_config.get("sample_config") or {}
            if isinstance(sample_config, Mapping):
                scan_name = str(sample_config.get("scan_name") or "").strip()
        self._file_ops = FileOperationService.start(
            mode=file_ops_mode,
            delete_method=self._delete_method,
            scan_name=scan_name or None,
        )
        if "handoff_to_staging" in self.worker_config:
            self._handoff_to_staging = bool(self.worker_config.get("handoff_to_staging"))
        else:
            cleanup = self.worker_config.get("cleanup_config") or {}
            strategy = str(cleanup.get("strategy", "direct")).strip().lower()
            self._handoff_to_staging = strategy == "mv_to_staging"
            archiver_cfg = self.worker_config.get("archiver_config") or {}
            if str(archiver_cfg.get("handoff", "direct")).strip().lower() == "staging":
                self._handoff_to_staging = True
        # Never mkdir staging/ for direct handoff (default).
        self._staging_dir = str(self.worker_config.get("staging_dir") or "").strip()
        if self._handoff_to_staging and not self._staging_dir:
            sample_config = self.worker_config.get("sample_config") or {}
            task_result_dir = str(sample_config.get("task_result_dir") or "").strip()
            if task_result_dir:
                self._staging_dir = resolve_staging_dir(task_result_dir, create=True)
        sample_dir_cfg = self.worker_config.get("sample_directory") or {}
        if isinstance(sample_dir_cfg, dict):
            self._sample_buckets_enabled = bool(sample_dir_cfg.get("enabled", True))
        self._mapper = build_mapper(self.worker_config.get("mapper"))
        expression_payload = {
            "opera_modules": self.worker_config.get("opera_modules"),
            "calculator_modules": self.worker_config.get("calculator_modules"),
            "likelihood_expressions": self.worker_config.get("likelihood_expressions"),
        }
        if operas_expression_functions_required(expression_payload):
            expression_context = build_operas_expression_context(required=True)
        else:
            expression_context = ExpressionContext()
        opera_configs = self.worker_config.get("opera_modules") or {}
        if isinstance(opera_configs, list):
            opera_configs = {
                str(item.get("name", f"Operas{index}")): item
                for index, item in enumerate(opera_configs)
            }
        self._operas = preload_operas(
            opera_configs,
            expression_context=expression_context,
        )
        calculator_configs = self.worker_config.get("calculator_modules") or []
        if isinstance(calculator_configs, dict):
            calculator_configs = list(calculator_configs.values())
        self._calculators = CalculatorModule.from_config_list(calculator_configs)
        self._shared_mode_sets = shared_mode_groups(calculator_configs)
        layer_width = max(
            1,
            int(self.worker_config.get("subprocess_max_concurrency", 0) or 0),
            max_layer_width(resolve_module_layers(calculator_configs)),
        )
        self._scheduler = AsyncSubprocessScheduler(
            SubprocessRuntimeConfig(
                max_concurrency=layer_width,
                log_policy="quiet",
                progress_interval_sec=3600.0,
            ),
            logger=get_jarvis_logger(
                "worker",
                worker_id=self.worker_id,
            ),
        )
        self._scheduler.start()
        self._command_parser = CommandParser.from_picklable(self.worker_config.get("command_parser"))
        for module in self._calculators.values():
            module.attach_scheduler(self._scheduler)
            module.attach_command_parser(self._command_parser)
            module.attach_expression_context(expression_context)
            module.attach_file_ops(self._file_ops)
            sources = resolve_env_setup_sources(
                module.env_setup,
                command_parser=self._command_parser,
            )
            if sources:
                module.bind_env(EnvCapture.merged_env(sources))
        likelihood_exprs = self.worker_config.get("likelihood_expressions") or []
        self._likelihood = LogLikelihoodEvaluator(
            likelihood_exprs,
            expression_context=expression_context,
        )
        self._expression_context = expression_context
        nuisance_cfg = self.worker_config.get("nuisance_config")
        if isinstance(nuisance_cfg, Mapping) and nuisance_cfg:
            from jarvishep2.Module.profile1d import Profile1DProfiler

            # Synthetic top-level config so from_config sees Nuisance block.
            self._nuisance_profiler = Profile1DProfiler.from_config(
                {"Nuisance": dict(nuisance_cfg)},
                context=expression_context,
            )
        else:
            self._nuisance_profiler = None

    def _heartbeat(self, status: str) -> None:
        if self._redis is None:
            return
        now = time.time()
        # Snapshot mutable fields under a lock: layer threads mutate
        # _held_calc_packs concurrently and the periodic heartbeat thread
        # reads them.
        with self._hb_lock():
            self._last_status = status
            held_packs = dict(self._held_calc_packs)
            current_task_ref = self._current_task
        current_task = ""
        if current_task_ref is not None:
            current_task = self._redis.encode_task_for_heartbeat(current_task_ref)
        active_pids: list[int] = []
        if self._scheduler is not None:
            active_pids = self._scheduler.active_subprocess_pids()
        file_operation_pid = self._file_ops.pid if self._file_ops is not None else None
        if file_operation_pid is not None:
            active_pids.append(file_operation_pid)
        self._redis.heartbeat(
            str(self.worker_id),
            status=status,
            pid=self.pid,
            current_sample=self._current_sample_uuid,
            last_heartbeat=now,
            ts=now,
            held_calc_packs=json.dumps(held_packs),
            active_subprocess_pids=json.dumps(active_pids),
            file_operation_pid=file_operation_pid,
            current_task=current_task,
        )

    def _heartbeat_loop(self, stop: threading.Event, interval_sec: float) -> None:
        """Refresh the heartbeat while long calculator commands run.

        Without this, a Worker only beats at acquire/release/pull, so any
        ``module.execute()`` longer than ``Watchdog.stale_sec`` is falsely
        recovered — fatal now that PackID slots alias shadow directories.
        """
        while not stop.wait(interval_sec):
            try:
                self._heartbeat(self._last_status)
                expected_owner = str(
                    self.worker_config.get("control_lock_owner") or ""
                ).strip()
                if (
                    expected_owner
                    and self._redis is not None
                    and self._redis.get_control_lock_owner() != expected_owner
                ):
                    self._is_running = False
                    return
            except Exception:  # pragma: no cover - heartbeat must never kill the Worker
                continue

    def _start_heartbeat_thread(self) -> None:
        interval = float(self.worker_config.get("heartbeat_interval_sec", 5.0) or 5.0)
        if interval <= 0:
            return
        self._heartbeat_stop = threading.Event()
        self._heartbeat_thread = threading.Thread(
            target=self._heartbeat_loop,
            args=(self._heartbeat_stop, max(0.1, interval)),
            name="JarvisWorkerHeartbeat",
            daemon=True,
        )
        self._heartbeat_thread.start()

    def _stop_heartbeat_thread(self) -> None:
        if self._heartbeat_stop is not None:
            self._heartbeat_stop.set()
        if self._heartbeat_thread is not None:
            self._heartbeat_thread.join(timeout=2.0)
        self._heartbeat_stop = None
        self._heartbeat_thread = None

    def _shutdown_runtime(self, worker_log: Any) -> None:
        """Best-effort child cleanup; one subsystem must never skip the next."""
        try:
            self._stop_heartbeat_thread()
        except Exception as exc:
            worker_log.warning("heartbeat shutdown failed -> %s", exc)
        try:
            self._heartbeat("stopped")
        except Exception as exc:
            worker_log.warning("final heartbeat failed -> %s", exc)
        if self._scheduler is not None:
            try:
                self._scheduler.shutdown(wait=True)
            except Exception as exc:
                # Cleanup is deliberately independent: a scheduler timeout
                # must never skip FileOperation termination.
                worker_log.warning("scheduler shutdown failed -> %s", exc)
            finally:
                self._scheduler = None
        if self._file_ops is not None:
            try:
                self._file_ops.shutdown()
            except Exception as exc:
                worker_log.warning("FileOperation shutdown failed -> %s", exc)
            finally:
                self._file_ops = None
        if self._redis is not None:
            try:
                self._redis.close()
            except Exception as exc:
                worker_log.warning("Worker Redis close failed -> %s", exc)
            finally:
                self._redis = None

    def _run_calculator_step(
        self,
        step_name: str,
        sample: Sample,
        *,
        selection_checked: bool = False,
    ) -> None:
        module = self._calculators.get(step_name)
        if module is None:
            raise KeyError(f"unknown calculator module '{step_name}'")
        if self._redis is None:
            raise RuntimeError("redis is not initialized in worker process")
        # Selection must precede slot acquisition and runtime preparation. In
        # shared mode, preparing a skipped step would rebuild and relabel a
        # PackID for a mode which never actually ran.
        if not selection_checked and not module.should_run(sample.info):
            return

        timeout = int(self.worker_config.get("calc_acquire_timeout", 30))
        info = mode_info(module.config)
        shared_parent = info[0] if info is not None else None
        target_mode = info[1] if info is not None else None
        held_key = SHARED_HELD_PREFIX + shared_parent if shared_parent else step_name
        if shared_parent and target_mode:
            acquired = self._redis.acquire_shared_calc(
                shared_parent,
                target_mode,
                modes=self._shared_mode_sets.get(shared_parent, [target_mode]),
                timeout=timeout,
                affinity_wait_sec=float(
                    self.worker_config.get("shared_mode_affinity_wait_sec", 3.0)
                ),
            )
            pack_id = acquired[0] if acquired is not None else None
        else:
            pack_id = self._redis.acquire_calc(step_name, timeout=timeout)
        if pack_id is None:
            raise TimeoutError(f"timed out acquiring calculator slot for '{step_name}'")
        with self._hb_lock():
            self._held_calc_packs[held_key] = pack_id
        self._heartbeat("busy")
        runtime_ready = False
        try:
            module.acquire_pack_id(pack_id)
            module.prepare_runtime(sample.info)
            # A successful preparation wrote a matching mode stamp (or reused
            # one).  Execution failure afterwards still leaves a valid pack.
            runtime_ready = True
            updated = module.execute(sample.info, runtime_prepared=True)
            lock = self._observables_lock
            if lock is not None:
                with lock:
                    self._merge_calculator_observables(sample, step_name, pack_id, updated)
            else:
                self._merge_calculator_observables(sample, step_name, pack_id, updated)
        finally:
            # Hard release: Redis socket errors must not leave the slot busy.
            self._force_release_pack(
                held_key,
                pack_id,
                shared_parent=shared_parent,
                ready_mode=target_mode if runtime_ready else None,
            )
            self._heartbeat("busy")

    def _force_release_pack(
        self,
        step_name: str,
        pack_id: str | None = None,
        *,
        shared_parent: str | None = None,
        ready_mode: str | None = None,
    ) -> None:
        """Return a calculator PackID, retaining local ownership on Redis errors."""
        if self._redis is None:
            return
        with self._hb_lock():
            held_pack = pack_id or self._held_calc_packs.get(step_name)
        if not held_pack:
            return
        try:
            # False is the benign already-free/double-release result; an
            # exception means Redis could not confirm the transition.
            if shared_parent:
                self._redis.release_shared_calc(
                    shared_parent, str(held_pack), ready_mode,
                )
            elif str(step_name).startswith(SHARED_HELD_PREFIX):
                self._redis.force_release_shared_calc(
                    str(step_name).removeprefix(SHARED_HELD_PREFIX), str(held_pack),
                )
            else:
                self._redis.force_release_calc(step_name, str(held_pack))
        except Exception as exc:
            try:
                get_jarvis_logger("worker", worker_id=self.worker_id).error(
                    "calculator slot release failed; retaining held pack %s/%s for watchdog sweep -> %s",
                    step_name,
                    held_pack,
                    exc,
                )
            except Exception:
                pass
            return
        with self._hb_lock():
            if self._held_calc_packs.get(step_name) == str(held_pack):
                self._held_calc_packs.pop(step_name, None)

    def _force_release_all_held_packs(self, logger: Any | None = None) -> None:
        """Safety net for process_task finally — release every slot this worker holds."""
        if self._redis is None:
            return
        with self._hb_lock():
            held = dict(self._held_calc_packs)
        for step_name, pack_id in held.items():
            self._force_release_pack(str(step_name), str(pack_id))

    def _merge_calculator_observables(
        self,
        sample: Sample,
        step_name: str,
        pack_id: str,
        updated: Mapping[str, Any],
    ) -> None:
        sample.merge_observables(updated)
        sample.record_pack_id(step_name, pack_id)

    def _force_serial_layers(self) -> bool:
        """Rollback switch: run same-layer calculators one after another."""
        return bool(self.worker_config.get("force_serial_layers", False))

    def _run_shared_mode_group(self, steps: list[ExecutionStep], sample: Sample) -> None:
        """Run one parent's modes serially, greedily preferring a warm mode."""
        pending = [
            step
            for step in steps
            if step.name in self._calculators
            and self._calculators[step.name].should_run(sample.info)
        ]
        while pending:
            first = self._calculators.get(pending[0].name)
            info = mode_info(first.config) if first is not None else None
            if info is None or self._redis is None:
                selected = pending.pop(0)
            else:
                parent, _ = info
                pending_modes = [
                    mode_info(self._calculators[step.name].config)[1]
                    for step in pending
                    if step.name in self._calculators
                    and mode_info(self._calculators[step.name].config) is not None
                ]
                counts = self._redis.shared_mode_free_counts(
                    parent,
                    pending_modes,
                )
                if any(counts.values()):
                    selected = max(
                        pending,
                        key=lambda step: (
                            counts.get(
                                (mode_info(self._calculators[step.name].config) or ("", ""))[1],
                                0,
                            ),
                            -pending.index(step),
                        ),
                    )
                else:
                    # No warm pack exists yet. Spread the initial builds across
                    # modes already claimed by peer Workers rather than sending
                    # every cold Worker into the first YAML mode.
                    busy = self._redis.shared_mode_busy_counts(parent, pending_modes)
                    selected = min(
                        pending,
                        key=lambda step: (
                            busy.get(
                                (mode_info(self._calculators[step.name].config) or ("", ""))[1],
                                0,
                            ),
                            pending.index(step),
                        ),
                    )
                pending.remove(selected)
            self._run_calculator_step(
                selected.name, sample, selection_checked=True,
            )

    def _run_calculator_steps(self, calc_steps: list[ExecutionStep], sample: Sample) -> None:
        """Fan out independent calculators, but serialize shared-pack siblings."""
        shared: dict[str, list[ExecutionStep]] = {}
        units: list[tuple[str, list[ExecutionStep]]] = []
        for step in calc_steps:
            module = self._calculators.get(step.name)
            info = mode_info(module.config) if module is not None else None
            if info is None:
                units.append((step.name, [step]))
            else:
                shared.setdefault(info[0], []).append(step)
        units.extend((f"{SHARED_HELD_PREFIX}{parent}", steps) for parent, steps in shared.items())
        if len(units) > 1 and not self._force_serial_layers():
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(units)) as pool:
                futures = [
                    pool.submit(
                        self._run_shared_mode_group if key.startswith(SHARED_HELD_PREFIX) else self._run_calculator_step,
                        steps if key.startswith(SHARED_HELD_PREFIX) else steps[0].name,
                        sample,
                    )
                    for key, steps in units
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()
            return
        for key, steps in units:
            if key.startswith(SHARED_HELD_PREFIX):
                self._run_shared_mode_group(steps, sample)
            else:
                self._run_calculator_step(steps[0].name, sample)

    def _run_layer(self, layer: list[ExecutionStep], sample: Sample) -> None:
        """Run one execution-plan layer; fan out same-layer calculators concurrently."""
        calc_steps = [step for step in layer if step.type == "calculator"]
        inline_steps = [step for step in layer if step.type != "calculator"]
        self._run_calculator_steps(calc_steps, sample)

        for step in inline_steps:
            if step.type == "opera":
                self._run_opera_step(step.name, sample)
            elif step.type == "nuisance_optimize":
                self._run_nuisance_optimize(sample)
            elif step.type == "likelihood":
                self._run_likelihood(sample)

    def _rerun_physics_pipeline(self, sample: Sample) -> None:
        """Re-run calculator + opera steps (nuisance profile attempts)."""
        plan = list(sample.execution_plan or [])
        layers = group_by_layer(plan)
        for layer in layers:
            calc_steps = [step for step in layer if step.type == "calculator"]
            opera_steps = [step for step in layer if step.type == "opera"]
            if calc_steps:
                self._run_calculator_steps(calc_steps, sample)
            for step in opera_steps:
                self._run_opera_step(step.name, sample)

    def _run_nuisance_optimize(self, sample: Sample) -> None:
        """Profile nuisance parameters (D13.4 Profile1D) and merge results."""
        profiler = self._nuisance_profiler
        if profiler is None:
            return
        from jarvishep2.Module.profile1d import ProfileProbeResult

        logger = sample.child_logger(
            module=f"{sample.info.get('logger_name', f'Sample@{sample.uuid}')} "
            f"(Nuisance:Profile1D)"
        )
        var_name = profiler.var_name

        def evaluate(z: float) -> ProfileProbeResult:
            # Inject free nuisance parameter into dual-truth fields.
            sample.params[var_name] = float(z)
            sample.observables[var_name] = float(z)
            sample._mirror_fields_to_info()
            if profiler.re_run_physics:
                self._rerun_physics_pipeline(sample)
            payload = dict(sample.params)
            payload.update(sample.observables)
            payload[var_name] = float(z)
            terms = profiler.logl_registry.evaluate_all(payload)
            total = float(sum(terms.values())) if terms else 0.0
            # Pass conditions may reference nuisance LogL term names.
            pass_payload = dict(payload)
            pass_payload.update(terms)
            pass_terms = profiler.pass_registry.evaluate_all(pass_payload)
            pass_ok = all(pass_terms.values()) if pass_terms else True
            if logger is not None:
                try:
                    logger.info(
                        "nuisance probe %s=%.6g LogL=%.6g pass=%s",
                        var_name,
                        z,
                        total,
                        pass_ok,
                    )
                except Exception:
                    pass
            return ProfileProbeResult(
                z=float(z),
                logl=total,
                terms=dict(terms),
                pass_ok=pass_ok,
                pass_terms=dict(pass_terms),
                observables=dict(sample.observables),
            )

        result = profiler.optimize(evaluate)
        # Apply best nuisance state onto the sample.
        sample.params[var_name] = float(result.best_z)
        sample.observables[var_name] = float(result.best_z)
        sample.observables.update({k: float(v) for k, v in result.best_terms.items()})
        sample.observables["NuisanceLogL"] = float(result.best_logl)
        sample.observables["nuisance_pass"] = bool(result.pass_ok)
        if not isinstance(sample.info, dict):
            sample.info = {}
        sample.info["nuisance"] = {
            "method": "Profile1D",
            "var": var_name,
            "best": float(result.best_z),
            "best_logl": float(result.best_logl),
            "pass": bool(result.pass_ok),
            "pass_terms": dict(result.pass_terms),
            "n_attempts": int(result.n_attempts),
            "status": result.status,
            "mode": result.mode,
            "history_z": list(result.history_z),
            "history_logl": list(result.history_logl),
        }
        sample._mirror_fields_to_info()
        if not result.pass_ok:
            # Soft-fail: mark Failed when pass conditions reject the profiled point.
            sample.status = "Failed"
            if isinstance(sample.info, dict):
                sample.info["error"] = (
                    f"nuisance pass conditions failed for {var_name}={result.best_z}"
                )
                sample.info["error_type"] = "NuisancePassCondition"
                sample.info["failed_module"] = "NuisanceOptimize"
        if logger is not None:
            try:
                logger.info(
                    "nuisance profile done %s=%.6g LogL=%.6g status=%s attempts=%d",
                    var_name,
                    result.best_z,
                    result.best_logl,
                    result.status,
                    result.n_attempts,
                )
            except Exception:
                pass

    def _run_opera_step(self, step_name: str, sample: Sample) -> None:
        module = self._operas.get(step_name)
        if module is None:
            raise KeyError(f"unknown opera module '{step_name}'")
        updated = module.execute(
            sample.observables,
            sample.info,
            sample_logger=sample.child_logger(
                module=f"{sample.info.get('logger_name', f'Sample@{sample.uuid}')} "
                f"(Operas:{step_name})"
            ),
        )
        sample.merge_observables(updated)

    def _run_likelihood(self, sample: Sample) -> None:
        if self._likelihood is None:
            return
        if not isinstance(sample.info, dict):
            sample.info = {}
        sample._mirror_fields_to_info()
        value = self._likelihood.calculate(sample.info)
        # Likelihood may rebind info["observables"] / params; adopt then re-share.
        sample.pull_dual_truth_from_info()
        sample._likelihood = value

    def _collect_cleanup_paths(self, sample: Sample) -> list[str]:
        if not isinstance(sample.info, dict):
            return []
        paths: list[str] = []
        for key in ("cleanup_paths", "staging_paths"):
            raw = sample.info.get(key)
            if isinstance(raw, str) and raw.strip():
                paths.append(raw.strip())
            elif isinstance(raw, list):
                paths.extend(str(item).strip() for item in raw if str(item).strip())
        return paths

    def _cleanup_transient_paths(self, sample: Sample) -> None:
        paths = self._collect_cleanup_paths(sample)
        if not paths:
            return
        if self._file_ops is not None:
            self._file_ops.delete(paths, missing_ok=True)
            return
        from jarvishep2.file_ops import delete_paths

        delete_paths(paths, method=self._delete_method, missing_ok=True)

    def _handoff_sample_to_staging(self, sample: Sample) -> None:
        """Fast metadata-only move of materialized work dirs into staging (WP-D4.1)."""
        if not self._handoff_to_staging or not self._staging_dir:
            return
        if not isinstance(sample.info, dict):
            return
        save_dir = str(sample.info.get("save_dir") or "").strip()
        if not save_dir or not os.path.isdir(save_dir):
            return

        # Logger must already be closed (sample.close writes SUMMARY first).
        sample.close_logger()
        staging_path = stage_sample_dir(save_dir, self._staging_dir, sample.uuid)
        sample.record_handoff(
            staging_path=staging_path,
            product_list=list_product_names(staging_path),
            clear_save_dir=True,
        )

    def _stage_and_submit(self, sample: Sample) -> None:
        if self._redis is None:
            return
        info = sample.to_info_dict()
        self._redis.submit_result(info)
        if bool(self.worker_config.get("publish_feedback", False)):
            from jarvishep2.feedback_return import build_feedback_record

            # Ensure unusable samples still carry LogL=-inf for flat feedback.
            if "LogL" not in sample.observables:
                sample.observables["LogL"] = float("-inf")
                sample._mirror_fields_to_info()
            spec = self.worker_config.get("feedback_return") or {
                "mode": "minimal",
                "include_logl": True,
                "fields": [],
            }
            self._redis.publish_feedback(
                build_feedback_record(sample, spec=spec)
            )

    def process_task(self, task: Mapping[str, Any]) -> None:
        """Core pipeline: rebuild Sample, execute workflow, submit result."""
        payload = dict(task)
        # Assign SAMPLE/<bucket>/<uuid> before materialize (Redis-owned numbering).
        if self._sample_buckets_enabled and self._redis is not None and not payload.get("bucket_dir"):
            allocation = None
            last_exc: Exception | None = None
            # High-worker Dynesty can stampede the bucket lock; retry before dying.
            for attempt in range(6):
                try:
                    allocation = self._redis.allocate_sample_bucket()
                    last_exc = None
                    break
                except RuntimeError as exc:
                    last_exc = exc
                    if "bucket lock" not in str(exc).lower() or attempt >= 5:
                        raise
                    time.sleep(0.02 * (attempt + 1))
            if last_exc is not None:
                raise last_exc
            if allocation is not None:
                payload["bucket_id"] = allocation["bucket_id"]
                payload["bucket_dir"] = allocation["bucket_dir"]
                payload["bucket_name"] = allocation["bucket_name"]
                payload["_bucket_parent"] = allocation["bucket_dir"]

        self._current_task = payload
        sample = Sample.from_task_dict(payload)
        self._current_sample_uuid = sample.uuid
        self._heartbeat("busy")
        top = get_jarvis_logger(
            "worker",
            worker_id=self.worker_id,
        ).bind(sample_uuid=sample.uuid)
        sample_config = dict(self.worker_config.get("sample_config") or {})
        # Stamp bucket parent before set_config so materialize lands under SAMPLE/<bucket>/<uuid>.
        if payload.get("bucket_dir"):
            sample_config["_bucket_parent"] = payload["bucket_dir"]
            sample_config["bucket_dir"] = payload["bucket_dir"]
        if payload.get("bucket_id") is not None:
            sample_config["bucket_id"] = payload["bucket_id"]
        if payload.get("bucket_name"):
            sample_config["bucket_name"] = payload["bucket_name"]
        sample.set_config(sample_config)
        sample.start()
        try:
            if self._mapper is not None:
                sample.bind_params(self._mapper)
            elif sample.opera_params:
                sample.adopt_params(sample.opera_params, as_observables=True)
            sample.materialize(worker_id=str(self.worker_id))
            delay_sec = float(self.worker_config.get("test_process_delay_sec", 0) or 0)
            if delay_sec > 0:
                time.sleep(delay_sec)
            layers = group_by_layer(sample.execution_plan)
            _profile_path = os.environ.get("JARVIS2_WORKER_PROFILE")
            _profile_started = time.monotonic()
            for layer in layers:
                self._run_layer(layer, sample)
            if _profile_path:
                with open(_profile_path, "a", encoding="utf-8") as _profile_handle:
                    _profile_handle.write(
                        json.dumps(
                            {
                                "worker_id": self.worker_id,
                                "layers": [
                                    {
                                        "calc_steps": [
                                            step.name
                                            for step in layer
                                            if step.type == "calculator"
                                        ],
                                        "width": len(layer),
                                    }
                                    for layer in layers
                                ],
                                "elapsed_sec": time.monotonic() - _profile_started,
                                "scheduler": (
                                    self._scheduler.snapshot()
                                    if self._scheduler is not None
                                    else {}
                                ),
                            }
                        )
                        + "\n"
                    )
            sample.set_status("Completed")
        except Exception as exc:
            sample.record_failure(exc)
            materialize_failure_artifacts(sample.info, error=exc)
            # Multi-line command failures embed cwd/cmd/stderr; keep them intact.
            top.error("sample failed; see sample log ->\n%s", exc)
        finally:
            # Always return calculator PackIDs first — even if handoff/archive fails.
            self._force_release_all_held_packs(logger=top)
            try:
                # Write Sample SUMMARY while the file logger is still open, then stage.
                sample.close()
                self._handoff_sample_to_staging(sample)
                self._stage_and_submit(sample)
                # Last Redis lifecycle touch for this sample's bucket (active→completed).
                # When sealed + idle, Redis enqueues the bucket for Archiver tar packing.
                bucket_id = None
                if isinstance(sample.info, dict):
                    bucket_id = sample.info.get("bucket_id")
                if bucket_id is None and isinstance(self._current_task, dict):
                    bucket_id = self._current_task.get("bucket_id")
                if self._redis is not None and bucket_id is not None:
                    try:
                        self._redis.finish_sample_bucket(bucket_id)
                    except Exception as bucket_exc:
                        top.error("bucket finish failed for %s -> %s", bucket_id, bucket_exc)
                self._cleanup_transient_paths(sample)
            except Exception as cleanup_exc:
                top.error("sample cleanup failed after release -> %s", cleanup_exc)
            self._current_sample_uuid = None
            with self._hb_lock():
                self._current_task = None

    def _main_loop(self) -> None:
        assert self._redis is not None
        pull_timeout = int(self.worker_config.get("pull_timeout", 5))
        while self._is_running:
            task = self._redis.pull_task(timeout=pull_timeout)
            if task is None:
                self._heartbeat("idle")
                continue
            self._heartbeat("busy")
            self.process_task(task)
            self._heartbeat("idle")

    def run(self) -> None:
        from jarvishep2.proc_title import set_process_title, worker_title

        # Prefer explicit worker_config.scan_name; fall back to sample_config.
        scan_name = str(self.worker_config.get("scan_name") or "").strip()
        if not scan_name:
            sample_config = self.worker_config.get("sample_config") or {}
            if isinstance(sample_config, dict):
                scan_name = str(sample_config.get("scan_name") or "").strip()
        set_process_title(worker_title(self.worker_id, scan_name=scan_name or None))
        # Always logs/<scan>/worker-NN.log under the project (never cwd jarvis_worker_PID).
        from jarvishep2.logging import component_log_path, scan_logs_dir

        logs_dir = str(self.worker_config.get("logs_dir") or "").strip()
        if not logs_dir:
            sample_config = self.worker_config.get("sample_config") or {}
            if not isinstance(sample_config, dict):
                sample_config = {}
            task_root = str(
                self.worker_config.get("task_root")
                or sample_config.get("task_root")
                or ""
            ).strip()
            task_result_dir = str(sample_config.get("task_result_dir") or "").strip()
            if not task_root and task_result_dir:
                trd = os.path.abspath(task_result_dir)
                parent = os.path.dirname(trd)
                if os.path.basename(parent) == "outputs":
                    task_root = os.path.dirname(parent)
                else:
                    task_root = parent
            if not task_root:
                task_root = os.getcwd()
            scan_for_logs = scan_name or "scan"
            logs_dir = scan_logs_dir(task_root, scan_for_logs)
        log_path = component_log_path(logs_dir, "worker", worker_id=self.worker_id)
        silence = bool(self.worker_config.get("log_silence", False))
        console_level = str(self.worker_config.get("console_level") or "WARNING")
        setup_jarvis_logging(
            role="worker",
            component="worker",
            worker_id=self.worker_id,
            console=not silence,
            console_level=console_level,
            silence=silence,
            use_queue=True,
            scan_logs_dir=logs_dir,
            log_path=log_path,
        )
        worker_log = get_jarvis_logger(
            "worker",
            worker_id=self.worker_id,
        )
        try:
            worker_log.info(
                "Worker process started pid=%s scan=%s log_path=%s",
                self.pid,
                scan_name or "?",
                log_path,
            )
        except Exception:
            pass
        signal.signal(signal.SIGTERM, self._handle_signal)
        signal.signal(signal.SIGINT, self._handle_signal)
        self._init_redis()
        try:
            # Keep initialization inside the cleanup boundary: FileOperation
            # starts early in _init_runtime, and any later setup exception must
            # still shut it down.
            self._init_runtime()
            self._heartbeat("starting")
            self._start_heartbeat_thread()
            self._main_loop()
        finally:
            self._shutdown_runtime(worker_log)


__all__ = ["Worker"]
