#!/usr/bin/env python3
"""Physics pipeline for one Sample, extracted from Worker (D25.5)."""

from __future__ import annotations

import concurrent.futures
import json
import os
import time
from collections.abc import Mapping
from typing import Any

from jarvishep2.calculator_modes import mode_info
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import SHARED_HELD_PREFIX
from jarvishep2.sample import Sample
from jarvishep2.sample import ExecutionStep
from jarvishep2.workflow import group_by_layer


def stamp_task_sample_bucket(
    redis: Any,
    payload: dict[str, Any],
    *,
    enabled: bool,
) -> dict[str, Any]:
    """Allocate SAMPLE/<bucket> and stamp the task payload (D25.5 layout helper)."""
    if not enabled or redis is None or payload.get("bucket_dir"):
        return payload
    allocation = None
    last_exc: Exception | None = None
    for attempt in range(6):
        try:
            allocation = redis.allocate_sample_bucket()
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
    return payload


class SampleExecutor:
    """Mapper + layered calculators / operas / likelihood / nuisance."""

    def __init__(self, worker: Any) -> None:
        self._worker = worker

    def process(self, sample: Sample) -> None:
        """Run mapper, layers, calculators, operas, likelihood, and nuisance."""
        worker = self._worker
        if worker._mapper is not None:
            sample.bind_params(worker._mapper)
        elif sample.opera_params:
            sample.adopt_params(sample.opera_params, as_observables=True)
        sample.apply_mcmc_identity()
        sample.materialize(worker_id=str(worker.worker_id))
        delay_sec = float(worker.worker_config.get("test_process_delay_sec", 0) or 0)
        if delay_sec > 0:
            time.sleep(delay_sec)
        layers = group_by_layer(sample.execution_plan)
        profile_path = os.environ.get("JARVIS2_WORKER_PROFILE")
        profile_started = time.monotonic()
        for layer in layers:
            worker._run_layer(layer, sample)
        if profile_path:
            with open(profile_path, "a", encoding="utf-8") as profile_handle:
                profile_handle.write(
                    json.dumps(
                        {
                            "worker_id": worker.worker_id,
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
                            "elapsed_sec": time.monotonic() - profile_started,
                            "scheduler": (
                                worker._scheduler.snapshot()
                                if worker._scheduler is not None
                                else {}
                            ),
                        }
                    )
                    + "\n"
                )
        sample.set_status("Completed")

    def _run_calculator_step(
        self,
        step_name: str,
        sample: Sample,
        *,
        selection_checked: bool = False,
    ) -> None:
        worker = self._worker
    
        module = worker._calculators.get(step_name)
        if module is None:
            raise KeyError(f"unknown calculator module '{step_name}'")
        if worker._redis is None:
            raise RuntimeError("redis is not initialized in worker process")
        # Selection must precede slot acquisition and runtime preparation. In
        # shared mode, preparing a skipped step would rebuild and relabel a
        # PackID for a mode which never actually ran.
        if not selection_checked and not module.should_run(sample.info):
            return
    
        timeout = int(worker.worker_config.get("calc_acquire_timeout", 30))
        info = mode_info(module.config)
        shared_parent = info[0] if info is not None else None
        target_mode = info[1] if info is not None else None
        held_key = SHARED_HELD_PREFIX + shared_parent if shared_parent else step_name
        if shared_parent and target_mode:
            acquired = worker._redis.acquire_shared_calc(
                shared_parent,
                target_mode,
                modes=worker._shared_mode_sets.get(shared_parent, [target_mode]),
                timeout=timeout,
                affinity_wait_sec=float(
                    worker.worker_config.get("shared_mode_affinity_wait_sec", 3.0)
                ),
            )
            pack_id = acquired[0] if acquired is not None else None
        else:
            pack_id = worker._redis.acquire_calc(step_name, timeout=timeout)
        if pack_id is None:
            raise TimeoutError(f"timed out acquiring calculator slot for '{step_name}'")
        with worker._hb_lock():
            worker._held_calc_packs[held_key] = pack_id
        worker._heartbeat("busy")
        runtime_ready = False
        try:
            module.acquire_pack_id(pack_id)
            module.prepare_runtime(sample.info)
            # A successful preparation wrote a matching mode stamp (or reused
            # one).  Execution failure afterwards still leaves a valid pack.
            runtime_ready = True
            updated = module.execute(sample.info, runtime_prepared=True)
            lock = worker._observables_lock
            if lock is not None:
                with lock:
                    worker._merge_calculator_observables(sample, step_name, pack_id, updated)
            else:
                worker._merge_calculator_observables(sample, step_name, pack_id, updated)
        finally:
            # Hard release: Redis socket errors must not leave the slot busy.
            worker._force_release_pack(
                held_key,
                pack_id,
                shared_parent=shared_parent,
                ready_mode=target_mode if runtime_ready else None,
            )
            worker._heartbeat("busy")

    def _force_release_pack(
        self,
        step_name: str,
        pack_id: str | None = None,
        *,
        shared_parent: str | None = None,
        ready_mode: str | None = None,
    ) -> None:
        """Return a calculator PackID, retaining local ownership on Redis errors."""
        worker = self._worker
    
        if worker._redis is None:
            return
        with worker._hb_lock():
            held_pack = pack_id or worker._held_calc_packs.get(step_name)
        if not held_pack:
            return
        try:
            # False is the benign already-free/double-release result; an
            # exception means Redis could not confirm the transition.
            if shared_parent:
                worker._redis.release_shared_calc(
                    shared_parent, str(held_pack), ready_mode,
                )
            elif str(step_name).startswith(SHARED_HELD_PREFIX):
                worker._redis.force_release_shared_calc(
                    str(step_name).removeprefix(SHARED_HELD_PREFIX), str(held_pack),
                )
            else:
                worker._redis.force_release_calc(step_name, str(held_pack))
        except Exception as exc:
            try:
                get_jarvis_logger("worker", worker_id=worker.worker_id).error(
                    "calculator slot release failed; retaining held pack %s/%s for watchdog sweep -> %s",
                    step_name,
                    held_pack,
                    exc,
                )
            except Exception:
                pass
            return
        with worker._hb_lock():
            if worker._held_calc_packs.get(step_name) == str(held_pack):
                worker._held_calc_packs.pop(step_name, None)

    def _force_release_all_held_packs(self, logger: Any | None = None) -> None:
        """Safety net for process_task finally — release every slot this worker holds."""
        worker = self._worker
    
        if worker._redis is None:
            return
        with worker._hb_lock():
            held = dict(worker._held_calc_packs)
        for step_name, pack_id in held.items():
            worker._force_release_pack(str(step_name), str(pack_id))

    def _merge_calculator_observables(
        self,
        sample: Sample,
        step_name: str,
        pack_id: str,
        updated: Mapping[str, Any],
    ) -> None:
        worker = self._worker
    
        sample.merge_observables(updated)
        sample.record_pack_id(step_name, pack_id)

    def _force_serial_layers(self) -> bool:
        """Rollback switch: run same-layer calculators one after another."""
        worker = self._worker
    
        return bool(worker.worker_config.get("force_serial_layers", False))

    def _run_shared_mode_group(self, steps: list[ExecutionStep], sample: Sample) -> None:
        """Run one parent's modes serially, greedily preferring a warm mode."""
        worker = self._worker
    
        pending = [
            step
            for step in steps
            if step.name in worker._calculators
            and worker._calculators[step.name].should_run(sample.info)
        ]
        while pending:
            first = worker._calculators.get(pending[0].name)
            info = mode_info(first.config) if first is not None else None
            if info is None or worker._redis is None:
                selected = pending.pop(0)
            else:
                parent, _ = info
                pending_modes = [
                    mode_info(worker._calculators[step.name].config)[1]
                    for step in pending
                    if step.name in worker._calculators
                    and mode_info(worker._calculators[step.name].config) is not None
                ]
                counts = worker._redis.shared_mode_free_counts(
                    parent,
                    pending_modes,
                )
                if any(counts.values()):
                    selected = max(
                        pending,
                        key=lambda step: (
                            counts.get(
                                (mode_info(worker._calculators[step.name].config) or ("", ""))[1],
                                0,
                            ),
                            -pending.index(step),
                        ),
                    )
                else:
                    # No warm pack exists yet. Spread the initial builds across
                    # modes already claimed by peer Workers rather than sending
                    # every cold Worker into the first YAML mode.
                    busy = worker._redis.shared_mode_busy_counts(parent, pending_modes)
                    selected = min(
                        pending,
                        key=lambda step: (
                            busy.get(
                                (mode_info(worker._calculators[step.name].config) or ("", ""))[1],
                                0,
                            ),
                            pending.index(step),
                        ),
                    )
                pending.remove(selected)
            worker._run_calculator_step(
                selected.name, sample, selection_checked=True,
            )

    def _run_calculator_steps(self, calc_steps: list[ExecutionStep], sample: Sample) -> None:
        """Fan out independent calculators, but serialize shared-pack siblings."""
        worker = self._worker
    
        shared: dict[str, list[ExecutionStep]] = {}
        units: list[tuple[str, list[ExecutionStep]]] = []
        for step in calc_steps:
            module = worker._calculators.get(step.name)
            info = mode_info(module.config) if module is not None else None
            if info is None:
                units.append((step.name, [step]))
            else:
                shared.setdefault(info[0], []).append(step)
        units.extend((f"{SHARED_HELD_PREFIX}{parent}", steps) for parent, steps in shared.items())
        if len(units) > 1 and not worker._force_serial_layers():
            with concurrent.futures.ThreadPoolExecutor(max_workers=len(units)) as pool:
                futures = [
                    pool.submit(
                        worker._run_shared_mode_group if key.startswith(SHARED_HELD_PREFIX) else worker._run_calculator_step,
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
                worker._run_shared_mode_group(steps, sample)
            else:
                worker._run_calculator_step(steps[0].name, sample)

    def _run_layer(self, layer: list[ExecutionStep], sample: Sample) -> None:
        """Run one execution-plan layer; fan out same-layer calculators concurrently."""
        worker = self._worker
    
        calc_steps = [step for step in layer if step.type == "calculator"]
        inline_steps = [step for step in layer if step.type != "calculator"]
        worker._run_calculator_steps(calc_steps, sample)
    
        for step in inline_steps:
            if step.type == "opera":
                worker._run_opera_step(step.name, sample)
            elif step.type == "nuisance_optimize":
                worker._run_nuisance_optimize(sample)
            elif step.type == "likelihood":
                worker._run_likelihood(sample)

    def _rerun_physics_pipeline(self, sample: Sample) -> None:
        """Re-run calculator + opera steps (nuisance profile attempts)."""
        worker = self._worker
    
        plan = list(sample.execution_plan or [])
        layers = group_by_layer(plan)
        for layer in layers:
            calc_steps = [step for step in layer if step.type == "calculator"]
            opera_steps = [step for step in layer if step.type == "opera"]
            if calc_steps:
                worker._run_calculator_steps(calc_steps, sample)
            for step in opera_steps:
                worker._run_opera_step(step.name, sample)

    def _run_nuisance_optimize(self, sample: Sample) -> None:
        """Profile nuisance parameters (D13.4 Profile1D) and merge results."""
        worker = self._worker
    
        profiler = worker._nuisance_profiler
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
                worker._rerun_physics_pipeline(sample)
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
        worker = self._worker
    
        module = worker._operas.get(step_name)
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
        worker = self._worker
    
        if worker._likelihood is None:
            return
        if not isinstance(sample.info, dict):
            sample.info = {}
        sample._mirror_fields_to_info()
        value = worker._likelihood.calculate(sample.info)
        # Likelihood may rebind info["observables"] / params; adopt then re-share.
        sample.pull_dual_truth_from_info()
        sample._likelihood = value
