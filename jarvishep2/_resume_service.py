#!/usr/bin/env python3
"""ResumeService collaborator for Jarvis2Core (D25.3)."""

from __future__ import annotations

import os
import sys
import threading
from collections.abc import Mapping
from typing import Any

from jarvishep2.database import (
    read_persisted_outcome_counts,
    read_persisted_sample_index_state,
)
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
from jarvishep2.task_config import sampling_method


class _ResumeService:
    """Private Jarvis2Core collaborator (D25.3)."""

    def __init__(self, core: Any) -> None:
        self._core = core

    def prepare_resume(self, *, resume: bool = False, fresh: bool = False) -> None:
        """Preload checkpoint state before sampler bring-up (CLI ``--resume``)."""
        core = self._core

        core._preload_resume_checkpoint(resume=resume, fresh=fresh)
        core._resume_checkpoint_preloaded = True

    def _resolve_checkpoint_sampler_name(self) -> str:
        """Directory name under ``checkpoints/<scan>/``.

        Must match what a running scan writes (``sampler.method``, e.g. ``Dynesty``).
        ``prepare_resume`` runs **before** ``init_sampler_from_config``, so
        ``info["sampler_name"]`` is still the placeholder ``SamplingVirtial`` —
        fall through to YAML ``Sampling.Method`` instead of looking for a
        non-existent ``…/SamplingVirtial/state.pkl`` (which silently skipped
        nested engine restore and made Dynesty restart from scratch).
        """
        core = self._core

        placeholders = {"", "sampler", "samplingvirtial", "samplingvirtual"}
        candidates: list[str] = []
        if core.sampler is not None:
            method = getattr(core.sampler, "method", None)
            if method:
                candidates.append(str(method))
            candidates.append(type(core.sampler).__name__)
        info_name = str(core.info.get("sampler_name") or "")
        if info_name:
            candidates.append(info_name)
        try:
            method_from_yaml = sampling_method(core.config)
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
        core = self._core

        task_root = str(
            core.info.get("task_root")
            or core.config.get("task_root")
            or os.getcwd()
        )
        scan_name = str(core.info.get("scan_name") or core.config.get("scan_name") or "scan")
        sampler_name = core._resolve_checkpoint_sampler_name()
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
        core = self._core

        if fresh:
            core._resume_policy = "fresh"
            core._resume_checkpoint_payload = None
            return

        path = core.checkpoint_file()
        if not os.path.exists(path):
            # Explicit --resume also supports DATABASE-only recovery when the
            # checkpoint is missing or deliberately removed.
            core._resume_policy = "resume" if resume else "fresh"
            core._resume_checkpoint_payload = None
            return

        if not resume:
            if not core.prompt_resume_from_checkpoint(path):
                core._resume_policy = "fresh"
                core._resume_checkpoint_payload = None
                try:
                    os.remove(path)
                except OSError:
                    pass
                core._logger.warning(
                    "Starting a fresh run from user confirmation; existing checkpoint was discarded."
                )
                return
        core._resume_policy = "resume"
        try:
            payload = load_checkpoint(path)
        except ValueError as exc:
            core._logger.warning("Checkpoint rejected -> %s", exc)
            core._resume_checkpoint_payload = None
            core._resume_policy = "fresh"
            return
        # D22.5: refuse resume when Sampling.Mapper / Variables fingerprint drifts.
        ok_mapper, mapper_reason = check_mapper_fingerprint(payload, core.config)
        if not ok_mapper:
            core._logger.error("Checkpoint rejected (starting fresh) -> %s", mapper_reason)
            core._resume_checkpoint_payload = None
            core._resume_policy = "fresh"
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
            core._logger.error("Checkpoint rejected (starting fresh) -> %s", consts_reason)
            core._resume_checkpoint_payload = None
            core._resume_policy = "fresh"
            try:
                os.remove(path)
            except OSError:
                pass
            return
        if str(consts_reason).startswith("skip:"):
            core._logger.warning("%s", consts_reason)
        core._resume_checkpoint_payload = payload

    def _load_persisted_database_state(self, database_dir: str) -> None:
        """Refresh index prefix, UUID set, and completed/failed counts from DATABASE.

        Index prefix drives resume fast-forward; outcome counts seed Redis
        SAMPLE_STATS so a resume does not re-label durable Failed rows as
        successful (D21.13).
        """
        core = self._core

        (
            core._persisted_index_prefix,
            pending_indices,
            core._persisted_uuids,
        ) = read_persisted_sample_index_state(database_dir)
        completed, failed = read_persisted_outcome_counts(database_dir)
        total_from_status = completed + failed
        # Prefer status counts when any durable row exists; fall back to the
        # index/legacy UUID reconstruction for pre-status shards.
        if total_from_status > 0:
            core._persisted_records_count = total_from_status
            core._persisted_failed_count = failed
        else:
            core._persisted_records_count = (
                core._persisted_index_prefix
                + len(pending_indices)
                + len(core._persisted_uuids)
            )
            core._persisted_failed_count = 0

    def apply_resume_checkpoint(self, worker_config: Mapping[str, Any] | None = None) -> int:
        core = self._core

        if core._resume_policy != "resume":
            return 0
        if core.redis is None:
            raise RuntimeError("init_redis() must run before apply_resume_checkpoint()")
        drained = prepare_resume(
            core.redis,
            worker_config=worker_config,
            persisted_count=core._persisted_records_count,
            persisted_failed=core._persisted_failed_count,
        )
        imported = core._import_sampler_checkpoint_state(core.sampler)
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_uuids"):
            core.sampler.set_persisted_uuids(core._persisted_uuids)
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_prefix"):
            core.sampler.set_persisted_prefix(core._persisted_index_prefix)
        if (
            bool(getattr(core.sampler, "uses_indexed_resume", False))
            and hasattr(core.sampler, "advance_to_persisted_prefix")
        ):
            core.sampler.advance_to_persisted_prefix(core._persisted_index_prefix)
        if (
            imported
            and core.sampler is not None
            and not bool(getattr(core.sampler, "uses_indexed_resume", False))
            and hasattr(core.sampler, "set_resume_repropose_hint")
        ):
            core.sampler.set_resume_repropose_hint(True)
        core._logger.info(
            "Resumed from checkpoint; drained %d stale task(s) from hep:task_queue",
            drained,
        )
        return drained

    def save_runtime_checkpoint(self, *, reason: str = "") -> str | None:
        core = self._core

        if core.sampler is None or not hasattr(core.sampler, "export_runtime_state"):
            return None
        persistence: dict[str, Any] = {}
        if core.archiver is not None and hasattr(core.archiver, "persistence_state"):
            persistence = dict(core.archiver.persistence_state())
        prefix = max(
            core._persisted_index_prefix,
            int(persistence.get("persisted_prefix", 0) or 0),
        )
        core._persisted_index_prefix = prefix
        if hasattr(core.sampler, "set_persisted_prefix"):
            core.sampler.set_persisted_prefix(prefix)
        task_root = str(core.info.get("task_root") or core.config.get("task_root") or os.getcwd())
        task_result_dir = str(
            core.info.get("task_result_dir")
            or core.config.get("task_result_dir")
            or task_root
        )
        scan_name = str(core.info.get("scan_name") or core.config.get("scan_name") or "scan")
        sampler_name = str(core.info.get("sampler_name") or type(core.sampler).__name__)
        run_spec = build_run_spec(
            config=core.config,
            scan_name=scan_name,
            task_root=task_root,
            task_result_dir=task_result_dir,
            sampler_name=sampler_name,
            worker_parallel=int(core.runtime.get("workers", 0) or 0),
        )
        indexed = bool(getattr(core.sampler, "uses_indexed_resume", False))
        if indexed:
            # ``prefix`` is the sole hot-path completion fact.  It arrives
            # from the Archiver after fsync and is independent of scan size.
            # Do not let ``checkpoint_runtime_state`` capture the current
            # generator state here: only a batch-boundary snapshot may become
            # a resume point.
            safe = False
        else:
            submitted = frozenset(getattr(core.sampler, "submitted_uuids", frozenset()))
            at_barrier = bool(getattr(core.sampler, "at_safe_barrier", lambda: False)())
            safe = at_barrier and submitted <= frozenset(
                getattr(core.sampler, "_persisted_uuids", set())
            )
        # Nested samplers override checkpoint_runtime_state to always re-export
        # the live dynesty engine (ignore stale FeedbackSampler barrier copies).
        if hasattr(core.sampler, "checkpoint_runtime_state"):
            # Interrupt / dynesty_engine_checkpoint: always force a live export
            # for methods that do not use generation barriers.
            force_live = str(reason or "") in {
                "interrupt",
                "dynesty_engine_checkpoint",
                "dynesty_finished",
            } or bool(getattr(core.sampler, "method", "") in {"Dynesty", "MultiNest"})
            sampler_state = core.sampler.checkpoint_runtime_state(
                safe=True if force_live else safe
            )
        else:
            sampler_state = core.sampler.export_runtime_state()
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
        path = core.checkpoint_file()
        save_checkpoint(path, payload)
        return path

    def _save_interrupt_checkpoint(self) -> str | None:
        """Force a resume checkpoint on SIGINT/SIGTERM (D8.3 human remainder).

        Agent-facing pieces (run_state interrupted, stop-ack) stay parked with D8.
        Failures are logged and never re-raised — teardown must still run.
        """
        core = self._core

        try:
            path = core.save_runtime_checkpoint(reason="interrupt")
            if path:
                try:
                    core._logger.warning(
                        "interrupt checkpoint written → %s "
                        "(resume: Jarvis run <task.yaml> --resume)",
                        path,
                    )
                except Exception:
                    pass
            return path
        except Exception as exc:
            try:
                core._logger.warning("interrupt checkpoint failed → %s", exc)
            except Exception:
                pass
            return None

    def _discard_resume_checkpoint(
        self,
        *,
        reason: str,
        level: str = "error",
    ) -> None:
        """Drop a rejected resume snapshot and force a fresh sampling start."""
        core = self._core

        log = core._logger.error if level == "error" else core._logger.warning
        log("Checkpoint sampler state rejected (starting fresh) -> %s", reason)
        core._resume_policy = "fresh"
        core._resume_checkpoint_payload = None
        path = core.checkpoint_file()
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
        core = self._core

        target = sampler if sampler is not None else core.sampler
        if core._resume_policy != "resume":
            return False
        if core._resume_checkpoint_payload is None or target is None:
            return False
        if not hasattr(target, "import_runtime_state"):
            return False
        sampler_state = core._resume_checkpoint_payload.get("sampler_state") or {}
        if not isinstance(sampler_state, Mapping):
            sampler_state = {}
        try:
            target.import_runtime_state(sampler_state)
        except ValueError as exc:
            # MCMC Bounds drift (num_chains / num_iters / scale / seed) and
            # similar sampler-shape mismatches: refuse the stale snapshot.
            core._discard_resume_checkpoint(reason=str(exc), level="error")
            return False
        if (
            not bool(getattr(target, "uses_indexed_resume", False))
            and hasattr(target, "set_resume_repropose_hint")
        ):
            target.set_resume_repropose_hint(True)
        return True

    def start_runtime_checkpoint(self) -> None:
        """Enable the configured checkpoint heartbeat once sampler and archiver exist."""
        core = self._core

        if core.sampler is None or not hasattr(core.sampler, "configure_checkpoint"):
            return
        core.sampler.configure_checkpoint(
            checkpoint_file=core.checkpoint_file(),
            save_callback=lambda reason="": core.save_runtime_checkpoint(reason=reason),
        )
