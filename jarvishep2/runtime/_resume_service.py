#!/usr/bin/env python3
"""ResumeService collaborator for Jarvis2Core (D25.3)."""

from __future__ import annotations

import os
import pickle
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
    check_sampling_fingerprint,
    uses_checkpoint_resume_protocol,
    check_operas_constants_fingerprint,
    checkpoint_path,
    load_checkpoint,
    save_checkpoint,
)
from jarvishep2.task_config import sampling_method
from jarvishep2.runtime_config import get_checkpoint_config
from jarvishep2.sampling.sample_checkpoint import (
    FORMAT as SAMPLE_CHECKPOINT_FORMAT,
    load_runtime,
    uses_mcmc_checkpoint,
    uses_sample_checkpoint,
)


class _ResumeService:
    """Private Jarvis2Core collaborator (D25.3)."""

    def __init__(self, core: Any) -> None:
        self._core = core

    def prepare_resume(self, *, resume: bool = False, fresh: bool = False) -> None:
        """Preload checkpoint state before sampler bring-up (CLI ``--resume``)."""
        core = self._core
        core._resume_checkpoint_source = None

        from jarvishep2.sampler_catalog import is_plugin
        if not fresh and is_plugin(sampling_method(core.config)):
            self._prepare_plugin_resume(resume=resume)
        else:
            core._preload_resume_checkpoint(resume=resume, fresh=fresh)
        core._resume_checkpoint_preloaded = True

    def _prepare_plugin_resume(self, *, resume: bool) -> None:
        """Never downgrade a refused plug-in resume to fresh sampling."""
        from jarvis_hep_sampler.errors import PluginResumeMismatch
        from jarvishep2.unit_runtime import preflight_units

        core = self._core
        path = core.checkpoint_file()
        exists = os.path.isfile(path)
        if exists and not resume and not core.prompt_resume_from_checkpoint(path):
            raise PluginResumeMismatch("existing plug-in scan cannot be restarted in place; use a new scan name")
        core._resume_policy = "resume" if resume or exists else "fresh"
        core._resume_checkpoint_payload = None
        if core._resume_policy == "resume":
            from jarvishep2.sampler_catalog import plugin_record

            record = plugin_record(sampling_method(core.config))
            if record is not None and not {"checkpoint", "resume"}.issubset(record.hooks):
                raise PluginResumeMismatch("sampler does not support resume", method=record.method)
            if not get_checkpoint_config(core.config)["enabled"]:
                raise PluginResumeMismatch("cannot resume with checkpoint.enabled=false")
            if not exists:
                legacy = os.path.join(os.path.dirname(path), "state.pkl")
                detail = f"legacy checkpoint {legacy} is unsupported" if os.path.isfile(legacy) else "checkpoint missing"
                raise PluginResumeMismatch(f"{detail}: expected {path}; use a new scan name")
        preflight_units(core, resume=core._resume_policy == "resume",
                        check=getattr(core, "_unit_check_run", False) is True)
        if core._resume_policy == "resume":
            try:
                payload = load_runtime(path, core.config)
                ok, reason = check_operas_constants_fingerprint(payload)
                if not ok:
                    raise ValueError(reason)
                core._resume_checkpoint_payload = payload
                core._resume_checkpoint_source = path
            except (ValueError, OSError, EOFError, pickle.UnpicklingError) as exc:
                raise PluginResumeMismatch(f"checkpoint rejected at {path}: {exc}; use a new scan name") from exc
        # This validates both source/configuration and the checkpoint before
        # DATABASE recovery, old cursors, or resource bootstrap can write.
        core.init_logger()
        core.preflight_sampler_plugin()

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
            check=bool(core.config.get("_check_modules_sample_layout")),
            sampler_pickle=getattr(core.sampler, "uses_sample_checkpoint", False) is True,
        )

    def _preload_sample_checkpoint(self, *, resume: bool, fresh: bool) -> None:
        core = self._core
        core._resume_policy = "fresh"
        core._resume_checkpoint_payload = None
        core._resume_checkpoint_source = None
        if fresh:
            return
        if not get_checkpoint_config(core.config)["enabled"]:
            if resume:
                raise ValueError("cannot resume with checkpoint.enabled=false")
            return
        path = core.checkpoint_file()
        legacy = False
        if (uses_mcmc_checkpoint(core.config) and not core.config.get("_unit_execution")
                and not os.path.exists(path)):
            legacy_path = os.path.join(os.path.dirname(path), "state.pkl")
            if os.path.isfile(legacy_path):
                path = legacy_path
                legacy = True
        if not os.path.isfile(path):
            if resume:
                raise FileNotFoundError(f"sampler checkpoint missing: {path}")
            return
        if not resume and not core.prompt_resume_from_checkpoint(path):
            os.remove(path)
            return
        if legacy:
            payload = load_checkpoint(path)
            if (payload.get("integrity") or {}).get("unit_execution"):
                raise ValueError("unit-aware MCMC resume requires runtime.pickle")
            for check in (check_mapper_fingerprint, check_sampling_fingerprint):
                ok, reason = check(payload, core.config)
                if not ok:
                    raise ValueError(reason)
        else:
            payload = load_runtime(path, core.config)
        if uses_mcmc_checkpoint(core.config):
            ok, reason = check_operas_constants_fingerprint(payload)
            if not ok:
                raise ValueError(reason)
            if reason.startswith("skip:"):
                core._logger.warning("%s", reason)
        core._resume_checkpoint_payload = payload
        core._resume_checkpoint_source = path
        core._resume_policy = "resume"

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

        from jarvishep2.unit_runtime import preflight_units
        preflight_units(core, resume=resume or (not fresh and os.path.isfile(core.checkpoint_file())),
                        check=getattr(core, "_unit_check_run", False) is True)

        if uses_sample_checkpoint(core.config):
            self._preload_sample_checkpoint(resume=resume, fresh=fresh)
            return

        if fresh:
            core._resume_policy = "fresh"
            core._resume_checkpoint_payload = None
            return

        protocol = uses_checkpoint_resume_protocol(core.config)
        if protocol and not get_checkpoint_config(core.config)["enabled"]:
            if resume:
                core._logger.warning("Resume requested with checkpoint.enabled: false; starting fresh")
            core._resume_policy = "fresh"
            core._resume_checkpoint_payload = None
            return

        path = core.checkpoint_file()
        if not os.path.exists(path):
            if resume and protocol:
                core._logger.warning("Resume requested without state.pkl; starting fresh")
            core._resume_policy = "fresh" if protocol else ("resume" if resume else "fresh")
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
        except (ValueError, OSError, EOFError, pickle.UnpicklingError) as exc:
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
        ok_sampling, sampling_reason = check_sampling_fingerprint(payload, core.config)
        if protocol and not ok_sampling:
            core._logger.error("Checkpoint rejected (starting fresh) -> %s", sampling_reason)
            core._resume_checkpoint_payload = None
            core._resume_policy = "fresh"
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
        """Set the resume cursor, UUID set, and completed/failed counts.

        Runtime-pickle samplers restore progress without scanning DATABASE.
        Their completion counts describe sampler feedback. Legacy and check
        paths can still read DATABASE to seed the Redis SAMPLE_STATS counters.
        """
        core = self._core

        payload = core._resume_checkpoint_payload
        if uses_sample_checkpoint(core.config):
            # Seed run statistics from sampler-confirmed completions only.
            # These are not claims about which output rows are on disk.
            state = (payload or {}).get("sampler_state") or {}
            completed = int(state.get("completed_count", 0))
            failed = int(state.get("failed_count", 0))
            core._persisted_index_prefix = 0
            core._persisted_pending_indices = set()
            core._persisted_uuids = set()
            core._persisted_records_count = completed + failed
            core._persisted_failed_count = failed
            return
        if (
            payload is not None
            and uses_checkpoint_resume_protocol(core.config)
            and not core.config.get("_check_modules_sample_layout")
        ):
            state = payload.get("sampler_state") or {}
            cursor = max(0, int(state.get("checkpoint_index", state.get("accepted_index", 0)) or 0))
            counts = (payload.get("persistence") or {}).get("outcome_counts") or {}
            completed = max(0, int(counts.get("completed", 0) or 0))
            failed = max(0, int(counts.get("failed", 0) or 0))
            core._persisted_index_prefix = cursor
            core._persisted_pending_indices = set()
            core._persisted_uuids = set()
            core._persisted_records_count = completed + failed
            core._persisted_failed_count = failed
            return

        (
            core._persisted_index_prefix,
            core._persisted_pending_indices,
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
                + len(core._persisted_pending_indices)
                + len(core._persisted_uuids)
            )
            core._persisted_failed_count = 0

    def save_runtime_checkpoint(self, *, reason: str = "") -> str | None:
        core = self._core

        if (getattr(core.sampler, "uses_sample_checkpoint", False) is True
                or uses_sample_checkpoint(core.config)):
            # Signals may interrupt a proposal halfway through an in-memory
            # transition. Keep the last complete periodic snapshot on interrupt.
            if core.sampler is None or reason == "interrupt":
                return None
            written = core.sampler.save_checkpoint(force=True)
            return core.checkpoint_file() if written else None

        if core.sampler is None or not hasattr(core.sampler, "export_runtime_state"):
            return None
        from jarvishep2.runtime_config import get_checkpoint_config

        if not get_checkpoint_config(getattr(core, "config", None))["enabled"]:
            return None
        # Counters include evaluations above the saved cursor; a resume
        # evaluates those again and DATABASE keeps both rows.
        counters = core._live_sample_counters()
        persistence: dict[str, Any] = {
            "outcome_counts": {
                "completed": int(counters.get("ok", 0) or 0),
                "failed": int(counters.get("failed", 0) or 0),
            }
        }
        prefix = core._persisted_index_prefix
        confirm = getattr(core.sampler, "confirm_archived_prefix", None)
        confirmed = confirm() if callable(confirm) else None
        if confirmed is not None:
            prefix = max(prefix, int(confirmed))
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
            # ``prefix`` is the sole hot-path completion fact.  Archive
            # receipts confirm it after fsync; the cost is the in-flight window.
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
        if hasattr(core.sampler, "checkpoint_runtime_state"):
            # Feedback engines need their live state on interrupt.
            # Indexed samplers must retain the newest batch snapshot proven
            # durable; their live cursor can include unarchived proposals.
            force_live = not indexed and str(reason or "") == "interrupt"
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
        sampler_log = getattr(core, "_sampler_log", None)
        if sampler_log is not None:
            sampler_log.checkpoint_saved(path, reason)
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
        core._persisted_index_prefix = 0
        core._persisted_pending_indices = set()
        core._persisted_uuids = set()
        core._persisted_records_count = 0
        core._persisted_failed_count = 0
        if core.sampler is not None and hasattr(core.sampler, "set_resume_mode"):
            core.sampler.set_resume_mode(False)
        path = core.checkpoint_file()
        try:
            if os.path.isfile(path):
                os.remove(path)
        except OSError:
            pass

    def _import_sampler_checkpoint_state(self, sampler: Any | None = None) -> bool:
        """Import sampler_state from the preloaded checkpoint, or start fresh.

        Shared by ``set_sampler``, late Method resolution, and bootstrap.
        On Bounds/shape mismatch (``ValueError``)
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
        if getattr(target, "uses_sample_checkpoint", False) is True:
            state = core._resume_checkpoint_payload["sampler_state"]
            from jarvishep2.sampler_catalog import is_plugin

            if is_plugin(sampling_method(core.config)):
                target.restore_checkpoint(state)
                return True
            if uses_mcmc_checkpoint(core.config):
                # A rejected MCMC snapshot must not start a new chain under
                # an existing scan name or remove the file being diagnosed.
                target.restore_checkpoint(
                    state, legacy=core._resume_checkpoint_payload.get("format") != SAMPLE_CHECKPOINT_FORMAT,
                )
                return True
            prepare_points = getattr(target, "prepare_resume_points", None)
            try:
                if callable(prepare_points):
                    prepare_points(core.checkpoint_file(), state)
                target.restore_checkpoint(state)
            except (ValueError, OSError) as exc:
                core._discard_resume_checkpoint(reason=str(exc), level="error")
                return False
            return True
        if not hasattr(target, "import_runtime_state"):
            return False
        sampler_state = core._resume_checkpoint_payload.get("sampler_state") or {}
        if not isinstance(sampler_state, Mapping):
            sampler_state = {}
        try:
            prepare_points = getattr(target, "prepare_resume_points", None)
            if callable(prepare_points):
                prepare_points(core.checkpoint_file(), sampler_state)
            target.import_runtime_state(sampler_state)
        except (ValueError, OSError) as exc:
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

        if uses_sample_checkpoint(core.config):
            policy = get_checkpoint_config(core.config)
            if getattr(core, "checkpoint_process", None) is not None:
                return
            if (core.sampler is not None and policy["enabled"]
                    and not core.config.get("_check_modules_sample_layout")
                    and getattr(core.sampler, "supports_checkpoint", True) is not False):
                from jarvishep2.runtime.checkpoint_process import CheckpointProcess

                try:
                    sampler_log = getattr(core, "_sampler_log", None)
                    writer = CheckpointProcess(
                        core.checkpoint_file(), core._logger, interval=policy["heartbeat"],
                        scan_name=str(core.info.get("scan_name") or core.config.get("scan_name") or ""),
                        redis_config=core.redis.connection_config() if core.redis is not None else None,
                        control_owner=str(core._control_lock_owner or ""),
                        logs_dir=core.info.get("logs_dir"),
                        on_saved=(sampler_log.checkpoint_saved if sampler_log is not None
                                  and getattr(core.sampler, "is_generational_plugin", False) else None),
                    )
                except Exception as exc:
                    core._logger.warning("Checkpoint process unavailable -> %s", exc)
                    return
                core.checkpoint_process = writer
                writer.start()
                core.sampler.configure_checkpoint(checkpoint_file=core.checkpoint_file(), writer=writer)
            return

        if core.sampler is None or not hasattr(core.sampler, "configure_checkpoint"):
            return
        core.sampler.configure_checkpoint(
            checkpoint_file=core.checkpoint_file(),
            save_callback=lambda reason="": core.save_runtime_checkpoint(reason=reason),
        )
        if (core._resume_policy == "fresh" and get_checkpoint_config(core.config)["enabled"]
                and not core.config.get("_check_modules_sample_layout")
                and uses_checkpoint_resume_protocol(core.config)):
            prepare_points = getattr(core.sampler, "prepare_fresh_points", None)
            if callable(prepare_points):
                prepare_points(core.checkpoint_file())
            core.save_runtime_checkpoint(reason="initial")
