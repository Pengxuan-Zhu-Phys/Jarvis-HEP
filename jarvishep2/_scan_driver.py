#!/usr/bin/env python3
"""ScanDriver collaborator for Jarvis2Core (D25.3)."""

from __future__ import annotations

import csv
import os
import time
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from jarvishep2.check_modules import (
    build_check_module_samples,
    verify_check_modules_golden,
)
from jarvishep2.log_kv import PermilleProgress, format_duration
from jarvishep2.runtime_config import get_runtime_block, get_sample_directory_config
from jarvishep2.sample import Sample
from jarvishep2.task_config import (
    check_modules_n_samples,
    check_modules_timeout_sec,
    resolve_check_modules_csv,
    sampling_method,
)


class _ScanDriver:
    """Private Jarvis2Core collaborator (D25.3)."""

    def __init__(self, core: Any) -> None:
        self._core = core

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
        core = self._core

        if core.sampler is None:
            raise RuntimeError("sampler is not configured")
        csv_path, raw_spec = resolve_check_modules_csv(core.config)
        if csv_path is not None:
            core._logger.warning(
                "Jarvis check: using fixed points from CSV -> %s (configured as %r)",
                csv_path,
                raw_spec,
            )
            rows, fieldnames = core._load_check_module_rows(csv_path)
            samples = build_check_module_samples(
                sampler=core.sampler,
                config=core.config,
                rows=rows,
                csv_fieldnames=fieldnames,
            )
            core._logger.warning(
                "Jarvis check: loaded %d fixed point(s) from CSV", len(samples)
            )
            return samples

        n_samples = check_modules_n_samples(core.config)
        method = sampling_method(core.config) or type(core.sampler).__name__
        core._logger.warning(
            "Jarvis check: CSV not found (configured as %r); "
            "drawing %d smoke point(s) from Sampling.Method=%r",
            raw_spec or "(none)",
            n_samples,
            method,
        )
        samples = core._sample_check_module_points_from_sampler(n_samples)
        core._logger.warning(
            "Jarvis check: prepared %d sampler-drawn smoke point(s)", len(samples)
        )
        return samples

    def _sample_check_module_points_from_sampler(self, n_samples: int) -> list[Sample]:
        """Draw up to *n_samples* smoke points from the active sampler (V1-like)."""
        core = self._core

        n_samples = max(1, int(n_samples))
        samples: list[Sample] = []

        # Prefer propose_next for Bridson/Random/Grid/CSV when available.
        propose = getattr(core.sampler, "propose_next", None)
        if callable(propose):
            for _ in range(n_samples):
                try:
                    sample = propose()
                except Exception as exc:
                    core._logger.warning(
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

            variables = load_variables(core.config)
            ndim = max(1, len(variables))
        except Exception:
            ndim = 2
        remaining = n_samples - len(samples)
        seed = 0
        try:
            sampling = dict(core.config.get("Sampling") or {})
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
            samples.append(core.sampler._build_sample(u_coords))

        if not samples:
            raise RuntimeError(
                "Jarvis check could not build any smoke samples: "
                "CSV missing and sampler produced no points"
            )
        return samples

    def run_distributed_scan(self, *, timeout: float = 3600.0) -> int:
        """Drive a stateless sampler through propose → Redis → Archiver."""
        core = self._core

        if core.sampler is None:
            raise RuntimeError("sampler is not configured")
        start_records = core._archiver_records_written()
        requeued = 0
        if core._resume_policy == "resume" and hasattr(core.sampler, "repropose_unfinished"):
            requeued = len(core.sampler.repropose_unfinished())
        if hasattr(core.sampler, "run_distributed"):
            pushed = int(core.sampler.run_distributed())
        elif hasattr(core.sampler, "submit_all_remaining"):
            pushed = len(core.sampler.submit_all_remaining())
        else:
            raise RuntimeError(
                f"Sampler {type(core.sampler).__name__} does not implement a distributed run API"
            )
        expected_records = start_records + requeued + pushed
        total_work = requeued + pushed
        if expected_records > start_records:
            core.wait_for_results(
                expected_records,
                timeout=timeout,
                progress_total=total_work,
                progress_base=start_records,
            )
        core._finalize_sample_buckets()
        return requeued + pushed

    def run_adaptive_scan(
        self,
        *,
        generation_timeout: float = 3600.0,
        timeout: float | None = None,
    ) -> int:
        """Drive feedback samplers; timeout is a per-generation wait, not a run budget."""
        core = self._core

        if timeout is not None:
            generation_timeout = timeout
        if core.sampler is None:
            raise RuntimeError("sampler is not configured")
        if core.redis is None:
            raise RuntimeError("adaptive scan requires redis")
        if hasattr(core.sampler, "set_redis"):
            core.sampler.set_redis(core.redis)
        # Stale feedback from a prior run must not poison barriers.
        drained_fb = int(core.redis.drain_feedback_queue())
        if drained_fb:
            core._logger.info("drained %d stale feedback record(s)", drained_fb)
        requeued = 0
        if core._resume_policy == "resume" and hasattr(core.sampler, "repropose_unfinished"):
            requeued = len(core.sampler.repropose_unfinished())
        if hasattr(core.sampler, "run_adaptive"):
            pushed = int(
                core.sampler.run_adaptive(generation_timeout=generation_timeout)
            )
        elif hasattr(core.sampler, "run_distributed"):
            pushed = int(core.sampler.run_distributed())
        else:
            raise RuntimeError(
                f"Sampler {type(core.sampler).__name__} does not implement run_adaptive"
            )
        # Nested/adaptive samplers barrier on *feedback* first; Archiver lags
        # behind on the archive queue. Wait until DATABASE rows catch Redis
        # completed counts so samples.csv export is not a partial HDF5 snapshot.
        core._wait_for_archive_caught_up(timeout=generation_timeout)
        core._finalize_sample_buckets()
        return requeued + pushed

    def _wait_for_archive_caught_up(self, *, timeout: float = 3600.0) -> None:
        """Block until Archiver has persisted every Redis-completed sample.

        Process-mode :class:`ArchiverProcess.drain` is a parent-side no-op; the
        child keeps consuming ``hep:archive``. Plot/CSV export must wait until
        ``records_written >= ok + failed`` (and workers are idle), otherwise
        ``DATABASE/samples.csv`` misses the last archive batches.
        """
        core = self._core

        if core.archiver is None:
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
            counters = core._live_sample_counters()
            expected = int(counters["ok"]) + int(counters["failed"])
            inflight = int(counters["running"]) + int(counters["queued"])
            if inflight == 0:
                break
            time.sleep(0.05)
        submitted_uuids = {
            str(item)
            for item in getattr(core.sampler, "submitted_uuids", frozenset())
            if str(item)
        }
        if bool(getattr(core.sampler, "uses_indexed_resume", False)):
            # Redis sample counters are run-local and can include work from a
            # killed control process.  Indexed samplers instead have one
            # durable logical target: the highest proposed sample index.
            expected = max(0, int(getattr(core.sampler, "_accepted_index", expected) or expected))
        elif submitted_uuids or core._persisted_uuids:
            expected = len(submitted_uuids | core._persisted_uuids)
        # Parent drain is meaningful only for in-process SimpleArchiver.
        try:
            core.archiver.drain(idle_timeout=2.0)
        except Exception:
            pass
        if expected <= 0:
            return
        remaining = max(1.0, deadline - time.monotonic())
        core.wait_for_results(
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
        core = self._core

        wait_timeout = (
            float(timeout)
            if timeout is not None
            else check_modules_timeout_sec(core.config)
        )
        if wait_timeout <= 0:
            wait_timeout = check_modules_timeout_sec(core.config)
        sample_root = core._resolve_sample_root()
        core._logger.warning(
            "Start Jarvis check smoke (assembly-line test) "
            "workers=1 sample_root=%s layout=flat-uuid pack=off timeout_sec=%.1f",
            sample_root,
            wait_timeout,
        )
        samples = core._build_check_module_samples()
        if not samples:
            raise RuntimeError("Jarvis check produced an empty sample list")
        # ``records_written`` is cumulative for an existing DATABASE.  A
        # check may reuse the same output directory, so wait for the new
        # archive delta instead of treating old rows as this run's results.
        start_records = core._archiver_records_written()
        core.submit_samples(samples)
        core.wait_for_results(
            start_records + len(samples),
            timeout=wait_timeout,
            progress_total=len(samples),
            progress_base=start_records,
            require_worker_completion=True,
        )
        core._finalize_sample_buckets()
        if verify_golden is not None:
            task_result_dir = str(core.info.get("task_result_dir") or os.getcwd())
            verify_check_modules_golden(task_result_dir=task_result_dir, golden=verify_golden)
        core._logger.warning(
            "Jarvis check finished: submitted %d sample(s) "
            "(pipeline smoke only — not a full MultiNest/Dynesty evidence run; "
            "inspect artifacts under %s)",
            len(samples),
            core._resolve_sample_root(),
        )
        return len(samples)

    def _apply_check_modules_runtime_policy(self) -> None:
        """Force smoke-friendly layout: 1 worker, SAMPLE/test, no tar pack.

        Applied before bootstrap so Factory/Archiver/Workers all see the policy.
        """
        core = self._core

        # --- single worker ---
        runtime = dict(get_runtime_block(core.config))
        runtime["workers"] = 1
        core.config["Runtime"] = runtime
        core.runtime = runtime

        envreqs = core.config.get("EnvReqs")
        envreqs = dict(envreqs) if isinstance(envreqs, Mapping) else {}
        v2 = envreqs.get("V2")
        v2 = dict(v2) if isinstance(v2, Mapping) else {}
        v2["workers"] = 1
        # SAMPLE/test layout: no numbered buckets (flat uuid dirs), no tar pack.
        sample_dir = dict(get_sample_directory_config(core.config))
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
        core.config["EnvReqs"] = envreqs

        # --- One calculator PackID only ---
        calculators = core.config.get("Calculators")
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
        core.config["Calculators"] = calculators

        # Layout flag: SAMPLE/test instead of SAMPLE/
        core.config["_check_modules_sample_layout"] = True
        # Eager path so logging / helpers work even before _init_sample_buckets.
        root = core._resolve_sample_root()
        os.makedirs(root, exist_ok=True)
        core.info["sample_root"] = root
        core.info["check_modules"] = True

    def _resolve_sample_root(self) -> str:
        """Return SAMPLE root; ``Jarvis check`` uses ``SAMPLE/test`` (no tar pack)."""
        core = self._core

        task_result_dir = str(
            core.info.get("task_result_dir")
            or core.config.get("task_result_dir")
            or os.getcwd()
        )
        base = os.path.join(task_result_dir, "SAMPLE")
        if bool(core.config.get("_check_modules_sample_layout")) or bool(
            core.info.get("check_modules")
        ):
            return os.path.join(base, "test")
        cached = core.info.get("sample_root")
        if isinstance(cached, str) and cached.strip():
            return os.path.abspath(cached)
        return base

    def submit_samples(self, samples: Sequence[Sample]) -> None:
        core = self._core

        if core.sampler is None:
            raise RuntimeError("sampler is not configured")
        if hasattr(core.sampler, "_submit_group"):
            core.sampler._submit_group(list(samples))
        else:
            for sample in samples:
                core.sampler._submit(sample)

    def _live_sample_counters(self) -> dict[str, int]:
        """Snapshot Redis sample/queue counters for completion progress lines."""
        core = self._core

        ok = failed = running = queued = archive_q = 0
        if core.redis is not None:
            try:
                stats = core.redis.fetch_sample_stats()
                lengths = core.redis.get_queue_lengths()
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
        require_worker_completion: bool = False,
    ) -> None:
        """Block until this batch is complete and Archiver has written its rows.

        When ``progress_total`` is set, emit V1-style ‰ completion heartbeats
        using archived count relative to ``progress_base`` (default 0).
        """
        core = self._core

        if core.archiver is None:
            raise RuntimeError("archiver is not configured")
        deadline = time.monotonic() + max(0.1, float(timeout))
        base = int(progress_base or 0)
        total = int(progress_total) if progress_total is not None else max(0, int(expected) - base)
        progress: PermilleProgress | None = None
        if total > 0:
            # Control-plane archive ‰ heartbeats are noisy vs DataRecorder/Archiver;
            # keep them at DEBUG. Final "sample drain complete" stays INFO below.
            progress = PermilleProgress(
                core._logger,
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
            written = core._archiver_records_written()
            counters = core._live_sample_counters()
            if progress is not None:
                done = max(0, written - base)
                extra = (
                    f"ok={counters['ok']} failed={counters['failed']} "
                    f"queued={counters['queued']} running={counters['running']} "
                    f"archive_q={counters['archive_q']} "
                    f"archived={written}/{expected}"
                )
                progress.update(done, extra=extra)
            # DATABASE rows are cumulative.  Check-mode additionally asks for
            # a terminal Worker count so an old row count cannot produce a
            # false outcome with completed=0.  Keep the legacy archive-only
            # barrier for normal scans, whose callers already establish their
            # own worker/feedback barrier.
            batch_workers_done = (
                not require_worker_completion
                or core.redis is None
                or counters["ok"] + counters["failed"] >= total
            )
            if written >= expected and batch_workers_done:
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
                    core._logger.info(
                        "sample drain complete: %d/%d archived in %s",
                        total,
                        total,
                        format_duration(time.time() - progress.t0),
                    )
                return
            # Detect permanent stall: workers done, archive queue empty, count frozen.
            workers_done = (
                core.redis is None
                or (
                    counters["running"] == 0
                    and counters["queued"] == 0
                    and counters["ok"] + counters["failed"] >= total
                )
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
            f"got {core._archiver_records_written()}"
        )

    def _finalize_sample_buckets(self) -> None:
        """Seal the open bucket after all samples finish; Archiver packs when idle."""
        core = self._core

        if core.redis is None:
            return
        try:
            core.redis.seal_current_sample_bucket()
        except Exception as exc:
            core._logger.warning("seal current SAMPLE bucket failed -> %s", exc)
        if core.archiver is not None:
            try:
                core.archiver.drain(idle_timeout=2.0)
            except Exception as exc:
                core._logger.warning("archiver bucket drain failed -> %s", exc)
