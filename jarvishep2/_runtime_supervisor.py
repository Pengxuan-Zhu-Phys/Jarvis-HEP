#!/usr/bin/env python3
"""RuntimeSupervisor collaborator for Jarvis2Core (D25.3)."""

from __future__ import annotations

import os
import signal
import threading
from collections.abc import Mapping
from typing import Any

from jarvishep2.archiver import ArchiverProcess, SimpleArchiver
from jarvishep2.command_parser import CommandParser, prepare_calculator_modules
from jarvishep2.calculator_modes import expand_calculator_modes, mode_info
from jarvishep2.library import LibraryInstaller
from jarvishep2.factory import TaskFactory
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import CONTROL_LOCK_TTL_SEC, RedisQueue
from jarvishep2.runtime_metadata import write_scan_metadata
from jarvishep2.runtime_config import (
    get_archiver_config,
    get_delete_method,
    get_redis_config,
    get_sample_directory_config,
    get_watchdog_config,
    pack_buckets_enabled,
)
from jarvishep2.worker_config import build_command_parser, build_worker_config
from jarvishep2.Module.runtime_preparer import (
    prepare_install_controls,
    refresh_install_control_summaries,
)
from jarvishep2.database import prepare_hdf5_database_for_writer
from jarvishep2.Sampling.runtime_checkpoint import prepare_resume


class _RuntimeSupervisor:
    """Private Jarvis2Core collaborator (D25.3)."""

    def __init__(self, core: Any) -> None:
        self._core = core

    def bootstrap_distributed_runtime(self) -> None:
        """Bring up Redis, Workers, and Archiver for a distributed run."""
        core = self._core

        if not core.is_redis_runtime():
            raise RuntimeError("distributed runtime requires the internal Redis runtime")
        from jarvishep2.process_cleanup import ensure_scan_name_available

        scan_name = str(core.info.get("scan_name") or core.config.get("scan_name") or "scan")
        ensure_scan_name_available(scan_name, cleanup_stale=core._resume_policy == "resume")
        db_path = os.path.join(
            str(core.info.get("task_result_dir") or os.getcwd()),
            "DATABASE",
            "samples.hdf5",
        )
        if core._resume_policy == "resume":
            try:
                # Control must not open samples.hdf5 for append (write lock would
                # block Archiver). Only clear SWMR flags / reclaim stale Archivers.
                # Logs a single summary line internally when recovery runs.
                prepare_hdf5_database_for_writer(
                    db_path,
                    logger=core._logger,
                    scan_name=scan_name,
                    probe_append=False,
                )
            except Exception as exc:
                # Non-fatal at this early stage — Archiver will retry recover on open.
                core._logger.warning(
                    "Early HDF5 prepare failed (will retry at Archiver start): %s",
                    exc,
                )
        if core._resume_policy == "resume":
            core._load_persisted_database_state(os.path.dirname(db_path))
        core.init_logger()
        # Preflight remains before command setup, Redis, Workers, and Archiver,
        # while init_logger above ensures all outcomes reach core.log.
        core.check_environment_requirements()
        # LibDeps is globally shared.  Resolve its V1 tokens with no registered
        # executable side effects, install/reuse it once, then register tools.
        root_path = core._environment_summary.get("ROOT path")
        library_parser = build_command_parser(
            core.config,
            root_path=str(root_path) if root_path else None,
            register_executables=False,
        )
        core.prepare_libraries(parser=library_parser)
        core.init_command_parser()
        core.init_redis()
        core._claim_redis_control_lock()
        core._publish_runtime_metadata()
        if core._resume_policy != "resume":
            core._reset_redis_for_fresh_run()
        core.init_sampler_from_config()
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_uuids"):
            core.sampler.set_persisted_uuids(core._persisted_uuids)
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_prefix"):
            core.sampler.set_persisted_prefix(core._persisted_index_prefix)
        if (
            core._resume_policy == "resume"
            and core.sampler is not None
            and hasattr(core.sampler, "advance_to_persisted_prefix")
        ):
            core.sampler.advance_to_persisted_prefix(core._persisted_index_prefix)
        core._init_sample_buckets()
        core.init_factory()
        # Immediately before Archiver spawn: reclaim orphan Archivers only —
        # still no append-open in Control (probe_append=False).
        if core._resume_policy == "resume" and os.path.isfile(db_path):
            try:
                # Only re-runs work if the file is still dirty; one summary log max.
                prepare_hdf5_database_for_writer(
                    db_path,
                    logger=core._logger,
                    scan_name=scan_name,
                    probe_append=False,
                )
            except Exception as exc:
                core._logger.warning(
                    "Pre-Archiver HDF5 prepare failed (Archiver will retry): %s",
                    exc,
                )
        core.init_archiver(db_path)
        core._load_persisted_database_state(os.path.dirname(db_path))
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_uuids"):
            core.sampler.set_persisted_uuids(core._persisted_uuids)
        if core.sampler is not None and hasattr(core.sampler, "set_persisted_prefix"):
            core.sampler.set_persisted_prefix(core._persisted_index_prefix)
        if (
            core._resume_policy == "resume"
            and core.sampler is not None
            and hasattr(core.sampler, "advance_to_persisted_prefix")
        ):
            core.sampler.advance_to_persisted_prefix(core._persisted_index_prefix)
        if core._resume_policy == "resume":
            core._logger.warning(
                "Resume DATABASE reconciliation found prefix=%d, "
                "completed=%d failed=%d, %d legacy UUID(s)",
                core._persisted_index_prefix,
                max(0, core._persisted_records_count - core._persisted_failed_count),
                core._persisted_failed_count,
                len(core._persisted_uuids),
            )
            # Nested / feedback engines: ensure runtime state was applied once Redis
            # and DATABASE facts exist (idempotent if set_sampler already imported).
            core._import_sampler_checkpoint_state(core.sampler)
            # Always arm nested engine resume under --resume, even if state.pkl
            # is missing: nested_engine.pkl alone is enough to continue.
            if core.sampler is not None and hasattr(core.sampler, "arm_engine_resume"):
                core.sampler.arm_engine_resume()
            ckpt = core.checkpoint_file()
            core._logger.warning(
                "Resume checkpoint path → %s (payload=%s)",
                ckpt,
                "yes" if core._resume_checkpoint_payload is not None else "MISSING",
            )
        core.start_runtime_checkpoint()
        core._export_workflow_flowchart()

    def init_redis(self, *, client: Any = None) -> RedisQueue:
        # EnvReqs.V2.redis (optional) overlays INTERNAL_REDIS_CONFIG defaults.
        core = self._core
        from jarvishep2.core import RedisQueue

        redis_config = get_redis_config(core.config)
        raw_runtime = core.config.get("Runtime")
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
            core._ensure_managed_redis(redis_config)
            redis_config = get_redis_config(core.config)
        if client is not None:
            core.redis = RedisQueue(redis_config, client=client)
        else:
            core.redis = RedisQueue(redis_config)
        core.redis.connect()
        try:
            core.redis.ping()
        except Exception as exc:
            host = redis_config.get("host", "127.0.0.1")
            port = redis_config.get("port", 6379)
            raise RuntimeError(
                f"internal Redis is unavailable at {host}:{port}; "
                "install redis-server or start a local Redis service"
            ) from exc
        return core.redis

    def _ensure_managed_redis(self, redis_config: Mapping[str, Any]) -> None:
        """Start redis-server as Jarvis-Redis:<scan> when the port is free."""
        core = self._core
        from jarvishep2.core import (
            ManagedRedisServer,
            find_available_redis_port,
            redis_port_open,
            update_default_redis_port,
        )

        if core._managed_redis is not None and core._managed_redis.started_by_us:
            return
        redis_config = dict(redis_config)
        reassigned = False
        host = str(redis_config.get("host") or "127.0.0.1")
        requested_port = int(redis_config.get("port") or 6379)
        if redis_port_open(host, requested_port):
            defaults_path = core.config.get("environment_defaults_path")
            if not defaults_path or core.config.get("redis_port_task_override"):
                raise RuntimeError(
                    f"Redis port {host}:{requested_port} is already in use. "
                    "Set EnvReqs.V2.redis.port in the project default environment YAML "
                    "(not in the task YAML) to enable automatic port reassignment."
                )
            selected_port = find_available_redis_port(host, requested_port + 1)
            update_default_redis_port(str(defaults_path), selected_port)
            redis_config["port"] = selected_port
            core._set_runtime_redis_config(redis_config)
            reassigned = True
            core._logger.warning(
                "Redis port %s:%s is occupied; reassigned this project to %s and updated %s",
                host,
                requested_port,
                selected_port,
                defaults_path,
            )
        scan_name = str(
            core.info.get("scan_name") or core.config.get("scan_name") or "scan"
        ).strip()
        task_result_dir = str(
            core.info.get("task_result_dir")
            or core.config.get("task_result_dir")
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
            core._logger.warning("managed redis-server ensure failed -> %s", exc)
            return
        core._managed_redis = managed
        if started:
            core._logger.info(
                "started managed redis-server as %s on %s:%s",
                managed.title,
                managed.host,
                managed.port,
            )
        else:
            core._logger.info(
                "using existing Redis at %s:%s "
                "(process title is only set when Jarvis starts redis-server)",
                managed.host,
                managed.port,
            )

    def _set_runtime_redis_config(self, redis_config: Mapping[str, Any]) -> None:
        """Keep the loaded task in sync with an automatically reassigned port."""
        core = self._core

        normalized = dict(redis_config)
        envreqs = dict(core.config.get("EnvReqs") or {})
        v2 = dict(envreqs.get("V2") or {})
        v2["redis"] = normalized
        envreqs["V2"] = v2
        core.config["EnvReqs"] = envreqs
        runtime = dict(core.config.get("Runtime") or {})
        runtime["redis"] = normalized
        core.config["Runtime"] = runtime
        core.runtime = runtime

    def _control_lock_owner_id(self) -> str:
        core = self._core

        run_id = str(core.info.get("run_id") or "jarvis2-run")
        return f"{os.uname().nodename}:{os.getpid()}:{run_id}"

    def _calculator_module_names(self) -> list[str]:
        core = self._core

        modules = (core.config.get("Calculators") or {}).get("Modules") or []
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
        core = self._core

        if core.redis is None:
            raise RuntimeError("init_redis() must run before claiming the control lock")
        owner = core._control_lock_owner_id()
        if core.redis.claim_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
            core._control_lock_owner = owner
            core._start_control_lease_refresh()
            core._logger.info("claimed Redis control lock (%s)", owner)
            return
        current = core.redis.get_control_lock_owner()
        if core._resume_policy == "resume" and current:
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
                core.redis.release_control_lock(str(current))
                if core.redis.claim_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
                    core._control_lock_owner = owner
                    core._start_control_lease_refresh()
                    core._logger.warning(
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
        core = self._core

        interval = max(1.0, CONTROL_LOCK_TTL_SEC / 3.0)
        while not stop.wait(interval):
            redis = core.redis
            if redis is None:
                return
            try:
                if not redis.refresh_control_lock(owner, ttl_sec=CONTROL_LOCK_TTL_SEC):
                    core._logger.error("lost Redis control lease; requesting shutdown")
                    core._interrupt_requested = True
                    return
            except Exception as exc:
                core._logger.warning("control lease refresh failed -> %s", exc)

    def _start_control_lease_refresh(self) -> None:
        core = self._core

        owner = str(core._control_lock_owner or "")
        if not owner:
            return
        core._stop_control_lease_refresh()
        stop = threading.Event()
        thread = threading.Thread(
            target=core._control_lease_loop,
            args=(stop, owner),
            name="Jarvis-ControlLease",
            daemon=True,
        )
        core._control_lease_stop = stop
        core._control_lease_thread = thread
        thread.start()

    def _stop_control_lease_refresh(self) -> None:
        core = self._core

        if core._control_lease_stop is not None:
            core._control_lease_stop.set()
        if core._control_lease_thread is not None:
            core._control_lease_thread.join(timeout=2.0)
        core._control_lease_stop = None
        core._control_lease_thread = None

    def _reset_redis_for_fresh_run(self) -> None:
        """Drop ephemeral queues/stats/calc pools before Workers are spawned."""
        core = self._core

        if core.redis is None:
            return
        workers = int(core.runtime.get("workers", 0) or 0)
        result = core.redis.reset_run_ephemeral_keys(
            calculator_names=core._calculator_module_names(),
            worker_ids=max(workers, 32),
        )
        core.redis.clear_archived_uuids(str(core.info.get("scan_name") or "scan"))
        core.redis.clear_archived_index_prefix(str(core.info.get("scan_name") or "scan"))
        core._logger.info(
            "reset Redis run keys (deleted=%s sample_stats_reset=%s)",
            result.get("deleted_keys"),
            result.get("sample_stats_reset"),
        )

    def _release_redis_control_lock(self) -> None:
        core = self._core

        core._stop_control_lease_refresh()
        owner = getattr(core, "_control_lock_owner", None)
        if not owner or core.redis is None:
            return
        try:
            core.redis.release_control_lock(str(owner))
            core._logger.info("released Redis control lock (%s)", owner)
        except Exception as exc:
            core._logger.warning("failed to release Redis control lock -> %s", exc)
        core._control_lock_owner = None

    def _publish_runtime_metadata(self) -> None:
        core = self._core

        if core.redis is None:
            return
        redis_config = core.redis.connection_config()
        path = write_scan_metadata(config=core.config, info=core.info, redis=redis_config)
        core.redis.set_runtime_metadata_path(path)
        core.info["runtime_metadata_path"] = path

    def init_command_parser(self) -> CommandParser:
        """Run Phase-1 static command resolution for the loaded task config."""
        core = self._core

        root_path = core._environment_summary.get("ROOT path")
        core.command_parser = build_command_parser(
            core.config,
            root_path=str(root_path) if root_path else None,
        )
        return core.command_parser

    def prepare_libraries(self, *, parser: CommandParser | None = None) -> dict[str, dict[str, Any]]:
        """Install or reuse LibDeps before Redis, Workers, or an Archiver start."""
        core = self._core

        root_path = core._environment_summary.get("ROOT path")
        library_parser = parser or build_command_parser(
            core.config,
            root_path=str(root_path) if root_path else None,
            register_executables=False,
        )
        logs_dir = str(
            core.info.get("logs_dir")
            or os.path.join(
                str(core.info.get("task_root") or os.getcwd()),
                "logs",
                str(core.info.get("scan_name") or "scan"),
            )
        )
        result = LibraryInstaller(
            core.config,
            parser=library_parser,
            logs_dir=logs_dir,
            logger=core._logger,
            skip_installation=core._skip_library_installation,
        ).prepare()
        if result:
            core.info["library_installations"] = result
        return result

    def _apply_command_parser_to_worker_config(self, worker_config: dict[str, Any]) -> dict[str, Any]:
        core = self._core

        if core.command_parser is None:
            core.init_command_parser()
        assert core.command_parser is not None
        merged = dict(worker_config)
        calculator_modules = merged.get("calculator_modules")
        if calculator_modules:
            merged["calculator_modules"] = prepare_calculator_modules(
                calculator_modules,
                core.command_parser,
            )
        merged["command_parser"] = core.command_parser.to_picklable()
        return merged

    def build_worker_config(self, **overrides: Any) -> dict[str, Any]:
        """Build a picklable Worker blueprint with Phase-1 command resolution applied."""
        core = self._core

        task_result_dir = str(
            overrides.pop("task_result_dir", None)
            or core.info.get("task_result_dir")
            or core.config.get("task_result_dir")
            or os.getcwd()
        )
        # Propagate CLI console policy to every Worker process.
        overrides.setdefault("log_silence", bool(getattr(core, "_log_silence", False)))
        console_level = str(getattr(core, "_console_level", "WARNING")).strip().upper()
        overrides.setdefault("console_level", console_level)
        sample_config = dict(overrides.get("sample_config") or {})
        sample_config.setdefault(
            "sample_console_level",
            "DEBUG" if console_level == "DEBUG" else "ERROR",
        )
        sample_config.setdefault(
            "sample_log_silence",
            bool(getattr(core, "_log_silence", False)),
        )
        overrides["sample_config"] = sample_config
        # Prefer control-resolved scan logs dir so Workers never land in cwd/logs/jarvis_worker_*.
        if core.info.get("logs_dir"):
            overrides.setdefault("logs_dir", str(core.info["logs_dir"]))
        if core.info.get("scan_name"):
            overrides.setdefault("scan_name", str(core.info["scan_name"]))
        if core.info.get("task_root"):
            overrides.setdefault("task_root", str(core.info["task_root"]))
        if core._control_lock_owner:
            overrides.setdefault("control_lock_owner", str(core._control_lock_owner))
        extra = dict(overrides or {})
        # Stamp active sampler into feedback_return resolution (D13.8).
        if "sampler" not in extra and core.sampler is not None:
            extra["sampler"] = core.sampler
        # Prefer Jarvis check SAMPLE/test (or any resolved sample_root) over SAMPLE/.
        sample_dirs = extra.pop("sample_dirs", None)
        if sample_dirs is None:
            sample_dirs = core._resolve_sample_root()
        return build_worker_config(
            core.config,
            task_result_dir=task_result_dir,
            parser=core.command_parser,
            calculator_modules=extra.pop("calculator_modules", None),
            likelihood_expressions=extra.pop("likelihood_expressions", None),
            opera_modules=extra.pop("opera_modules", None),
            sample_dirs=sample_dirs,
            extra=extra or None,
        )

    def init_factory(self, worker_config: Mapping[str, Any] | None = None) -> TaskFactory | None:
        core = self._core

        if not core.is_redis_runtime():
            core._logger.info("Runtime.mode != redis; skipping TaskFactory bring-up")
            return None

        redis_config = get_redis_config(core.config)
        # Explicit instance ownership (D9.4) — not the deprecated singleton shell.
        if core.factory is None:
            core.factory = TaskFactory(redis_config)
        elif redis_config and not core.factory.redis_config:
            core.factory.redis_config.update(redis_config)
        if core.redis is not None:
            core.factory.redis = core.redis
        else:
            core.factory.init_redis()

        workers = int(core.runtime.get("workers", 1) or 1)
        if workers <= 0:
            workers = 1

        if worker_config is None:
            if core.command_parser is None:
                core.init_command_parser()
            merged_config = core.build_worker_config()
        else:
            merged_config = dict(worker_config)
            if "command_parser" not in merged_config:
                merged_config = core._apply_command_parser_to_worker_config(merged_config)
        if core._control_lock_owner:
            merged_config.setdefault("control_lock_owner", str(core._control_lock_owner))
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
            logger=core._logger,
        )
        core._install_control_modules = merged_config.get("calculator_modules") or []
        if core._resume_policy == "resume":
            prepare_resume(
                core.redis,
                worker_config=merged_config,
                persisted_count=core._persisted_records_count,
                persisted_failed=core._persisted_failed_count,
            )
        core.factory.start_workers(workers, **merged_config)
        core.factory.start_monitor()
        watchdog = get_watchdog_config(core.config)
        core.factory.start_watchdog(**watchdog)
        core._logger.info("TaskFactory started with %d worker(s)", workers)
        return core.factory

    def init_archiver(self, db_path: str | None = None) -> SimpleArchiver | ArchiverProcess:
        core = self._core

        if core.redis is None:
            raise RuntimeError("init_redis() must run before init_archiver()")
        task_result_dir = str(
            core.info.get("task_result_dir")
            or core.config.get("task_result_dir")
            or os.getcwd()
        )
        database_dir = os.path.join(task_result_dir, "DATABASE")
        sample_root = core._resolve_sample_root()
        os.makedirs(database_dir, exist_ok=True)
        os.makedirs(sample_root, exist_ok=True)
        resolved_db_path = db_path or os.path.join(database_dir, "samples.hdf5")
        archiver_config = dict(get_archiver_config(core.config))
        if bool(core.config.get("_check_modules_sample_layout")):
            # Belt-and-suspenders: never tar SAMPLE/test during check smoke.
            archiver_config["pack_buckets"] = False
        # Same CLI console policy as Workers.
        archiver_config.setdefault("log_silence", bool(getattr(core, "_log_silence", False)))
        archiver_config.setdefault("console_level", str(getattr(core, "_console_level", "WARNING")))
        if core._control_lock_owner:
            archiver_config["control_lock_owner"] = str(core._control_lock_owner)
        delete_method = get_delete_method(core.config)
        redis_config = get_redis_config(core.config)
        if core.redis is not None:
            # Prefer the live connection if init_redis already bound a client.
            try:
                redis_config = dict(core.redis.connection_config())
            except Exception:
                pass

        from jarvishep2.logging import component_log_path

        scan_name = str(core.info.get("scan_name") or "").strip() or None
        logs_dir = str(core.info.get("logs_dir") or "").strip() or None
        if str(archiver_config.get("mode", "process")).strip().lower() == "process":
            archiver_log = component_log_path(logs_dir, "archiver") if logs_dir else None
            core.archiver = ArchiverProcess(
                redis_config,
                db_path=resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
                scan_name=scan_name,
                log_dir=logs_dir,
                log_path=archiver_log,
            )
            core.archiver.start()
            core._logger.info(
                "Archiver process started (own logger → %s)",
                archiver_log or "logs/<scan>/archiver.log",
            )
        else:
            # Thread mode: same process as control; Archiver logger → archiver.log.
            core.archiver = SimpleArchiver(
                core.redis,
                resolved_db_path,
                sample_root=sample_root,
                delete_method=delete_method,
                archiver_config=archiver_config,
                scan_name=scan_name,
                logger=get_jarvis_logger("archiver"),
            )
            core.archiver.start()
            core._logger.info("Archiver thread started (Jarvis-HEP.Archiver)")
        core._restore_archiver_persistence()
        return core.archiver

    def _restore_archiver_persistence(self) -> None:
        """Compatibility hook; DATABASE seeding now happens in SimpleArchiver.

        Never promote checkpoint/Redis acknowledgements into completion truth:
        only UUIDs actually read from HDF5 may seed the deduplication set.
        """
        core = self._core

        return

    def _archiver_records_written(self) -> int:
        core = self._core

        archiver = core.archiver
        if archiver is None:
            return 0
        counter = getattr(archiver, "records_written", 0)
        if hasattr(counter, "value"):
            return int(counter.value)
        return int(counter)

    def _init_sample_buckets(self) -> None:
        """Register Redis SAMPLE bucket meta (numbering + active/completed state)."""
        core = self._core

        if core.redis is None:
            return
        sample_dir = get_sample_directory_config(core.config)
        sample_root = core._resolve_sample_root()
        os.makedirs(sample_root, exist_ok=True)
        meta = core.redis.init_sample_buckets(
            sample_root=sample_root,
            limit=int(sample_dir.get("limit", 200)),
            width=int(sample_dir.get("width", 6)),
            start_bucket=int(sample_dir.get("start_bucket", 1)),
            pack=pack_buckets_enabled(core.config),
            enabled=bool(sample_dir.get("enabled", True)),
        )
        core.info["sample_root"] = sample_root
        core.info["sample_directory"] = dict(sample_dir)
        enabled = bool(sample_dir.get("enabled", True))
        core._logger.info(
            "SAMPLE ready -> root=%s enabled=%s limit=%s width=%s pack=%s",
            sample_root,
            enabled,
            meta.get("limit"),
            meta.get("width"),
            bool(int(meta.get("pack", 0))),
        )

    def _install_control_signal_handlers(self) -> None:
        """SIGINT/SIGTERM → clean shutdown; refuse Ctrl+Z suspend."""
        core = self._core

        if core._signal_handlers_installed:
            return

        def _stop_handler(signum: int, _frame: Any) -> None:
            try:
                sig_name = signal.Signals(signum).name
            except Exception:
                sig_name = str(signum)
            core._interrupt_requested = True
            try:
                core._logger.warning(
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
                core._logger.warning(
                    "Ctrl+Z / SIGTSTP ignored — process suspend leaves Workers and "
                    "Redis half-alive. Use Ctrl+C to stop the scan cleanly."
                )
            except Exception:
                pass

        for sig in (signal.SIGINT, signal.SIGTERM):
            try:
                core._previous_signal_handlers[sig] = signal.getsignal(sig)
                signal.signal(sig, _stop_handler)
            except (ValueError, OSError):
                # Not in main thread or signal unsupported.
                pass
        if hasattr(signal, "SIGTSTP"):
            try:
                core._previous_signal_handlers[signal.SIGTSTP] = signal.getsignal(
                    signal.SIGTSTP
                )
                signal.signal(signal.SIGTSTP, _refuse_suspend)
            except (ValueError, OSError):
                pass
        core._signal_handlers_installed = True

    def _restore_control_signal_handlers(self) -> None:
        core = self._core

        if not core._signal_handlers_installed:
            return
        for sig, previous in list(core._previous_signal_handlers.items()):
            try:
                signal.signal(sig, previous)
            except (ValueError, OSError):
                pass
        core._previous_signal_handlers.clear()
        core._signal_handlers_installed = False

    def _stop_managed_redis(self) -> None:
        core = self._core

        managed = core._managed_redis
        if managed is None:
            return
        core._managed_redis = None
        try:
            if managed.started_by_us:
                try:
                    core._logger.info(
                        "stopping managed redis-server (%s)",
                        managed.title or "Jarvis-Redis",
                    )
                except Exception:
                    pass
            managed.stop()
        except Exception as exc:
            try:
                core._logger.warning("managed redis stop failed -> %s", exc)
            except Exception:
                pass

    def _atexit_cleanup(self) -> None:
        """Best-effort stop of managed Redis if the process dies unexpectedly."""
        core = self._core

        if core._shutdown_done:
            return
        try:
            core._stop_managed_redis()
        except Exception:
            pass

    def shutdown(self, *, wait: bool = True, write_run_summary: bool = False) -> None:
        """Stop Archiver, Workers, release the control lock, and managed Redis.

        Idempotent: safe to call from ``finally``, signal paths, and atexit.
        Managed Redis is always stopped last in an inner ``finally`` so a failure
        earlier in the chain cannot leave ``Jarvis-Redis:*`` orphaned.
        """
        core = self._core

        if core._shutdown_done:
            # Still try Redis in case a previous partial shutdown left it.
            core._stop_managed_redis()
            return
        core._shutdown_done = True
        try:
            if core.sampler is not None and hasattr(core.sampler, "shutdown_checkpointing"):
                try:
                    core.sampler.shutdown_checkpointing()
                except Exception as exc:
                    core._logger.warning("sampler checkpoint shutdown failed -> %s", exc)
            if not core._interrupt_requested:
                try:
                    core._finalize_sample_buckets()
                except Exception:
                    pass
            if core.archiver is not None:
                try:
                    if isinstance(core.archiver, ArchiverProcess):
                        # Interrupt path: don't hang forever on tar packing.
                        timeout = 5.0 if core._interrupt_requested else 30.0
                        core.archiver.stop(wait=wait, timeout=timeout)
                    else:
                        core.archiver.stop(
                            wait=wait,
                            drain=not core._interrupt_requested,
                        )
                except Exception as exc:
                    core._logger.warning("archiver stop failed -> %s", exc)
                try:
                    core._final_archived_records = core._archiver_records_written()
                except Exception:
                    core._final_archived_records = None
                core.archiver = None
            if write_run_summary and core.factory is not None:
                try:
                    core.write_run_summary()
                except Exception as exc:
                    core._logger.warning("run_summary write failed -> %s", exc)
            try:
                if core.redis is not None:
                    core.redis.clear_runtime_metadata_path()
                core._release_redis_control_lock()
            except Exception:
                pass
            if core.factory is not None:
                try:
                    core.factory.shutdown(wait=wait)
                except Exception as exc:
                    core._logger.warning("factory shutdown failed -> %s", exc)
                core.factory = None
            elif core.redis is not None:
                try:
                    core.redis.close()
                except Exception:
                    pass
            core.redis = None
            try:
                refresh_install_control_summaries(getattr(core, "_install_control_modules", []))
            except Exception as exc:
                core._logger.debug("install-control summary refresh skipped -> %s", exc)
        finally:
            core._stop_managed_redis()
