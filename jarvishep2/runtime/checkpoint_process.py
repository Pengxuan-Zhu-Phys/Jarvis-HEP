"""Control-owned spawn process for small sampler snapshots.

The local IPC window contains at most one detached snapshot. Redis carries
health only: this component never inspects task/archive queues or HDF5 files.
"""

from __future__ import annotations

import os
import signal
import socket
import time
from collections.abc import Callable, Mapping
from queue import Empty, Full
from typing import Any

from jarvishep2.mp_context import get_spawn_context
from jarvishep2.sampling.sample_checkpoint import atomic_pickle, atomic_write_bytes


def _payload_field(payload: Mapping[str, Any], field: str) -> Any:
    current: Any = payload
    for part in field.split("."):
        if not part or not isinstance(current, Mapping) or part not in current:
            return None
        current = current[part]
    return current


def _write_payload_sidecars(checkpoint_path: str, payload: Mapping[str, Any], logger: Any) -> None:
    """Write sibling files from bytes already stored in the snapshot.

    A sidecar is ``{"name": basename, "field": "sampler_state.some_blob"}``.
    The scan thread does not copy those bytes a second time, and a sidecar
    failure does not discard the resume pickle.
    """
    specs = payload.get("sidecars")
    if not isinstance(specs, (list, tuple)):
        return
    directory = os.path.dirname(os.path.abspath(checkpoint_path))
    for spec in specs:
        if not isinstance(spec, Mapping):
            continue
        name = spec.get("name")
        field = spec.get("field")
        if not isinstance(name, str) or not isinstance(field, str):
            continue
        if name != os.path.basename(name) or name in {"", ".", ".."}:
            logger.warning("Checkpoint sidecar name rejected: %s", name)
            continue
        blob = _payload_field(payload, field)
        if not isinstance(blob, (bytes, bytearray)) or not blob:
            continue
        try:
            atomic_write_bytes(os.path.join(directory, name), bytes(blob))
        except Exception as exc:
            logger.warning("Checkpoint sidecar %s was not written -> %s", name, exc)


def _run_checkpoint(
    path, inbox, idle, requested, stop, written, failures, *,
    interval, parent_pid, scan_name, redis_config, control_owner, logs_dir,
) -> None:
    from jarvishep2.logging import get_jarvis_logger, setup_jarvis_logging, shutdown_jarvis_logging
    from jarvishep2.proc_title import checkpoint_title, set_process_title
    from jarvishep2.redis_queue import (
        RedisQueue, classify_control_lock, control_lock_missing_grace_sec,
        next_control_lock_watch,
    )

    set_process_title(checkpoint_title(scan_name=scan_name))
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    if logs_dir:
        setup_jarvis_logging(role="checkpoint", component="checkpoint",
                             scan_logs_dir=logs_dir, console=False, use_queue=False)
    logger = get_jarvis_logger("checkpoint")
    redis = None
    if redis_config is not None:
        try:
            redis = RedisQueue(redis_config)
            redis.connect()
        except Exception as exc:
            logger.warning("Checkpoint telemetry unavailable -> %s", exc)
            if redis is not None:
                redis.close()
            redis = None

    last_saved = 0.0
    last_error = ""
    last_index = 0
    next_request = time.monotonic() + interval
    next_board = 0.0
    missing_since = None
    owner_lost = False

    def publish(status):
        if redis is None:
            return
        try:
            redis.publish_proc_board(
                "checkpoint", role="checkpoint", status=status,
                pid=os.getpid(), parent_pid=parent_pid, pgid=os.getpgid(0),
                host=socket.gethostname(),
                owner=control_owner, scan_name=scan_name or "",
                last_saved=last_saved, last_index=last_index,
                last_error=last_error, snapshots_written=written.value,
                write_failures=failures.value, ts=time.time(),
            )
        except Exception as exc:
            logger.debug("Checkpoint telemetry update skipped -> %s", exc)

    try:
        logger.info("Checkpoint process started pid=%d parent=%d", os.getpid(), parent_pid)
        while os.getppid() == parent_pid:
            now = time.monotonic()
            if now >= next_board:
                if redis is not None and control_owner:
                    try:
                        decision, missing_since = next_control_lock_watch(
                            classify_control_lock(redis.get_control_lock_owner(), control_owner),
                            missing_since=missing_since, now=now,
                            grace_sec=control_lock_missing_grace_sec(),
                        )
                        if decision in {"shutdown_stolen", "shutdown_missing"}:
                            owner_lost = True
                            break
                    except Exception as exc:
                        logger.debug("Checkpoint control lease check skipped -> %s", exc)
                publish("error" if last_error else "idle")
                next_board = now + 1.0
            if now >= next_request:
                requested.value = 1
                next_request = now + interval
            try:
                payload = inbox.get(timeout=0.1)
            except Empty:
                if stop.value:
                    break
                continue
            publish("writing")
            try:
                _write_payload_sidecars(path, payload, logger)
                atomic_pickle(path, payload)
                last_saved = time.time()
                last_index = int((payload.get("sampler_state") or {}).get("next_sample_index", 0))
                last_error = ""
                written.value += 1
                logger.info("Checkpoint saved to %s (next sample index=%d)", path, last_index)
            except Exception as exc:
                last_error = str(exc)[:512]
                failures.value += 1
                logger.warning("Checkpoint write failed -> %s", exc)
            finally:
                idle.value = 1
            next_request = time.monotonic() + interval
            publish("error" if last_error else "idle")
            if stop.value:
                break
    finally:
        if not owner_lost and os.getppid() == parent_pid:
            publish("stopped")
        if redis is not None:
            redis.close()
        inbox.close()
        shutdown_jarvis_logging()


class CheckpointProcess:
    """Small nonblocking producer handle; all file writes run in the child."""

    def __init__(
        self, path: str, logger: Any, *, interval: float = 30.0,
        scan_name: str = "", redis_config: Mapping[str, Any] | None = None,
        control_owner: str = "", logs_dir: str | None = None,
        on_saved: Callable[..., None] | None = None,
    ) -> None:
        self.path = str(path)
        self._logger = logger
        self._on_saved = on_saved
        self._pending_log: tuple[int, int, str] | None = None
        self._last_saved_index: int | None = None
        ctx = get_spawn_context()
        self._inbox = ctx.Queue(maxsize=1)
        # The parent must not join an IPC feeder blocked by a dead child.
        self._inbox.cancel_join_thread()
        # Single-writer scalar flags avoid locks left held by a killed child.
        self._idle = ctx.RawValue("b", 1)
        self._requested = ctx.RawValue("b", 0)
        self._stop = ctx.RawValue("b", 0)
        self._written = ctx.RawValue("Q", 0)
        self._failures = ctx.RawValue("Q", 0)
        self._process = ctx.Process(
            name=f"Jarvis-Checkpoint:{scan_name}", target=_run_checkpoint, daemon=False,
            args=(self.path, self._inbox, self._idle, self._requested, self._stop,
                  self._written, self._failures),
            kwargs=dict(interval=max(0.01, float(interval)), parent_pid=os.getpid(),
                        scan_name=scan_name, redis_config=(dict(redis_config) if redis_config is not None else None),
                        control_owner=str(control_owner or ""), logs_dir=logs_dir),
        )
        self._started = False
        self._closed = False
        self._reported_dead = False
        self.pid: int | None = None

    def start(self) -> bool:
        try:
            self._process.start()
            self._started = True
            self.pid = self._process.pid
            return True
        except Exception as exc:
            self._logger.warning("Checkpoint process could not start -> %s", exc)
            self.finish(timeout=0)
            return False

    def is_alive(self) -> bool:
        return self._started and not self._closed and self._process.is_alive()

    def requested(self) -> bool:
        return bool(self._requested.value)

    def request_snapshot(self) -> None:
        self._requested.value = 1

    @property
    def last_saved_index(self) -> int | None:
        self._report_save()
        return self._last_saved_index

    def _report_save(self) -> None:
        if self._pending_log is None or self._on_saved is None:
            return
        previous, index, progress = self._pending_log
        if self._written.value == previous:
            return
        self._pending_log = None
        self._last_saved_index = index
        self._on_saved(self.path, "generation barrier", progress=progress)

    def available(self) -> bool:
        self._report_save()
        if not self.is_alive():
            if self._started and not self._closed and not self._reported_dead:
                self._logger.warning("Checkpoint process exited; sampling continues (pid=%s)", self.pid)
                self._reported_dead = True
            return False
        return bool(self._idle.value) and not self._stop.value

    def submit(self, payload: dict[str, Any]) -> bool:
        if not self.available():
            return False
        pending_log = None
        if self._on_saved is not None:
            state = payload["sampler_state"]
            progress = f"generation={state['generation']} proposals={state['index']}"
            retained = state["session"]["retained"]
            if retained:
                progress += f" retained={retained} chain_bytes={state['trajectory_offset']}"
            pending_log = (self._written.value, state["next_sample_index"], progress)
        self._idle.value = 0
        try:
            self._inbox.put_nowait(payload)
        except (Full, OSError, ValueError) as exc:
            self._idle.value = 1
            self._logger.warning("Checkpoint snapshot was not submitted -> %s", exc)
            return False
        self._pending_log = pending_log
        self._requested.value = 0
        return True

    def wait(self, timeout: float = 10.0) -> bool:
        deadline = time.monotonic() + max(0.0, timeout)
        while not self._idle.value and self.is_alive():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                return False
            time.sleep(min(0.01, remaining))
        self._report_save()
        return bool(self._idle.value)

    def finish(self, timeout: float = 5.0) -> bool:
        if self._closed:
            return True
        self.wait(timeout)
        self._stop.value = 1
        if self._started:
            self._process.join(timeout=1.0)
            if self._process.is_alive():
                self._logger.warning("Stopping stalled checkpoint process pid=%s", self.pid)
                self._process.terminate()
                self._process.join(timeout=1.0)
            if self._process.is_alive():
                self._process.kill()
                self._process.join(timeout=1.0)
            if self._process.is_alive():
                self._logger.warning("Checkpoint process did not stop pid=%s", self.pid)
                return False
        self._inbox.close()
        self._process.close()
        self._closed = True
        return True
