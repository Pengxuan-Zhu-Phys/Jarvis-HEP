#!/usr/bin/env python3
"""Dedicated FileOperation service (HEP-owned save/copy/delete).

Runs as a **separate process** (default) inside each Worker so calculator
threads never block on NAS-style copies on the calculator critical path
more than a queue hop. Inline mode is available for tests.

YAML is unchanged: ``save: true`` still means SAMPLE copy; execution is here,
not in Jarvis-Portal adapters.
"""

from __future__ import annotations

import os
import queue
import signal
import threading
import traceback
import uuid
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from jarvishep2.file_ops import (
    DEFAULT_DELETE_METHOD,
    apply_io_save_policy,
    delete_paths,
    normalize_delete_method,
)
from jarvishep2.mp_context import get_spawn_context


@dataclass
class FileOperationService:
    """Client side of the FileOperation process (or inline backend)."""

    mode: str = "process"  # process | inline
    delete_method: str = DEFAULT_DELETE_METHOD
    scan_name: str = ""
    _request: Any = None
    _response: Any = None
    _process: Any = None
    _started: bool = False

    @property
    def pid(self) -> int | None:
        """PID of the dedicated child, if process mode has started it."""
        process = self._process
        raw = getattr(process, "pid", None) if process is not None else None
        return int(raw) if raw is not None else None

    @classmethod
    def start(
        cls,
        *,
        mode: str = "process",
        delete_method: str = DEFAULT_DELETE_METHOD,
        scan_name: str | None = None,
    ) -> FileOperationService:
        service = cls(
            mode=str(mode or "process").strip().lower() or "process",
            delete_method=normalize_delete_method(delete_method),
            scan_name=str(scan_name or "").strip(),
        )
        service._start()
        return service

    def _start(self) -> None:
        if self._started:
            return
        if self.mode == "inline":
            self._started = True
            return
        ctx = get_spawn_context()
        self._request = ctx.Queue()
        self._response = ctx.Queue()
        self._process = ctx.Process(
            target=_file_operation_main,
            args=(
                self._request,
                self._response,
                self.delete_method,
                self.scan_name,
                os.getpid(),
            ),
            name="Jarvis2-FileOperation",
            daemon=True,
        )
        try:
            self._process.start()
        except Exception:
            process = self._process
            if process is not None and getattr(process, "pid", None) is not None:
                try:
                    if process.is_alive():
                        _signal_process_tree(process, signal.SIGKILL)
                        process.join(timeout=1.0)
                except Exception:
                    pass
            self._close_queues()
            self._process = None
            raise
        self._started = True

    def shutdown(self, *, timeout: float = 5.0) -> None:
        if not self._started:
            return
        if self.mode == "inline":
            self._started = False
            return
        try:
            if self._request is not None:
                self._request.put({"op": "shutdown"})
        except Exception:
            pass
        process = self._process
        try:
            if process is not None:
                process.join(timeout=timeout)
                if process.is_alive():
                    _signal_process_tree(process, signal.SIGTERM)
                    process.join(timeout=1.0)
                if process.is_alive():
                    # File operations can be blocked in filesystem I/O.  A
                    # second hard-stop tier prevents a half-shutdown Worker
                    # from abandoning the service when SIGTERM is ineffective.
                    _signal_process_tree(process, signal.SIGKILL)
                    process.join(timeout=1.0)
        finally:
            self._close_queues()
            if process is not None:
                try:
                    if not process.is_alive():
                        process.close()
                except (AttributeError, ValueError):
                    pass
            self._process = None
            self._started = False

    def _close_queues(self) -> None:
        """Close parent-side queue handles without waiting on feeder threads."""
        for channel in (self._request, self._response):
            if channel is None:
                continue
            try:
                channel.close()
            except Exception:
                pass
            try:
                channel.cancel_join_thread()
            except Exception:
                pass
        self._request = None
        self._response = None

    def save_io(
        self,
        *,
        source_path: str,
        sample_save_dir: str | None,
        module: str | None,
        spec: Mapping[str, Any],
        direction: str,
    ) -> str | None:
        """Apply SAMPLE save/.temp policy; returns destination path or None."""
        payload = {
            "op": "save_io",
            "source_path": str(source_path),
            "sample_save_dir": sample_save_dir,
            "module": module,
            "spec": dict(spec),
            "direction": str(direction),
        }
        return self._call(payload)

    def delete(self, paths: list[str] | tuple[str, ...], *, missing_ok: bool = True) -> None:
        self._call(
            {
                "op": "delete",
                "paths": [str(p) for p in paths],
                "missing_ok": bool(missing_ok),
                "method": self.delete_method,
            }
        )

    def _call(self, payload: dict[str, Any]) -> Any:
        if not self._started:
            self._start()
        if self.mode == "inline":
            return _execute_job(payload, delete_method=self.delete_method)

        job_id = str(uuid.uuid4())
        payload = dict(payload)
        payload["job_id"] = job_id
        assert self._request is not None and self._response is not None
        self._request.put(payload)
        # Wait for matching job_id (single-threaded client per worker is fine).
        while True:
            try:
                reply = self._response.get(timeout=0.5)
            except queue.Empty:
                process = self._process
                if process is None or not process.is_alive():
                    raise RuntimeError(
                        "FileOperation process exited before returning a result"
                    )
                continue
            if not isinstance(reply, dict):
                continue
            if reply.get("job_id") != job_id:
                # Unexpected; re-queue is not supported — drop.
                continue
            if not reply.get("ok", False):
                raise RuntimeError(str(reply.get("error") or "FileOperation failed"))
            return reply.get("result")


def _execute_job(payload: Mapping[str, Any], *, delete_method: str) -> Any:
    op = str(payload.get("op") or "").strip()
    if op == "shutdown":
        return None
    if op == "save_io":
        return apply_io_save_policy(
            source_path=str(payload.get("source_path") or ""),
            sample_save_dir=(
                str(payload["sample_save_dir"])
                if payload.get("sample_save_dir") is not None
                else None
            ),
            module=str(payload["module"]) if payload.get("module") is not None else None,
            spec=dict(payload.get("spec") or {}),
            direction=str(payload.get("direction") or "output"),
        )
    if op == "delete":
        delete_paths(
            list(payload.get("paths") or []),
            method=str(payload.get("method") or delete_method),
            missing_ok=bool(payload.get("missing_ok", True)),
        )
        return None
    raise ValueError(f"unknown FileOperation op: {op!r}")


def _signal_process_tree(process: Any, sig: int) -> None:
    """Signal the FileOperation session, falling back to its direct PID."""
    raw_pid = getattr(process, "pid", None)
    if raw_pid is None:
        return
    pid = int(raw_pid)
    try:
        if hasattr(os, "getpgid") and os.getpgid(pid) == pid:
            os.killpg(pid, sig)
            return
    except ProcessLookupError:
        return
    except (PermissionError, OSError):
        pass
    try:
        os.kill(pid, sig)
    except (ProcessLookupError, PermissionError, OSError):
        pass


def _file_operation_main(
    request_queue: Any,
    response_queue: Any,
    delete_method: str,
    scan_name: str = "",
    owner_pid: int | None = None,
) -> None:
    """Child process entry — process title + job loop."""
    if owner_pid is not None and os.getppid() != int(owner_pid):
        # The Worker died during spawn/bootstrap; never enter the queue loop.
        return
    owns_process_group = False
    try:
        # Make the service a process-group leader.  The Factory watchdog can
        # then safely SIGKILL this exact child from a recorded heartbeat PID
        # without risking an unrelated recycled PID.
        os.setsid()
        owns_process_group = True
    except (AttributeError, OSError):
        pass
    try:
        from jarvishep2.proc_title import file_operation_title, set_process_title

        set_process_title(file_operation_title(scan_name=scan_name or None))
    except Exception:
        pass
    owner_watch_stop = threading.Event()
    owner_watch: threading.Thread | None = None
    if owner_pid is not None:
        owner_watch = threading.Thread(
            target=_watch_owner_process,
            args=(int(owner_pid), owner_watch_stop, owns_process_group),
            name="Jarvis2-FileOperationOwnerWatch",
            daemon=True,
        )
        owner_watch.start()
    method = normalize_delete_method(delete_method)
    try:
        while True:
            try:
                job = request_queue.get(timeout=2.0)
            except queue.Empty:
                if owner_pid and os.getppid() != int(owner_pid):
                    break
                continue
            except Exception:
                break
            if not isinstance(job, dict):
                continue
            job_id = job.get("job_id")
            if str(job.get("op") or "") == "shutdown":
                if job_id is not None:
                    response_queue.put({"job_id": job_id, "ok": True, "result": None})
                break
            try:
                result = _execute_job(job, delete_method=method)
                response_queue.put({"job_id": job_id, "ok": True, "result": result})
            except Exception as exc:
                response_queue.put(
                    {
                        "job_id": job_id,
                        "ok": False,
                        "error": f"{exc}\n{traceback.format_exc()}",
                    }
                )
    finally:
        owner_watch_stop.set()
        if owner_watch is not None:
            owner_watch.join(timeout=0.5)


def _watch_owner_process(
    owner_pid: int,
    stop: threading.Event,
    owns_process_group: bool,
    *,
    poll_interval_sec: float = 0.25,
) -> None:
    """Hard-exit the service when its Worker parent disappears.

    This runs independently of the job loop, so it still fires while the main
    thread is blocked in a long NAS copy, delete, or response-queue write.
    """
    expected = int(owner_pid)
    interval = max(0.05, float(poll_interval_sec))
    while not stop.wait(interval):
        if os.getppid() != expected:
            if owns_process_group:
                try:
                    os.killpg(os.getpid(), signal.SIGKILL)
                except (ProcessLookupError, PermissionError, OSError):
                    pass
            os._exit(0)


__all__ = ["FileOperationService"]
