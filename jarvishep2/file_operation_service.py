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
    _request: Any = None
    _response: Any = None
    _process: Any = None
    _started: bool = False

    @classmethod
    def start(
        cls,
        *,
        mode: str = "process",
        delete_method: str = DEFAULT_DELETE_METHOD,
    ) -> FileOperationService:
        service = cls(
            mode=str(mode or "process").strip().lower() or "process",
            delete_method=normalize_delete_method(delete_method),
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
            args=(self._request, self._response, self.delete_method),
            name="Jarvis2-FileOperation",
            daemon=True,
        )
        self._process.start()
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
        if self._process is not None:
            self._process.join(timeout=timeout)
            if self._process.is_alive():
                self._process.terminate()
                self._process.join(timeout=1.0)
        self._process = None
        self._request = None
        self._response = None
        self._started = False

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
            reply = self._response.get()
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


def _file_operation_main(request_queue: Any, response_queue: Any, delete_method: str) -> None:
    """Child process entry — process title + job loop."""
    try:
        from jarvishep2.proc_title import set_process_title

        set_process_title("Jarvis2-FileOperation")
    except Exception:
        pass
    method = normalize_delete_method(delete_method)
    while True:
        try:
            job = request_queue.get()
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


__all__ = ["FileOperationService"]
