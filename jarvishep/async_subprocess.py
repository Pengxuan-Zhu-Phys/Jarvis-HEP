#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import concurrent.futures
import json
import os
import signal
import subprocess
import threading
import time
import traceback
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping
from uuid import uuid4

from jarvishep.log_kv import format_two_column_log


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _safe_json(obj: Any) -> Any:
    try:
        json.dumps(obj)
        return obj
    except Exception:
        return str(obj)


def _emit_log_line(logger_obj: Any | None, message: str, *, raw: bool = False) -> None:
    if logger_obj is None:
        return

    payload = str(message)
    target = logger_obj
    if raw:
        binder = getattr(logger_obj, "bind", None)
        if callable(binder):
            try:
                target = binder(raw=True)
            except Exception:
                target = logger_obj

    try:
        target.info(payload)
    except Exception:
        return


def _normalize_cmd_for_shell(raw: str | Iterable[str]) -> str:
    if isinstance(raw, str):
        return raw
    return " ".join(str(x) for x in raw)


def _normalize_cmd_for_exec(raw: str | Iterable[str]) -> list[str]:
    if isinstance(raw, str):
        return [raw]
    return [str(x) for x in raw]


def _best_effort_fd_count() -> int | None:
    # Linux
    for p in ("/proc/self/fd", "/dev/fd"):
        try:
            return len(os.listdir(p))
        except Exception:
            continue
    return None


def _best_effort_rss_mb() -> float | None:
    try:
        import resource

        rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        # macOS reports bytes, Linux reports KB.
        if rss > 10_000_000:
            return rss / (1024.0 * 1024.0)
        return rss / 1024.0
    except Exception:
        return None


@dataclass
class SubprocessRuntimeConfig:
    max_concurrency: int = 8
    max_pending: int = 256
    queue_put_timeout_sec: float = 30.0
    per_task_timeout_sec: float | None = None
    terminate_grace_sec: float = 5.0
    chunk_size: int = 65536
    log_policy: str = "logger"  # logger | file | quiet | tee-limited
    progress_interval_sec: float = 5.0
    diagnostics_enabled: bool = False
    diagnostics_interval_sec: float = 10.0

    def __post_init__(self):
        self.max_concurrency = max(1, int(self.max_concurrency))
        self.max_pending = max(1, int(self.max_pending))
        self.queue_put_timeout_sec = max(0.1, float(self.queue_put_timeout_sec))
        if self.per_task_timeout_sec is not None:
            self.per_task_timeout_sec = max(0.1, float(self.per_task_timeout_sec))
        self.terminate_grace_sec = max(0.1, float(self.terminate_grace_sec))
        self.chunk_size = max(1024, int(self.chunk_size))
        self.progress_interval_sec = max(0.2, float(self.progress_interval_sec))
        self.diagnostics_interval_sec = max(0.2, float(self.diagnostics_interval_sec))
        mode = str(self.log_policy).strip().lower()
        if mode not in {"logger", "file", "quiet", "tee-limited"}:
            mode = "logger"
        self.log_policy = mode


@dataclass
class SubprocessJob:
    cmd: str | list[str]
    cwd: str | None = None
    env: Mapping[str, str] | None = None
    shell: bool = True
    timeout_sec: float | None = None
    log_dir: str | None = None
    log_policy: str | None = None
    task_id: str | None = None
    stream_logger: Any | None = None
    meta: Dict[str, Any] | None = None

    def __post_init__(self):
        if self.task_id is None:
            self.task_id = str(uuid4())


@dataclass
class SubprocessResult:
    task_id: str
    returncode: int
    started_at: str
    finished_at: str
    duration_sec: float
    timed_out: bool
    stdout_path: str | None
    stderr_path: str | None
    stdout_bytes: int
    stderr_bytes: int
    cmd_display: str

    @property
    def ok(self) -> bool:
        return int(self.returncode) == 0 and not bool(self.timed_out)


class SubprocessExecutionError(RuntimeError):
    def __init__(self, message: str, result: SubprocessResult | None = None):
        super().__init__(message)
        self.result = result


class AsyncSubprocessScheduler:
    """Single-loop bounded scheduler for subprocess workloads.

    - Fixed worker count (== max_concurrency)
    - Bounded pending queue for producer backpressure
    - Streaming stdout/stderr draining to avoid pipe deadlocks
    """

    def __init__(
        self,
        config: SubprocessRuntimeConfig | None = None,
        logger=None,
        status_path: str | None = None,
    ) -> None:
        self.config = config or SubprocessRuntimeConfig()
        self.logger = logger
        self.status_path = status_path

        self._loop: asyncio.AbstractEventLoop | None = None
        self._thread: threading.Thread | None = None
        self._ready = threading.Event()
        self._stopped = threading.Event()
        self._shutdown_lock = threading.Lock()

        self._queue: asyncio.Queue | None = None
        self._semaphore: asyncio.Semaphore | None = None
        self._workers: list[asyncio.Task] = []
        self._diagnostics_task: asyncio.Task | None = None
        self._progress_task: asyncio.Task | None = None
        self._shutdown_started = False

        self._submitted = 0
        self._completed = 0
        self._failed = 0
        self._running = 0
        self._peak_running = 0
        self._timed_out = 0

    def start(self) -> None:
        if self._thread and self._thread.is_alive():
            return
        self._thread = threading.Thread(
            target=self._thread_main,
            name="JarvisSubprocessScheduler",
            daemon=True,
        )
        self._thread.start()
        if not self._ready.wait(timeout=10.0):
            raise RuntimeError("AsyncSubprocessScheduler failed to start event loop thread.")

    def snapshot(self) -> Dict[str, Any]:
        pending = 0
        if self._queue is not None:
            try:
                pending = int(self._queue.qsize())
            except Exception:
                pending = 0
        return {
            "submitted": int(self._submitted),
            "completed": int(self._completed),
            "failed": int(self._failed),
            "running": int(self._running),
            "peak_running": int(self._peak_running),
            "timed_out": int(self._timed_out),
            "pending": int(pending),
            "fd_count": _best_effort_fd_count(),
            "rss_mb": _best_effort_rss_mb(),
        }

    def submit(self, job: SubprocessJob) -> concurrent.futures.Future:
        self.start()
        if self._loop is None or self._queue is None:
            raise RuntimeError("Subprocess scheduler not initialized.")
        if self._shutdown_started:
            raise RuntimeError("Subprocess scheduler is shutting down; reject new jobs.")

        out_future: concurrent.futures.Future = concurrent.futures.Future()
        entry = (job, out_future)

        async def _enqueue() -> None:
            assert self._queue is not None
            await asyncio.wait_for(self._queue.put(entry), timeout=self.config.queue_put_timeout_sec)

        put_future = asyncio.run_coroutine_threadsafe(_enqueue(), self._loop)
        try:
            put_future.result(timeout=self.config.queue_put_timeout_sec + 0.5)
        except concurrent.futures.TimeoutError as exc:
            put_future.cancel()
            out_future.set_exception(
                TimeoutError(
                    "Scheduler pending queue is full; put timed out. "
                    "Consider reducing producer speed or increasing max_pending."
                )
            )
            raise exc
        except TimeoutError as exc:
            put_future.cancel()
            out_future.set_exception(
                TimeoutError(
                    "Scheduler pending queue is full; put timed out. "
                    "Consider reducing producer speed or increasing max_pending."
                )
            )
            raise exc
        except Exception:
            if not out_future.done():
                out_future.set_exception(RuntimeError("Failed to enqueue subprocess job."))
            raise
        return out_future

    async def arun(self, job: SubprocessJob) -> SubprocessResult:
        fut = self.submit(job)
        return await asyncio.wrap_future(fut)

    def run(self, job: SubprocessJob, timeout: float | None = None) -> SubprocessResult:
        fut = self.submit(job)
        return fut.result(timeout=timeout)

    def shutdown(self, wait: bool = True, timeout: float = 30.0) -> None:
        with self._shutdown_lock:
            if self._shutdown_started:
                if wait:
                    self._stopped.wait(timeout=max(0.1, timeout))
                return
            self._shutdown_started = True
        if self._loop is None:
            self._stopped.set()
            return

        close_future = asyncio.run_coroutine_threadsafe(self._async_shutdown(), self._loop)

        def _stop_loop(_fut):
            if self._loop is not None:
                self._loop.call_soon_threadsafe(self._loop.stop)

        close_future.add_done_callback(_stop_loop)
        if wait:
            close_future.result(timeout=max(0.1, timeout))
            self._stopped.wait(timeout=max(0.1, timeout))

    def _thread_main(self) -> None:
        self._loop = asyncio.new_event_loop()
        asyncio.set_event_loop(self._loop)
        self._queue = asyncio.Queue(maxsize=self.config.max_pending)
        self._semaphore = asyncio.Semaphore(self.config.max_concurrency)

        for ii in range(self.config.max_concurrency):
            self._workers.append(self._loop.create_task(self._worker(ii)))

        if self.config.diagnostics_enabled and self.status_path:
            self._diagnostics_task = self._loop.create_task(self._diagnostics_loop())
        self._progress_task = self._loop.create_task(self._progress_loop())

        self._ready.set()
        try:
            self._loop.run_forever()
        finally:
            pending = [t for t in asyncio.all_tasks(loop=self._loop) if not t.done()]
            for task in pending:
                task.cancel()
            if pending:
                self._loop.run_until_complete(asyncio.gather(*pending, return_exceptions=True))
            self._loop.close()
            self._stopped.set()

    async def _async_shutdown(self) -> None:
        if self._queue is None:
            return

        await self._queue.join()
        for _ in self._workers:
            await self._queue.put(None)
        if self._workers:
            await asyncio.gather(*self._workers, return_exceptions=True)
        self._workers = []

        if self._diagnostics_task is not None:
            self._diagnostics_task.cancel()
            await asyncio.gather(self._diagnostics_task, return_exceptions=True)
            self._diagnostics_task = None
        if self.config.diagnostics_enabled and self.status_path:
            self._write_status_snapshot(reason="shutdown-final")
        if self._progress_task is not None:
            self._progress_task.cancel()
            await asyncio.gather(self._progress_task, return_exceptions=True)
            self._progress_task = None

        loop = asyncio.get_running_loop()
        if loop.get_debug() and self.logger is not None:
            dangling = [
                t
                for t in asyncio.all_tasks(loop=loop)
                if t is not asyncio.current_task(loop=loop) and not t.done()
            ]
            if dangling:
                self.logger.warning(
                    f"AsyncSubprocessScheduler shutdown with {len(dangling)} dangling task(s)."
                )

    async def _worker(self, worker_id: int):
        assert self._queue is not None
        while True:
            item = await self._queue.get()
            if item is None:
                self._queue.task_done()
                return

            job, out_future = item
            if out_future.cancelled():
                self._queue.task_done()
                continue

            self._submitted += 1
            self._running += 1
            if self._running > self._peak_running:
                self._peak_running = self._running

            try:
                if self._semaphore is None:
                    result = await self._run_job(job)
                else:
                    async with self._semaphore:
                        result = await self._run_job(job)
                if not out_future.done():
                    out_future.set_result(result)
                self._completed += 1
                if result.timed_out:
                    self._timed_out += 1
            except Exception as exc:
                self._failed += 1
                if not out_future.done():
                    out_future.set_exception(exc)
            finally:
                self._running = max(0, self._running - 1)
                self._queue.task_done()

    async def _run_job(self, job: SubprocessJob) -> SubprocessResult:
        mode = (job.log_policy or self.config.log_policy or "logger").strip().lower()
        if mode not in {"logger", "file", "quiet", "tee-limited"}:
            mode = "logger"

        stdout_path = None
        stderr_path = None
        stdout_fh = None
        stderr_fh = None

        task_id = str(job.task_id or uuid4())
        run_meta = dict(job.meta or {})
        stream_logger = getattr(job, "stream_logger", None)

        if mode in {"file", "tee-limited"} and job.log_dir:
            task_dir = Path(job.log_dir).expanduser().resolve()
            task_dir.mkdir(parents=True, exist_ok=True)
            stdout_path = str(task_dir / "stdout.log")
            stderr_path = str(task_dir / "stderr.log")
            stdout_fh = open(stdout_path, "wb")
            stderr_fh = open(stderr_path, "wb")

        started_ts = time.monotonic()
        started_at = _utc_now_iso()
        timeout = job.timeout_sec
        if timeout is None:
            timeout = self.config.per_task_timeout_sec
        timed_out = False
        stdout_bytes = 0
        stderr_bytes = 0

        if job.shell:
            cmd_display = _normalize_cmd_for_shell(job.cmd)
        else:
            cmd_display = " ".join(_normalize_cmd_for_exec(job.cmd))

        job_logger = self.logger
        if job_logger is not None:
            module_name = str(run_meta.get("module") or "Subprocess").strip() or "Subprocess"
            pack_id = run_meta.get("pack_id")
            if pack_id not in (None, "", "None"):
                module_name = f"{module_name}-{pack_id}"
            binder = getattr(job_logger, "bind", None)
            if callable(binder):
                try:
                    job_logger = binder(module=module_name)
                except Exception:
                    job_logger = self.logger
            if mode in {"logger", "tee-limited"}:
                job_logger.info(
                    f" Run command -> \n\t{cmd_display} \n in path -> \n\t{job.cwd or '.'} \n Screen output -> "
                )

        env = None
        if isinstance(job.env, Mapping):
            env = dict(os.environ)
            env.update({str(k): str(v) for k, v in job.env.items()})

        try:
            if job.shell:
                process = await asyncio.create_subprocess_shell(
                    _normalize_cmd_for_shell(job.cmd),
                    cwd=job.cwd,
                    env=env,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                    start_new_session=True,
                )
            else:
                exec_cmd = _normalize_cmd_for_exec(job.cmd)
                process = await asyncio.create_subprocess_exec(
                    *exec_cmd,
                    cwd=job.cwd,
                    env=env,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                    start_new_session=True,
                )
        except OSError as exc:
            fd_hint = (
                "Likely file-descriptor pressure. Check 'ulimit -n' and lower max_concurrency."
                if exc.errno in {23, 24}
                else ""
            )
            raise SubprocessExecutionError(
                f"Failed to spawn subprocess: {exc}. {fd_hint}".strip()
            ) from exc

        async def _drain_stream(stream, sink_fh, stream_name: str) -> int:
            total = 0
            emitted = 0
            text_buffer = ""

            def _emit_text(text: str) -> None:
                _emit_log_line(job_logger, text, raw=True)
                if stream_logger is not None and stream_logger is not job_logger:
                    _emit_log_line(stream_logger, text, raw=True)

            while True:
                chunk = await stream.read(self.config.chunk_size)
                if not chunk:
                    break
                total += len(chunk)
                if sink_fh is not None:
                    sink_fh.write(chunk)
                if mode == "logger" and job_logger is not None:
                    text_buffer += chunk.decode(errors="replace")
                    while "\n" in text_buffer:
                        line, text_buffer = text_buffer.split("\n", 1)
                        line = line.rstrip("\r")
                        if line:
                            _emit_text(line)
                if mode == "tee-limited" and job_logger is not None and emitted < 3:
                    text = chunk.decode(errors="ignore")
                    if text.strip():
                        emitted += 1
                        _emit_text(text.strip()[:200])
            if mode == "logger" and job_logger is not None:
                tail = text_buffer.rstrip("\r")
                if tail:
                    _emit_text(tail)
            if sink_fh is not None:
                sink_fh.flush()
            return total

        drain_out = asyncio.create_task(_drain_stream(process.stdout, stdout_fh, "stdout"))
        drain_err = asyncio.create_task(_drain_stream(process.stderr, stderr_fh, "stderr"))

        try:
            if timeout is None:
                rc = await process.wait()
            else:
                rc = await asyncio.wait_for(process.wait(), timeout=float(timeout))
        except asyncio.TimeoutError:
            timed_out = True
            await self._terminate_with_fallback(process)
            rc = await process.wait()
        finally:
            try:
                stdout_bytes, stderr_bytes = await asyncio.gather(
                    drain_out, drain_err, return_exceptions=False
                )
            finally:
                if stdout_fh is not None:
                    stdout_fh.close()
                if stderr_fh is not None:
                    stderr_fh.close()

        finished_at = _utc_now_iso()
        duration = max(0.0, time.monotonic() - started_ts)
        result = SubprocessResult(
            task_id=task_id,
            returncode=int(rc),
            started_at=started_at,
            finished_at=finished_at,
            duration_sec=duration,
            timed_out=timed_out,
            stdout_path=stdout_path,
            stderr_path=stderr_path,
            stdout_bytes=int(stdout_bytes),
            stderr_bytes=int(stderr_bytes),
            cmd_display=cmd_display,
        )

        if job.log_dir:
            meta_path = str(Path(job.log_dir).expanduser().resolve() / "meta.json")
            payload = {
                "task_id": result.task_id,
                "cmd": result.cmd_display,
                "cwd": job.cwd,
                "returncode": result.returncode,
                "timed_out": result.timed_out,
                "started_at": result.started_at,
                "finished_at": result.finished_at,
                "duration_sec": result.duration_sec,
                "stdout_path": result.stdout_path,
                "stderr_path": result.stderr_path,
                "stdout_bytes": result.stdout_bytes,
                "stderr_bytes": result.stderr_bytes,
                "meta": _safe_json(run_meta),
            }
            try:
                tmp_meta = meta_path + ".tmp"
                with open(tmp_meta, "w", encoding="utf-8") as f:
                    json.dump(payload, f, indent=2, ensure_ascii=False)
                os.replace(tmp_meta, meta_path)
            except Exception:
                if self.logger is not None:
                    self.logger.warning(
                        format_two_column_log(
                            "Failed to write subprocess meta file",
                            [
                                ("meta_path", meta_path),
                                ("error", traceback.format_exc().strip()),
                            ],
                        )
                    )

        return result

    async def _terminate_with_fallback(self, process: asyncio.subprocess.Process) -> None:
        pid = int(getattr(process, "pid", 0) or 0)
        try:
            if pid > 0:
                os.killpg(pid, signal.SIGTERM)
            else:
                process.terminate()
        except Exception:
            try:
                process.terminate()
            except Exception:
                pass

        try:
            await asyncio.wait_for(process.wait(), timeout=self.config.terminate_grace_sec)
            return
        except Exception:
            pass

        try:
            if pid > 0:
                os.killpg(pid, signal.SIGKILL)
            else:
                process.kill()
        except Exception:
            try:
                process.kill()
            except Exception:
                pass

    async def _progress_loop(self):
        while True:
            await asyncio.sleep(self.config.progress_interval_sec)
            if self.logger is None:
                continue
            snap = self.snapshot()
            if snap["submitted"] <= 0:
                continue
            self.logger.warning(
                format_two_column_log(
                    "Subprocess scheduler progress",
                    [
                        ("submitted", snap["submitted"]),
                        ("completed", snap["completed"]),
                        ("failed", snap["failed"]),
                        ("running", snap["running"]),
                        ("pending", snap["pending"]),
                        ("peak_running", snap["peak_running"]),
                    ],
                )
            )

    async def _diagnostics_loop(self):
        while True:
            await asyncio.sleep(self.config.diagnostics_interval_sec)
            self._write_status_snapshot(reason="periodic")

    def _write_status_snapshot(self, reason: str) -> None:
        if not self.status_path:
            return
        status_path = Path(self.status_path).expanduser().resolve()
        status_path.parent.mkdir(parents=True, exist_ok=True)
        snap = self.snapshot()
        payload = {
            "ts": _utc_now_iso(),
            "component": "subprocess_scheduler",
            "reason": str(reason),
            **snap,
        }
        try:
            with open(status_path, "a", encoding="utf-8") as f:
                f.write(json.dumps(payload, ensure_ascii=False) + "\n")
        except Exception:
            if self.logger is not None:
                self.logger.warning(
                    format_two_column_log(
                        "Failed to write diagnostics status file",
                        [
                            ("status_path", status_path),
                            ("error", traceback.format_exc().strip()),
                        ],
                    )
                )
