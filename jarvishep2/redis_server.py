#!/usr/bin/env python3
"""Manage a local redis-server with an OS-visible process title.

External ``brew services`` Redis cannot be renamed from outside. When Jarvis
starts the server itself, it launches via ``exec -a`` so ``ps`` / Activity
Monitor show::

    Jarvis-Redis:<scan-name>

Lifecycle is optional: if something already answers on the target host:port,
we leave it alone (title is not changed).
"""

from __future__ import annotations

import os
import shlex
import shutil
import signal
import socket
import subprocess
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Any

from jarvishep2.proc_title import redis_title


def redis_port_open(host: str, port: int, *, timeout: float = 0.35) -> bool:
    """Return True when a TCP connect to host:port succeeds."""
    try:
        with socket.create_connection((str(host), int(port)), timeout=float(timeout)):
            return True
    except OSError:
        return False


def find_available_redis_port(host: str, start_port: int, *, attempts: int = 1000) -> int:
    """Find a locally bindable TCP port at or above ``start_port``.

    The bind probe is more reliable than a connect-only check for choosing a
    replacement port.  Redis still owns the final bind, so callers must retain
    normal startup error handling for the small unavoidable race window.
    """
    bind_host = "127.0.0.1" if str(host) in {"localhost", "127.0.0.1", "::1"} else str(host)
    first = max(1, int(start_port))
    last = min(65535, first + max(1, int(attempts)) - 1)
    for port in range(first, last + 1):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
            probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            try:
                probe.bind((bind_host, port))
            except OSError:
                continue
            return port
    raise RuntimeError(f"no free Redis port in range {first}..{last} on {host}")


def find_redis_server_binary() -> str | None:
    path = shutil.which("redis-server")
    return path if path else None


def write_redis_conf(
    conf_path: str,
    *,
    host: str,
    port: int,
    data_dir: str,
    pidfile: str,
    logfile: str,
) -> str:
    """Write a minimal non-daemon redis conf for a scan-local instance."""
    os.makedirs(os.path.dirname(os.path.abspath(conf_path)) or ".", exist_ok=True)
    os.makedirs(data_dir, exist_ok=True)
    bind = "127.0.0.1" if host in {"localhost", "127.0.0.1", "::1"} else str(host)
    body = "\n".join(
        [
            f"bind {bind}",
            "protected-mode yes",
            f"port {int(port)}",
            "daemonize no",
            f"dir {data_dir}",
            f"pidfile {pidfile}",
            f"logfile {logfile}",
            "save \"\"",
            "appendonly no",
            "databases 16",
            "timeout 0",
            "tcp-backlog 511",
            "",
        ]
    )
    with open(conf_path, "w", encoding="utf-8") as handle:
        handle.write(body)
    return conf_path


@dataclass
class ManagedRedisServer:
    """Owns an optional redis-server child started by Jarvis."""

    host: str = "127.0.0.1"
    port: int = 6379
    db: int = 0
    scan_name: str = ""
    work_dir: str = ""
    process: subprocess.Popen[Any] | None = field(default=None, repr=False)
    started_by_us: bool = False
    title: str = ""
    conf_path: str = ""
    pidfile: str = ""
    binary: str = ""

    @classmethod
    def from_redis_config(
        cls,
        redis_config: Mapping[str, Any] | None,
        *,
        scan_name: str | None = None,
        work_dir: str | None = None,
    ) -> ManagedRedisServer:
        cfg = dict(redis_config or {})
        return cls(
            host=str(cfg.get("host") or "127.0.0.1"),
            port=int(cfg.get("port") or 6379),
            db=int(cfg.get("db") or 0),
            scan_name=str(scan_name or "").strip(),
            work_dir=str(work_dir or "").strip(),
        )

    def is_reachable(self) -> bool:
        return redis_port_open(self.host, self.port)

    def ensure(
        self,
        *,
        scan_name: str | None = None,
        work_dir: str | None = None,
        force_start: bool = False,
        ready_timeout: float = 8.0,
    ) -> bool:
        """Ensure a Redis endpoint is available.

        Returns True when this call started a managed redis-server process.
        """
        if scan_name is not None:
            self.scan_name = str(scan_name or "").strip()
        if work_dir is not None:
            self.work_dir = str(work_dir or "").strip()

        if self.is_reachable() and not force_start:
            # Pre-existing service (e.g. brew redis) — process title stays whatever
            # the OS already shows; we cannot rename another process.
            self.started_by_us = False
            return False

        if self.is_reachable() and force_start:
            raise RuntimeError(
                f"cannot force-start managed Redis: {self.host}:{self.port} already in use"
            )

        binary = find_redis_server_binary()
        if not binary:
            raise RuntimeError(
                "redis-server not found on PATH; install Redis or start a local service"
            )
        self.binary = binary

        base = self.work_dir or os.path.join(os.getcwd(), ".jarvis2-redis")
        base = os.path.abspath(base)
        os.makedirs(base, exist_ok=True)
        data_dir = os.path.join(base, "data")
        conf_path = os.path.join(base, "redis.conf")
        pidfile = os.path.join(base, "redis.pid")
        logfile = os.path.join(base, "redis.log")
        write_redis_conf(
            conf_path,
            host=self.host,
            port=self.port,
            data_dir=data_dir,
            pidfile=pidfile,
            logfile=logfile,
        )
        self.conf_path = conf_path
        self.pidfile = pidfile
        self.title = redis_title(scan_name=self.scan_name, port=self.port, db=self.db)

        # argv0 rename so `ps` shows Jarvis-Redis:<scan> instead of redis-server.
        # shell=True + exec -a is the portable macOS/Linux approach for non-Python bins.
        launch = (
            f"exec -a {shlex.quote(self.title)} "
            f"{shlex.quote(binary)} {shlex.quote(conf_path)}"
        )
        self.process = subprocess.Popen(
            ["bash", "-lc", launch],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        self.started_by_us = True

        deadline = time.monotonic() + max(0.5, float(ready_timeout))
        while time.monotonic() < deadline:
            if self.process.poll() is not None:
                detail = ""
                try:
                    if os.path.isfile(logfile):
                        with open(logfile, "r", encoding="utf-8", errors="replace") as handle:
                            detail = handle.read()[-2000:]
                except OSError:
                    pass
                raise RuntimeError(
                    f"managed redis-server exited early (code={self.process.returncode})"
                    + (f":\n{detail}" if detail else "")
                )
            if self.is_reachable():
                return True
            time.sleep(0.05)

        self.stop()
        raise RuntimeError(
            f"managed redis-server did not become ready on {self.host}:{self.port} "
            f"within {ready_timeout:.1f}s"
        )

    def stop(self, *, timeout: float = 3.0) -> None:
        """Stop only a redis-server that this instance started."""
        if not self.started_by_us:
            return
        proc = self.process
        self.started_by_us = False
        self.process = None
        if proc is None:
            # Try pidfile if the Popen handle was lost.
            pid = self._read_pidfile()
            if pid is not None:
                self._terminate_pid(pid, timeout=timeout)
            return
        if proc.poll() is not None:
            return
        try:
            os.killpg(proc.pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                proc.terminate()
            except Exception:
                pass
        try:
            proc.wait(timeout=max(0.1, float(timeout)))
            return
        except subprocess.TimeoutExpired:
            pass
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            try:
                proc.kill()
            except Exception:
                pass
        try:
            proc.wait(timeout=1.0)
        except Exception:
            pass

    def _read_pidfile(self) -> int | None:
        path = str(self.pidfile or "").strip()
        if not path or not os.path.isfile(path):
            return None
        try:
            with open(path, "r", encoding="utf-8") as handle:
                return int(handle.read().strip())
        except (OSError, TypeError, ValueError):
            return None

    @staticmethod
    def _terminate_pid(pid: int, *, timeout: float = 3.0) -> None:
        try:
            os.kill(pid, signal.SIGTERM)
        except (ProcessLookupError, PermissionError, OSError):
            return
        deadline = time.monotonic() + max(0.1, float(timeout))
        while time.monotonic() < deadline:
            try:
                os.kill(pid, 0)
            except OSError:
                return
            time.sleep(0.05)
        try:
            os.kill(pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError, OSError):
            return


def ensure_local_redis(
    *,
    redis_config: Mapping[str, Any] | None = None,
    scan_name: str | None = None,
    work_dir: str | None = None,
    force_start: bool = False,
) -> ManagedRedisServer:
    """Ensure Redis is up; start a titled redis-server when the port is free."""
    managed = ManagedRedisServer.from_redis_config(
        redis_config,
        scan_name=scan_name,
        work_dir=work_dir,
    )
    managed.ensure(
        scan_name=scan_name,
        work_dir=work_dir,
        force_start=force_start,
    )
    return managed


__all__ = [
    "ManagedRedisServer",
    "ensure_local_redis",
    "find_redis_server_binary",
    "redis_port_open",
    "find_available_redis_port",
    "write_redis_conf",
]
