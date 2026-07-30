#!/usr/bin/env python3
"""Set OS-visible process titles (Activity Monitor / ps / top).

``multiprocessing.Process(name=...)`` only changes the Python-side name.
Resource managers usually show ``python3.x`` unless the process title is
updated via ``setproctitle`` (same approach as V1).
"""

from __future__ import annotations


def set_process_title(title: str) -> bool:
    """Set the process title. Returns True if the OS title was updated."""
    text = str(title or "").strip()
    if not text:
        return False
    try:
        import setproctitle  # type: ignore import-not-found

        setproctitle.setproctitle(text)
        return True
    except Exception:
        return False


def control_title(*, scan_name: str | None = None) -> str:
    scan = str(scan_name or "").strip()
    return f"Jarvis2:{scan}" if scan else "Jarvis2"


def worker_title(worker_id: int, *, scan_name: str | None = None) -> str:
    scan = str(scan_name or "").strip()
    base = f"Jarvis2-Worker-{int(worker_id)}"
    return f"{base}:{scan}" if scan else base


def archiver_title(*, scan_name: str | None = None) -> str:
    scan = str(scan_name or "").strip()
    return f"Jarvis2-Archiver:{scan}" if scan else "Jarvis2-Archiver"


def file_operation_title(*, scan_name: str | None = None) -> str:
    """OS title for a Worker-owned file-operation child process."""
    scan = str(scan_name or "").strip()
    return f"Jarvis2-FileOperation:{scan}" if scan else "Jarvis2-FileOperation"


def redis_title(*, scan_name: str | None = None, port: int | None = None, db: int = 0) -> str:
    """OS process title for a Jarvis-managed redis-server."""
    scan = str(scan_name or "").strip()
    base = f"Jarvis-Redis:{scan}" if scan else "Jarvis-Redis"
    return f"{base}@{int(port)}/{int(db)}" if port is not None else base


__all__ = [
    "archiver_title",
    "control_title",
    "file_operation_title",
    "redis_title",
    "set_process_title",
    "worker_title",
]
