"""Runtime scan metadata published through the scan's own Redis instance."""

from __future__ import annotations

import json
import os
import tempfile
import time
from collections.abc import Mapping
from typing import Any


RUNTIME_METADATA_KEY = "hep:runtime:metadata_path"


def write_scan_metadata(*, config: Mapping[str, Any], info: Mapping[str, Any], redis: Mapping[str, Any]) -> str:
    """Write the non-secret scan identity record and return its absolute path."""
    task_result_dir = os.path.abspath(str(info.get("task_result_dir") or config.get("task_result_dir") or os.getcwd()))
    directory = os.path.join(task_result_dir, ".jarvis2")
    os.makedirs(directory, exist_ok=True)
    path = os.path.join(directory, "runtime.json")
    payload: dict[str, Any] = {
        "schema": 1,
        "written_at": time.time(),
        "scan_name": str(info.get("scan_name") or config.get("scan_name") or "scan"),
        "task_yaml": str(config.get("task_yaml") or ""),
        "project_root": str(config.get("project_root") or config.get("task_root") or ""),
        "task_result_dir": task_result_dir,
        "redis": {"host": str(redis["host"]), "port": int(redis["port"]), "db": int(redis["db"])},
        "control_pid": os.getpid(),
    }
    fd, temporary = tempfile.mkstemp(prefix=".runtime.", suffix=".json", dir=directory)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise
    return path


def read_scan_metadata(path: str, *, redis: Mapping[str, Any], expected_scan: str) -> dict[str, Any] | None:
    """Read and validate the record advertised by Redis before using it."""
    try:
        with open(path, encoding="utf-8") as handle:
            payload = json.load(handle)
    except (OSError, ValueError, TypeError):
        return None
    if not isinstance(payload, dict) or payload.get("schema") != 1:
        return None
    if str(payload.get("scan_name") or "") != str(expected_scan):
        return None
    recorded = payload.get("redis")
    if not isinstance(recorded, Mapping):
        return None
    if (str(recorded.get("host")), int(recorded.get("port", -1)), int(recorded.get("db", -1))) != (
        str(redis.get("host")), int(redis.get("port", -2)), int(redis.get("db", -2))
    ):
        return None
    return payload
