"""Minimal HDF5 DATABASE writer for the D1.1 Archiver MVP."""

from __future__ import annotations

import csv
import errno
import glob
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import threading
import time
from collections.abc import Mapping, Sequence
from typing import Any

import h5py
import numpy as np


def make_json_compatible(value: Any) -> Any:
    """Recursively convert values into JSON-serializable Python objects."""
    if isinstance(value, Mapping):
        return {str(k): make_json_compatible(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [make_json_compatible(v) for v in value]
    if isinstance(value, np.ndarray):
        return make_json_compatible(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def _csv_cell(value: Any) -> Any:
    """Flatten a record value for CSV (scalars as-is; nested → JSON text)."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return value
    if isinstance(value, (str, int, float)):
        return value
    if hasattr(value, "item") and callable(value.item):
        try:
            return value.item()
        except Exception:
            pass
    try:
        return json.dumps(value, ensure_ascii=False, default=str)
    except Exception:
        return str(value)


def _fieldnames_for_records(
    records: Sequence[Mapping[str, Any]],
    *,
    preferred: Sequence[str] | None = None,
) -> list[str]:
    preferred_keys = list(preferred or ("x", "y", "LogL", "uuid"))
    fieldnames: list[str] = []
    seen: set[str] = set()
    for key in preferred_keys:
        if any(key in row for row in records):
            fieldnames.append(str(key))
            seen.add(str(key))
    for row in records:
        for key in row.keys():
            text = str(key)
            if text not in seen:
                seen.add(text)
                fieldnames.append(text)
    return fieldnames


def _file_md5(path: str) -> str:
    """Return a content digest for one HDF5 export source."""
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _read_export_md5(path: str) -> str | None:
    try:
        with open(path, encoding="ascii") as handle:
            value = handle.read().strip().lower()
    except OSError:
        return None
    return value if len(value) == 32 and all(char in "0123456789abcdef" for char in value) else None


def _write_export_md5(path: str, digest: str) -> None:
    tmp = path + ".tmp"
    try:
        with open(tmp, "w", encoding="ascii") as handle:
            handle.write(digest + "\n")
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def convert_hdf5_to_csv(
    hdf5_path: str,
    csv_path: str | None = None,
    *,
    preferred_columns: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Convert one V2 DATABASE HDF5 (JSON-string rows) into a flat CSV.

    Returns a status dict with keys:
      ``hdf5``, ``csv``, ``status`` (``converted`` | ``skipped_unchanged`` |
      ``empty`` | ``missing``), and ``rows`` when converted.
    """
    source = os.path.abspath(str(hdf5_path))
    target = os.path.abspath(str(csv_path or (os.path.splitext(source)[0] + ".csv")))
    result: dict[str, Any] = {"hdf5": source, "csv": target, "status": "missing", "rows": 0}

    if not os.path.isfile(source):
        result["status"] = "missing"
        return result

    source_md5 = _file_md5(source)
    md5_path = target + ".md5"
    if os.path.isfile(target) and _read_export_md5(md5_path) == source_md5:
        result["status"] = "skipped_unchanged"
        result["md5"] = source_md5
        return result

    writer = SimpleHDF5Writer(source)
    records = writer.read_records()
    if not records:
        result["status"] = "empty"
        return result

    fieldnames = _fieldnames_for_records(records, preferred=preferred_columns)
    if not fieldnames:
        result["status"] = "empty"
        return result

    parent = os.path.dirname(target)
    if parent:
        os.makedirs(parent, exist_ok=True)
    tmp = target + ".tmp"
    try:
        with open(tmp, "w", encoding="utf-8", newline="") as handle:
            csv_writer = csv.DictWriter(handle, fieldnames=fieldnames)
            csv_writer.writeheader()
            for row in records:
                csv_writer.writerow({key: _csv_cell(row.get(key)) for key in fieldnames})
        os.replace(tmp, target)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise

    result["status"] = "converted"
    result["rows"] = len(records)
    if _file_md5(source) == source_md5:
        _write_export_md5(md5_path, source_md5)
        result["md5"] = source_md5
    else:
        # The CSV is still a valid point-in-time snapshot, but must be
        # regenerated next time so it is never marked as current incorrectly.
        result["source_changed_during_export"] = True
    return result


def discover_database_hdf5(database_dir: str) -> list[str]:
    """List convert targets under ``DATABASE/`` (``*.hdf5``, no ``.snap``)."""
    root = os.path.abspath(str(database_dir))
    if not os.path.isdir(root):
        return []
    paths: list[str] = []
    for path in sorted(glob.glob(os.path.join(root, "*.hdf5"))):
        if path.endswith(".snap") or path.endswith(".hdf5.snap"):
            continue
        if os.path.isfile(path):
            paths.append(os.path.abspath(path))
    return paths


def discover_samples_hdf5_shards(database_dir: str) -> list[str]:
    """Return scan-record shards in chronological order, then the live file.

    A completed shard is named ``samples.0001.hdf5`` (and so on); the current
    writer always keeps the compatibility name ``samples.hdf5``.  Keeping the
    active file separate means readers never need to open an unbounded HDF5
    archive.
    """
    root = os.path.abspath(str(database_dir))
    if not os.path.isdir(root):
        return []
    shards: list[tuple[int, str]] = []
    for path in glob.glob(os.path.join(root, "samples.*.hdf5")):
        stem = os.path.basename(path)[len("samples.") : -len(".hdf5")]
        if stem.isdigit() and os.path.isfile(path):
            shards.append((int(stem), os.path.abspath(path)))
    ordered = [path for _, path in sorted(shards)]
    active = os.path.join(root, "samples.hdf5")
    if os.path.isfile(active):
        ordered.append(os.path.abspath(active))
    return ordered


def _open_hdf5_readonly(path: str) -> h5py.File:
    """Open DATABASE for control-side reads without blocking the Archiver writer.

    Cross-process note (macOS/HDF5): a reader opened with default file locking
    holds an advisory lock that makes Archiver's append open fail with
    ``errno 35 unable to lock file``.  Prefer ``locking=False``.

    Same-process note: if a writer already has the file open with default
    locking, HDF5 rejects a ``locking=False`` reader (``locking flag values
    don't match``).  Fall back to the writer's locking mode in that case.
    """
    attempts: list[dict[str, Any]] = [
        {"libver": "latest", "swmr": True, "locking": False},
        {"libver": "latest", "swmr": True},
        {"locking": False},
        {},
    ]
    last_error: BaseException | None = None
    for kwargs in attempts:
        try:
            return h5py.File(path, "r", **kwargs)
        except TypeError as exc:
            # Older h5py may not accept locking= / swmr=.
            last_error = exc
            continue
        except OSError as exc:
            last_error = exc
            text = str(exc).lower()
            if "locking flag" in text or "swmr" in text:
                continue
            # Fall through other OSErrors to try a simpler open once more.
            continue
    if last_error is not None:
        raise last_error
    return h5py.File(path, "r")

def read_hdf5_uuids(hdf5_path: str) -> set[str]:
    """Read only the UUID field from one JSON-row HDF5 database.

    DATABASE is the resume authority.  Keep this reader deliberately
    independent of Redis/checkpoint state so it also works after a power loss
    or when the checkpoint was deleted.
    """
    path = os.path.abspath(str(hdf5_path))
    if not os.path.isfile(path):
        return set()
    uuids: set[str] = set()
    with _open_hdf5_readonly(path) as handle:
        if "records" not in handle:
            return uuids
        for item in handle["records"]:
            if isinstance(item, bytes):
                item = item.decode("utf-8")
            try:
                record = json.loads(item)
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            uuid = str(record.get("uuid") or "").strip() if isinstance(record, Mapping) else ""
            if uuid:
                uuids.add(uuid)
    return uuids


def read_persisted_uuids(database_dir: str) -> set[str]:
    """Return the union of durable UUIDs across sealed and live scan shards."""
    persisted: set[str] = set()
    for path in discover_samples_hdf5_shards(database_dir):
        persisted.update(read_hdf5_uuids(path))
    return persisted


def _record_is_failed(record: Mapping[str, Any]) -> bool:
    """Whether a durable DATABASE row should count as a failed sample."""
    return str(record.get("status") or "Completed").strip() == "Failed"


def read_hdf5_outcome_counts(hdf5_path: str) -> tuple[int, int]:
    """Return ``(completed, failed)`` counts from one JSON-row HDF5 shard.

    Missing ``status`` is treated as ``Completed`` (historical rows).  Only the
    explicit ``Failed`` status increments the failure counter — used so resume
    can seed Redis SAMPLE_STATS honestly (D21.13).
    """
    path = os.path.abspath(str(hdf5_path))
    if not os.path.isfile(path):
        return 0, 0
    completed = 0
    failed = 0
    with _open_hdf5_readonly(path) as handle:
        if "records" not in handle:
            return 0, 0
        for item in handle["records"]:
            if isinstance(item, bytes):
                item = item.decode("utf-8")
            try:
                record = json.loads(item)
            except (TypeError, ValueError, json.JSONDecodeError):
                continue
            if not isinstance(record, Mapping):
                continue
            if _record_is_failed(record):
                failed += 1
            else:
                completed += 1
    return completed, failed


def read_persisted_outcome_counts(database_dir: str) -> tuple[int, int]:
    """Return ``(completed, failed)`` across sealed and live DATABASE shards."""
    completed = 0
    failed = 0
    for path in discover_samples_hdf5_shards(database_dir):
        ok, bad = read_hdf5_outcome_counts(path)
        completed += ok
        failed += bad
    return completed, failed


def read_persisted_sample_index_state(
    database_dir: str,
) -> tuple[int, set[int], set[str]]:
    """Return ``(prefix, out_of_order, legacy_uuids)`` from DATABASE.

    ``prefix`` is the largest durable contiguous logical index range
    ``[0, prefix)``.  The trailing index set is normally bounded by the
    runtime in-flight window; UUIDs are retained only for records written by
    pre-index checkpoint formats. DATABASE is still the sole authority.
    """
    indices: set[int] = set()
    legacy_uuids: set[str] = set()
    for path in discover_samples_hdf5_shards(database_dir):
        with _open_hdf5_readonly(path) as handle:
            if "records" not in handle:
                continue
            for item in handle["records"]:
                if isinstance(item, bytes):
                    item = item.decode("utf-8")
                try:
                    record = json.loads(item)
                except (TypeError, ValueError, json.JSONDecodeError):
                    continue
                if not isinstance(record, Mapping):
                    continue
                raw_index = record.get("sample_index")
                try:
                    index = int(raw_index) if raw_index is not None else -1
                except (TypeError, ValueError):
                    index = -1
                if index >= 0:
                    indices.add(index)
                    continue
                uuid = str(record.get("uuid") or "").strip()
                if uuid:
                    legacy_uuids.add(uuid)
    prefix = 0
    while prefix in indices:
        indices.remove(prefix)
        prefix += 1
    return prefix, indices, legacy_uuids


def ensure_samples_csv_shards(
    database_dir: str,
    *,
    preferred_columns: Sequence[str] | None = None,
) -> list[str]:
    """Return CSV sources for scan shards, refreshing only the live shard.

    Sealed HDF5 shards are converted at rollover time and their CSV is
    immutable.  Re-hashing each multi-gigabyte historical shard when plotting
    would defeat sharding, so only a missing sealed CSV is regenerated.  The
    current ``samples.hdf5`` remains mutable and uses the normal MD5 refresh.
    """
    paths = discover_samples_hdf5_shards(database_dir)
    active = os.path.abspath(os.path.join(str(database_dir), "samples.hdf5"))
    csv_paths: list[str] = []
    for hdf5_path in paths:
        csv_path = os.path.splitext(hdf5_path)[0] + ".csv"
        if os.path.abspath(hdf5_path) == active or not os.path.isfile(csv_path):
            convert_hdf5_to_csv(
                hdf5_path,
                csv_path,
                preferred_columns=preferred_columns,
            )
        if os.path.isfile(csv_path):
            csv_paths.append(os.path.abspath(csv_path))
    return csv_paths


def convert_database_dir(
    database_dir: str,
    *,
    preferred_columns: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """Convert every ``*.hdf5`` under a DATABASE directory to sibling ``*.csv``."""
    results: list[dict[str, Any]] = []
    for hdf5_path in discover_database_hdf5(database_dir):
        results.append(
            convert_hdf5_to_csv(
                hdf5_path,
                preferred_columns=preferred_columns,
            )
        )
    return results


class SimpleHDF5Writer:
    """One-shot HDF5 reader/writer used by utilities and compatibility paths."""

    def __init__(self, db_path: str) -> None:
        self.db_path = str(db_path)
        self._lock = threading.Lock()
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)

    @staticmethod
    def _ensure_records_dataset(handle: h5py.File) -> h5py.Dataset:
        dtype = h5py.string_dtype(encoding="utf-8")
        if "records" not in handle:
            return handle.create_dataset(
                "records",
                shape=(0,),
                maxshape=(None,),
                chunks=True,
                dtype=dtype,
            )
        return handle["records"]

    @staticmethod
    def _encode_records(records: Sequence[Mapping[str, Any]]) -> list[str]:
        return [
            json.dumps(make_json_compatible(dict(record)), ensure_ascii=False)
            for record in records
        ]

    @classmethod
    def _append_records(cls, handle: h5py.File, records: Sequence[Mapping[str, Any]]) -> int:
        payloads = cls._encode_records(records)
        if not payloads:
            return 0
        dataset = cls._ensure_records_dataset(handle)
        start = int(dataset.shape[0])
        end = start + len(payloads)
        dataset.resize((end,))
        dataset[start:end] = payloads
        return len(payloads)

    def add_records(self, records: Sequence[Mapping[str, Any]]) -> int:
        """Append one complete batch while opening the HDF5 file only once."""
        if not records:
            return 0
        with self._lock:
            with h5py.File(self.db_path, "a", libver="latest") as handle:
                return self._append_records(handle, records)

    def add_record(self, observables: Mapping[str, Any]) -> None:
        self.add_records([observables])

    def close(self) -> None:
        """Compatibility no-op: one-shot writes close their own HDF5 handle."""

    def read_records(self, *, limit: int | None = None) -> list[dict[str, Any]]:
        if not os.path.exists(self.db_path):
            return []
        rows: list[dict[str, Any]] = []
        with _open_hdf5_readonly(self.db_path) as handle:
            if "records" not in handle:
                return []
            dataset = handle["records"]
            size = int(dataset.shape[0])
            if limit is not None:
                size = min(size, max(0, int(limit)))
            for item in dataset[:size]:
                if isinstance(item, bytes):
                    item = item.decode("utf-8")
                rows.append(json.loads(item))
        return rows


def _is_stale_swmr_write_flag(error: BaseException) -> bool:
    """Whether HDF5 rejected a file due to a stale superblock write/SWMR flag."""
    text = str(error).lower()
    if "already open for write" in text:
        return True
    if "file consistency" in text or "status_flags" in text:
        return True
    if "may use" in text and "h5clear" in text:
        return True
    return False


def _is_hdf5_lock_error(error: BaseException) -> bool:
    """Whether open failed because the OS/HDF5 file lock is held or stale.

    Seen after crash/Ctrl-C as::

        BlockingIOError: [Errno 35] Unable to … open file
        (unable to lock file, errno = 35, … 'Resource temporarily unavailable')

    Distinct from the SWMR *write consistency flag* message handled by
    :func:`_is_stale_swmr_write_flag`.
    """
    if isinstance(error, BlockingIOError):
        return True
    err_no = getattr(error, "errno", None)
    # Linux EAGAIN/EWOULDBLOCK=11; macOS often reports 35 for the same case.
    if err_no in {errno.EAGAIN, errno.EWOULDBLOCK, 11, 35}:
        text = str(error).lower()
        if "lock" in text or "resource temporarily unavailable" in text:
            return True
    text = str(error).lower()
    return (
        "unable to lock file" in text
        or "resource temporarily unavailable" in text
        or ("lock file" in text and "unavailable" in text)
    )


def _is_recoverable_hdf5_open_error(error: BaseException) -> bool:
    return _is_stale_swmr_write_flag(error) or _is_hdf5_lock_error(error)


def _needs_swmr_status_flag_clear(db_path: str) -> bool:
    """True when append/r+ open fails solely due to a stale write/SWMR flag."""
    path = os.path.abspath(str(db_path))
    if not os.path.isfile(path):
        return False
    open_attempts: list[dict[str, Any]] = [
        {"mode": "r+", "libver": "latest", "locking": False},
        {"mode": "r+", "libver": "latest"},
        {"mode": "a", "libver": "latest"},
    ]
    for kwargs in open_attempts:
        mode = str(kwargs.pop("mode"))
        try:
            with h5py.File(path, mode, **kwargs):
                return False
        except TypeError:
            # Older h5py without locking=
            try:
                with h5py.File(path, mode, libver="latest"):
                    return False
            except OSError as exc:
                return _is_stale_swmr_write_flag(exc)
        except OSError as exc:
            if _is_stale_swmr_write_flag(exc):
                return True
            if _is_hdf5_lock_error(exc):
                # Locked by a live process — not a superblock flag problem.
                return False
            # Try next open style.
            continue
    return False


def _clear_stale_swmr_write_flag_python(db_path: str) -> None:
    """Clear stale write/SWMR superblock flags using pure h5py (no ``h5clear``).

    Strategy: SWMR **readers** may still open a file that has the write
    consistency flag set.  Copy every root attribute and object into a fresh
    file (clean superblock), then atomically replace the original.  Content is
    preserved; only the dirty superblock status bits are dropped.

    This intentionally avoids the external ``h5clear`` CLI (not shipped with
    the h5py wheel) so resume works in a stock ``pip install h5py`` environment.
    """
    path = os.path.abspath(str(db_path))
    parent = os.path.dirname(path) or "."
    os.makedirs(parent, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".samples_swmr_clear_", suffix=".hdf5", dir=parent)
    os.close(fd)
    try:
        # Prefer SWMR read + no file locking — allowed while write flag is set.
        try:
            src = h5py.File(path, "r", libver="latest", swmr=True, locking=False)
        except TypeError:
            src = h5py.File(path, "r", libver="latest", swmr=True)
        except OSError:
            try:
                src = h5py.File(path, "r", locking=False)
            except TypeError:
                src = h5py.File(path, "r")
        try:
            with h5py.File(tmp, "w", libver="latest") as dst:
                for key, value in src.attrs.items():
                    dst.attrs[key] = value
                for name in src.keys():
                    src.copy(name, dst, name=name)
        finally:
            src.close()
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _clear_stale_swmr_write_flag(db_path: str) -> str:
    """Clear stale SWMR/write superblock flags; return action label.

    Primary path is pure Python (:func:`_clear_stale_swmr_write_flag_python`).
    Optional fallback to the system ``h5clear -s`` CLI when the rewrite path
    fails (e.g. exotic superblock / truncated file).
    """
    path = os.path.abspath(str(db_path))
    try:
        _clear_stale_swmr_write_flag_python(path)
        return "swmr_flag_rewrite"
    except Exception as py_exc:
        h5clear = shutil.which("h5clear")
        if not h5clear:
            raise RuntimeError(
                "failed to clear stale HDF5 write/SWMR flag with pure Python "
                f"({path}): {py_exc}"
            ) from py_exc
        completed = subprocess.run(
            [h5clear, "-s", path], capture_output=True, text=True, check=False
        )
        if completed.returncode != 0:
            detail = (
                completed.stderr or completed.stdout or f"exit {completed.returncode}"
            ).strip()
            raise RuntimeError(
                f"pure-Python SWMR clear failed ({py_exc}); "
                f"h5clear -s also failed for {path}: {detail}"
            ) from py_exc
        return "swmr_flag_h5clear_fallback"


def _list_file_holder_pids(path: str) -> list[int]:
    """Best-effort PIDs that currently have *path* open (via ``lsof``)."""
    if not shutil.which("lsof"):
        return []
    try:
        completed = subprocess.run(
            ["lsof", "-t", "--", path],
            capture_output=True,
            text=True,
            check=False,
            timeout=5.0,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []
    pids: list[int] = []
    for line in (completed.stdout or "").splitlines():
        text = line.strip()
        if not text:
            continue
        try:
            pids.append(int(text))
        except ValueError:
            continue
    return sorted(set(pids))


def _pid_is_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        # Exists but we cannot signal it — treat as live.
        return True
    except OSError:
        return False
    return True


def _remove_sidecar_lock_files(db_path: str) -> list[str]:
    """Delete known HDF5 sidecar lock files if present (best effort)."""
    removed: list[str] = []
    for candidate in (
        f"{db_path}.lock",
        f"{db_path}.lck",
        f"{db_path}__lock",
    ):
        try:
            if os.path.isfile(candidate):
                os.unlink(candidate)
                removed.append(candidate)
        except OSError:
            pass
    return removed


def _break_stale_file_lock_by_inode_replace(db_path: str) -> None:
    """Copy+atomic-replace so a new open is not blocked by a stale flock.

    POSIX advisory locks are bound to the open file description / inode.  When
    a crashed writer left the kernel thinking the file is locked (or a
    leftover process died without an unlock the VFD still believes), rewriting
    the path onto a fresh inode lets a new Archiver open it.  **Only** call
    this after confirming no live process holds the path.
    """
    path = os.path.abspath(str(db_path))
    parent = os.path.dirname(path) or "."
    os.makedirs(parent, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".samples_relock_", suffix=".hdf5", dir=parent)
    os.close(fd)
    try:
        shutil.copy2(path, tmp)
        os.replace(tmp, path)
    except Exception:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def _probe_hdf5_append_open(db_path: str) -> None:
    """Open for append briefly and close; raises OSError on failure."""
    with h5py.File(db_path, "a", libver="latest") as handle:
        # Touch the handle so lazy open failures surface here.
        _ = handle.id


def _process_command(pid: int) -> str:
    try:
        completed = subprocess.run(
            ["ps", "-p", str(int(pid)), "-o", "command="],
            capture_output=True,
            text=True,
            check=False,
            timeout=3.0,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return (completed.stdout or "").strip()


def _is_reclaimable_hdf5_holder(pid: int, *, scan_name: str | None) -> bool:
    """True when *pid* is a stale Jarvis Archiver (safe to SIGTERM on resume)."""
    cmd = _process_command(pid)
    if not cmd:
        return False
    first = cmd.split(None, 1)[0]
    # Only reclaim Archiver holders automatically — never the live Control.
    if not first.startswith("Jarvis2-Archiver"):
        return False
    if not scan_name:
        return True
    return f":{scan_name}" in cmd or cmd.endswith(f":{scan_name}")


def _terminate_pids(pids: Sequence[int], *, logger: Any | None = None) -> list[int]:
    """SIGTERM then SIGKILL reclaimable holders; return pids that exited."""
    import signal

    gone: list[int] = []
    for pid in pids:
        try:
            os.kill(int(pid), signal.SIGTERM)
        except ProcessLookupError:
            gone.append(int(pid))
        except OSError:
            continue
    deadline = time.monotonic() + 2.0
    remaining = {int(p) for p in pids if int(p) not in gone}
    while remaining and time.monotonic() < deadline:
        time.sleep(0.05)
        for pid in list(remaining):
            if not _pid_is_alive(pid):
                remaining.discard(pid)
                gone.append(pid)
    for pid in list(remaining):
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            gone.append(pid)
        except OSError:
            pass
        if not _pid_is_alive(pid):
            gone.append(pid)
    if logger is not None and gone:
        try:
            logger.warning("HDF5 recovery: terminated lock holder pid(s) %s", sorted(set(gone)))
        except Exception:
            pass
    return sorted(set(gone))


def prepare_hdf5_database_for_writer(
    db_path: str,
    *,
    logger: Any | None = None,
    retries: int = 8,
    retry_sleep_sec: float = 0.12,
    scan_name: str | None = None,
    probe_append: bool = True,
    allow_holder_pids: Sequence[int] | None = None,
) -> list[str]:
    """Make ``samples.hdf5`` openable for a new Archiver writer.

    Handles both:

    * SWMR / write *file consistency flags* (pure-Python rewrite; optional
      ``h5clear -s`` fallback only if rewrite fails)
    * OS/HDF5 *file lock* (``errno 35`` / ``unable to lock file``)

    Callers must first stop any live Archiver/control fleet for this scan
    (``ensure_scan_name_available(..., cleanup_stale=True)``).  Returns a list
    of recovery action labels that were applied (empty if already openable).
    Emits **at most one** summary WARNING when recovery actually ran.

    Parameters
    ----------
    probe_append:
        When False (control bootstrap), only run flag/lock recovery steps and
        **do not** open the file for append in this process — so the Control
        never holds a write lock that would block the Archiver child.
    allow_holder_pids:
        PIDs that may keep the file open (e.g. parent Control while Archiver
        starts).  They are not auto-killed; we wait/retry instead.
    """
    path = os.path.abspath(str(db_path))
    if not os.path.isfile(path):
        return []

    actions: list[str] = []
    last_error: BaseException | None = None
    protected = {int(p) for p in (allow_holder_pids or []) if int(p) > 0}
    protected.add(os.getpid())

    def _log_summary() -> None:
        if logger is None or not actions:
            return
        try:
            logger.warning(
                "Recovered HDF5 database (%s): %s",
                ", ".join(actions),
                path,
            )
        except Exception:
            pass

    def _try_clear_swmr_flag() -> bool:
        """Clear superblock write/SWMR flags when present. Return True if ran."""
        if not _needs_swmr_status_flag_clear(path):
            return False
        label = _clear_stale_swmr_write_flag(path)
        actions.append(label)
        return True

    for attempt in range(max(1, int(retries))):
        if probe_append:
            try:
                _probe_hdf5_append_open(path)
                _log_summary()
                return actions
            except OSError as exc:
                last_error = exc
                if not _is_recoverable_hdf5_open_error(exc):
                    raise
        else:
            # Control bootstrap: never append-open here. Recover only when the
            # superblock flag is actually dirty, reclaim orphan Archivers, done.
            if attempt == 0 and not actions:
                try:
                    _try_clear_swmr_flag()
                except RuntimeError as clear_exc:
                    # Non-fatal for control: Archiver will retry with probe_append.
                    if logger is not None:
                        try:
                            logger.warning(
                                "HDF5 SWMR flag clear deferred to Archiver: %s",
                                clear_exc,
                            )
                        except Exception:
                            pass
                removed = _remove_sidecar_lock_files(path)
                if removed:
                    actions.append("removed_sidecar_locks")
                holders = _list_file_holder_pids(path)
                reclaimable = [
                    pid
                    for pid in holders
                    if pid not in protected
                    and _pid_is_alive(pid)
                    and _is_reclaimable_hdf5_holder(pid, scan_name=scan_name)
                ]
                if reclaimable:
                    killed = _terminate_pids(reclaimable, logger=None)
                    if killed:
                        actions.append(
                            f"killed_holders:{','.join(str(p) for p in killed)}"
                        )
                _log_summary()
                return actions
            _log_summary()
            return actions

        # --- recovery steps (ordered, cheapest first) ---
        if last_error is not None and _is_stale_swmr_write_flag(last_error):
            try:
                # Force clear when the open error itself reported the flag.
                label = _clear_stale_swmr_write_flag(path)
                if label not in actions:
                    actions.append(label)
                continue
            except RuntimeError as clear_exc:
                if not _is_hdf5_lock_error(last_error):
                    raise clear_exc
        elif attempt == 0:
            try:
                if _try_clear_swmr_flag():
                    continue
            except RuntimeError:
                pass

        removed = _remove_sidecar_lock_files(path)
        if removed:
            actions.append("removed_sidecar_locks")
            continue

        holders = _list_file_holder_pids(path)
        live = [
            pid
            for pid in holders
            if _pid_is_alive(pid) and pid not in protected
        ]
        reclaimable = [
            pid for pid in live if _is_reclaimable_hdf5_holder(pid, scan_name=scan_name)
        ]
        if reclaimable:
            killed = _terminate_pids(reclaimable, logger=None)
            if killed:
                actions.append(f"killed_holders:{','.join(str(p) for p in killed)}")
                continue

        blocking = [pid for pid in live if pid not in reclaimable]
        if blocking and attempt >= max(1, int(retries) - 2):
            raise RuntimeError(
                f"HDF5 database is locked by live process(es) {blocking}: {path}. "
                "Stop the owning process (e.g. `Jarvis2 kill ZP -y`) then retry. "
                "Control-side DATABASE reads must use locking=False so they do not "
                f"block Archiver. Last open error: {last_error}"
            ) from last_error

        if (
            not blocking
            and (attempt >= 1 or (last_error is not None and _is_hdf5_lock_error(last_error)))
        ):
            try:
                _break_stale_file_lock_by_inode_replace(path)
                actions.append("inode_replace_break_lock")
                continue
            except OSError as replace_exc:
                last_error = replace_exc

        time.sleep(max(0.0, float(retry_sleep_sec)) * (1 + attempt))

    if not probe_append:
        _log_summary()
        return actions

    try:
        _probe_hdf5_append_open(path)
        _log_summary()
        return actions
    except OSError as exc:
        raise RuntimeError(
            f"HDF5 database still not openable for append after recovery: {path}. "
            f"Last error: {exc}. Tried actions: {actions or ['none']}."
        ) from exc


def recover_stale_hdf5_write_flag(db_path: str) -> bool:
    """Backward-compatible resume probe: recover SWMR flag **and** stale locks.

    Returns True when any recovery action was applied.  Prefer
    :func:`prepare_hdf5_database_for_writer` for new call sites.
    """
    actions = prepare_hdf5_database_for_writer(db_path, probe_append=True)
    return bool(actions)


class StreamingHDF5Writer:
    """Single-writer HDF5 stream with batch publication and explicit fsync.

    The Archiver owns this object for its full lifetime. New files use HDF5
    SWMR so a converter can read batches already published by ``flush``.
    """

    def __init__(
        self,
        db_path: str,
        *,
        logger: Any | None = None,
        recover_stale: bool = False,
    ) -> None:
        self.db_path = str(db_path)
        self._logger = logger
        self._lock = threading.Lock()
        self._pending: list[Mapping[str, Any]] = []
        self._batch_open = False
        parent = os.path.dirname(self.db_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        self._handle: h5py.File | None = None
        self._fallback: SimpleHDF5Writer | None = None
        self._open_streaming_handle(recover_stale=recover_stale)
        self.records_persisted = (
            int(self._handle["records"].shape[0]) if self._handle is not None else 0
        )

    def _try_open_swmr_writer(self) -> None:
        self._handle = h5py.File(self.db_path, "a", libver="latest")
        SimpleHDF5Writer._ensure_records_dataset(self._handle)
        self._handle.flush()
        self._handle.swmr_mode = True

    def _open_streaming_handle(self, *, recover_stale: bool) -> None:
        try:
            self._try_open_swmr_writer()
            return
        except OSError as exc:
            if _is_recoverable_hdf5_open_error(exc):
                if not recover_stale:
                    raise RuntimeError(
                        "stale HDF5 writer state or file lock detected; cleanup "
                        "the owning runtime before requesting recovery"
                    ) from exc
                # prepare_hdf5_database_for_writer emits the single recovery
                # summary when actions run; do not log again here.
                prepare_hdf5_database_for_writer(
                    self.db_path, logger=self._logger
                )
                try:
                    self._try_open_swmr_writer()
                    return
                except OSError:
                    raise
                except RuntimeError:
                    # Pre-SWMR superblock — fall through to one-shot writer.
                    pass
            else:
                # Unknown OSError — do not silently swallow.
                raise
        except RuntimeError:
            # Existing pre-SWMR files cannot change their superblock in place.
            pass
        if self._handle is not None:
            try:
                self._handle.close()
            except Exception:
                pass
        self._handle = None
        self._fallback = SimpleHDF5Writer(self.db_path)

    @staticmethod
    def _fsync_descriptor(vfd_handle: Any) -> None:
        if isinstance(vfd_handle, tuple):
            vfd_handle = vfd_handle[0]
        if not isinstance(vfd_handle, int):
            raise OSError("HDF5 driver does not expose a file descriptor for fsync")
        os.fsync(vfd_handle)

    def add_records(self, records: Sequence[Mapping[str, Any]]) -> int:
        """Append, publish, and physically sync one durable Archiver batch."""
        if not records:
            return 0
        with self._lock:
            if self._handle is None:
                if self._fallback is None:
                    raise RuntimeError("HDF5 writer is closed")
                written = self._fallback.add_records(records)
                with open(self.db_path, "rb", buffering=0) as handle:
                    os.fsync(handle.fileno())
                self.records_persisted += written
                return written
            written = SimpleHDF5Writer._append_records(self._handle, records)
            self._handle["records"].flush()
            self._handle.flush()
            self._fsync_descriptor(self._handle.id.get_vfd_handle())
            self.records_persisted += written
            return written

    def begin_batch(self) -> None:
        """Start collecting one Archiver batch without touching the file."""
        with self._lock:
            if self._batch_open:
                raise RuntimeError("HDF5 writer batch is already open")
            self._pending = []
            self._batch_open = True

    def add_record(self, observables: Mapping[str, Any]) -> None:
        """Queue a record in an open batch, or immediately persist it otherwise."""
        with self._lock:
            if self._batch_open:
                self._pending.append(dict(observables))
                return
        self.add_records([observables])

    def commit_batch(self) -> int:
        """Persist all queued records as one flush + fsync transaction."""
        with self._lock:
            if not self._batch_open:
                raise RuntimeError("HDF5 writer batch is not open")
            pending = self._pending
            self._pending = []
            self._batch_open = False
        return self.add_records(pending)

    def abort_batch(self) -> None:
        """Discard records that have not yet been written to HDF5."""
        with self._lock:
            self._pending = []
            self._batch_open = False

    def close(self) -> None:
        with self._lock:
            if self._handle is not None:
                self._handle.close()
                self._handle = None
            self._fallback = None


class RollingHDF5Writer:
    """Bounded ``samples.hdf5`` writer with immutable HDF5/CSV shards.

    After a durable batch pushes the live HDF5 above ``max_bytes``, it is
    closed, renamed to ``samples.NNNN.hdf5``, and converted to a sibling CSV
    before a fresh live ``samples.hdf5`` writer is opened.  The rollover point
    is always a batch boundary, so a shard is self-contained and safe for
    conversion.
    """

    def __init__(
        self,
        db_path: str,
        *,
        max_bytes: int = 1024 * 1024 * 1024,
        logger: Any | None = None,
        recover_stale: bool = True,
    ) -> None:
        self.db_path = os.path.abspath(str(db_path))
        self.max_bytes = max(1, int(max_bytes))
        self._lock = threading.Lock()
        self._logger = logger
        self._recover_stale = bool(recover_stale)
        self._writer = StreamingHDF5Writer(
            self.db_path, logger=logger, recover_stale=self._recover_stale
        )
        self._records_persisted = 0

    @property
    def records_persisted(self) -> int:
        return self._records_persisted

    def _next_shard_path(self) -> str:
        directory = os.path.dirname(self.db_path)
        highest = 0
        for path in glob.glob(os.path.join(directory, "samples.*.hdf5")):
            stem = os.path.basename(path)[len("samples.") : -len(".hdf5")]
            if stem.isdigit():
                highest = max(highest, int(stem))
        return os.path.join(directory, f"samples.{highest + 1:04d}.hdf5")

    def _rollover_if_needed_locked(self) -> None:
        if not os.path.isfile(self.db_path) or os.path.getsize(self.db_path) <= self.max_bytes:
            return
        self._writer.close()
        shard_path = self._next_shard_path()
        os.replace(self.db_path, shard_path)
        try:
            convert_hdf5_to_csv(shard_path)
        finally:
            # The former live CSV describes the just-sealed file; never let it
            # masquerade as data from the newly-created live shard.
            for suffix in (".csv", ".csv.md5"):
                try:
                    os.unlink(os.path.splitext(self.db_path)[0] + suffix)
                except FileNotFoundError:
                    pass
            self._writer = StreamingHDF5Writer(
                self.db_path,
                logger=self._logger,
                recover_stale=self._recover_stale,
            )

    def add_records(self, records: Sequence[Mapping[str, Any]]) -> int:
        with self._lock:
            written = self._writer.add_records(records)
            self._records_persisted += written
            self._rollover_if_needed_locked()
            return written

    def begin_batch(self) -> None:
        self._writer.begin_batch()

    def add_record(self, observables: Mapping[str, Any]) -> None:
        self._writer.add_record(observables)

    def commit_batch(self) -> int:
        with self._lock:
            written = self._writer.commit_batch()
            self._records_persisted += written
            self._rollover_if_needed_locked()
            return written

    def abort_batch(self) -> None:
        self._writer.abort_batch()

    def close(self) -> None:
        with self._lock:
            self._writer.close()
            # Final live shard stays under the compatibility name but gets its
            # CSV export as soon as the scan has finished.
            if os.path.isfile(self.db_path):
                convert_hdf5_to_csv(self.db_path)


__all__ = [
    "SimpleHDF5Writer",
    "StreamingHDF5Writer",
    "RollingHDF5Writer",
    "make_json_compatible",
    "convert_hdf5_to_csv",
    "convert_database_dir",
    "discover_database_hdf5",
    "discover_samples_hdf5_shards",
    "read_hdf5_uuids",
    "read_hdf5_outcome_counts",
    "read_persisted_outcome_counts",
    "read_persisted_sample_index_state",
    "read_persisted_uuids",
    "prepare_hdf5_database_for_writer",
    "recover_stale_hdf5_write_flag",
    "ensure_samples_csv_shards",
]
