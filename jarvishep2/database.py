"""Minimal HDF5 DATABASE writer for the D1.1 Archiver MVP."""

from __future__ import annotations

import csv
import glob
import hashlib
import json
import os
import shutil
import subprocess
import threading
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
    try:
        handle = h5py.File(path, "r", libver="latest", swmr=True)
    except OSError:
        handle = h5py.File(path, "r")
    with handle:
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
    try:
        handle = h5py.File(path, "r", libver="latest", swmr=True)
    except OSError:
        handle = h5py.File(path, "r")
    with handle:
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
        try:
            handle = h5py.File(path, "r", libver="latest", swmr=True)
        except OSError:
            handle = h5py.File(path, "r")
        with handle:
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
        try:
            handle = h5py.File(self.db_path, "r", libver="latest", swmr=True)
        except OSError:
            # Compatibility for HDF5 files created before the SWMR writer.
            handle = h5py.File(self.db_path, "r")
        with handle:
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


def _is_stale_swmr_write_flag(error: OSError) -> bool:
    """Whether HDF5 rejected a file solely because a prior SWMR writer crashed."""
    text = str(error).lower()
    return "already open for write" in text and "swmr write" in text


def _clear_stale_swmr_write_flag(db_path: str) -> None:
    """Clear HDF5's stale write-status flag with the HDF5 utility.

    ``h5clear`` ships with the system/HDF5 **command-line tools**, not the
    ``h5py`` pip wheel.  Install the distro package (e.g. ``hdf5-tools`` on
    Debian/Ubuntu, ``hdf5`` on Homebrew) so resume can recover a crashed SWMR
    writer flag.
    """
    h5clear = shutil.which("h5clear")
    if not h5clear:
        raise RuntimeError(
            "h5clear is required to recover a stale HDF5 SWMR write flag "
            f"({db_path}). It is an HDF5 CLI tool, not part of the h5py wheel — "
            "install hdf5-tools (apt) / hdf5 (brew), confirm no live Archiver "
            f"owns the file, then run `h5clear -s {db_path}`"
        )
    completed = subprocess.run(
        [h5clear, "-s", db_path], capture_output=True, text=True, check=False
    )
    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or f"exit {completed.returncode}").strip()
        raise RuntimeError(f"h5clear could not recover {db_path}: {detail}")


def recover_stale_hdf5_write_flag(db_path: str) -> bool:
    """Clear a crashed SWMR writer flag, returning whether recovery occurred.

    Callers must first establish that no live Archiver owns the file.  The
    resume bootstrap does that by cleaning an exact-name stale runtime before
    touching DATABASE.
    """
    path = os.path.abspath(str(db_path))
    if not os.path.isfile(path):
        return False
    try:
        # A SWMR reader is deliberately allowed to open a file whose writer
        # consistency flag is set, so read-only probing cannot distinguish a
        # clean file from a crashed writer.  The resume caller has already
        # removed the exact-name runtime fleet; an append-capable probe is
        # therefore both safe and the only reliable way to test this flag.
        with h5py.File(path, "r+", libver="latest"):
            return False
    except OSError as exc:
        if not _is_stale_swmr_write_flag(exc):
            raise
    _clear_stale_swmr_write_flag(path)
    return True


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
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)
        self._handle: h5py.File | None = None
        self._fallback: SimpleHDF5Writer | None = None
        try:
            self._handle = h5py.File(self.db_path, "a", libver="latest")
            SimpleHDF5Writer._ensure_records_dataset(self._handle)
            self._handle.flush()
            self._handle.swmr_mode = True
        except OSError as exc:
            if not _is_stale_swmr_write_flag(exc):
                raise
            if not recover_stale:
                raise RuntimeError(
                    "stale HDF5 writer state detected; cleanup the owning runtime "
                    "before requesting recovery"
                ) from exc
            _clear_stale_swmr_write_flag(self.db_path)
            if self._logger is not None:
                self._logger.warning(
                    "Recovered stale HDF5 SWMR write flag; resuming append to %s",
                    self.db_path,
                )
            self._handle = h5py.File(self.db_path, "a", libver="latest")
            SimpleHDF5Writer._ensure_records_dataset(self._handle)
            self._handle.flush()
            self._handle.swmr_mode = True
        except RuntimeError:
            # Existing pre-SWMR files cannot change their superblock in place.
            # Retain compatibility while still writing an entire persistence
            # batch per open/close cycle instead of one cycle per record.
            if self._handle is not None:
                self._handle.close()
            self._handle = None
            self._fallback = SimpleHDF5Writer(self.db_path)
        self.records_persisted = (
            int(self._handle["records"].shape[0]) if self._handle is not None else 0
        )

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
    ) -> None:
        self.db_path = os.path.abspath(str(db_path))
        self.max_bytes = max(1, int(max_bytes))
        self._lock = threading.Lock()
        self._logger = logger
        self._writer = StreamingHDF5Writer(self.db_path, logger=logger)
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
            self._writer = StreamingHDF5Writer(self.db_path, logger=self._logger)

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
    "recover_stale_hdf5_write_flag",
    "ensure_samples_csv_shards",
]
