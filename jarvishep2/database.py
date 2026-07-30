"""Minimal HDF5 DATABASE writer for the D1.1 Archiver MVP."""

from __future__ import annotations

import csv
import glob
import json
import os
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


def convert_hdf5_to_csv(
    hdf5_path: str,
    csv_path: str | None = None,
    *,
    force: bool = False,
    preferred_columns: Sequence[str] | None = None,
) -> dict[str, Any]:
    """Convert one V2 DATABASE HDF5 (JSON-string rows) into a flat CSV.

    Returns a status dict with keys:
      ``hdf5``, ``csv``, ``status`` (``converted`` | ``skipped_exists`` |
      ``empty`` | ``missing``), and ``rows`` when converted.
    """
    source = os.path.abspath(str(hdf5_path))
    target = os.path.abspath(str(csv_path or (os.path.splitext(source)[0] + ".csv")))
    result: dict[str, Any] = {"hdf5": source, "csv": target, "status": "missing", "rows": 0}

    if not os.path.isfile(source):
        result["status"] = "missing"
        return result

    if os.path.isfile(target) and not force:
        result["status"] = "skipped_exists"
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
    with open(target, "w", encoding="utf-8", newline="") as handle:
        csv_writer = csv.DictWriter(handle, fieldnames=fieldnames)
        csv_writer.writeheader()
        for row in records:
            csv_writer.writerow({key: _csv_cell(row.get(key)) for key in fieldnames})

    result["status"] = "converted"
    result["rows"] = len(records)
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


def convert_database_dir(
    database_dir: str,
    *,
    force: bool = False,
    preferred_columns: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    """Convert every ``*.hdf5`` under a DATABASE directory to sibling ``*.csv``."""
    results: list[dict[str, Any]] = []
    for hdf5_path in discover_database_hdf5(database_dir):
        results.append(
            convert_hdf5_to_csv(
                hdf5_path,
                force=force,
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

    def read_records(self) -> list[dict[str, Any]]:
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
            for item in handle["records"][()]:
                if isinstance(item, bytes):
                    item = item.decode("utf-8")
                rows.append(json.loads(item))
        return rows


class StreamingHDF5Writer:
    """Single-writer HDF5 stream with batch publication and explicit fsync.

    The Archiver owns this object for its full lifetime. New files use HDF5
    SWMR so a converter can read batches already published by ``flush``.
    """

    def __init__(self, db_path: str) -> None:
        self.db_path = str(db_path)
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


__all__ = [
    "SimpleHDF5Writer",
    "StreamingHDF5Writer",
    "make_json_compatible",
    "convert_hdf5_to_csv",
    "convert_database_dir",
    "discover_database_hdf5",
]
