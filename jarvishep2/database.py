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
    """Append JSON-encoded observables rows to samples.hdf5."""

    def __init__(self, db_path: str) -> None:
        self.db_path = str(db_path)
        self._lock = threading.Lock()
        os.makedirs(os.path.dirname(self.db_path), exist_ok=True)

    def add_record(self, observables: Mapping[str, Any]) -> None:
        payload = json.dumps(make_json_compatible(dict(observables)), ensure_ascii=False)
        with self._lock:
            with h5py.File(self.db_path, "a") as handle:
                dtype = h5py.string_dtype(encoding="utf-8")
                if "records" not in handle:
                    handle.create_dataset(
                        "records",
                        shape=(0,),
                        maxshape=(None,),
                        chunks=True,
                        dtype=dtype,
                    )
                records = handle["records"]
                start = int(records.shape[0])
                end = start + 1
                records.resize((end,))
                records[start:end] = payload

    def read_records(self) -> list[dict[str, Any]]:
        if not os.path.exists(self.db_path):
            return []
        rows: list[dict[str, Any]] = []
        with h5py.File(self.db_path, "r") as handle:
            if "records" not in handle:
                return []
            for item in handle["records"][()]:
                if isinstance(item, bytes):
                    item = item.decode("utf-8")
                rows.append(json.loads(item))
        return rows


__all__ = [
    "SimpleHDF5Writer",
    "make_json_compatible",
    "convert_hdf5_to_csv",
    "convert_database_dir",
    "discover_database_hdf5",
]
