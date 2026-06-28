#!/usr/bin/env python3
"""Minimal HDF5 DATABASE writer for the D1.1 Archiver MVP."""

from __future__ import annotations

import json
import os
import threading
from collections.abc import Mapping
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


__all__ = ["SimpleHDF5Writer", "make_json_compatible"]