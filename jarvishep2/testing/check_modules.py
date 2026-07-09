#!/usr/bin/env python3
"""check-modules helpers (golden verify + CSV → Sample coords)."""

from __future__ import annotations

import glob
import os
from collections.abc import Mapping, Sequence
from typing import Any

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.sample import Sample
from jarvishep2.Sampling.variables import load_variables


def sample_tree_file_sets(sample_root: str) -> list[list[str]]:
    manifests: list[list[str]] = []
    if not os.path.isdir(sample_root):
        return manifests
    for child in sorted(os.listdir(sample_root)):
        child_path = os.path.join(sample_root, child)
        if not os.path.isdir(child_path):
            continue
        files = sorted(
            os.path.relpath(path, child_path)
            for path in glob.glob(os.path.join(child_path, "**", "*"), recursive=True)
            if os.path.isfile(path)
        )
        manifests.append(files)
    return sorted(manifests)


def normalize_check_module_records(
    records: list[dict[str, Any]],
    *,
    keys: Sequence[str] | None = None,
) -> list[dict[str, float]]:
    key_list = list(keys) if keys else ("x", "y", "z", "LogL")
    normalized: list[dict[str, float]] = []
    for row in records:
        item = {key: round(float(row[key]), 12) for key in key_list if key in row}
        normalized.append(item)
    sort_keys = [key for key in key_list if key in (normalized[0] if normalized else {})]
    if not sort_keys:
        return normalized
    return sorted(normalized, key=lambda item: tuple(item.get(key, 0.0) for key in sort_keys))


def verify_check_modules_golden(
    *,
    task_result_dir: str,
    golden: Mapping[str, Any],
    record_keys: Sequence[str] | None = None,
) -> None:
    db_path = os.path.join(task_result_dir, "DATABASE", "samples.hdf5")
    records = SimpleHDF5Writer(db_path).read_records()
    expected_records = list(golden.get("records") or [])
    expected_files = golden.get("sample_files")
    if expected_records:
        normalized = normalize_check_module_records(records, keys=record_keys)
        expected_normalized = normalize_check_module_records(expected_records, keys=record_keys)
        if normalized != expected_normalized:
            raise RuntimeError("check-modules DATABASE parity mismatch")
    if expected_files:
        tree = sample_tree_file_sets(os.path.join(task_result_dir, "SAMPLE"))
        for files in tree:
            if list(files) != list(expected_files):
                raise RuntimeError("check-modules SAMPLE parity mismatch")


def coordinate_columns_for_check_modules(
    config: Mapping[str, Any],
    csv_fieldnames: Sequence[str] | None,
) -> list[str]:
    """Infer u-coord CSV columns from Sampling.Variables, else CSV headers minus uuid."""
    try:
        variables = load_variables(config)
        names = [str(var.name) for var in variables if str(getattr(var, "name", "")).strip()]
        if names:
            return names
    except Exception:
        pass
    if csv_fieldnames:
        return [name for name in csv_fieldnames if str(name).strip() and str(name) != "uuid"]
    return ["x", "y"]


def build_check_module_u_coords(
    row: Mapping[str, Any],
    columns: Sequence[str],
) -> list[float]:
    missing = [name for name in columns if name not in row]
    if missing:
        raise KeyError(f"check-modules CSV row missing columns: {missing}")
    return [float(row[name]) for name in columns]


def build_check_module_samples(
    *,
    sampler: Any,
    config: Mapping[str, Any],
    rows: Sequence[Mapping[str, Any]],
    csv_fieldnames: Sequence[str] | None = None,
) -> list[Sample]:
    import numpy as np

    columns = coordinate_columns_for_check_modules(config, csv_fieldnames)
    samples: list[Sample] = []
    for row in rows:
        coords = build_check_module_u_coords(row, columns)
        sample = sampler._build_sample(np.array(coords, dtype=np.float64))
        if "uuid" in row:
            sample.uuid = str(row["uuid"])
        samples.append(sample)
    return samples


__all__ = [
    "build_check_module_samples",
    "build_check_module_u_coords",
    "coordinate_columns_for_check_modules",
    "normalize_check_module_records",
    "sample_tree_file_sets",
    "verify_check_modules_golden",
]
