#!/usr/bin/env python3
"""CSV replay sampler for Jarvis-HEP V2."""

from __future__ import annotations

import csv
import re
from collections.abc import Iterator, Mapping, Sequence
from copy import deepcopy
from typing import Any
from uuid import uuid4

from jarvishep2.Sampling.checkpointed_sampler import CheckpointedSampler
from jarvishep2.Sampling.sampling_utils import evaluate_selection
from jarvishep2.Sampling.stateless_batch import run_stateless_distributed
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.runtime_config import get_runtime_block
from jarvishep2.sample import Sample
from jarvishep2.task_config import resolve_sampling_path

_INT_PATTERN = re.compile(r"^[+-]?\d+$")


class CSVSampler(CheckpointedSampler):
    method = "CSV"

    def __init__(self) -> None:
        super().__init__()
        self._logger = get_jarvis_logger("sampler.csv")
        self._csv_path: str | None = None
        self._csv_delimiter = ","
        self._csv_encoding = "utf-8"
        self._uuid_column_requested = "uuid"
        self._uuid_column_resolved: str | None = None
        self._selected_variables: list[str] | None = None
        self._selectionexp: str | None = None
        self._batch_size = 16
        self._runtime_csv_cursor = 0
        self._runtime_seen_source_uuid: set[str] = set()
        self._records_exhausted = False
        self._params_by_uuid: dict[str, dict[str, Any]] = {}

    def set_config(self, config_info: Mapping[str, Any]) -> None:
        super().set_config(config_info)
        sampling = dict(self.config.get("Sampling") or {})
        runtime = get_runtime_block(self.config)
        csv_cfg = dict(sampling.get("CSV") or {})
        if not csv_cfg:
            raise ValueError("CSV sampler requires Sampling.CSV configuration")
        raw_path = str(csv_cfg.get("path") or "").strip()
        if not raw_path:
            raise ValueError("Sampling.CSV.path is required")
        self._csv_path = resolve_sampling_path(self.config, raw_path)
        self._uuid_column_requested = str(csv_cfg.get("uuid_column", "uuid")).strip() or "uuid"
        self._selected_variables = self._normalize_selected_variables(csv_cfg.get("variables"))
        self._selectionexp = sampling.get("selection")
        self._csv_delimiter = str(csv_cfg.get("delimiter", ","))
        self._csv_encoding = str(csv_cfg.get("encoding", "utf-8"))
        workers = int(runtime.get("workers", 1) or 1)
        self._batch_size = max(1, int(runtime.get("batch_size", workers) or workers))
        self._uuid_column_resolved = self._resolve_uuid_column_from_file()

    @staticmethod
    def _normalize_selected_variables(raw_variables: Any) -> list[str] | None:
        if raw_variables is None:
            return None
        if not isinstance(raw_variables, Sequence) or isinstance(raw_variables, (str, bytes)):
            raise ValueError("Sampling.CSV.variables must be a list of column names")
        cols: list[str] = []
        seen: set[str] = set()
        for item in raw_variables:
            name = str(item).strip()
            if not name or name in seen:
                continue
            seen.add(name)
            cols.append(name)
        return cols or None

    def _resolve_uuid_column_from_file(self) -> str | None:
        assert self._csv_path is not None
        with open(self._csv_path, "r", encoding=self._csv_encoding, newline="") as handle:
            reader = csv.DictReader(handle, delimiter=self._csv_delimiter)
            fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError(f"CSV file has no header row: {self._csv_path}")
        return self._resolve_uuid_column(fieldnames)

    def _resolve_uuid_column(self, fieldnames: Sequence[str]) -> str | None:
        req = self._uuid_column_requested.lower()
        if not req:
            return None
        for name in fieldnames:
            if str(name).strip().lower() == req:
                return str(name).strip()
        return None

    @staticmethod
    def _coerce_cell(value: Any) -> Any:
        if value is None:
            return None
        if isinstance(value, (bool, int, float)):
            return value
        text = str(value).strip()
        if not text:
            return None
        lower = text.lower()
        if lower == "true":
            return True
        if lower == "false":
            return False
        if _INT_PATTERN.match(text):
            try:
                return int(text)
            except ValueError:
                pass
        try:
            return float(text)
        except ValueError:
            return text

    def _iter_records(self) -> Iterator[dict[str, Any]]:
        assert self._csv_path is not None
        with open(self._csv_path, "r", encoding=self._csv_encoding, newline="") as handle:
            reader = csv.reader(handle, delimiter=self._csv_delimiter)
            header_row = next(reader, None)
            fieldnames = [str(name).strip() for name in (header_row or [])]
            if not fieldnames:
                raise ValueError(f"CSV file has no header row: {self._csv_path}")

            name_to_idx = {name: idx for idx, name in enumerate(fieldnames) if name}
            uuid_col = self._resolve_uuid_column(fieldnames)
            uuid_idx = name_to_idx.get(uuid_col) if uuid_col else None

            if self._selected_variables is None:
                selected_variables = [name for name in fieldnames if name != uuid_col]
            else:
                selected_variables = list(self._selected_variables)

            selected_idx: list[tuple[str, int]] = []
            for col in selected_variables:
                if col == uuid_col:
                    continue
                if col not in name_to_idx:
                    raise ValueError(f"CSV required column not found: {col}")
                selected_idx.append((col, int(name_to_idx[col])))

            for row_index, row in enumerate(reader, start=1):
                if row_index <= int(self._runtime_csv_cursor):
                    continue
                params: dict[str, Any] = {}
                for key, col_idx in selected_idx:
                    raw_value = row[col_idx] if col_idx < len(row) else ""
                    params[key] = self._coerce_cell(raw_value)
                if self._selectionexp and not evaluate_selection(self._selectionexp, params):
                    continue
                source_uuid = None
                if uuid_idx is not None and uuid_idx < len(row):
                    source_uuid = str(row[uuid_idx]).strip() or None
                    if source_uuid is not None:
                        if source_uuid in self._runtime_seen_source_uuid:
                            raise ValueError(
                                f"Duplicate CSV uuid found at row {row_index}: {source_uuid}"
                            )
                        self._runtime_seen_source_uuid.add(source_uuid)
                effective_uuid = source_uuid or str(uuid4())
                yield {
                    "row_index": row_index,
                    "source_uuid": source_uuid,
                    "uuid": effective_uuid,
                    "params": params,
                }

    def propose_next(self) -> Sample | None:
        if self._records_exhausted:
            return None
        if not hasattr(self, "_record_iter"):
            self._record_iter = iter(self._iter_records())
        try:
            record = next(self._record_iter)
        except StopIteration:
            self._records_exhausted = True
            return None
        self._runtime_csv_cursor = int(record["row_index"])
        sample = self._build_sample([])
        sample.uuid = str(record["uuid"])
        sample.opera_params = dict(record["params"])
        self._params_by_uuid[sample.uuid] = dict(record["params"])
        return sample

    def repropose_unfinished(self) -> list[str]:
        if not self._repropose_after_resume:
            return []
        pending = [uuid for uuid in self._submitted_uuids if uuid not in self._completed_uuids]
        requeued: list[str] = []
        for uuid in pending:
            params = self._params_by_uuid.get(uuid)
            if params is None:
                continue
            sample = self._build_sample([])
            sample.uuid = uuid
            sample.opera_params = dict(params)
            self._submit(sample)
            requeued.append(uuid)
        return requeued

    def run_distributed(self) -> int:
        return run_stateless_distributed(self, propose_next=self.propose_next)

    def at_safe_barrier(self) -> bool:
        if not self._records_exhausted:
            return False
        if not self._submitted_uuids:
            return True
        return set(self._submitted_uuids) <= self._completed_uuids

    def export_runtime_state(self) -> dict[str, Any]:
        return {
            "csv_path": self._csv_path,
            "csv_delimiter": self._csv_delimiter,
            "csv_encoding": self._csv_encoding,
            "uuid_column_requested": self._uuid_column_requested,
            "uuid_column_resolved": self._uuid_column_resolved,
            "selected_variables": deepcopy(self._selected_variables),
            "selectionexp": self._selectionexp,
            "cursor_row_index": int(self._runtime_csv_cursor),
            "seen_source_uuid": sorted(self._runtime_seen_source_uuid),
            "params_by_uuid": deepcopy(self._params_by_uuid),
            "submitted_uuids": list(self._submitted_uuids),
            "completed_uuids": sorted(self._completed_uuids),
            "chains": [],
            "ready_queue": [],
            "control_state": self._checkpoint_control_state(),
            "numpy_random_state": None,
        }

    def import_runtime_state(self, state: Mapping[str, Any]) -> None:
        self._csv_path = state.get("csv_path", self._csv_path)
        self._csv_delimiter = state.get("csv_delimiter", self._csv_delimiter)
        self._csv_encoding = state.get("csv_encoding", self._csv_encoding)
        self._uuid_column_requested = state.get(
            "uuid_column_requested",
            self._uuid_column_requested,
        )
        self._uuid_column_resolved = state.get(
            "uuid_column_resolved",
            self._uuid_column_resolved,
        )
        self._selected_variables = deepcopy(state.get("selected_variables", self._selected_variables))
        self._selectionexp = state.get("selectionexp", self._selectionexp)
        self._runtime_csv_cursor = int(state.get("cursor_row_index", self._runtime_csv_cursor) or 0)
        self._runtime_seen_source_uuid = {
            str(item) for item in state.get("seen_source_uuid") or []
        }
        self._params_by_uuid = {
            str(key): dict(value)
            for key, value in dict(state.get("params_by_uuid") or {}).items()
        }
        self._records_exhausted = False
        if hasattr(self, "_record_iter"):
            delattr(self, "_record_iter")
        self._import_checkpoint_control_state(state)


__all__ = ["CSVSampler"]