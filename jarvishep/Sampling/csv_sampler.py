#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import csv
import os
import re
from typing import Any, Dict, Iterator, Sequence
from uuid import uuid4

from jarvishep.log_kv import format_two_column_log
from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.sample import Sample


_INT_PATTERN = re.compile(r"^[+-]?\d+$")


class CSVSampler(SamplingVirtial):
    """Sampler that replays sample points from a CSV table."""

    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "CSV"
        self.tasks = set()
        self.future_to_sample = {}
        self.config = {}
        self.info = {}
        self._selectionexp = None
        self._csv_path = None
        self._csv_delimiter = ","
        self._csv_encoding = "utf-8"
        self._uuid_column_requested = "uuid"
        self._uuid_column_resolved = None
        self._selected_variables = None
        self._records_iter = None

    def load_schema_file(self):
        self.schema = self.path["CSVSchema"]

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self) -> None:
        sampling_cfg = self.config.get("Sampling", {}) if isinstance(self.config, dict) else {}
        csv_cfg = sampling_cfg.get("CSV", {})
        if not isinstance(csv_cfg, dict):
            raise ValueError("Sampling.CSV must be an object.")

        raw_path = csv_cfg.get("path")
        if not isinstance(raw_path, str) or not raw_path.strip():
            raise ValueError("Sampling.CSV.path is required.")
        self._csv_path = self.decode_path(raw_path.strip())

        if not os.path.exists(self._csv_path):
            raise FileNotFoundError(f"CSV source not found: {self._csv_path}")

        self._uuid_column_requested = str(csv_cfg.get("uuid_column", "uuid")).strip() or "uuid"
        self._selected_variables = self._normalize_selected_variables(csv_cfg.get("variables"))
        self._selectionexp = sampling_cfg.get("selection")

        # CSV sampler consumes direct values from table rows; no distribution mapping required.
        self.vars = tuple()
        self._records_iter = None
        self._uuid_column_resolved = self._resolve_uuid_column_from_file()

    def initialize(self):
        self.logger.warning("Initializing the CSV sampler")
        self._records_iter = None
        if not self._selectionexp:
            return
        try:
            records = self._iter_records()
            first = next(records, None)
            if first is not None:
                self.evaluate_selection(self._selectionexp, first["params"])
        except Exception as exc:
            self.logger.error(f"CSV Sampler meets error when checking selection expression: {exc}")
            raise

    def _resolve_uuid_column_from_file(self) -> str | None:
        assert self._csv_path is not None
        with open(self._csv_path, "r", encoding=self._csv_encoding, newline="") as f1:
            reader = csv.DictReader(f1, delimiter=self._csv_delimiter)
            fieldnames = reader.fieldnames or []
        if not fieldnames:
            raise ValueError(f"CSV file has no header row: {self._csv_path}")
        return self._resolve_uuid_column(fieldnames)

    def _resolve_uuid_column(self, fieldnames) -> str | None:
        req = str(self._uuid_column_requested).strip().lower()
        if not req:
            return None
        normalized = [str(name).strip() for name in fieldnames if name is not None]
        for name in normalized:
            if name.lower() == req:
                return name
        return None

    @staticmethod
    def _normalize_selected_variables(raw_variables) -> list[str] | None:
        if raw_variables is None:
            return None
        if not isinstance(raw_variables, Sequence) or isinstance(raw_variables, (str, bytes)):
            raise ValueError("Sampling.CSV.variables must be a list of column names.")
        cols = []
        seen = set()
        for item in raw_variables:
            if not isinstance(item, str):
                raise ValueError("Sampling.CSV.variables must contain only string names.")
            name = item.strip()
            if not name:
                continue
            if name in seen:
                continue
            seen.add(name)
            cols.append(name)
        return cols if cols else []

    def __iter__(self):
        self._records_iter = self._iter_records()
        return self

    def __next__(self):
        if self._records_iter is None:
            self._records_iter = self._iter_records()
        record = next(self._records_iter)
        return record["params"]

    def next_sample(self):
        return self.__next__()

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
            except Exception:
                pass
        try:
            return float(text)
        except Exception:
            return text

    def _iter_records(self) -> Iterator[Dict[str, Any]]:
        assert self._csv_path is not None
        seen_source_uuid = set()
        with open(self._csv_path, "r", encoding=self._csv_encoding, newline="") as f1:
            reader = csv.reader(f1, delimiter=self._csv_delimiter)
            header_row = next(reader, None)
            fieldnames = [str(name).strip() for name in (header_row or [])]
            if not fieldnames:
                raise ValueError(f"CSV file has no header row: {self._csv_path}")

            name_to_idx = {}
            for idx, name in enumerate(fieldnames):
                if name and name not in name_to_idx:
                    name_to_idx[name] = idx
            uuid_col = self._resolve_uuid_column(fieldnames)
            uuid_idx = name_to_idx.get(uuid_col) if uuid_col else None

            if self._selected_variables is None:
                selected_variables = [name for name in fieldnames if name != uuid_col]
            else:
                selected_variables = list(self._selected_variables)

            selected_idx = []
            for col in selected_variables:
                if col == uuid_col:
                    continue
                if col not in name_to_idx:
                    raise ValueError(f"CSV required column not found: {col}")
                selected_idx.append((col, int(name_to_idx[col])))

            for row_index, row in enumerate(reader, start=1):
                params: Dict[str, Any] = {}
                for key, col_idx in selected_idx:
                    raw_value = row[col_idx] if col_idx < len(row) else ""
                    params[key] = self._coerce_cell(raw_value)

                if self._selectionexp:
                    if not self.evaluate_selection(self._selectionexp, params):
                        continue

                source_uuid = None
                if uuid_idx is not None and uuid_idx < len(row):
                    raw_uuid = row[uuid_idx]
                    source_uuid = str(raw_uuid).strip() or None
                    if source_uuid is not None:
                        if source_uuid in seen_source_uuid:
                            raise ValueError(
                                f"Duplicate CSV uuid found at row {row_index}: {source_uuid}"
                            )
                        seen_source_uuid.add(source_uuid)

                effective_uuid = source_uuid or str(uuid4())
                yield {
                    "row_index": row_index,
                    "source_uuid": source_uuid,
                    "uuid": effective_uuid,
                    "params": params,
                }

    def _resolve_uuid_map_path(self) -> str | None:
        sample_cfg = self.info.get("sample", {}) if isinstance(self.info, dict) else {}
        task_result_dir = sample_cfg.get("task_result_dir")
        if not isinstance(task_result_dir, str) or not task_result_dir:
            return None
        return os.path.join(task_result_dir, "DATABASE", "csv_uuid_map.csv")

    def run_nested(self):
        if not self._csv_path:
            raise RuntimeError("CSV sampler is not configured.")
        if not hasattr(self, "factory") or self.factory is None:
            raise RuntimeError("CSV sampler has no WorkerFactory.")

        worker_slots = max(1, int(getattr(self, "max_workers", os.cpu_count() or 1)))
        self.tasks = set()
        self.future_to_sample = {}
        exhausted = False
        records = self._iter_records()
        base_sample_cfg = self.info["sample"]
        submitted = 0
        completed = 0

        map_path = self._resolve_uuid_map_path()
        map_fp = None
        map_writer = None
        if map_path is not None:
            os.makedirs(os.path.dirname(map_path), exist_ok=True)
            map_fp = open(map_path, "w", encoding="utf-8", newline="")
            map_writer = csv.DictWriter(
                map_fp,
                fieldnames=["row_index", "source_uuid", "effective_uuid"],
            )
            map_writer.writeheader()

        try:
            while (not exhausted) or self.tasks:
                while not exhausted and len(self.tasks) < worker_slots:
                    try:
                        record = next(records)
                    except StopIteration:
                        exhausted = True
                        break

                    sample = Sample(record["params"])
                    sample.update_uuid(record["uuid"])
                    save_dir = self._next_bucket_dir_for_sample()
                    sample_cfg = self.build_sample_config(base_sample_cfg, save_dir=save_dir)
                    sample.set_config(sample_cfg)

                    if map_writer is not None:
                        map_writer.writerow(
                            {
                                "row_index": int(record["row_index"]),
                                "source_uuid": record["source_uuid"] or "",
                                "effective_uuid": sample.uuid,
                            }
                        )
                    submitted += 1

                    try:
                        future = self.factory.submit_task(sample.info)
                    except Exception:
                        self._on_sample_completed(sample.info)
                        sample.close()
                        raise

                    self.tasks.add(future)
                    self.future_to_sample[future] = sample

                if not self.tasks:
                    continue

                done, _ = concurrent.futures.wait(
                    self.tasks,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                self.tasks.difference_update(done)

                for future in done:
                    sample = self.future_to_sample.pop(future, None)
                    try:
                        future.result()
                    except Exception as exc:
                        suuid = sample.uuid if sample else "UNKNOWN"
                        self.logger.error(
                            format_two_column_log(
                                "[WorkerFactory] future exception consumed",
                                [("uuid", suuid), ("error", exc)],
                            )
                        )
                        raise
                    finally:
                        if sample is not None:
                            self._on_sample_completed(sample.info)
                            sample.close()
                            completed += 1
        finally:
            if map_fp is not None:
                map_fp.close()

        self.logger.warning(
            format_two_column_log(
                "CSV sampler finished",
                [
                    ("submitted", submitted),
                    ("completed", completed),
                    ("map", map_path or "DISABLED"),
                ],
            )
        )

    def finalize(self):
        pass

    def combine_data(self, df_full) -> None:
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.set_bucket_alloc()
        self.logger.warning("WorkerFactory is ready for CSV sampler")
