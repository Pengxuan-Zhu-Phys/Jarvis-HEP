from __future__ import annotations

import csv
import json
import os
import threading
import time
import traceback
from queue import Empty, Full, Queue

import h5py
import numpy as np
from loguru import logger

from jarvishep.log_kv import format_two_column_log
from jarvishep.observable_io import (
    flatten_records_for_csv,
    load_schema,
    make_json_compatible,
    resolve_schema_path,
    save_schema,
    update_schema_with_record,
)


class GlobalHDF5Writer:
    def __init__(
        self,
        dbinfo,
        write_interval=15,
        queue_size=20000,
        batch_size=512,
        enqueue_timeout=1.0,
        max_size_mb=100,
    ):
        super().__init__()
        self.logger = logger.bind(
            module="Jarvis-HEP.hdf5-Writter",
            to_console=True,
            Jarvis=True,
            _log_domain="jarvis_hep",
        )
        self._info_lock = threading.Lock()
        self._schema_lock = threading.Lock()
        self._metric_lock = threading.Lock()
        self.initialize(dbinfo)

        self.write_interval = max(0.05, float(write_interval))
        self.batch_size = max(1, int(batch_size))
        self.enqueue_timeout = max(0.05, float(enqueue_timeout))
        self.max_size = max(1, int(max_size_mb)) * 1024 * 1024

        self.data_queue = Queue(maxsize=max(1, int(queue_size)))
        self.shutdown_event = threading.Event()
        self._stop_token = object()
        self._backpressure_log_ts = 0.0
        self._writer_error = None
        self._writer_trace = ""

        self.enqueued_count = 0
        self.flushed_count = 0
        self.max_queue_depth = 0
        self.flush_count = 0

        self.writer_thread = threading.Thread(
            target=self._writer_loop, name="GlobalHDF5Writer", daemon=False
        )

    def initialize(self, dbinfo):
        pthroot, pthext = os.path.splitext(dbinfo["path"])
        self.schema_path = resolve_schema_path(dbinfo["path"])
        self.schema = load_schema(self.schema_path, pthroot)
        self._log_schema_warnings(self.schema, context="initialize")
        self._save_schema()
        self.infos = {
            "info": dbinfo["info"],
            "schema": self.schema_path,
            "activeNO": 0,
            "totalDB": 1,
            "converted": [],
            "pending_converted": [],
            "errors": [],
            "paths": [],
            "pathroot": pthroot,
            "pathext": pthext,
            "active path": "".join([pthroot, ".0", pthext]),
        }
        self._save_infos()

    def _atomic_write_json(self, path: str, data: dict) -> None:
        tmp_path = path + ".tmp"
        with open(tmp_path, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)
            f.flush()
        os.replace(tmp_path, path)

    def _save_infos(self) -> None:
        with self._info_lock:
            self._atomic_write_json(self.infos["info"], self.infos)

    def _save_schema(self) -> None:
        with self._schema_lock:
            save_schema(self.schema_path, self.schema)

    def _log_schema_warnings(self, schema: dict, context: str) -> None:
        warn_msgs = schema.get("_warnings", [])
        if not warn_msgs:
            return
        for msg in warn_msgs:
            self.logger.warning(f"[schema:{context}] {msg}")

    def _update_schema(self, record) -> None:
        if not isinstance(record, dict):
            return
        with self._schema_lock:
            changed = update_schema_with_record(self.schema, record)
            if changed:
                save_schema(self.schema_path, self.schema)

    def _raise_if_writer_failed(self):
        if self._writer_error is None:
            return
        raise RuntimeError("Global HDF5 writer thread failed") from self._writer_error

    def _set_writer_failed(self, exc):
        self._writer_error = exc
        self._writer_trace = traceback.format_exc()
        self.logger.error(f"Global HDF5 writer thread failed -> {exc}")
        if self._writer_trace:
            self.logger.error(self._writer_trace)

    def _record_enqueue(self):
        with self._metric_lock:
            self.enqueued_count += 1
            depth = self.data_queue.qsize()
            if depth > self.max_queue_depth:
                self.max_queue_depth = depth

    def _record_flush(self, nrows):
        with self._metric_lock:
            self.flushed_count += int(nrows)
            self.flush_count += 1

    def start(self):
        """Start the writer thread."""
        if not self.writer_thread.is_alive():
            self.writer_thread.start()

    def add_data(self, data):
        """Add one record to the write queue with bounded backpressure."""
        self._raise_if_writer_failed()
        self._update_schema(data)
        payload = convert(data)

        while not self.shutdown_event.is_set():
            self._raise_if_writer_failed()
            try:
                self.data_queue.put(payload, timeout=self.enqueue_timeout)
                self._record_enqueue()
                return
            except Full:
                now = time.monotonic()
                if now - self._backpressure_log_ts >= 5.0:
                    self.logger.warning(
                        format_two_column_log(
                            "Global writer queue full; backpressure active",
                            [
                                ("queue_depth", self.data_queue.qsize()),
                                ("queue_capacity", self.data_queue.maxsize),
                            ],
                        )
                    )
                    self._backpressure_log_ts = now

        raise RuntimeError("Global HDF5 writer has been stopped; rejecting new data.")

    def _writer_loop(self):
        """Event-driven writer loop: flush on batch/fullness or timeout."""
        pending = []
        last_flush = time.monotonic()

        while True:
            timeout = max(0.0, self.write_interval - (time.monotonic() - last_flush))
            try:
                item = self.data_queue.get(timeout=timeout)
                if item is self._stop_token:
                    pending.extend(self._take_pending_from_queue())
                    if pending:
                        try:
                            self._flush_batch_to_hdf5(pending)
                        except Exception as exc:
                            self._set_writer_failed(exc)
                    break

                pending.append(item)
                if len(pending) >= self.batch_size:
                    try:
                        self._flush_batch_to_hdf5(pending)
                    except Exception as exc:
                        self._set_writer_failed(exc)
                        break
                    pending = []
                    last_flush = time.monotonic()

            except Empty:
                if pending:
                    try:
                        self._flush_batch_to_hdf5(pending)
                    except Exception as exc:
                        self._set_writer_failed(exc)
                        break
                    pending = []
                    last_flush = time.monotonic()
                elif self.shutdown_event.is_set():
                    # Keep waiting for stop token while allowing fast wake-up.
                    continue
            except Exception as exc:
                self._set_writer_failed(exc)
                break

    def _take_pending_from_queue(self):
        pending = []
        while True:
            try:
                item = self.data_queue.get_nowait()
            except Empty:
                break
            if item is self._stop_token:
                continue
            pending.append(item)
        return pending

    def _rotate_if_needed(self):
        active_path = self.infos["active path"]
        current_file_size = os.path.getsize(active_path) if os.path.exists(active_path) else 0
        if current_file_size < self.max_size:
            return

        try:
            self.hdf5_to_csv()
        except Exception as exc:
            self.logger.warning(f"hdf5_to_csv during rotation failed -> {exc}")

        self.infos["paths"].append(active_path)
        self.infos["activeNO"] += 1
        self.infos["totalDB"] += 1
        self.infos["active path"] = "{}.{}{}".format(
            self.infos["pathroot"],
            self.infos["activeNO"],
            self.infos["pathext"],
        )
        self.logger.warning(f"Creating new HDF5 file due to size limit -> {self.infos['active path']}")
        self._save_infos()

    def _flush_batch_to_hdf5(self, batch):
        if not batch:
            return

        self._rotate_if_needed()
        payload = [json.dumps(item, ensure_ascii=False) for item in batch]

        with h5py.File(self.infos["active path"], "a") as h5f:
            dtype = h5py.string_dtype(encoding="utf-8")
            if "records" not in h5f:
                h5f.create_dataset(
                    "records",
                    shape=(0,),
                    maxshape=(None,),
                    chunks=True,
                    dtype=dtype,
                )
            records = h5f["records"]
            start = int(records.shape[0])
            end = start + len(payload)
            records.resize((end,))
            records[start:end] = payload

        self._record_flush(len(payload))
        self.logger.info(
            "Global writer flushed {} rows -> {} (queue depth {})".format(
                len(payload),
                self.infos["active path"],
                self.data_queue.qsize(),
            )
        )

    def _write_data_to_hdf5(self):
        """Immediate drain helper used by tests and shutdown fallback."""
        pending = self._take_pending_from_queue()
        if pending:
            self._flush_batch_to_hdf5(pending)

    def stop(self):
        """Signal writer to stop, drain queue, and finalize CSV conversion."""
        self.shutdown_event.set()

        if self.writer_thread.is_alive():
            stop_enqueued = False
            while True:
                try:
                    self.data_queue.put(self._stop_token, timeout=self.enqueue_timeout)
                    stop_enqueued = True
                    break
                except Full:
                    # Keep single-writer semantics: do not flush from main thread
                    # while writer thread is active.
                    if not self.writer_thread.is_alive():
                        break
                    continue
            if not stop_enqueued and not self.writer_thread.is_alive():
                # Writer exited before stop token could be queued. Drain
                # synchronously as a fallback (single-writer remains true).
                self._write_data_to_hdf5()
            self.writer_thread.join()
        else:
            # If writer was never started, still flush synchronously.
            self._write_data_to_hdf5()

        self.logger.warning("Global HDF5 writer stopped")
        with self._metric_lock:
            self.logger.warning(
                format_two_column_log(
                    "Global writer summary",
                    [
                        ("enqueued", self.enqueued_count),
                        ("flushed", self.flushed_count),
                        ("flush_count", self.flush_count),
                        ("max_queue_depth", self.max_queue_depth),
                    ],
                )
            )
            pending_gap = self.enqueued_count - self.flushed_count
        if pending_gap != 0:
            msg = format_two_column_log(
                "Global writer data mismatch",
                [
                    ("enqueued", self.enqueued_count),
                    ("flushed", self.flushed_count),
                ],
            )
            self.logger.error(msg)
            if self._writer_error is None:
                raise RuntimeError(msg)

        try:
            self.hdf5_to_csv()
        except Exception as exc:
            self.logger.warning(f"Final hdf5_to_csv failed (ignored) -> {exc}")

        self._raise_if_writer_failed()

    @staticmethod
    def _decode_one_record(raw, dataset_name):
        if isinstance(raw, (bytes, bytearray, np.bytes_)):
            text = bytes(raw).decode("utf-8")
        elif isinstance(raw, str):
            text = raw
        else:
            text = str(raw)

        try:
            return json.loads(text)
        except Exception as exc:
            raise ValueError(f"Failed to decode JSON record in dataset '{dataset_name}': {exc}") from exc

    @classmethod
    def _iter_dataset_records(cls, raw, dataset_name):
        if isinstance(raw, np.ndarray):
            if raw.ndim == 0:
                yield cls._decode_one_record(raw.item(), dataset_name)
                return
            for item in raw.reshape(-1):
                yield cls._decode_one_record(item, dataset_name)
            return
        yield cls._decode_one_record(raw, dataset_name)

    def _iter_records_from_hdf5(self):
        with h5py.File(self.infos["active path"], "r") as hdf5_file:
            for dataset_name in hdf5_file:
                raw = hdf5_file[dataset_name][()]
                for record in self._iter_dataset_records(raw, dataset_name):
                    if isinstance(record, dict):
                        yield record

    def hdf5_to_csv(self) -> bool:
        """Convert active HDF5 file to CSV (best-effort)."""
        csv_path = "{}.{}.csv".format(self.infos["pathroot"], self.infos["activeNO"])

        if csv_path in self.infos.get("converted", []):
            return True

        try:
            schema_for_csv = load_schema(self.schema_path, self.infos["pathroot"])
            self._log_schema_warnings(schema_for_csv, context="hdf5_to_csv")
            if schema_for_csv.get("_warnings"):
                save_schema(self.schema_path, schema_for_csv)

            schema_changed = False
            fieldnames = []
            seen = set()
            row_count = 0
            for record in self._iter_records_from_hdf5():
                csv_rows, row_changed = flatten_records_for_csv(
                    records=[record],
                    schema=schema_for_csv,
                    populate_name_map=True,
                )
                if not csv_rows:
                    continue
                row = csv_rows[0]
                row_count += 1
                schema_changed = schema_changed or row_changed
                for key in row.keys():
                    if key not in seen:
                        seen.add(key)
                        fieldnames.append(key)

            if schema_changed:
                save_schema(self.schema_path, schema_for_csv)

            if row_count == 0:
                self.infos.setdefault("converted", []).append(csv_path)
                self._save_infos()
                return True

            if not fieldnames:
                self.logger.warning(
                    f"CSV conversion produced no columns under schema rules; skip writing -> {csv_path}"
                )
                self.infos.setdefault("converted", []).append(csv_path)
                self._save_infos()
                return True

            with open(csv_path, "w", newline="", encoding="utf-8") as csv_file:
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
                writer.writeheader()
                for record in self._iter_records_from_hdf5():
                    csv_rows, _ = flatten_records_for_csv(
                        records=[record],
                        schema=schema_for_csv,
                        populate_name_map=False,
                    )
                    if not csv_rows:
                        continue
                    writer.writerow(csv_rows[0])

            self.logger.warning(f"Converted HDF5 data to CSV at -> {csv_path}.")
            self.infos.setdefault("converted", []).append(csv_path)
            if csv_path in self.infos.get("pending_converted", []):
                self.infos["pending_converted"].remove(csv_path)
            self._save_infos()
            return True

        except (BlockingIOError, OSError) as exc:
            msg = (
                "CSV conversion skipped (file busy/locked) ->\n"
                "\tactive -> {}\n"
                "\terror  -> {}"
            ).format(self.infos["active path"], exc)
            self.logger.warning(msg)
            self.infos.setdefault("pending_converted", [])
            if csv_path not in self.infos["pending_converted"]:
                self.infos["pending_converted"].append(csv_path)
            self.infos.setdefault("errors", []).append(msg)
            self.infos["errors"] = self.infos["errors"][-50:]
            self._save_infos()
            return False

        except Exception as exc:
            msg = (
                "CSV conversion failed (ignored) ->\n"
                "\tactive -> {}\n"
                "\terror  -> {}"
            ).format(self.infos["active path"], exc)
            self.logger.warning(msg)
            self.infos.setdefault("pending_converted", [])
            if csv_path not in self.infos["pending_converted"]:
                self.infos["pending_converted"].append(csv_path)
            self.infos.setdefault("errors", []).append(msg)
            self.infos["errors"] = self.infos["errors"][-50:]
            self._save_infos()
            return False


def convert(data):
    return make_json_compatible(data)
