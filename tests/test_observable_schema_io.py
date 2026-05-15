#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import os
import shutil
import sys
import tempfile
import time
import unittest
import warnings

import h5py
import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.core import Core  # noqa: E402
from jarvishep.hdf5writer import GlobalHDF5Writer  # noqa: E402
from jarvishep.observable_io import (  # noqa: E402
    create_default_schema,
    load_schema,
    resolve_schema_path,
    save_schema,
    update_schema_with_record,
    validate_csv_export_schema,
)
from jarvishep.utils import convert_hdf5_to_csv  # noqa: E402


class ObservableSchemaIOTests(unittest.TestCase):
    def _build_writer(self, tmpdir: str) -> tuple[GlobalHDF5Writer, str, str]:
        h5_path = os.path.join(tmpdir, "samples.hdf5")
        info_path = os.path.join(tmpdir, "running.json")
        writer = GlobalHDF5Writer({"path": h5_path, "info": info_path}, write_interval=999)
        return writer, h5_path, resolve_schema_path(h5_path)

    def test_add_data_supports_array_and_list_and_creates_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data(
                {
                    "uuid": "u1",
                    "LogL": 1.25,
                    "arr": np.array([1.0, 2.0], dtype=np.float64),
                    "lst": [np.float32(3.0), np.float64(4.0)],
                }
            )
            writer._write_data_to_hdf5()

            self.assertTrue(os.path.exists(schema_path))
            with open(schema_path, "r") as f:
                schema = json.load(f)

            self.assertIn("arr", schema["columns"])
            self.assertIn("lst", schema["columns"])
            self.assertEqual(schema["columns"]["arr"]["kind"], "ndarray")
            self.assertEqual(schema["columns"]["lst"]["kind"], "list")
            self.assertNotIn("flatten_defaults", schema)
            self.assertNotIn("flatten", schema["columns"]["arr"])
            self.assertNotIn("csv_export", schema["columns"]["arr"])
            self.assertTrue(schema["columns"]["arr"]["keep"])
            self.assertEqual(schema["columns"]["arr"]["name"]["0"]["name"], "arr_0")
            self.assertEqual(schema["columns"]["arr"]["name"]["1"]["name"], "arr_1")
            self.assertTrue(schema["columns"]["lst"]["keep"])
            self.assertEqual(schema["columns"]["lst"]["name"]["0"]["name"], "lst_0")
            self.assertEqual(schema["columns"]["lst"]["name"]["1"]["name"], "lst_1")

    def test_hdf5_to_csv_uses_csv_export_from_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0]), "meta": {"x": 1}})
            writer.add_data({"uuid": "u2", "LogL": 2.0, "arr": np.array([5.0, 6.0]), "meta": {"x": 2}})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema["columns"]["arr"]["name"]["0"]["name"] = "arr_left"
            schema["columns"]["arr"]["name"]["1"]["name"] = "arr_right"
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            ok = writer.hdf5_to_csv()
            self.assertTrue(ok)

            csv_path = os.path.join(tmpdir, "samples.0.csv")
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames or []
                rows = list(reader)

            self.assertIn("arr_left", fieldnames)
            self.assertIn("arr_right", fieldnames)
            self.assertEqual(rows[0]["arr_left"], "1.0")
            self.assertEqual(rows[0]["arr_right"], "2.0")

            # User switches a dropped dict to a one-column JSON export without rescanning.
            schema["columns"]["meta"]["keep"] = True
            schema["columns"]["meta"]["name"] = "meta"
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            snapshot = os.path.join(tmpdir, "samples.0.hdf5.snap")
            shutil.copy2(os.path.join(tmpdir, "samples.0.hdf5"), snapshot)
            csv_custom = os.path.join(tmpdir, "samples.0.custom.csv")
            convert_hdf5_to_csv(snapshot, csv_custom, schema_path=schema_path)

            with open(csv_custom, newline="") as f:
                reader = csv.DictReader(f)
                fieldnames_custom = reader.fieldnames or []
                rows_custom = list(reader)

            self.assertIn("arr_left", fieldnames_custom)
            self.assertIn("arr_right", fieldnames_custom)
            self.assertIn("meta", fieldnames_custom)
            self.assertEqual(rows_custom[0]["arr_left"], "1.0")
            self.assertEqual(rows_custom[0]["arr_right"], "2.0")
            self.assertEqual(rows_custom[0]["meta"], '{"x":1}')

    def test_utils_convert_requires_valid_csv_export_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1, 2, 3])})
            writer._write_data_to_hdf5()

            with open(schema_path, "w") as f:
                f.write("{ invalid json")

            snapshot = os.path.join(tmpdir, "samples.0.hdf5.snap")
            shutil.copy2(os.path.join(tmpdir, "samples.0.hdf5"), snapshot)
            out_csv = os.path.join(tmpdir, "out.csv")

            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                with self.assertRaises(ValueError) as ctx:
                    convert_hdf5_to_csv(snapshot, out_csv)

            self.assertFalse(os.path.exists(out_csv))
            self.assertIn("must exist for CSV export", str(ctx.exception))
            self.assertTrue(
                any("fallback to default schema" in str(item.message) for item in caught),
                msg="Expected fallback warning when schema is broken.",
            )

            # Broken schema should be repaired and rewritten.
            with open(schema_path, "r") as f:
                repaired = json.load(f)
            self.assertIn("columns", repaired)

    def test_invalid_simplified_export_schema_raises_clear_error(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0])})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema["columns"]["arr"]["keep"] = "yes"
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            with self.assertRaises(ValueError):
                writer.hdf5_to_csv()

            self.assertTrue(
                any(".keep must be boolean" in msg for msg in writer.infos.get("errors", []))
            )

    def test_schema_version_missing_is_set_to_current(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, h5_path, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0])})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema.pop("version", None)
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            loaded = load_schema(schema_path, os.path.splitext(h5_path)[0])

            self.assertEqual(loaded["version"], 3)
            self.assertTrue(loaded["columns"]["arr"]["keep"])
            self.assertTrue(
                any("schema.version missing; set to 3" in msg for msg in loaded.get("_warnings", []))
            )

    def test_future_schema_version_is_preserved_with_compatibility_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, h5_path, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0])})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema["version"] = 99
            schema["custom_meta"] = {"owner": "user"}
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            loaded = load_schema(schema_path, os.path.splitext(h5_path)[0])

            self.assertEqual(loaded["version"], 99)
            self.assertEqual(loaded["custom_meta"]["owner"], "user")
            self.assertTrue(loaded["columns"]["arr"]["keep"])
            self.assertTrue(
                any("newer than supported 3" in msg for msg in loaded.get("_warnings", []))
            )

    def test_large_and_nested_dicts_drop_by_default(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data(
                {
                    "uuid": "u1",
                    "rge_g83_gut": {f"k{i}": i for i in range(83)},
                    "vacuum_fields": {"a": {"nested": 1}, "b": 2},
                    "rge_output": "/tmp/rge.out",
                }
            )
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)

            validate_csv_export_schema(schema)
            self.assertFalse(schema["columns"]["rge_g83_gut"]["keep"])
            self.assertEqual(schema["columns"]["rge_g83_gut"]["num_keys"], 83)
            self.assertFalse(schema["columns"]["vacuum_fields"]["keep"])
            self.assertTrue(schema["columns"]["rge_output"]["keep"])
            self.assertEqual(schema["columns"]["rge_output"]["name"], "rge_output")

            self.assertTrue(writer.hdf5_to_csv())
            csv_path = os.path.join(tmpdir, "samples.0.csv")
            with open(csv_path, newline="", encoding="utf-8") as f:
                fieldnames = csv.DictReader(f).fieldnames or []
            self.assertIn("rge_output", fieldnames)
            self.assertNotIn("rge_g83_gut", fieldnames)
            self.assertNotIn("vacuum_fields", fieldnames)

    def test_old_csv_export_blocks_are_simplified_using_current_policy(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pathroot = os.path.join(tmpdir, "samples")
            schema_path = f"{pathroot}.schema.json"
            with open(schema_path, "w", encoding="utf-8") as f1:
                json.dump(
                    {
                        "version": 2,
                        "pathroot": pathroot,
                        "csv_export_policy": {
                            "scalar": "keep",
                            "path_like_string": "drop",
                            "small_list": "split",
                            "small_flat_dict": "split",
                            "large_list": "drop",
                            "large_dict": "drop",
                            "nested_dict": "drop",
                            "mixed": "drop",
                            "unsupported": "drop",
                        },
                        "columns": {
                            "rge_output": {
                                "kind": "str",
                                "dtype": "str",
                                "csv_export": {
                                    "enabled": True,
                                    "action": "keep",
                                    "mode": "path",
                                    "default_columns": ["rge_output"],
                                    "reason": "path_like_string",
                                },
                            }
                        },
                    },
                    f1,
                )

            loaded = load_schema(schema_path, pathroot)

            self.assertFalse(loaded["columns"]["rge_output"]["keep"])
            self.assertNotIn("csv_export", loaded["columns"]["rge_output"])
            self.assertEqual(loaded["columns"]["rge_output"]["name"], "rge_output")

    def test_writer_stores_records_in_single_appendable_dataset(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, h5_path, _ = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0})
            writer.add_data({"uuid": "u2", "LogL": 2.0})
            writer._write_data_to_hdf5()

            with h5py.File(os.path.join(tmpdir, "samples.0.hdf5"), "r") as h5f:
                self.assertIn("records", h5f)
                self.assertEqual(int(h5f["records"].shape[0]), 2)

            self.assertTrue(os.path.exists(h5_path.replace(".hdf5", ".0.hdf5")))

    def test_writer_stop_without_start_still_flushes_and_converts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, _ = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0})
            writer.stop()

            csv_path = os.path.join(tmpdir, "samples.0.csv")
            self.assertTrue(os.path.exists(csv_path))
            with open(csv_path, newline="", encoding="utf-8") as f1:
                rows = list(csv.DictReader(f1))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["uuid"], "u1")

    def test_core_convert_converts_missing_csv_for_all_hdf5_shards_and_skips_existing_csv(self):
        class _Logger:
            def __init__(self):
                self.messages = []

            def warning(self, message):
                self.messages.append(str(message))

            def error(self, message):
                self.messages.append(str(message))

        def _write_hdf5(path, records):
            with h5py.File(path, "w") as h5f:
                dtype = h5py.string_dtype(encoding="utf-8")
                payload = [json.dumps(record, ensure_ascii=False) for record in records]
                h5f.create_dataset("records", data=payload, dtype=dtype)

        with tempfile.TemporaryDirectory() as tmpdir:
            pathroot = os.path.join(tmpdir, "samples")
            schema_path = f"{pathroot}.schema.json"
            record0 = {"uuid": "u0", "LogL": 0.0}
            record1 = {"uuid": "u1", "LogL": 1.0}
            schema = create_default_schema(pathroot)
            update_schema_with_record(schema, record0)
            save_schema(schema_path, schema)

            hdf5_0 = f"{pathroot}.0.hdf5"
            hdf5_1 = f"{pathroot}.1.hdf5"
            csv_0 = f"{pathroot}.0.csv"
            csv_1 = f"{pathroot}.1.csv"
            _write_hdf5(hdf5_0, [record0])
            _write_hdf5(hdf5_1, [record1])
            with open(csv_0, "w", encoding="utf-8") as f1:
                f1.write("existing\n")

            running_json = os.path.join(tmpdir, "running.json")
            with open(running_json, "w", encoding="utf-8") as f1:
                json.dump(
                    {
                        "pathroot": pathroot,
                        "pathext": ".hdf5",
                        "activeNO": 1,
                        "active path": hdf5_1,
                        "paths": [hdf5_0],
                        "converted": [],
                        "pending_converted": [],
                    },
                    f1,
                )

            core = Core.__new__(Core)
            core.info = {"db": {"info": running_json, "path": f"{pathroot}.hdf5"}}
            core.logger = _Logger()
            core.convert()

            with open(csv_0, encoding="utf-8") as f1:
                self.assertEqual(f1.read(), "existing\n")
            with open(csv_1, newline="", encoding="utf-8") as f1:
                rows = list(csv.DictReader(f1))
            self.assertEqual(rows[0]["uuid"], "u1")
            self.assertTrue(any("CSV file already exists" in msg for msg in core.logger.messages))

            with open(running_json, encoding="utf-8") as f1:
                dbinfo = json.load(f1)
            self.assertIn(csv_0, dbinfo["converted"])
            self.assertIn(csv_1, dbinfo["converted"])

    def test_writer_failure_is_exposed_to_producer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            h5_path = os.path.join(tmpdir, "samples.hdf5")
            info_path = os.path.join(tmpdir, "running.json")
            writer = GlobalHDF5Writer({"path": h5_path, "info": info_path}, write_interval=0.05)

            def _boom(_batch):
                raise RuntimeError("forced flush failure")

            writer._flush_batch_to_hdf5 = _boom
            writer.start()
            writer.add_data({"uuid": "u1", "LogL": 1.0})

            deadline = time.time() + 2.0
            while writer._writer_error is None and time.time() < deadline:
                time.sleep(0.02)

            self.assertIsNotNone(writer._writer_error)
            with self.assertRaises(RuntimeError):
                writer.add_data({"uuid": "u2", "LogL": 2.0})
            with self.assertRaises(RuntimeError):
                writer.stop()


if __name__ == "__main__":
    unittest.main()
