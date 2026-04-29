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

from jarvishep.hdf5writer import GlobalHDF5Writer  # noqa: E402
from jarvishep.observable_io import load_schema, resolve_schema_path  # noqa: E402
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
            self.assertEqual(schema["columns"]["arr"]["flatten"]["mode"], "split")
            self.assertEqual(schema["columns"]["lst"]["flatten"]["mode"], "split")

    def test_hdf5_to_csv_uses_split_mode_from_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0]), "meta": {"x": 1}})
            writer.add_data({"uuid": "u2", "LogL": 2.0, "arr": np.array([5.0, 6.0]), "meta": {"x": 2}})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema["columns"]["arr"]["flatten"] = {"mode": "split"}
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            ok = writer.hdf5_to_csv()
            self.assertTrue(ok)

            csv_path = os.path.join(tmpdir, "samples.0.csv")
            with open(csv_path, newline="") as f:
                reader = csv.DictReader(f)
                fieldnames = reader.fieldnames or []
                rows = list(reader)

            self.assertIn("arr[0]", fieldnames)
            self.assertIn("arr[1]", fieldnames)
            self.assertEqual(rows[0]["arr[0]"], "1.0")
            self.assertEqual(rows[0]["arr[1]"], "2.0")

            with open(schema_path, "r") as f:
                schema = json.load(f)
            name_map = schema["columns"]["arr"]["flatten"].get("name_map", {})
            self.assertEqual(name_map.get("arr[0]"), "arr[0]")
            self.assertEqual(name_map.get("arr[1]"), "arr[1]")

            # User customizes output names via name_map.
            schema["columns"]["arr"]["flatten"]["name_map"]["arr[0]"] = "arr_left"
            schema["columns"]["arr"]["flatten"]["name_map"]["arr[1]"] = "arr_right"
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
            self.assertEqual(rows_custom[0]["arr_left"], "1.0")
            self.assertEqual(rows_custom[0]["arr_right"], "2.0")

    def test_utils_convert_warns_and_falls_back_on_broken_schema(self):
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
                convert_hdf5_to_csv(snapshot, out_csv)

            self.assertTrue(os.path.exists(out_csv))
            self.assertTrue(
                any("fallback to default schema" in str(item.message) for item in caught),
                msg="Expected fallback warning when schema is broken.",
            )

            # Broken schema should be repaired and rewritten.
            with open(schema_path, "r") as f:
                repaired = json.load(f)
            self.assertIn("columns", repaired)

    def test_invalid_modes_are_normalized_and_persisted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, _, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0])})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema["flatten_defaults"]["array"] = "???"
            schema["columns"]["arr"]["flatten"] = {"mode": "explode"}
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            ok = writer.hdf5_to_csv()
            self.assertTrue(ok)

            with open(schema_path, "r") as f:
                repaired = json.load(f)

            self.assertEqual(repaired["flatten_defaults"]["array"], "split")
            self.assertEqual(repaired["columns"]["arr"]["flatten"]["mode"], "split")

    def test_schema_version_missing_is_migrated_with_backward_compatible_normalization(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            writer, h5_path, schema_path = self._build_writer(tmpdir)
            writer.add_data({"uuid": "u1", "LogL": 1.0, "arr": np.array([1.0, 2.0])})
            writer._write_data_to_hdf5()

            with open(schema_path, "r") as f:
                schema = json.load(f)
            schema.pop("version", None)
            schema["columns"]["arr"]["flatten"] = {"mode": "split"}
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            loaded = load_schema(schema_path, os.path.splitext(h5_path)[0])

            self.assertEqual(loaded["version"], 1)
            self.assertEqual(loaded["columns"]["arr"]["flatten"]["mode"], "split")
            self.assertTrue(
                any("schema.version missing; migrated to 1" in msg for msg in loaded.get("_warnings", []))
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
            schema["columns"]["arr"]["flatten"] = {"mode": "explode"}
            with open(schema_path, "w") as f:
                json.dump(schema, f, indent=2)

            loaded = load_schema(schema_path, os.path.splitext(h5_path)[0])

            self.assertEqual(loaded["version"], 99)
            self.assertEqual(loaded["custom_meta"]["owner"], "user")
            self.assertEqual(loaded["columns"]["arr"]["flatten"]["mode"], "split")
            self.assertTrue(
                any("newer than supported 1" in msg for msg in loaded.get("_warnings", []))
            )

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
