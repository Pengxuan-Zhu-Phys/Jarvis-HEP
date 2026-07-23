#!/usr/bin/env python3
"""HDF5 → CSV convert helpers and CLI dispatch."""

from __future__ import annotations

import csv
import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.client import (
    build_parser,
    dispatch,
    dispatch_convert,
    normalize_argv,
)
from jarvishep2.database import (
    SimpleHDF5Writer,
    convert_database_dir,
    convert_hdf5_to_csv,
    discover_database_hdf5,
)
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE


def _write_records(path: str, records: list[dict]) -> None:
    writer = SimpleHDF5Writer(path)
    for row in records:
        writer.add_record(row)


class ConvertHdf5HelpersTests(unittest.TestCase):
    def test_convert_hdf5_to_csv_writes_full_columns(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hdf5 = os.path.join(tmp, "samples.hdf5")
            csv_path = os.path.join(tmp, "samples.csv")
            _write_records(
                hdf5,
                [
                    {"uuid": "u0", "LogL": -1.5, "MChi": 100.0, "nested": {"a": 1}},
                    {"uuid": "u1", "LogL": -2.0, "MChi": 200.0, "extra": 3},
                ],
            )
            result = convert_hdf5_to_csv(hdf5, csv_path)
            self.assertEqual(result["status"], "converted")
            self.assertEqual(result["rows"], 2)
            self.assertTrue(os.path.isfile(csv_path))
            with open(csv_path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["uuid"], "u0")
            self.assertEqual(rows[0]["MChi"], "100.0")
            self.assertIn("nested", rows[0])
            self.assertIn("extra", rows[1])

    def test_skip_existing_unless_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hdf5 = os.path.join(tmp, "samples.hdf5")
            csv_path = os.path.join(tmp, "samples.csv")
            _write_records(hdf5, [{"uuid": "u0", "LogL": 1.0}])
            with open(csv_path, "w", encoding="utf-8") as handle:
                handle.write("stale\n")
            skipped = convert_hdf5_to_csv(hdf5, csv_path)
            self.assertEqual(skipped["status"], "skipped_exists")
            with open(csv_path, encoding="utf-8") as handle:
                self.assertEqual(handle.read(), "stale\n")
            forced = convert_hdf5_to_csv(hdf5, csv_path, force=True)
            self.assertEqual(forced["status"], "converted")
            with open(csv_path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["uuid"], "u0")

    def test_discover_and_convert_database_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            db = os.path.join(tmp, "DATABASE")
            os.makedirs(db)
            path = os.path.join(db, "samples.hdf5")
            _write_records(path, [{"uuid": "a", "x": 1.0, "y": 2.0, "LogL": 0.0}])
            found = discover_database_hdf5(db)
            self.assertEqual(found, [os.path.abspath(path)])
            results = convert_database_dir(db)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["status"], "converted")
            self.assertTrue(os.path.isfile(os.path.join(db, "samples.csv")))


class ConvertCliTests(unittest.TestCase):
    def test_normalize_legacy_convert_flag(self) -> None:
        self.assertEqual(
            normalize_argv(["task.yaml", "--convert", "--force"]),
            ["convert", "task.yaml", "--force"],
        )

    def test_parse_convert_subcommand(self) -> None:
        args = build_parser().parse_args(["convert", "task.yaml", "--force"])
        self.assertEqual(args.command, "convert")
        self.assertEqual(args.task_yaml, "task.yaml")
        self.assertTrue(args.force)

    def test_dispatch_convert_routes_to_core(self) -> None:
        args = build_parser().parse_args(
            normalize_argv(["/tmp/task.yaml", "--convert"])
        )
        core = mock.Mock()
        core.convert.return_value = [
            {
                "hdf5": "/tmp/samples.hdf5",
                "csv": "/tmp/samples.csv",
                "status": "converted",
                "rows": 3,
            }
        ]
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch(args)
        self.assertEqual(code, EXIT_OK)
        core.load_task_yaml.assert_called_once_with("/tmp/task.yaml", validate=False)
        core.convert.assert_called_once_with(force=False)

    def test_dispatch_convert_missing_db_fails(self) -> None:
        core = mock.Mock()
        core.convert.return_value = []
        with mock.patch("jarvishep2.client.Jarvis2Core", return_value=core):
            code = dispatch_convert("task.yaml")
        self.assertEqual(code, EXIT_RUN_FAILED)

    def test_dispatch_convert_requires_yaml(self) -> None:
        self.assertEqual(dispatch_convert(""), EXIT_USAGE)


if __name__ == "__main__":
    unittest.main()
