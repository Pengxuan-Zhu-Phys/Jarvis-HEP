#!/usr/bin/env python3
from __future__ import annotations

import concurrent.futures
import csv
import os
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.csv_sampler import CSVSampler  # noqa: E402
from jarvishep.distributor import Distributor  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _ImmediateFactory:
    def __init__(self):
        self.submitted = []

    def submit_task(self, sample_info):
        self.submitted.append(dict(sample_info))
        fut = concurrent.futures.Future()
        fut.set_result(0.0)
        return fut


def _build_base_config(csv_path: str, with_uuid: bool = True, variables=None) -> dict:
    sampling_cfg = {
        "Method": "CSV",
        "CSV": {
            "path": csv_path,
        },
        "LogLikelihood": [
            {"name": "LogL", "expression": "0.0"},
        ],
    }
    if with_uuid:
        sampling_cfg["CSV"]["uuid_column"] = "uuid"
    if variables is not None:
        sampling_cfg["CSV"]["variables"] = list(variables)
    return {
        "Scan": {"name": "csv-test", "save_dir": os.path.dirname(csv_path)},
        "Sampling": sampling_cfg,
        "EnvReqs": {"OS": [], "Python": {"version": ">=3.10", "Dependencies": []}},
        "Calculators": {"make_paraller": 1, "Modules": []},
    }


class TestCSVSampler(unittest.TestCase):
    def test_distributor_registers_csv_sampler(self):
        sampler = Distributor.set_method("CSV")
        self.assertIsInstance(sampler, CSVSampler)

    def test_csv_uuid_column_is_used_and_missing_uuid_is_generated(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-sampler-") as tmp:
            csv_path = os.path.join(tmp, "points.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f1:
                writer = csv.DictWriter(f1, fieldnames=["uuid", "x", "y"])
                writer.writeheader()
                writer.writerow({"uuid": "csv-row-001", "x": "1.0", "y": "2"})
                writer.writerow({"uuid": "", "x": "3.5", "y": "4"})

            sampler = CSVSampler()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_logger(_NoopLogger())
            sampler.set_config(_build_base_config(csv_path, with_uuid=True))
            factory = _ImmediateFactory()
            sampler.set_factory(factory)
            sampler.run_nested()

            self.assertEqual(len(factory.submitted), 2)
            self.assertEqual(factory.submitted[0]["uuid"], "csv-row-001")
            self.assertTrue(factory.submitted[1]["uuid"])
            self.assertNotEqual(factory.submitted[1]["uuid"], "csv-row-001")
            self.assertEqual(factory.submitted[0]["params"]["x"], 1.0)
            self.assertEqual(factory.submitted[0]["params"]["y"], 2)

            map_path = os.path.join(tmp, "DATABASE", "csv_uuid_map.csv")
            self.assertTrue(os.path.exists(map_path))
            with open(map_path, "r", encoding="utf-8", newline="") as f2:
                rows = list(csv.DictReader(f2))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["source_uuid"], "csv-row-001")
            self.assertEqual(rows[0]["effective_uuid"], "csv-row-001")
            self.assertEqual(rows[1]["source_uuid"], "")
            self.assertTrue(rows[1]["effective_uuid"])

    def test_csv_without_uuid_column_generates_uuid_map(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-sampler-") as tmp:
            csv_path = os.path.join(tmp, "points_no_uuid.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f1:
                writer = csv.DictWriter(f1, fieldnames=["x", "flag"])
                writer.writeheader()
                writer.writerow({"x": "10", "flag": "true"})
                writer.writerow({"x": "20.5", "flag": "false"})

            sampler = CSVSampler()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_logger(_NoopLogger())
            sampler.set_config(_build_base_config(csv_path, with_uuid=False))
            factory = _ImmediateFactory()
            sampler.set_factory(factory)
            sampler.run_nested()

            self.assertEqual(len(factory.submitted), 2)
            self.assertEqual(factory.submitted[0]["params"]["x"], 10)
            self.assertEqual(factory.submitted[0]["params"]["flag"], True)
            self.assertEqual(factory.submitted[1]["params"]["x"], 20.5)
            self.assertEqual(factory.submitted[1]["params"]["flag"], False)
            self.assertNotEqual(factory.submitted[0]["uuid"], factory.submitted[1]["uuid"])

            map_path = os.path.join(tmp, "DATABASE", "csv_uuid_map.csv")
            with open(map_path, "r", encoding="utf-8", newline="") as f2:
                rows = list(csv.DictReader(f2))
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["source_uuid"], "")
            self.assertEqual(rows[1]["source_uuid"], "")

    def test_csv_variables_filter_only_selected_inputs(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-sampler-") as tmp:
            csv_path = os.path.join(tmp, "points_subset.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f1:
                writer = csv.DictWriter(f1, fieldnames=["uuid", "x", "y", "z"])
                writer.writeheader()
                writer.writerow({"uuid": "u1", "x": "1", "y": "2", "z": "3"})

            sampler = CSVSampler()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_logger(_NoopLogger())
            sampler.set_config(_build_base_config(csv_path, with_uuid=True, variables=["x", "z"]))
            factory = _ImmediateFactory()
            sampler.set_factory(factory)
            sampler.run_nested()

            self.assertEqual(len(factory.submitted), 1)
            params = factory.submitted[0]["params"]
            self.assertIn("x", params)
            self.assertIn("z", params)
            self.assertNotIn("y", params)
            self.assertEqual(params["x"], 1)
            self.assertEqual(params["z"], 3)

    def test_csv_variables_filter_raises_on_missing_column(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-sampler-") as tmp:
            csv_path = os.path.join(tmp, "points_missing_col.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f1:
                writer = csv.DictWriter(f1, fieldnames=["x", "y"])
                writer.writeheader()
                writer.writerow({"x": "1", "y": "2"})

            sampler = CSVSampler()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_logger(_NoopLogger())
            sampler.set_config(_build_base_config(csv_path, with_uuid=False, variables=["x", "w"]))
            sampler.set_factory(_ImmediateFactory())
            with self.assertRaises(ValueError):
                sampler.run_nested()

    def test_csv_sampler_supports_initialize_contract(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-csv-sampler-") as tmp:
            csv_path = os.path.join(tmp, "points_init.csv")
            with open(csv_path, "w", encoding="utf-8", newline="") as f1:
                writer = csv.DictWriter(f1, fieldnames=["uuid", "x"])
                writer.writeheader()
                writer.writerow({"uuid": "u1", "x": "1"})

            sampler = CSVSampler()
            sampler.info["sample"] = {
                "task_result_dir": tmp,
                "sample_dirs": os.path.join(tmp, "SAMPLE"),
                "archive_samples": False,
            }
            sampler.set_logger(_NoopLogger())
            sampler.set_config(_build_base_config(csv_path, with_uuid=True))
            sampler.initialize()

            self.assertIsNone(sampler._records_iter)


if __name__ == "__main__":
    unittest.main()
