from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.process_cleanup import JarvisProcess, JarvisScan, runtime_metadata_for_scan
from jarvishep2.runtime_metadata import (
    calculator_pools_from_metadata,
    read_scan_metadata,
    sampler_metadata_from_config,
    write_scan_metadata,
)


class RuntimeMetadataTests(unittest.TestCase):
    def test_sampler_metadata_is_allowlisted_and_hides_csv_path(self) -> None:
        metadata = sampler_metadata_from_config(
            {
                "Sampling": {
                    "Method": "CSV",
                    "Bounds": {
                        "path": "/secret/project/points.csv",
                        "uuid_column": "uuid",
                        "delimiter": ",",
                    },
                }
            }
        )
        self.assertEqual(metadata["method"], "CSV")
        self.assertEqual(metadata["family"], "simple")
        self.assertEqual(
            metadata["config"],
            {"basename": "points.csv", "uuid_column": "uuid"},
        )
        self.assertNotIn("/secret", str(metadata))

    def test_grid_sampler_metadata_derives_shape_and_total(self) -> None:
        metadata = sampler_metadata_from_config(
            {
                "Sampling": {
                    "Method": "Grid",
                    "Variables": [
                        {"distribution": {"parameters": {"num": 10}}},
                        {"distribution": {"parameters": {"num": 30}}},
                    ],
                }
            }
        )
        self.assertEqual(metadata["dimensions"], 2)
        self.assertEqual(metadata["config"]["shape"], [10, 30])
        self.assertEqual(metadata["config"]["total"], 300)

    def test_metadata_file_binds_scan_and_redis_endpoint(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            redis = {"host": "127.0.0.1", "port": 6381, "db": 2}
            path = write_scan_metadata(
                config={"scan_name": "alpha", "task_yaml": "/tmp/task.yaml"},
                info={"scan_name": "alpha", "task_result_dir": root},
                redis=redis,
            )
            payload = read_scan_metadata(path, redis=redis, expected_scan="alpha")

        self.assertIsNotNone(payload)
        assert payload is not None
        self.assertEqual(payload["redis"]["port"], 6381)
        self.assertFalse(os.path.exists(path))

    def test_discovery_reads_advertised_redis_metadata(self) -> None:
        scan = JarvisScan(
            reference="R1",
            name="alpha",
            processes=(JarvisProcess(10, "Jarvis-Redis:alpha@6381/2"),),
        )
        client = mock.Mock()
        client.get_runtime_metadata_path.return_value = "/tmp/runtime.json"
        with mock.patch("jarvishep2.process_cleanup.RedisQueue", return_value=client), mock.patch(
            "jarvishep2.process_cleanup.read_scan_metadata",
            return_value={"scan_name": "alpha"},
        ) as read:
            payload = runtime_metadata_for_scan(scan)

        self.assertEqual(payload, {"scan_name": "alpha"})
        read.assert_called_once_with(
            "/tmp/runtime.json",
            redis={"host": "127.0.0.1", "port": 6381, "db": 2},
            expected_scan="alpha",
        )

    def test_metadata_records_exclusive_and_shared_calculator_pools(self) -> None:
        with tempfile.TemporaryDirectory() as root:
            redis = {"host": "127.0.0.1", "port": 6381, "db": 2}
            path = write_scan_metadata(
                config={
                    "scan_name": "alpha",
                    "Calculators": {
                        "Pools": {"Slow": 3, "Prep": 2},
                        "Modules": [
                            {"name": "Slow"},
                            {
                                "name": "Prep",
                                "modes": [{"name": "fast"}, {"name": "full"}],
                            },
                        ],
                    },
                },
                info={"scan_name": "alpha", "task_result_dir": root},
                redis=redis,
            )
            payload = read_scan_metadata(path, redis=redis, expected_scan="alpha")

        assert payload is not None
        self.assertEqual(payload["calculator_pools"], {"Slow": 3})
        self.assertEqual(
            payload["calculator_shared"],
            {"Prep": {"modes": ["fast", "full"], "n": 2}},
        )
        self.assertEqual(
            calculator_pools_from_metadata(payload),
            (["Slow", "Prep"], {"Slow": 3, "Prep": 2}),
        )
