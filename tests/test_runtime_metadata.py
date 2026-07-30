from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock

from jarvishep2.process_cleanup import JarvisProcess, JarvisScan, runtime_metadata_for_scan
from jarvishep2.runtime_metadata import read_scan_metadata, write_scan_metadata


class RuntimeMetadataTests(unittest.TestCase):
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
