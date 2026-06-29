#!/usr/bin/env python3
"""Jarvis2Core distributed run integration gates (check-modules e2e)."""

from __future__ import annotations

import json
import os
import tempfile
import unittest
from unittest import mock

from fakeredis import TcpFakeServer

from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.factory import TaskFactory
from test_worker_calculator import (
    _normalize_database_records,
    _sample_tree_file_sets,
    _start_tcp_fakeredis,
)


TESTS_ROOT = os.path.dirname(__file__)
PARITY_PROJECT = os.path.join(TESTS_ROOT, "parity_project")
FIXTURES = os.path.join(TESTS_ROOT, "fixtures", "parity_m1")
CHECK_MODULES_YAML = os.path.join(PARITY_PROJECT, "check_modules.yaml")


def _stop_factory_workers() -> None:
    factory = TaskFactory._instance
    if factory is not None:
        try:
            factory.shutdown(wait=True)
        except Exception:
            pass
    TaskFactory.reset_instance()


class CoreRunDistributedTests(unittest.TestCase):
    def setUp(self) -> None:
        _stop_factory_workers()

    def tearDown(self) -> None:
        _stop_factory_workers()

    def test_check_modules_yaml_end_to_end_golden_parity(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(
                    os.path.join(FIXTURES, "expected_calculator_records.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_records = json.load(handle)
                with open(
                    os.path.join(FIXTURES, "expected_sample_files.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_files = json.load(handle)

                core = Jarvis2Core()
                core.load_task_yaml(CHECK_MODULES_YAML)
                core.config["task_result_dir"] = tmpdir
                core.config["Runtime"]["redis"] = redis_config
                core.runtime = core.config["Runtime"]
                core._populate_info_from_config()

                golden = {
                    "records": expected_records,
                    "sample_files": expected_files,
                }
                count = core.run(check_modules=True, verify_golden=golden, write_run_summary=False)
                self.assertEqual(count, 10)

                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))

                tree = _sample_tree_file_sets(os.path.join(tmpdir, "SAMPLE"))
                self.assertEqual(len(tree), 10)
                for files in tree:
                    self.assertEqual(files, expected_files)
        finally:
            server.shutdown()
            server.server_close()

    def test_check_modules_via_cli_dispatch_mocked_runtime(self) -> None:
        from jarvishep2.client import main

        with mock.patch("jarvishep2.client.Jarvis2Core") as core_cls:
            core = core_cls.return_value
            core.check_modules.return_value = 10
            code = main([CHECK_MODULES_YAML, "--check-modules"])
        self.assertEqual(code, 0)
        core.load_task_yaml.assert_called_once_with(CHECK_MODULES_YAML)
        core.check_modules.assert_called_once_with()


if __name__ == "__main__":
    unittest.main()