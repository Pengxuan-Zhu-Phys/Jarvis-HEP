#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from datetime import datetime


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.sample import Sample  # noqa: E402
from jarvishep.sample_logger import SampleLogger  # noqa: E402
from jarvishep.Module.likelihood import LogLikelihood  # noqa: E402


class SampleSaveDirFallbackTests(unittest.TestCase):
    def test_set_config_uses_sample_dirs_when_save_dir_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config({"sample_dirs": tmpdir, "task_result_dir": tmpdir})
            try:
                self.assertTrue(sample.info["save_dir"].startswith(tmpdir))
                self.assertTrue(sample.info["run_log"].startswith(tmpdir))
                self.assertTrue(os.path.isdir(sample.info["save_dir"]))
            finally:
                sample.close()

    def test_sample_logger_preserves_legacy_format_and_raw_behavior(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
            logger = SampleLogger.open(
                log_path,
                module="Sample@test-uuid",
                time_provider=lambda: fixed_dt,
            )
            child = logger.bind(module="Sample@test-uuid (Calc-No.7)")

            child.info("hello world")
            child.bind(raw=True).info("\tstdout line\n")
            child.warning("warn")
            logger.close()

            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(),
                    (
                        "\n·•· Sample@test-uuid (Calc-No.7) \n"
                        "\t-> 01-02 03:04:05.678 - [INFO] >>> \n"
                        "hello world"
                        "\tstdout line\n"
                        "\n·•· Sample@test-uuid (Calc-No.7) \n"
                        "\t-> 01-02 03:04:05.678 - [WARNING] >>> \n"
                        "warn"
                    ),
                )

    def test_sample_logger_inserts_blank_line_between_formatted_blocks(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
            logger = SampleLogger.open(
                log_path,
                module="Sample@test-uuid",
                time_provider=lambda: fixed_dt,
            )

            logger.info("first block")
            logger.info("second block")
            logger.close()

            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(),
                    (
                        "\n·•· Sample@test-uuid \n"
                        "\t-> 01-02 03:04:05.678 - [INFO] >>> \n"
                        "first block"
                        "\n\n·•· Sample@test-uuid \n"
                        "\t-> 01-02 03:04:05.678 - [INFO] >>> \n"
                        "second block"
                    ),
                )

    def test_likelihood_logger_binds_to_sample_logger_backend(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
            logger = SampleLogger.open(
                log_path,
                module="Sample@test-uuid",
                time_provider=lambda: fixed_dt,
            )

            sample_info = {
                "logger": logger,
                "logger_name": "Sample@test-uuid",
            }
            LogLikelihood([]).update_logger(sample_info).info("likelihood path check")
            logger.close()

            with open(log_path, "r", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(),
                    (
                        "\n·•· Sample@test-uuid (Likelihood) \n"
                        "\t-> 01-02 03:04:05.678 - [INFO] >>> \n"
                        "likelihood path check"
                    ),
                )

    def test_sample_close_logs_summary_and_not_close_kv_footer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
            sample = Sample({"x": 1.0})
            sample.update_uuid("test-uuid")
            sample.logger = SampleLogger.open(
                log_path,
                module="Sample@test-uuid",
                time_provider=lambda: fixed_dt,
            )
            sample.info = {
                "observables": {
                    "x": 1.0,
                    "uuid": "test-uuid",
                    "LogL": -2.5,
                },
                "logger": sample.logger,
            }

            sample.close()

            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("Sample SUMMARY", text)
            self.assertIn("Sample closed", text)
            self.assertIn("uuid", text)
            self.assertNotIn("Field | Value", text)
            self.assertNotIn("[Sample] sample close", text)

    def test_set_config_uses_task_result_dir_when_both_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample = Sample({"x": 1.0})
            sample.set_config({"task_result_dir": tmpdir})
            try:
                expected_root = os.path.join(tmpdir, "SAMPLE")
                self.assertTrue(sample.info["save_dir"].startswith(expected_root))
                self.assertTrue(sample.info["run_log"].startswith(expected_root))
                self.assertTrue(os.path.isdir(sample.info["save_dir"]))
            finally:
                sample.close()


if __name__ == "__main__":
    unittest.main()
