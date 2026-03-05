#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.sample import Sample  # noqa: E402


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
