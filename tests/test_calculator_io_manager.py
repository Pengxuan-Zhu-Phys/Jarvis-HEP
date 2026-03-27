#!/usr/bin/env python3
from __future__ import annotations

import asyncio
import json
import os
import sys
import tempfile
import threading
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.IOs.Input import JsonInputFile  # noqa: E402
from jarvishep.IOs.Output import JsonOutputFile  # noqa: E402


class _NoopLogger:
    def __init__(self):
        self._lock = threading.Lock()
        self.messages = []

    def _record(self, level, message):
        with self._lock:
            self.messages.append((level, str(message)))

    def debug(self, message, *args, **kwargs):
        self._record("DEBUG", message)

    def info(self, message, *args, **kwargs):
        self._record("INFO", message)

    def warning(self, message, *args, **kwargs):
        self._record("WARNING", message)

    def error(self, message, *args, **kwargs):
        self._record("ERROR", message)


class _CountingIOManager:
    def __init__(self):
        self.run_blocking_calls = 0

    async def run_blocking(self, fn, *args, **kwargs):
        self.run_blocking_calls += 1
        return fn(*args, **kwargs)

    async def read_text(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def write_text(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def make_dirs(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")

    async def exists(self, *args, **kwargs):
        raise AssertionError("file handler should offload the full file task via run_blocking")


class CalculatorIOManagerTests(unittest.TestCase):
    def test_json_input_file_runs_whole_write_in_io_manager(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "input.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"alpha": 1}, handle)

            io_manager = _CountingIOManager()
            logger = _NoopLogger()
            handler = JsonInputFile(
                "input_cfg",
                path,
                "JSON",
                [{"type": "Dump", "variables": [{"name": "alpha"}]}],
                False,
                logger,
                None,
                tmpdir,
                "DemoModule",
                {},
                io_manager,
            )

            observables = asyncio.run(handler.write({"alpha": 7}))

            self.assertEqual(io_manager.run_blocking_calls, 1)
            self.assertEqual(observables, {})
            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            self.assertEqual(payload["alpha"], 7)

    def test_json_output_file_runs_whole_read_in_io_manager(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, "output.json")
            with open(path, "w", encoding="utf-8") as handle:
                json.dump({"alpha": 9}, handle)

            io_manager = _CountingIOManager()
            logger = _NoopLogger()
            handler = JsonOutputFile(
                "output_cfg",
                path,
                "JSON",
                [{"name": "alpha"}],
                True,
                logger,
                None,
                tmpdir,
                "DemoModule",
                {},
                io_manager,
            )

            observables = asyncio.run(handler.read())

            self.assertEqual(io_manager.run_blocking_calls, 1)
            self.assertEqual(observables["alpha"], 9)
            self.assertTrue(os.path.exists(observables["output_cfg"]))


if __name__ == "__main__":
    unittest.main()
