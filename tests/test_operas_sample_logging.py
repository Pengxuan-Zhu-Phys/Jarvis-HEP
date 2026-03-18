#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from datetime import datetime

import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Module.operas import OperasModule  # noqa: E402
from jarvishep.sample_logger import SampleLogger  # noqa: E402


class _FakeOperasRegistry:
    @staticmethod
    def resolve_name(name):
        return name

    @staticmethod
    def get(_name):
        def _operator(x=None, y=None, uuid=None, observables=None):
            return {"z": float(x) + float(y)}

        return _operator

    @staticmethod
    def call(_name, logger=None, **kwargs):
        if logger is not None:
            logger.debug("dispatched call")
        return {"z": float(kwargs["x"]) + float(kwargs["y"])}


class OperasSampleLoggingTests(unittest.TestCase):
    def test_operas_logs_use_likelihood_style_and_suppress_duplicate_kwargs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            log_path = os.path.join(tmpdir, "Sample_running.log")
            fixed_dt = datetime(2026, 1, 2, 3, 4, 5, 678000)
            slogger = SampleLogger.open(
                log_path,
                module="Sample@test-uuid",
                time_provider=lambda: fixed_dt,
            )
            module = OperasModule(
                "EggBox",
                {
                    "operator": "helper.eggbox2d",
                    "input": [],
                    "output": [{"name": "z"}],
                    "kwargs": {},
                    "call_mode": "call",
                },
            )
            module._registry = _FakeOperasRegistry()

            sample_info = {
                "logger": slogger,
                "logger_name": "Sample@test-uuid",
            }
            observables = {
                "x": np.float64(3.5),
                "y": np.float64(0.25),
                "uuid": "test-uuid",
            }
            result = module.execute(observables, sample_info)
            slogger.close()

            self.assertEqual(result, {"z": 3.75})
            with open(log_path, "r", encoding="utf-8") as handle:
                text = handle.read()
            self.assertIn("Operas input dispatch:", text)
            self.assertIn("   module \t-> EggBox", text)
            self.assertIn("   operator \t-> helper.eggbox2d", text)
            self.assertIn("Operas input observables:\n   with input \t-> [x : 3.5, y : 0.25, uuid : test-uuid]", text)
            self.assertNotIn("Operas input kwargs", text)
            self.assertNotIn("np.float64", text)
            self.assertNotIn("dispatched call", text)


if __name__ == "__main__":
    unittest.main()
