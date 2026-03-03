#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import threading
import tempfile
import unittest

import numpy as np


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.grid import grid_sampling  # noqa: E402
from jarvishep.Sampling.variables import Variable  # noqa: E402
from jarvishep.modulePool import ModulePool  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class _FakeModule:
    def __init__(self, base_dir):
        self.name = "FakeModule"
        self.type = "Calculator"
        self.config = {"path": os.path.join(base_dir, "@PackID")}


class _FakeInstance:
    def __init__(self):
        self.PackID = "001"
        self.is_busy = False
        self.is_installed = True
        self.installation_event = threading.Event()
        self.installation_event.set()
        self.calls = 0

    def execute(self, params, sample_info):
        self.calls += 1
        return {"ok": True, "params": params, "uuid": sample_info.get("uuid")}


class _PoisonExecutor:
    def submit(self, *_args, **_kwargs):
        raise AssertionError("Nested submit should not be used in ModulePool.execute hot path")


class TestV165RemainingP1(unittest.TestCase):
    def test_modulepool_execute_without_nested_submit(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-modulepool-") as tmp_dir:
            pool = ModulePool(_FakeModule(tmp_dir), max_workers=2)
            pool.logger = _NoopLogger()
            inst = _FakeInstance()
            pool.instances = [inst]
            pool.executor = _PoisonExecutor()

            out = pool.execute({"x": 1.0}, {"uuid": "sample-1"})
            self.assertTrue(out["ok"])
            self.assertEqual(inst.calls, 1)
            self.assertFalse(inst.is_busy)

    def test_grid_sampling_open_interval_endpoints(self):
        pts = grid_sampling(np.array([3, 4]))
        self.assertTrue(np.all(pts > 0.0))
        self.assertTrue(np.all(pts < 1.0))

        var_normal = Variable(
            name="x",
            description="",
            distribution="Normal",
            parameters={"mean": 0.0, "stddev": 1.0},
        )
        var_logit = Variable(
            name="y",
            description="",
            distribution="Logit",
            parameters={"location": 0.0, "scale": 1.0},
        )

        mapped_x = np.asarray([var_normal.map_standard_random_to_distribution(v) for v in pts[:, 0]], dtype=float)
        mapped_y = np.asarray([var_logit.map_standard_random_to_distribution(v) for v in pts[:, 1]], dtype=float)
        self.assertTrue(np.all(np.isfinite(mapped_x)))
        self.assertTrue(np.all(np.isfinite(mapped_y)))

    def test_grid_sampling_single_step_center(self):
        pts = grid_sampling(np.array([1]))
        self.assertEqual(pts.shape, (1, 1))
        self.assertAlmostEqual(float(pts[0, 0]), 0.5, places=15)


if __name__ == "__main__":
    unittest.main()
