#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import types
import unittest

import yaml


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.plot import JarvisPLOT  # noqa: E402


class _NoopLogger:
    def info(self, *_args, **_kwargs):
        return None

    def warning(self, *_args, **_kwargs):
        return None

    def error(self, *_args, **_kwargs):
        return None


class PlotMultiNestSupportTests(unittest.TestCase):
    def test_emit_jplot_uses_multinest_result_dataset_for_multinest_method(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            db_dir = os.path.join(tmpdir, "DATABASE")
            os.makedirs(db_dir, exist_ok=True)
            multinest_csv = os.path.join(db_dir, "multinest_result.csv")
            with open(multinest_csv, "w", encoding="utf-8") as f:
                f.write("uuid,log_Like\ns1,1.0\n")

            out_yaml = os.path.join(tmpdir, "jplot.yaml")
            scan_yaml = types.SimpleNamespace(
                config={
                    "Sampling": {
                        "Method": "MultiNest",
                        "Variables": [
                            {
                                "name": "x",
                                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
                            },
                            {
                                "name": "y",
                                "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}},
                            },
                        ],
                    }
                }
            )

            plotter = JarvisPLOT()
            plotter.logger = _NoopLogger()
            plotter.info = {
                "plot": {"config": out_yaml},
                "sample": {"task_result_dir": tmpdir},
                "db": {"out_csv": False, "info": os.path.join(tmpdir, "dbinfo.json")},
            }
            plotter.emit_jplot(scan_yaml)

            with open(out_yaml, "r", encoding="utf-8") as f:
                payload = yaml.safe_load(f)

            datasets = payload.get("DataSet", [])
            dynesty_like = [d for d in datasets if d.get("name") == "dynesty"]
            self.assertEqual(len(dynesty_like), 1)
            self.assertEqual(dynesty_like[0].get("path"), "DATABASE/multinest_result.csv")


if __name__ == "__main__":
    unittest.main()
