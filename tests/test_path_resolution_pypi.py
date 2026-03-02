#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.base import Base  # noqa: E402


class PathResolutionPyPITests(unittest.TestCase):
    def test_src_marker_points_to_package_root(self):
        base = Base()
        self.assertTrue(base.path["src_root"].endswith("jarvishep"))
        self.assertTrue(os.path.exists(base.path["args_info"]))
        self.assertTrue(os.path.exists(base.path["logger_config_path"]))
        self.assertTrue(os.path.exists(base.path["preference"]))

    def test_legacy_src_marker_is_not_remapped(self):
        base = Base()
        legacy = "&SRC/" + "src/card/preference.json"
        with self.assertRaises(ValueError):
            base.decode_path(legacy)

    def test_runtime_legacy_src_marker_is_not_remapped(self):
        base = Base()
        with tempfile.TemporaryDirectory() as tmpdir:
            base.configure_runtime_context(cwd=tmpdir)
            legacy = "&J/" + "src/card/preference.json"
            with self.assertRaises(ValueError):
                base.decode_path(legacy)


if __name__ == "__main__":
    unittest.main()
