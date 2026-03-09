#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import types
import unittest
from unittest import mock


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.workflow import Workflow  # noqa: E402


class _CacheClearProbe:
    def __init__(self):
        self.calls = 0

    def cache_clear(self):
        self.calls += 1


class WorkflowFontHandleCleanupTests(unittest.TestCase):
    def test_release_font_handles_calls_cache_clear_when_available(self):
        probe = _CacheClearProbe()
        fake_matplotlib = types.ModuleType("matplotlib")
        fake_matplotlib.font_manager = types.SimpleNamespace(_get_font=probe)

        with mock.patch.dict(sys.modules, {"matplotlib": fake_matplotlib}):
            Workflow._release_matplotlib_font_handles()

        self.assertEqual(probe.calls, 1)

    def test_release_font_handles_noop_when_cache_clear_missing(self):
        fake_matplotlib = types.ModuleType("matplotlib")
        fake_matplotlib.font_manager = types.SimpleNamespace(_get_font=object())

        with mock.patch.dict(sys.modules, {"matplotlib": fake_matplotlib}):
            Workflow._release_matplotlib_font_handles()


if __name__ == "__main__":
    unittest.main()
