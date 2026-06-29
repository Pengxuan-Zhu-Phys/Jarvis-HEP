#!/usr/bin/env python3
"""LibDeps / LibraryManager tests (WP-D2.3 design-doc parity)."""

from __future__ import annotations

import os
import tempfile
import threading
import unittest

from jarvishep2.library import LibraryManager

TESTS_ROOT = os.path.dirname(__file__)
SAFE_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "safe_calc", "source")
SHADOW_SOURCE = os.path.join(TESTS_ROOT, "fixtures", "shadow_calc", "source")


class LibraryManagerTests(unittest.TestCase):
    def test_requires_shadow_flag(self) -> None:
        self.assertTrue(LibraryManager.requires_shadow(True))
        self.assertFalse(LibraryManager.requires_shadow(False))

    def test_link_into_sample_creates_symlink(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "000001")
            link_path = LibraryManager().link_into_sample(SAFE_SOURCE, sample_dir, "SafeCalc")
            self.assertTrue(os.path.islink(link_path))
            self.assertEqual(os.path.realpath(link_path), os.path.realpath(SAFE_SOURCE))

    def test_safe_vs_shadow_isolation_policy(self) -> None:
        manager = LibraryManager()
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "safe-sample")
            link_path = manager.link_into_sample(SAFE_SOURCE, sample_dir, "SafeCalc")
            self.assertTrue(os.path.islink(link_path))
            self.assertFalse(manager.requires_shadow(False))
            self.assertTrue(manager.requires_shadow(True))

    def test_concurrent_link_into_sample_is_race_safe(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            errors: list[BaseException] = []
            lock = threading.Lock()

            def _link(sample_index: int) -> None:
                try:
                    sample_dir = os.path.join(tmpdir, "SAMPLE", f"{sample_index:06d}")
                    link_path = LibraryManager().link_into_sample(
                        SAFE_SOURCE,
                        sample_dir,
                        "SafeCalc",
                    )
                    if not os.path.islink(link_path):
                        raise AssertionError(f"expected symlink at {link_path}")
                except BaseException as exc:
                    with lock:
                        errors.append(exc)

            threads = [threading.Thread(target=_link, args=(index,)) for index in range(8)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join(timeout=5.0)
            self.assertEqual(errors, [])

    def test_missing_source_raises_clear_boot_error(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = os.path.join(tmpdir, "SAMPLE", "missing")
            missing = os.path.join(tmpdir, "does-not-exist", "tool")
            with self.assertRaises(FileNotFoundError) as ctx:
                LibraryManager().link_into_sample(missing, sample_dir, "MissingTool")
            self.assertIn("LibDeps source does not exist", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()