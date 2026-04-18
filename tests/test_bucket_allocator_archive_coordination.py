#!/usr/bin/env python3
from __future__ import annotations

import os
import queue
import sys
import tempfile
import unittest


PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

from jarvishep.Sampling.bucketallocator import BucketAllocator  # noqa: E402
from jarvishep.Sampling.sample_archive import SampleArchiveManager  # noqa: E402


class TestBucketAllocatorArchiveCoordination(unittest.TestCase):
    def test_sealed_bucket_dispatched_after_inflight_drops_to_zero(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-bucket-alloc-") as tmp:
            callback_paths = []
            alloc = BucketAllocator(
                base_path=tmp,
                limit=2,
                width=3,
                on_bucket_sealed=lambda p: callback_paths.append(os.path.basename(p)),
            )

            b1_a = alloc.next_bucket_dir()
            b1_b = alloc.next_bucket_dir()
            self.assertEqual(callback_paths, [])

            b2_a = alloc.next_bucket_dir()  # seals bucket 001 but should not dispatch yet
            self.assertEqual(os.path.basename(b2_a), "002")
            self.assertEqual(callback_paths, [])

            alloc.mark_sample_finished(b1_a)
            self.assertEqual(callback_paths, [])

            alloc.mark_sample_finished(b1_b)
            self.assertEqual(callback_paths, ["001"])

            # No duplicate dispatch.
            alloc.mark_sample_finished(b1_b)
            self.assertEqual(callback_paths, ["001"])

    def test_seal_current_bucket_dispatches_when_idle(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-bucket-alloc-") as tmp:
            callback_paths = []
            alloc = BucketAllocator(
                base_path=tmp,
                limit=10,
                width=3,
                on_bucket_sealed=lambda p: callback_paths.append(os.path.basename(p)),
            )

            b1 = alloc.next_bucket_dir()
            alloc.mark_sample_finished(b1)
            alloc.seal_current_bucket()
            self.assertEqual(callback_paths, ["001"])

    def test_archive_manager_force_reenqueue_existing_bucket(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-archive-force-") as tmp:
            sample_root = os.path.join(tmp, "SAMPLE")
            os.makedirs(sample_root, exist_ok=True)
            bucket_dir = os.path.join(sample_root, "000001")
            os.makedirs(bucket_dir, exist_ok=True)

            manager = SampleArchiveManager(sample_root=sample_root, enabled=True)
            manager.start = lambda: None
            manager._task_q = queue.Queue(maxsize=8)
            manager._scheduled.add(bucket_dir)

            self.assertFalse(manager.enqueue_bucket_dir(bucket_dir, force=False))
            self.assertTrue(manager.enqueue_bucket_dir(bucket_dir, force=True))
            queued = manager._task_q.get_nowait()
            self.assertEqual(os.path.abspath(bucket_dir), queued)

    def test_archive_manager_uses_project_database_for_check_modules_sample_root(self):
        with tempfile.TemporaryDirectory(prefix="jarvis-archive-opc-") as tmp:
            sample_root = os.path.join(tmp, "SAMPLE", "tests")
            manager = SampleArchiveManager(sample_root=sample_root, enabled=False)
            self.assertEqual(manager.database_root, os.path.join(tmp, "DATABASE"))
            self.assertEqual(
                manager.manifest_jsonl_path,
                os.path.join(tmp, "DATABASE", "archive_manifest.jsonl"),
            )


if __name__ == "__main__":
    unittest.main()
