#!/usr/bin/env python3
"""WP-D3.3 Runtime.FileOperation.delete_method tests."""

from __future__ import annotations

import os
import tempfile
import unittest

import numpy as np

from jarvishep2.archiver import SimpleArchiver
from jarvishep2.file_ops import (
    DEFAULT_DELETE_METHOD,
    delete_path,
    delete_paths,
    normalize_delete_method,
)
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.runtime_config import get_delete_method, normalize_file_operation, normalize_runtime_block
from jarvishep2.sample import Sample
from jarvishep2.worker import Worker


class DeletePathTests(unittest.TestCase):
    def test_default_method_is_shutil(self) -> None:
        self.assertEqual(DEFAULT_DELETE_METHOD, "shutil")

    def test_shutil_deletes_file(self) -> None:
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            path = handle.name
        try:
            delete_path(path, method="shutil")
            self.assertFalse(os.path.lexists(path))
        finally:
            if os.path.lexists(path):
                os.remove(path)

    def test_shutil_deletes_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            nested = os.path.join(tmpdir, "nested")
            os.makedirs(nested)
            marker = os.path.join(nested, "marker.txt")
            with open(marker, "w", encoding="utf-8") as handle:
                handle.write("x")
            delete_path(tmpdir, method="shutil")
            self.assertFalse(os.path.lexists(tmpdir))

    def test_rm_deletes_file(self) -> None:
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            path = handle.name
        try:
            delete_path(path, method="rm")
            self.assertFalse(os.path.lexists(path))
        finally:
            if os.path.lexists(path):
                os.remove(path)

    def test_rm_deletes_directory(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            nested = os.path.join(tmpdir, "nested")
            os.makedirs(nested)
            delete_path(tmpdir, method="rm")
            self.assertFalse(os.path.lexists(tmpdir))

    def test_bad_method_raises(self) -> None:
        with tempfile.NamedTemporaryFile(delete=False) as handle:
            path = handle.name
        try:
            with self.assertRaises(ValueError) as ctx:
                delete_path(path, method="rsync")
            self.assertIn("invalid delete_method", str(ctx.exception))
        finally:
            os.remove(path)

    def test_missing_ok_skips_absent_path(self) -> None:
        delete_path("/tmp/jarvis2-does-not-exist-delete-path", method="shutil", missing_ok=True)

    def test_delete_paths_ignores_blank_entries(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            first = os.path.join(tmpdir, "first")
            second = os.path.join(tmpdir, "second")
            os.makedirs(first)
            os.makedirs(second)
            delete_paths(["", first, "  ", second], method="shutil")
            self.assertFalse(os.path.lexists(first))
            self.assertFalse(os.path.lexists(second))


class RuntimeConfigDeleteMethodTests(unittest.TestCase):
    def test_normalize_delete_method_defaults_invalid_to_shutil(self) -> None:
        self.assertEqual(normalize_delete_method(None), "shutil")
        self.assertEqual(normalize_delete_method("RM"), "rm")
        self.assertEqual(normalize_delete_method("bogus"), "shutil")

    def test_normalize_file_operation_parses_runtime_block(self) -> None:
        parsed = normalize_file_operation({"delete_method": "rm"})
        self.assertEqual(parsed, {"delete_method": "rm"})

    def test_runtime_block_includes_file_operation(self) -> None:
        runtime = normalize_runtime_block(
            {"FileOperation": {"delete_method": "rm"}, "workers": 2}
        )
        self.assertEqual(runtime["FileOperation"], {"delete_method": "rm"})
        self.assertEqual(runtime["workers"], 2)

    def test_get_delete_method_from_config(self) -> None:
        cfg = {"Runtime": {"FileOperation": {"delete_method": "rm"}}}
        self.assertEqual(get_delete_method(cfg), "rm")
        self.assertEqual(get_delete_method({}), "shutil")


class WorkerCleanupTests(unittest.TestCase):
    def test_collect_cleanup_paths_reads_sample_info(self) -> None:
        worker = Worker(0, {"host": "127.0.0.1", "port": 6379, "db": 0}, {})
        sample = Sample(uuid="u1", observables={}, info={})
        sample.info["cleanup_paths"] = ["/tmp/a", "/tmp/b"]
        sample.info["staging_path"] = "/tmp/c"
        self.assertEqual(
            worker._collect_cleanup_paths(sample),
            ["/tmp/a", "/tmp/b", "/tmp/c"],
        )

    def test_cleanup_transient_paths_deletes_configured_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            target = os.path.join(tmpdir, "staging")
            os.makedirs(target)
            worker = Worker(
                0,
                {"host": "127.0.0.1", "port": 6379, "db": 0},
                {"delete_method": "rm"},
            )
            worker._delete_method = "rm"
            sample = Sample(uuid="u1", observables={}, info={"cleanup_paths": [target]})
            worker._cleanup_transient_paths(sample)
            self.assertFalse(os.path.lexists(target))

    def test_process_task_deletes_cleanup_paths_from_sample_config(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            staging = os.path.join(tmpdir, "staging", "sample-uuid")
            os.makedirs(staging)
            marker = os.path.join(staging, "done.txt")
            with open(marker, "w", encoding="utf-8") as handle:
                handle.write("staged")

            worker = Worker(
                0,
                {"host": "127.0.0.1", "port": 6379, "db": 0},
                {
                    "delete_method": "shutil",
                    "sample_config": {
                        "task_result_dir": tmpdir,
                        "sample_dirs": os.path.join(tmpdir, "SAMPLE"),
                        "sample_artifacts": "never",
                        "workflow_has_calculator": False,
                        "workflow_references_sdir": False,
                        "cleanup_paths": [staging],
                    },
                    "mapper": {"type": "identity", "keys": []},
                    "opera_modules": {},
                    "calculator_modules": [],
                    "likelihood_expressions": [],
                },
            )
            worker._init_runtime()
            try:
                worker.process_task(
                    {
                        "uuid": "sample-uuid",
                        "u_coords": np.array([], dtype=np.float64),
                        "execution_plan": [],
                        "opera_params": {},
                        "sample_artifacts": "never",
                    }
                )
                self.assertFalse(os.path.lexists(staging))
            finally:
                if worker._scheduler is not None:
                    worker._scheduler.shutdown(wait=True)


class ArchiverCleanupTests(unittest.TestCase):
    def test_cleanup_staging_uses_delete_method(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            target = os.path.join(tmpdir, "staging")
            os.makedirs(target)
            queue = make_fakeredis_queue()
            archiver = SimpleArchiver(queue, os.path.join(tmpdir, "samples.hdf5"), delete_method="rm")
            archiver.cleanup_staging([target])
            self.assertFalse(os.path.lexists(target))


if __name__ == "__main__":
    unittest.main()