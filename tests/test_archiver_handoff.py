#!/usr/bin/env python3
"""WP-D4.1 Worker staging mv + Archiver handoff tests."""

from __future__ import annotations

import json
import os
import tempfile
import threading
import time
import unittest
from typing import Any

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.archive_handoff import move_tree, stage_sample_dir
from jarvishep2.archiver import ArchiveProcessor, ArchiverProcess, SimpleArchiver
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.sample import Sample

from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    FIXTURES,
    LIKELIHOOD_EXPRESSIONS,
    _load_csv_points,
    _normalize_database_records,
    _sample_tree_file_sets,
    _start_tcp_fakeredis,
    _worker_config,
)


class ArchiveHandoffUnitTests(unittest.TestCase):
    def test_stage_sample_dir_moves_work_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = os.path.join(tmpdir, "SAMPLE", "uuid-1")
            os.makedirs(work_dir)
            marker = os.path.join(work_dir, "marker.txt")
            with open(marker, "w", encoding="utf-8") as handle:
                handle.write("staged")

            staging_root = os.path.join(tmpdir, "staging")
            staging_path = stage_sample_dir(work_dir, staging_root, "uuid-1")
            self.assertEqual(staging_path, os.path.join(staging_root, "uuid-1"))
            self.assertFalse(os.path.exists(work_dir))
            self.assertTrue(os.path.isfile(os.path.join(staging_path, "marker.txt")))

    def test_same_volume_move_preserves_file_inode(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            source = os.path.join(tmpdir, "src")
            os.makedirs(source)
            marker = os.path.join(source, "marker.txt")
            with open(marker, "w", encoding="utf-8") as handle:
                handle.write("inode")
            inode_before = os.stat(marker).st_ino
            destination = os.path.join(tmpdir, "dst")
            move_tree(source, destination, strategy="move")
            moved_marker = os.path.join(destination, "marker.txt")
            self.assertEqual(os.stat(moved_marker).st_ino, inode_before)


class ArchiveHandoffIntegrationTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_calculator_handoff_lands_in_sample_and_drains_queue(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                with open(
                    os.path.join(FIXTURES, "expected_calculator_records.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_records = json.load(handle)
                with open(
                    os.path.join(FIXTURES, "expected_sample_files.json"),
                    encoding="utf-8",
                ) as handle:
                    expected_files = json.load(handle)

                core = Jarvis2Core(
                    {
                        "Runtime": {
                            "mode": "redis",
                            "workers": 1,
                            "redis": redis_config,
                        },
                        "Calculators": {
                            "Cleanup": {"strategy": "mv_to_staging"},
                            "Archiver": {"mode": "thread", "batch_size": 1},
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                core.init_factory(_worker_config(tmpdir))
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                from jarvishep2.Sampling.sampler import SamplingVirtial

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[EGGBOX_CALC_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                samples = []
                for row in _load_csv_points():
                    sample = sampler._build_sample(
                        np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                    )
                    sample.uuid = str(row["uuid"])
                    samples.append(sample)

                core.submit_samples(samples)
                core.wait_for_results(len(samples), timeout=90.0)
                core.shutdown()

                staging_root = os.path.join(tmpdir, "staging")
                if os.path.isdir(staging_root):
                    self.assertEqual(os.listdir(staging_root), [])

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(records, _normalize_database_records(expected_records))

                sample_root = os.path.join(tmpdir, "SAMPLE")
                tree = _sample_tree_file_sets(sample_root)
                self.assertEqual(len(tree), 10)
                for files in tree:
                    self.assertEqual(files, expected_files)
        finally:
            server.shutdown()
            server.server_close()

    def test_archiver_process_mode_handoff(self) -> None:
        server, redis_config = _start_tcp_fakeredis()
        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                core = Jarvis2Core(
                    {
                        "Runtime": {"mode": "redis", "workers": 1, "redis": redis_config},
                        "Calculators": {
                            "Cleanup": {"strategy": "mv_to_staging"},
                            "Archiver": {"mode": "process", "batch_size": 1},
                        },
                        "task_result_dir": tmpdir,
                    }
                )
                core.init_redis()
                core.init_factory(_worker_config(tmpdir))
                db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
                core.init_archiver(db_path)

                from jarvishep2.Sampling.sampler import SamplingVirtial

                sampler = SamplingVirtial()
                sampler.set_config(core.config)
                sampler.set_execution_plan_template(
                    calculator_modules=[EGGBOX_CALC_MODULE],
                    include_likelihood=True,
                )
                core.set_sampler(sampler)

                row = next(iter(_load_csv_points()))
                sample = sampler._build_sample(
                    np.array([float(row["x"]), float(row["y"])], dtype=np.float64)
                )
                sample.uuid = str(row["uuid"])
                core.submit_samples([sample])
                core.wait_for_results(1, timeout=60.0)
                core.shutdown()

                sample_dir = os.path.join(tmpdir, "SAMPLE", sample.uuid)
                self.assertTrue(os.path.isdir(sample_dir))
                self.assertTrue(os.path.exists(os.path.join(sample_dir, "input.json")))
        finally:
            server.shutdown()
            server.server_close()

    def test_processor_idempotent_when_sample_dir_exists(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_root = os.path.join(tmpdir, "SAMPLE")
            staging_root = os.path.join(tmpdir, "staging")
            uuid = "repeat-sample"
            existing = os.path.join(sample_root, uuid)
            os.makedirs(existing)
            with open(os.path.join(existing, "keep.txt"), "w", encoding="utf-8") as handle:
                handle.write("keep")

            staging_path = os.path.join(staging_root, uuid)
            os.makedirs(staging_path)
            with open(os.path.join(staging_path, "new.txt"), "w", encoding="utf-8") as handle:
                handle.write("new")

            processor = ArchiveProcessor(
                SimpleHDF5Writer(os.path.join(tmpdir, "DATABASE", "samples.hdf5")),
                sample_root=sample_root,
                batch_size=1,
            )
            processor.ingest(
                {
                    "uuid": uuid,
                    "staging_path": staging_path,
                    "observables": {"x": 1.0, "y": 2.0, "z": 3.0},
                }
            )
            processor.flush_batch()

            self.assertTrue(os.path.isfile(os.path.join(existing, "keep.txt")))
            self.assertFalse(os.path.exists(os.path.join(existing, "new.txt")))


if __name__ == "__main__":
    unittest.main()