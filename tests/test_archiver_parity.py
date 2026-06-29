#!/usr/bin/env python3
"""WP-D4.2 Archiver batch persistence + parity gate tests."""

from __future__ import annotations

import json
import os
import tempfile
import threading
import unittest

import numpy as np
from fakeredis import TcpFakeServer

from jarvishep2.archiver import SimpleArchiver
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.factory import TaskFactory
from jarvishep2.redis_queue import make_fakeredis_queue

from test_worker_calculator import (
    EGGBOX_CALC_MODULE,
    FIXTURES,
    _load_csv_points,
    _normalize_database_records,
    _sample_tree_file_sets,
    _start_tcp_fakeredis,
    _worker_config,
)


class ArchiverBatchParityTests(unittest.TestCase):
    def setUp(self) -> None:
        TaskFactory.reset_instance()

    def tearDown(self) -> None:
        TaskFactory.reset_instance()

    def test_batch_size_three_preserves_database_and_sample_parity(self) -> None:
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
                        "Runtime": {"mode": "redis", "workers": 1, "redis": redis_config},
                        "Calculators": {
                            "Cleanup": {"strategy": "mv_to_staging"},
                            "Archiver": {
                                "mode": "thread",
                                "batch_size": 3,
                                "flush_interval_sec": 30.0,
                            },
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

                records = _normalize_database_records(SimpleHDF5Writer(db_path).read_records())
                self.assertEqual(len(records), 10)
                self.assertEqual(records, _normalize_database_records(expected_records))

                tree = _sample_tree_file_sets(os.path.join(tmpdir, "SAMPLE"))
                self.assertEqual(len(tree), 10)
                for files in tree:
                    self.assertEqual(files, expected_files)
        finally:
            server.shutdown()
            server.server_close()

    def test_flush_interval_drains_partial_batch(self) -> None:
        queue = make_fakeredis_queue()
        queue.connect()
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_root = os.path.join(tmpdir, "SAMPLE")
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            archiver = SimpleArchiver(
                queue,
                db_path,
                sample_root=sample_root,
                archiver_config={
                    "batch_size": 10,
                    "flush_interval_sec": 0.05,
                },
            )
            for index in range(2):
                staging = os.path.join(tmpdir, "staging", f"sample-{index}")
                os.makedirs(staging)
                queue.submit_result(
                    {
                        "uuid": f"sample-{index}",
                        "status": "Completed",
                        "staging_path": staging,
                        "observables": {"x": float(index), "z": float(index)},
                    }
                )
            archiver.drain(idle_timeout=0.2)
            self.assertEqual(archiver.records_written, 2)
            self.assertTrue(os.path.isdir(os.path.join(sample_root, "sample-0")))
            self.assertTrue(os.path.isdir(os.path.join(sample_root, "sample-1")))


if __name__ == "__main__":
    unittest.main()