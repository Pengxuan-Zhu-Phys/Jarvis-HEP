#!/usr/bin/env python3
"""WP-D4.1 Worker staging mv + Archiver handoff tests."""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
import unittest
from unittest import mock

import numpy as np

from jarvishep2.archive_handoff import move_tree, stage_sample_dir
from jarvishep2.archiver import ArchiveProcessor, ArchiverProcess, SimpleArchiver
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import (
    RollingHDF5Writer,
    SimpleHDF5Writer,
    StreamingHDF5Writer,
    read_persisted_outcome_counts,
    read_persisted_uuids,
    recover_stale_hdf5_write_flag,
)
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


class ArchiveHandoffUnitTests(unittest.TestCase):
    def test_resume_probe_uses_write_access_before_clearing_stale_flag(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "samples.hdf5")
            open(db_path, "wb").close()
            with (
                mock.patch(
                    "jarvishep2.database.h5py.File",
                    side_effect=OSError("file is already open for write/SWMR write"),
                ) as open_hdf5,
                mock.patch("jarvishep2.database._clear_stale_swmr_write_flag") as clear,
            ):
                self.assertTrue(recover_stale_hdf5_write_flag(db_path))
            open_hdf5.assert_called_once_with(db_path, "r+", libver="latest")
            clear.assert_called_once_with(db_path)

    def test_streaming_writer_recovers_stale_swmr_status_and_logs_warning(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            os.makedirs(os.path.dirname(db_path))
            open(db_path, "wb").close()
            handle = mock.MagicMock()
            handle.__contains__.return_value = True
            handle.__getitem__.return_value.shape = (0,)
            logger = mock.Mock()
            with (
                mock.patch(
                    "jarvishep2.database.h5py.File",
                    side_effect=[
                        OSError("file is already open for write/SWMR write"),
                        handle,
                    ],
                ) as open_hdf5,
                mock.patch("jarvishep2.database.shutil.which", return_value="/usr/bin/h5clear"),
                mock.patch("jarvishep2.database.subprocess.run") as clear,
            ):
                clear.return_value = subprocess.CompletedProcess([], 0, "", "")
                writer = StreamingHDF5Writer(db_path, logger=logger, recover_stale=True)
            self.assertEqual(open_hdf5.call_count, 2)
            clear.assert_called_once_with(
                ["/usr/bin/h5clear", "-s", db_path],
                capture_output=True,
                text=True,
                check=False,
            )
            logger.warning.assert_called_once()
            writer.close()

    def test_streaming_writer_publishes_one_fsynced_batch(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            with mock.patch("jarvishep2.database.os.fsync") as fsync:
                writer = StreamingHDF5Writer(db_path)
                writer.begin_batch()
                writer.add_record({"uuid": "u-1", "x": 1.0})
                writer.add_record({"uuid": "u-2", "x": 2.0})
                self.assertEqual(writer.commit_batch(), 2)
                self.assertEqual(writer.records_persisted, 2)
                fsync.assert_called_once()
                self.assertEqual(len(SimpleHDF5Writer(db_path).read_records()), 2)
                writer.close()

    def test_processor_persists_when_either_batch_or_interval_is_due(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            writer = StreamingHDF5Writer(db_path)
            processor = ArchiveProcessor(
                writer,
                sample_root=os.path.join(tmpdir, "SAMPLE"),
                batch_size=2,
                flush_interval_sec=60.0,
            )
            self.assertEqual(
                processor.ingest({"uuid": "u-0", "observables": {"x": 0}}), 0
            )
            # A full batch must flush without waiting for the interval.
            self.assertEqual(
                processor.ingest({"uuid": "u-1", "observables": {"x": 1}}), 2
            )
            self.assertEqual(writer.records_persisted, 2)

            # A partial tail must also flush when the interval expires.
            processor._last_flush -= 60.0
            self.assertEqual(
                processor.ingest({"uuid": "u-2", "observables": {"x": 2}}),
                1,
            )
            self.assertEqual(writer.records_persisted, 3)
            writer.close()

    def test_index_prefix_waits_for_delayed_record_and_deduplicates_replay(self) -> None:
        """D21.12 A8: a durable hole may not advance the resume watermark."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            writer = StreamingHDF5Writer(db_path)
            processor = ArchiveProcessor(
                writer,
                sample_root=os.path.join(tmpdir, "SAMPLE"),
                batch_size=1,
            )
            self.assertEqual(
                processor.ingest({"uuid": "u0", "sample_index": 0, "observables": {"x": 0}}),
                1,
            )
            self.assertEqual(
                processor.ingest({"uuid": "u2", "sample_index": 2, "observables": {"x": 2}}),
                1,
            )
            state = processor.persistence_state()
            self.assertEqual(state["persisted_prefix"], 1)
            self.assertEqual(state["pending_persisted_indices"], 1)
            # Replaying the ahead-of-watermark item must remain a no-op.
            self.assertEqual(
                processor.ingest({"uuid": "u2-replay", "sample_index": 2, "observables": {"x": 2}}),
                0,
            )
            self.assertEqual(
                processor.ingest({"uuid": "u1", "sample_index": 1, "observables": {"x": 1}}),
                1,
            )
            self.assertEqual(processor.persistence_state()["persisted_prefix"], 3)
            self.assertEqual(len(SimpleHDF5Writer(db_path).read_records()), 3)
            writer.close()

    def test_process_persistence_state_reads_only_the_o1_index_prefix(self) -> None:
        process = ArchiverProcess.__new__(ArchiverProcess)
        process.redis_config = {"host": "127.0.0.1", "port": 6379, "db": 0}
        process.scan_name = "watermark"
        process.drain = mock.Mock(return_value=17)
        redis = mock.Mock()
        redis.get_archived_index_prefix.return_value = 12
        with mock.patch("jarvishep2.archiver.RedisQueue", return_value=redis):
            state = process.persistence_state()
        self.assertEqual(state, {"persisted_prefix": 12, "records_written": 17})
        redis.get_archived_uuids.assert_not_called()
        redis.connect.assert_called_once()
        redis.close.assert_called_once()

    def test_thread_persistence_state_is_o1_highwater_only(self) -> None:
        """D21.14: do not sorted()-materialise the full UUID set on the hot path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            writer = StreamingHDF5Writer(db_path)
            processor = ArchiveProcessor(
                writer,
                sample_root=os.path.join(tmpdir, "SAMPLE"),
                batch_size=1,
                persisted_uuids={"legacy-a", "legacy-b"},
            )
            state = processor.persistence_state()
            self.assertNotIn("acked_uuids", state)
            self.assertEqual(state["acked_uuids_highwater"], 2)
            self.assertEqual(state["persisted_prefix"], 0)
            writer.close()

    def test_read_persisted_outcome_counts_split_failed_and_completed(self) -> None:
        """D21.13: DATABASE status is the source of truth for resume SAMPLE_STATS."""
        with tempfile.TemporaryDirectory() as tmpdir:
            database_dir = os.path.join(tmpdir, "DATABASE")
            os.makedirs(database_dir, exist_ok=True)
            db_path = os.path.join(database_dir, "samples.hdf5")
            writer = SimpleHDF5Writer(db_path)
            for index, status in enumerate(
                ["Completed", "Failed", "Completed", "Failed", "Completed"]
            ):
                writer.add_record(
                    {
                        "uuid": f"u-{index}",
                        "sample_index": index,
                        "status": status,
                        "x": float(index),
                    }
                )
            completed, failed = read_persisted_outcome_counts(database_dir)
            self.assertEqual((completed, failed), (3, 2))

    def test_database_uuid_seed_prevents_duplicate_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            database_dir = os.path.join(tmpdir, "DATABASE")
            db_path = os.path.join(database_dir, "samples.hdf5")
            SimpleHDF5Writer(db_path).add_record({"uuid": "u-1", "x": 1.0})
            redis = make_fakeredis_queue()
            redis.connect()
            archiver = SimpleArchiver(
                redis,
                db_path,
                sample_root=os.path.join(tmpdir, "SAMPLE"),
                archiver_config={"batch_size": 1},
                scan_name="resume",
            )
            try:
                self.assertEqual(
                    archiver._ingest_result(
                        {"uuid": "u-1", "observables": {"uuid": "u-1", "x": 99.0}}
                    ),
                    0,
                )
                self.assertEqual(len(SimpleHDF5Writer(db_path).read_records()), 1)
                self.assertEqual(read_persisted_uuids(database_dir), {"u-1"})
                self.assertEqual(redis.get_archived_uuids("resume"), {"u-1"})
            finally:
                archiver.stop(wait=False, drain=False)
                redis.close()

    def test_failed_artifact_only_sample_gets_database_status_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            sample_dir = os.path.join(tmpdir, "SAMPLE", "failed-1")
            os.makedirs(sample_dir)
            writer = StreamingHDF5Writer(db_path)
            processor = ArchiveProcessor(
                writer,
                sample_root=os.path.join(tmpdir, "SAMPLE"),
                batch_size=1,
            )
            self.assertEqual(
                processor.ingest(
                    {"uuid": "failed-1", "status": "Failed", "save_dir": sample_dir}
                ),
                1,
            )
            self.assertEqual(
                SimpleHDF5Writer(db_path).read_records(),
                [{"uuid": "failed-1", "status": "Failed"}],
            )
            writer.close()

    def test_rolling_writer_seals_hdf5_and_exports_csv(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = os.path.join(tmpdir, "DATABASE", "samples.hdf5")
            writer = RollingHDF5Writer(db_path, max_bytes=1)
            writer.begin_batch()
            writer.add_record({"uuid": "sealed", "x": 1.0, "LogL": -1.0})
            self.assertEqual(writer.commit_batch(), 1)

            sealed_hdf5 = os.path.join(tmpdir, "DATABASE", "samples.0001.hdf5")
            sealed_csv = os.path.join(tmpdir, "DATABASE", "samples.0001.csv")
            self.assertTrue(os.path.isfile(sealed_hdf5))
            self.assertTrue(os.path.isfile(sealed_csv))
            self.assertEqual(
                SimpleHDF5Writer(sealed_hdf5).read_records(),
                [{"uuid": "sealed", "x": 1.0, "LogL": -1.0}],
            )

            writer.begin_batch()
            writer.add_record({"uuid": "live", "x": 2.0, "LogL": -2.0})
            writer.abort_batch()
            writer.close()
            # A newly-created, empty live shard never creates a bogus CSV.
            self.assertFalse(os.path.isfile(os.path.join(tmpdir, "DATABASE", "samples.csv")))

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
