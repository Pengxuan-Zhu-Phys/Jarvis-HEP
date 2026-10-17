"""Failure boundaries and shard preservation for completed replacement scans."""
from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import pytest

from jarvishep2.io import database
from jarvishep2.io.archiver import ArchiverProcess, SimpleArchiver
from jarvishep2.runtime._runtime_supervisor import _RuntimeSupervisor


def _rows(path: Path, indices, *, tag="new"):
    database.SimpleHDF5Writer(str(path)).add_records([
        {"sample_index": index, "uuid": f"{tag}-{index}", "value": tag} for index in indices
    ])


def _records(root):
    return [row for path in database.discover_samples_hdf5_shards(str(root))
            for row in database.SimpleHDF5Writer(path).read_records()]


def test_trim_preserves_partitions_and_latest_rows(tmp_path):
    first = tmp_path / "samples.0001.hdf5"
    second = tmp_path / "samples.0002.hdf5"
    live = tmp_path / "samples.hdf5"
    _rows(first, [0, 1, 7], tag="old")
    _rows(second, [2, 3, 8])
    _rows(live, [0, 4, 9])
    assert database.trim_scan_records(str(tmp_path), 5) == 4
    assert [Path(p).name for p in database.discover_samples_hdf5_shards(str(tmp_path))] == [
        first.name, second.name, live.name,
    ]
    assert [r["sample_index"] for r in database.SimpleHDF5Writer(str(first)).read_records()] == [1]
    assert [r["sample_index"] for r in database.SimpleHDF5Writer(str(second)).read_records()] == [2, 3]
    rows = _records(tmp_path)
    assert sorted(r["sample_index"] for r in rows) == list(range(5))
    assert next(r for r in rows if r["sample_index"] == 0)["uuid"] == "new-0"


def test_clean_shards_are_not_rewritten_or_reexported(tmp_path):
    writer = database.RollingHDF5Writer(str(tmp_path / "samples.hdf5"), max_bytes=1)
    for i in range(3):
        writer.add_records([{"sample_index": i, "uuid": str(i)}])
    writer.close()
    def snapshot():
        return {p.name: (p.stat().st_ino, p.stat().st_mtime_ns, hashlib.sha256(p.read_bytes()).digest())
                for p in tmp_path.iterdir() if p.is_file()}
    before = snapshot()
    with mock.patch.object(database, "_trim_shard") as rewrite:
        assert database.trim_scan_records(str(tmp_path), 3) == 0
    rewrite.assert_not_called()
    assert snapshot() == before
    assert len(database.discover_samples_hdf5_shards(str(tmp_path))) == 4


def test_failed_shard_replace_keeps_original_data_and_is_retryable(tmp_path):
    first = tmp_path / "samples.0001.hdf5"
    live = tmp_path / "samples.hdf5"
    _rows(first, [0, 1, 5], tag="old")
    _rows(live, [0, 2, 6])
    before = {p: p.read_bytes() for p in (first, live)}
    with mock.patch.object(database.os, "replace", side_effect=OSError("disk error")):
        with pytest.raises(OSError, match="disk error"):
            database.trim_scan_records(str(tmp_path), 3)
    assert all(p.read_bytes() == contents for p, contents in before.items())
    assert database.trim_scan_records(str(tmp_path), 3) == 3
    assert sorted(r["sample_index"] for r in _records(tmp_path)) == [0, 1, 2]
    assert first.exists() and live.exists()


def test_partial_trim_failure_never_removes_latest_copy(tmp_path):
    first, live = tmp_path / "samples.0001.hdf5", tmp_path / "samples.hdf5"
    _rows(first, [0, 1, 5], tag="old")
    _rows(live, [0, 2, 6])
    original = database._trim_shard
    def fail_live(path, rows, scratch):
        if Path(path) == live:
            raise OSError("second shard")
        return original(path, rows, scratch)
    with mock.patch.object(database, "_trim_shard", side_effect=fail_live):
        with pytest.raises(OSError, match="second shard"):
            database.trim_scan_records(str(tmp_path), 3)
    assert next(r for r in _records(tmp_path) if r["sample_index"] == 0)["uuid"] == "new-0"
    database.trim_scan_records(str(tmp_path), 3)
    assert sorted(r["sample_index"] for r in _records(tmp_path)) == [0, 1, 2]


def _core():
    return SimpleNamespace(
        _shutdown_done=False, _interrupt_requested=False, _scan_completed=True,
        _replacing_existing=True, config={"Sampling": {"Method": "Grid"}},
        sampler=mock.Mock(uses_indexed_resume=True, _accepted_index=3),
        archiver=None, factory=mock.Mock(), redis=mock.Mock(), _logger=mock.Mock(),
        _finalize_sample_buckets=mock.Mock(), _archiver_records_written=mock.Mock(return_value=3),
        _resolve_database_dir=mock.Mock(return_value="unused"), write_run_summary=mock.Mock(),
        _release_redis_control_lock=mock.Mock(), _stop_managed_redis=mock.Mock(),
    )


def test_trim_failure_still_writes_summary_releases_lock_and_stops_workers():
    core = _core()
    factory, redis = core.factory, core.redis
    with mock.patch.object(database, "trim_scan_records", side_effect=OSError("disk full")):
        with pytest.raises(RuntimeError, match="DATABASE cleanup.*disk full"):
            _RuntimeSupervisor(core).shutdown(write_run_summary=True)
    core.write_run_summary.assert_called_once()
    core._release_redis_control_lock.assert_called_once()
    redis.clear_runtime_metadata_path.assert_called_once()
    factory.shutdown.assert_called_once_with(wait=True)
    core._stop_managed_redis.assert_called_once()
    assert core._shutdown_done is False
    with mock.patch.object(database, "trim_scan_records", return_value=0):
        _RuntimeSupervisor(core).shutdown()
    assert core._shutdown_done is True


def test_archiver_timeout_forces_cleanup_and_never_trims_even_on_retry():
    core = _core()
    factory = core.factory
    archiver = mock.Mock(spec=ArchiverProcess)
    archiver.stop.return_value = False
    archiver.is_alive.side_effect = [True, True, False, False]
    core.archiver = archiver
    with mock.patch.object(database, "trim_scan_records") as trim:
        with pytest.raises(RuntimeError, match="Archiver did not stop"):
            _RuntimeSupervisor(core).shutdown(write_run_summary=True)
        with pytest.raises(RuntimeError, match="previous Archiver exit was unclean"):
            _RuntimeSupervisor(core).shutdown()
    trim.assert_not_called()
    archiver.terminate.assert_called_once()
    factory.shutdown.assert_called_once()
    core.write_run_summary.assert_called_once()
    assert core._stop_managed_redis.call_count == 2


def test_stop_exception_skips_trim_and_still_finishes_other_cleanup():
    core = _core()
    factory = core.factory
    core.archiver = mock.Mock()
    core.archiver.stop.side_effect = OSError("writer close error")
    with mock.patch.object(database, "trim_scan_records") as trim:
        with pytest.raises(RuntimeError, match="writer close error"):
            _RuntimeSupervisor(core).shutdown(write_run_summary=True)
    trim.assert_not_called()
    factory.shutdown.assert_called_once()
    core._release_redis_control_lock.assert_called_once()
    core.write_run_summary.assert_called_once()


def test_thread_stop_retains_live_thread_and_does_not_close_its_writer():
    archiver = object.__new__(SimpleArchiver)
    thread = mock.Mock()
    thread.is_alive.return_value = True
    archiver._thread = thread
    archiver._stop_event = mock.Mock()
    archiver.processor = mock.Mock()
    with mock.patch.object(archiver, "drain") as drain:
        assert archiver.stop(wait=True) is False
    drain.assert_not_called()
    archiver.processor.writer.close.assert_not_called()
    assert archiver._thread is thread and archiver.is_alive()


def test_original_run_failure_is_not_masked_by_shutdown_error():
    from jarvishep2.runtime.core import Jarvis2Core
    core = Jarvis2Core({"Sampling": {"Method": "Grid"}})
    with mock.patch.object(core, "prepare_resume"), \
         mock.patch.object(core, "_install_control_signal_handlers"), \
         mock.patch.object(core, "_restore_control_signal_handlers"), \
         mock.patch.object(core, "bootstrap_distributed_runtime", side_effect=ValueError("original")), \
         mock.patch.object(core, "shutdown", side_effect=RuntimeError("cleanup")) as shutdown:
        with pytest.raises(ValueError, match="original"):
            core._run_configured()
    shutdown.assert_called_once()


@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_replacement_scan_only_emits_plots_after_successful_cleanup(cleanup_fails):
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.run_outcome import RunOutcome
    core = Jarvis2Core({"Sampling": {"Method": "Grid"}})
    core._replacing_existing = True
    calls = []
    def shutdown(**_kwargs):
        calls.append("cleanup")
        if cleanup_fails:
            raise RuntimeError("trim failed")
    with mock.patch.object(core, "prepare_resume"), \
         mock.patch.object(core, "_install_control_signal_handlers"), \
         mock.patch.object(core, "_restore_control_signal_handlers"), \
         mock.patch.object(core, "bootstrap_distributed_runtime"), \
         mock.patch.object(core, "run_distributed_scan", return_value=1), \
         mock.patch.object(core, "_wait_for_archive_caught_up"), \
         mock.patch.object(core, "_capture_run_outcome", return_value=RunOutcome(completed=1)), \
         mock.patch.object(core, "_finalize_nested_result_csv"), \
         mock.patch.object(core, "_emit_plot_scenes", side_effect=lambda: calls.append("plots")), \
         mock.patch.object(core, "shutdown", side_effect=shutdown):
        if cleanup_fails:
            with pytest.raises(RuntimeError, match="trim failed"):
                core._run_configured()
            assert calls == ["cleanup"]
        else:
            core._run_configured()
            assert calls == ["cleanup", "plots"]


@pytest.mark.parametrize("failure", ["thread-alive", "drain", "close", None])
def test_archiver_child_reports_unclean_finalization_and_closes_redis(failure):
    from contextlib import nullcontext
    from jarvishep2.io import archiver as module
    counter = SimpleNamespace(value=0, get_lock=lambda: nullcontext())
    process = SimpleNamespace(
        scan_name="test", archiver_config={}, log_dir=None, log_path=None, redis_config={},
        db_path="unused", sample_root="unused", poll_timeout=1, delete_method="unlink",
        pending_batch=None, records_written=counter, _set_loop_alive=mock.Mock(),
        _stop_event=mock.Mock(is_set=mock.Mock(return_value=True)),
    )
    child = mock.Mock(records_written=3, buckets_packed=0)
    child.stop.return_value = failure != "thread-alive"
    if failure == "drain":
        child.drain.side_effect = OSError("drain failed")
    if failure == "close":
        child.stop.side_effect = OSError("close failed")
    queue = mock.Mock()
    with mock.patch.object(module, "SimpleArchiver", return_value=child), \
         mock.patch.object(module, "RedisQueue", return_value=queue), \
         mock.patch("jarvishep2.proc_title.set_process_title"), \
         mock.patch("jarvishep2.logging.setup_jarvis_logging"), \
         mock.patch("jarvishep2.logging.get_jarvis_logger", return_value=mock.Mock()):
        if failure is not None:
            with pytest.raises(RuntimeError, match="Archiver finalization failed"):
                ArchiverProcess.run(process)
        else:
            ArchiverProcess.run(process)
    child.stop.assert_called_once_with(wait=True, drain=False)
    queue.close.assert_called_once()
    assert counter.value == 3


def test_trim_retry_repairs_sealed_csv_after_hdf5_was_already_replaced(tmp_path):
    import csv
    sealed, live = tmp_path / "samples.0001.hdf5", tmp_path / "samples.hdf5"
    _rows(sealed, [0, 1, 8], tag="old")
    _rows(live, [0, 2])
    database.convert_hdf5_to_csv(str(sealed))
    database.convert_hdf5_to_csv(str(live))
    with mock.patch.object(database, "convert_hdf5_to_csv", side_effect=OSError("CSV disk error")):
        with pytest.raises(OSError, match="CSV disk error"):
            database.trim_scan_records(str(tmp_path), 3)
    assert [r["sample_index"] for r in database.SimpleHDF5Writer(str(sealed)).read_records()] == [1]
    database.trim_scan_records(str(tmp_path), 3)
    with sealed.with_suffix(".csv").open() as handle:
        assert [int(row["sample_index"]) for row in csv.DictReader(handle)] == [1]


def test_trim_to_empty_live_removes_stale_csv_and_its_digest(tmp_path):
    live = tmp_path / "samples.hdf5"
    _rows(live, [0, 1])
    database.convert_hdf5_to_csv(str(live))
    assert live.with_suffix(".csv").exists()
    assert database.trim_scan_records(str(tmp_path), 0) == 2
    assert live.exists() and database.SimpleHDF5Writer(str(live)).read_records() == []
    assert not live.with_suffix(".csv").exists()
    assert not (tmp_path / "samples.csv.md5").exists()
