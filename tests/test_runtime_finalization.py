"""Shutdown failure boundaries and Archiver configuration."""
from __future__ import annotations

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


def _core():
    return SimpleNamespace(
        _shutdown_done=False, _interrupt_requested=False, _scan_completed=True,
        config={"Sampling": {"Method": "Grid"}},
        sampler=mock.Mock(uses_indexed_resume=True, _accepted_index=3),
        archiver=None, factory=mock.Mock(), redis=mock.Mock(), _logger=mock.Mock(),
        _finalize_sample_buckets=mock.Mock(), _archiver_records_written=mock.Mock(return_value=3),
        _resolve_database_dir=mock.Mock(return_value="unused"), write_run_summary=mock.Mock(),
        _release_redis_control_lock=mock.Mock(), _stop_managed_redis=mock.Mock(),
    )


def test_summary_failure_still_releases_lock_and_stops_workers():
    core = _core()
    factory, redis = core.factory, core.redis
    core.write_run_summary.side_effect = OSError("disk full")
    with pytest.raises(RuntimeError, match="run_summary write.*disk full"):
        _RuntimeSupervisor(core).shutdown(write_run_summary=True)
    core._release_redis_control_lock.assert_called_once()
    redis.clear_runtime_metadata_path.assert_called_once()
    factory.shutdown.assert_called_once_with(wait=True)
    core._stop_managed_redis.assert_called_once()
    assert core._shutdown_done is False
    _RuntimeSupervisor(core).shutdown()
    assert core._shutdown_done is True


def test_archiver_timeout_forces_cleanup_and_reports_unclean_exit_on_retry():
    core = _core()
    factory = core.factory
    archiver = mock.Mock(spec=ArchiverProcess)
    archiver.stop.return_value = False
    archiver.is_alive.side_effect = [True, True, False, False]
    core.archiver = archiver
    with pytest.raises(RuntimeError, match="Archiver did not stop"):
        _RuntimeSupervisor(core).shutdown(write_run_summary=True)
    with pytest.raises(RuntimeError, match="previous Archiver exit was unclean"):
        _RuntimeSupervisor(core).shutdown()
    archiver.terminate.assert_called_once()
    factory.shutdown.assert_called_once()
    core.write_run_summary.assert_called_once()
    assert core._stop_managed_redis.call_count == 2


def test_stop_exception_still_finishes_other_cleanup():
    core = _core()
    factory = core.factory
    core.archiver = mock.Mock()
    core.archiver.stop.side_effect = OSError("writer close error")
    with pytest.raises(RuntimeError, match="writer close error"):
        _RuntimeSupervisor(core).shutdown(write_run_summary=True)
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


def test_completed_scan_emits_plots_before_shutdown():
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.run_outcome import RunOutcome
    core = Jarvis2Core({"Sampling": {"Method": "Grid"}})
    calls = []
    with mock.patch.object(core, "prepare_resume"), \
         mock.patch.object(core, "_install_control_signal_handlers"), \
         mock.patch.object(core, "_restore_control_signal_handlers"), \
         mock.patch.object(core, "bootstrap_distributed_runtime"), \
         mock.patch.object(core, "run_distributed_scan", return_value=1), \
         mock.patch.object(core, "_wait_for_archive_caught_up"), \
         mock.patch.object(core, "_capture_run_outcome", return_value=RunOutcome(completed=1)), \
         mock.patch.object(core, "_finalize_nested_result_csv"), \
         mock.patch.object(core, "_emit_plot_scenes", side_effect=lambda: calls.append("plots")), \
         mock.patch.object(core, "shutdown", side_effect=lambda **_kw: calls.append("cleanup")):
        core._run_configured()
    assert calls == ["plots", "cleanup"]


@pytest.mark.parametrize("failure", ["thread-alive", "drain", "close", None])
def test_archiver_child_reports_unclean_finalization_and_closes_redis(failure):
    from jarvishep2.io import archiver as module
    process = ArchiverProcess({}, db_path="unused", sample_root="unused", scan_name="test")
    process._stop_event.set()
    child = mock.Mock(records_written=3, buckets_packed=0)
    child.stop.return_value = failure != "thread-alive"
    if failure in {"drain", "close"}:
        child.stop.side_effect = OSError(f"{failure} failed")
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
    child.stop.assert_called_once_with(wait=True)
    queue.close.assert_called_once()
    assert process.records_written.value == 3


def _archiver_config_for(sampler, method, tmp_path, policy="fresh"):
    core = SimpleNamespace(
        config={"Sampling": {"Method": method}, "EnvReqs": {"V2": {"archiver": {"mode": "thread"}}}},
        sampler=sampler, redis=mock.Mock(), _logger=mock.Mock(), info={},
        _resume_policy=policy, _resume_checkpoint_payload=None, _control_lock_owner=None,
        _resolve_database_dir=mock.Mock(return_value=str(tmp_path / "DATABASE")),
        _resolve_sample_root=mock.Mock(return_value=str(tmp_path / "SAMPLE")),
        _restore_archiver_persistence=mock.Mock(),
    )
    core.redis.connection_config.return_value = {}
    module = "jarvishep2.runtime._runtime_supervisor"
    with mock.patch(f"{module}.SimpleArchiver") as archiver:
        _RuntimeSupervisor(core).init_archiver()
    return archiver.call_args.kwargs["archiver_config"]


@pytest.mark.parametrize("policy", ["fresh", "resume"])
@pytest.mark.parametrize("method", ["Random", "Grid", "Dynesty", "MultiNest"])
def test_archiver_config_carries_no_resume_state(method, policy, tmp_path):
    # The Archiver only appends; it never learns the sampler or resume policy.
    sampler = SimpleNamespace(method=method)
    config = _archiver_config_for(sampler, method, tmp_path, policy)
    assert "replace_duplicates" not in config
    assert "initial_persisted_prefix" not in config
