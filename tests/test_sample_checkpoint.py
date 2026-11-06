"""Completion-based finite sampler snapshots, independent of HDF5/Archiver."""

from collections import Counter
from copy import deepcopy
from pathlib import Path
import pickle
from types import SimpleNamespace

import numpy as np
import pytest

from jarvishep2.distributor import Distributor
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.completion_sampler import PointSetSampler
from jarvishep2.sampling.sample_checkpoint import (
    atomic_pickle, load_runtime, read_pickle,
)

METHODS = ["Grid", "Random", "Bridson", "CSV"]


def config(tmp_path, method, *, enabled=True, selection=None, window=4):
    csv = tmp_path / "points.csv"
    csv.write_text("uuid,x,y\n" + "".join(f"row-{i},{i / 20},0.5\n" for i in range(20)))
    bounds = {
        "Grid": {"seed": 7},
        "Random": {"seed": 7, "point_number": 20},
        "Bridson": {"seed": 7, "radius": 0.2, "max_attempt": 20, "max_worker": window},
        "CSV": {"path": str(csv)},
    }[method]
    return {
        "task_root": str(tmp_path), "scan_name": "completion",
        "Sampling": {
            "Method": method, "Bounds": bounds, "selection": selection,
            "Variables": [
                {"name": name, "distribution": {"type": "Flat", "parameters": {
                    "min": 0, "max": 1, "length": 1, "num": 5,
                }}} for name in ("x", "y")
            ] if method != "CSV" else [],
        },
        "EnvReqs": {"V2": {
            "workers": window, "batch_size": 3,
            "checkpoint": {"enabled": enabled, "heartbeat": 30},
        }},
    }


def make_sampler(cfg):
    sampler = Distributor.set_method(cfg["Sampling"]["Method"])
    sampler.set_config(cfg)
    if sampler.method == "Bridson":
        # The checkpoint tests exercise a precomputed point cloud, not its
        # Poisson-disk construction or external calculator programs.
        sampler._P = np.arange(40, dtype=float).reshape(20, 2) / 40
        sampler.info["NSamples"] = 20
    return sampler


def terminal(sample, status="Completed"):
    return {"uuid": sample.uuid, "sample_index": sample.sample_index, "terminal_status": status}


def prepare_points(sampler, path):
    if isinstance(sampler, PointSetSampler):
        sampler._prepare_checkpoint_data(str(path))


@pytest.mark.parametrize("method", METHODS)
def test_out_of_order_completion_restores_holes_and_then_new_points(tmp_path, method):
    cfg = config(tmp_path, method)
    sampler = make_sampler(cfg)
    path = tmp_path / method / "runtime.pickle"
    prepare_points(sampler, path)
    samples = [sampler.propose_next() for _ in range(9)]
    for i in [8, 1, 7, 0, 4, 2, 5]:
        assert sampler.accept_completion(terminal(samples[i], "Failed" if i == 4 else "Completed"))
    assert not sampler.accept_completion(terminal(samples[8]))
    state = sampler.checkpoint_state()
    assert state["next_sample_index"] == 9
    assert set(state["unfinished"]) == {3, 6}
    assert (state["completed_count"], state["failed_count"]) == (6, 1)

    restored = make_sampler(cfg)
    if isinstance(restored, PointSetSampler):
        restored.prepare_resume_points(str(path), state)
    restored.restore_checkpoint(state)
    for i in (3, 6):
        sample = restored._rebuild_unfinished(i)
        assert sample.uuid == samples[i].uuid
        assert sample.sample_index == i
        np.testing.assert_array_equal(sample.u_coords, samples[i].u_coords)
        assert sample.opera_params == samples[i].opera_params
    expected, actual = sampler.propose_next(), restored.propose_next()
    assert expected.uuid == actual.uuid
    np.testing.assert_array_equal(expected.u_coords, actual.u_coords)
    assert expected.opera_params == actual.opera_params


class ReverseQueue:
    """Return newest results first; keep index zero outstanding until the end."""

    def __init__(self, total, sampler):
        self.total, self.sampler = total, sampler
        self.sent, self.pending = [], []
        self.peak = 0

    def read_proc_board(self, role):
        return {}

    def push_task(self, task):
        assert task["sample_index"] in self.sampler._unfinished
        self.sent.append(task)
        self.pending.append(task)
        self.peak = max(self.peak, len(self.pending))

    def push_many_tasks(self, tasks):
        for task in tasks:
            self.push_task(task)

    def poll_feedback_batch(self, limit):
        records = []
        for task in list(reversed(self.pending)):
            if task["sample_index"] == 0 and len(self.sent) < self.total:
                continue
            if len(records) >= limit:
                break
            self.pending.remove(task)
            records.append({
                "uuid": task["uuid"], "sample_index": task["sample_index"],
                "terminal_status": "Failed" if task["sample_index"] == 5 else "Completed",
            })
        return records


@pytest.mark.parametrize("method", METHODS)
def test_disabled_checkpoint_still_consumes_feedback_with_bounded_window(tmp_path, method):
    sampler = make_sampler(config(tmp_path, method, enabled=False))
    total = 25 if method == "Grid" else 20
    queue = ReverseQueue(total, sampler)
    sampler.set_redis(queue)
    sampler.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"))
    assert sampler.run_distributed(timeout=2) == total
    assert not sampler._unfinished
    assert sampler._completed_count == total - 1
    assert sampler._failed_count == 1
    assert queue.peak <= sampler._max_inflight
    assert Counter(t["sample_index"] for t in queue.sent) == Counter(range(total))
    assert not (tmp_path / "runtime.pickle").exists()


@pytest.mark.parametrize("method", METHODS)
def test_resume_drains_saved_window_even_with_fewer_workers(tmp_path, method):
    cfg = config(tmp_path, method, window=8)
    original = make_sampler(cfg)
    path = tmp_path / method / "runtime.pickle"
    prepare_points(original, path)
    for _ in range(9):
        original.propose_next()
    state = original.checkpoint_state()
    resumed_cfg = deepcopy(cfg)
    resumed_cfg["EnvReqs"]["V2"].update(workers=1, batch_size=1)
    if method == "Bridson":
        resumed_cfg["Sampling"]["Bounds"]["max_worker"] = 1
    restored = make_sampler(resumed_cfg)
    if isinstance(restored, PointSetSampler):
        restored.prepare_resume_points(str(path), state)
    restored.restore_checkpoint(state)
    restored._max_inflight = restored._batch_size = 1
    total = 25 if method == "Grid" else 20
    queue = ReverseQueue(total, restored)
    # Don't delay the replay's first entry: this test isolates the reduced cap.
    queue.total = 0
    restored.set_redis(queue)
    assert restored.run_distributed(timeout=2) == total
    assert queue.peak == 1
    assert not restored._unfinished
    assert Counter(t["sample_index"] for t in queue.sent) == Counter(range(total))


@pytest.mark.parametrize("method", METHODS)
def test_selection_does_not_confuse_source_positions_with_sample_indices(tmp_path, method):
    sampler = make_sampler(config(tmp_path, method, selection="x > 0.4"))
    sample = sampler.propose_next()
    assert sample.sample_index == 0
    assert set(sampler._unfinished) == {0}
    state = sampler.checkpoint_state()
    assert state["next_sample_index"] == 1
    assert sampler.accept_completion(terminal(sample))
    assert state["unfinished"]  # Snapshot is detached from the live state.


def test_random_resume_does_not_depend_on_global_numpy_state(tmp_path):
    sampler = make_sampler(config(tmp_path, "Random"))
    sampler.propose_next()
    state = sampler.checkpoint_state()
    expected = sampler.propose_next()
    np.random.seed(999)
    np.random.random(100)
    restored = make_sampler(sampler.config)
    restored.restore_checkpoint(state)
    np.testing.assert_array_equal(restored.propose_next().u_coords, expected.u_coords)


@pytest.mark.parametrize("method", ["Grid", "Bridson"])
def test_all_pickle_is_written_once_and_resume_never_generates_points(tmp_path, method, monkeypatch):
    sampler = make_sampler(config(tmp_path, method))
    path = tmp_path / method / "runtime.pickle"
    sampler.configure_checkpoint(checkpoint_file=str(path))
    try:
        all_path = path.with_name("all.pickle")
        original = all_path.read_bytes()
        for _ in range(4):
            sample = sampler.propose_next()
            sampler.accept_completion(terminal(sample))
        sampler.save_checkpoint(force=True)
    finally:
        sampler.shutdown_checkpointing()
    assert all_path.read_bytes() == original
    payload = load_runtime(str(path), sampler.config)
    restored = make_sampler(sampler.config)
    monkeypatch.setattr(restored, "initialize", lambda **kwargs: pytest.fail("regenerated point set"))
    restored._P = None
    restored.prepare_resume_points(str(path), payload["sampler_state"])
    restored.restore_checkpoint(payload["sampler_state"])
    assert restored.propose_next().sample_index == 4


def test_runtime_size_depends_on_unfinished_not_completed_history(tmp_path):
    sampler = make_sampler(config(tmp_path, "Random"))
    sample = sampler.propose_next()
    before = len(pickle.dumps(sampler.checkpoint_state()))
    sampler._accepted_index = sampler._index = 10_000_001
    sampler._completed_count = 10_000_000
    after = len(pickle.dumps(sampler.checkpoint_state()))
    assert after < before + 64
    assert set(sampler._unfinished) == {sample.sample_index}


def test_busy_writer_skips_capture_without_blocking_sampler(tmp_path, monkeypatch):
    sampler = make_sampler(config(tmp_path, "Random"))
    path = tmp_path / "runtime.pickle"
    sampler.configure_checkpoint(checkpoint_file=str(path))
    writer = sampler._pickle_writer
    try:
        assert writer.wait()
        writer._idle.value = 0
        writer.request_snapshot()
        monkeypatch.setattr(sampler, "checkpoint_state", lambda: pytest.fail("captured while writer busy"))
        assert not sampler.save_checkpoint()
        assert sampler.propose_next().sample_index == 0
    finally:
        writer._idle.value = 1
        sampler.shutdown_checkpointing()
    assert read_pickle(str(path))["sampler_state"]["next_sample_index"] == 0


def test_runtime_resume_refuses_missing_legacy_and_changed_configuration(tmp_path):
    cfg = config(tmp_path, "Random")
    core = Jarvis2Core(cfg)
    with pytest.raises(FileNotFoundError):
        core.prepare_resume(resume=True)
    path = Path(core.checkpoint_file())
    atomic_pickle(str(path), {"format": "jarvis-hep.v2-distributed"})
    with pytest.raises(ValueError, match="format"):
        core.prepare_resume(resume=True)
    sampler = make_sampler(cfg)
    sampler.configure_checkpoint(checkpoint_file=str(path))
    sampler.shutdown_checkpointing()
    core.prepare_resume(resume=True)
    core.config["Sampling"]["Bounds"]["seed"] = 8
    with pytest.raises(ValueError, match="configuration"):
        core.prepare_resume(resume=True)


def test_checkpoint_and_restore_never_read_database_state(tmp_path, monkeypatch):
    import jarvishep2.runtime._resume_service as resume

    def forbidden(*args, **kwargs):
        pytest.fail("finite sampler checkpoint touched the database")

    monkeypatch.setattr(resume, "read_persisted_sample_index_state", forbidden)
    monkeypatch.setattr(resume, "read_persisted_outcome_counts", forbidden)
    core = Jarvis2Core(config(tmp_path, "Random"))
    sampler = make_sampler(core.config)
    sample = sampler.propose_next()
    sampler.accept_completion(terminal(sample, "Failed"))
    core._resume_checkpoint_payload = {"sampler_state": sampler.checkpoint_state()}
    core._load_persisted_database_state(str(tmp_path / "missing-DATABASE"))
    assert core._persisted_records_count == core._persisted_failed_count == 1
    core.sampler = sampler
    core._resume_policy = "resume"
    assert core._resumed_record_base() == 0


def test_bad_counters_or_feedback_do_not_silently_advance(tmp_path):
    sampler = make_sampler(config(tmp_path, "Random"))
    sample = sampler.propose_next()
    with pytest.raises(ValueError, match="index"):
        sampler.accept_completion({**terminal(sample), "sample_index": 99})
    assert sample.sample_index in sampler._unfinished
    state = sampler.checkpoint_state()
    state["completed_count"] = 1
    with pytest.raises(ValueError, match="partition"):
        make_sampler(sampler.config).restore_checkpoint(state)


@pytest.mark.parametrize("method", METHODS)
def test_worker_blueprint_requires_completion_feedback_even_without_checkpoint(tmp_path, method):
    from jarvishep2.runtime.worker_config import build_worker_config

    blueprint = build_worker_config(
        config(tmp_path, method, enabled=False), task_result_dir=str(tmp_path),
        extra={"publish_feedback": False, "completion_feedback": False},
    )
    assert blueprint["publish_feedback"] is True
    assert blueprint["completion_feedback"] is True


@pytest.mark.parametrize("status", ["Completed", "Failed"])
def test_worker_reports_terminal_identity_independently_of_archive_ack(status):
    from jarvishep2.runtime.worker import Worker

    events = []
    queue = SimpleNamespace(
        submit_result=lambda info: events.append(("archive", info)),
        publish_feedback=lambda info: events.append(("feedback", info)),
    )
    worker = SimpleNamespace(_redis=queue, worker_config={"completion_feedback": True})
    sample = SimpleNamespace(
        uuid="point-3", sample_index=3, status=status,
        to_info_dict=lambda: {"uuid": "point-3", "status": status},
    )
    Worker._stage_and_submit(worker, sample)
    assert events[-1] == ("feedback", {
        "uuid": "point-3", "sample_index": 3, "terminal_status": status,
    })


@pytest.mark.parametrize("damage", ["missing", "corrupt", "different-dataset"])
def test_missing_or_damaged_all_pickle_refuses_resume(tmp_path, damage):
    sampler = make_sampler(config(tmp_path, "Grid"))
    path = tmp_path / "runtime.pickle"
    sampler.configure_checkpoint(checkpoint_file=str(path))
    sampler.shutdown_checkpointing()
    state = read_pickle(str(path))["sampler_state"]
    all_path = path.with_name("all.pickle")
    if damage == "missing":
        all_path.unlink()
    elif damage == "corrupt":
        all_path.write_bytes(b"not a pickle")
    else:
        payload = read_pickle(str(all_path))
        payload["dataset_id"] = "a-different-point-set"
        atomic_pickle(str(all_path), payload)
    with pytest.raises((OSError, ValueError, pickle.UnpicklingError)):
        make_sampler(sampler.config).prepare_resume_points(str(path), state)
    assert read_pickle(str(path))["sampler_state"] == state


def test_csv_changed_source_is_rejected_without_overwriting_runtime(tmp_path):
    cfg = config(tmp_path, "CSV")
    sampler = make_sampler(cfg)
    path = tmp_path / "runtime.pickle"
    sampler.configure_checkpoint(checkpoint_file=str(path))
    sampler.shutdown_checkpointing()
    original = path.read_bytes()
    source = Path(cfg["Sampling"]["Bounds"]["path"])
    source.write_text(source.read_text().replace("0.5", "0.6"))
    with pytest.raises(ValueError, match="CSV source file changed"):
        make_sampler(cfg).restore_checkpoint(read_pickle(str(path))["sampler_state"])
    assert path.read_bytes() == original


def test_snapshot_failure_keeps_the_previous_atomic_pickle(tmp_path, monkeypatch):
    import jarvishep2.sampling.sample_checkpoint as storage

    path = tmp_path / "runtime.pickle"
    atomic_pickle(str(path), {"previous": 1})
    original = path.read_bytes()

    def fail_dump(*args, **kwargs):
        raise OSError("injected write failure")

    monkeypatch.setattr(storage.pickle, "dump", fail_dump)
    with pytest.raises(OSError):
        atomic_pickle(str(path), {"next": 2})
    assert path.read_bytes() == original
    assert not list(tmp_path.glob(".pickle-*.tmp"))


def _adaptive_config(tmp_path):
    return {
        "task_root": str(tmp_path),
        "scan_name": "adaptive",
        "Sampling": {
            "Method": "AdaptiveBridson",
            "Bounds": {
                "target_expression": "r2",
                "target_value": 0.04,
                "initial_radius": 0.1,
                "seed": 3,
            },
            "Variables": [
                {"name": name, "distribution": {"type": "Flat", "parameters": {
                    "min": 0, "max": 1,
                }}} for name in ("x", "y")
            ],
        },
        "EnvReqs": {"V2": {"workers": 1, "checkpoint": {"enabled": True, "heartbeat": 30}}},
    }


class _Heartbeat:
    """In-process stand-in: the scan only observes a due flag and an idle slot."""

    def __init__(self):
        self.flag = False
        self.idle = True
        self.payloads = []
        self.waits = 0

    def requested(self):
        return self.flag

    def wait(self, timeout=10.0):
        self.waits += 1
        return True

    def available(self):
        return self.idle

    def submit(self, payload):
        assert self.idle
        self.payloads.append(payload)
        self.flag = False
        return True


def test_adaptive_bridson_heartbeat_skips_copy_and_ignores_redis(tmp_path):
    from jarvishep2.Sampling.adaptive_bridson import AdaptiveBridsonSampler, LevelSetPoint

    sampler = AdaptiveBridsonSampler()
    sampler.set_config(_adaptive_config(tmp_path))
    writer = _Heartbeat()
    sampler.configure_checkpoint(checkpoint_file=str(tmp_path / "runtime.pickle"), writer=writer)
    assert writer.waits == 1
    assert len(writer.payloads) == 1

    sampler._points = [LevelSetPoint(u=[0.2, 0.3], x={"x": 0.2, "y": 0.3}, f=0.04, uuid="p1")]
    sampler._uuid_to_index = {"p1": 0}
    sampler._pending_uuids = {"p1"}
    sampler._accepted_index = 1
    sampler._radius = 0.05
    sampler._last_densify_new = 4
    sampler._last_densify_attempted = True

    def forbidden_copy():
        pytest.fail("idle heartbeat copied AdaptiveBridson state")

    original = sampler.checkpoint_state
    sampler.checkpoint_state = forbidden_copy
    try:
        assert sampler.save_checkpoint() is False
        writer.flag = True
        writer.idle = False
        assert sampler.save_checkpoint() is False
        assert sampler.checkpoint_at_barrier(reason="als_generation_0") is False
    finally:
        sampler.checkpoint_state = original
    assert writer.waits == 1
    assert len(writer.payloads) == 1

    class BoomRedis:
        def __getattr__(self, name):
            raise AssertionError(f"checkpoint touched redis.{name}")

    sampler.redis = BoomRedis()
    writer.idle = True
    writer.flag = True
    assert sampler.save_checkpoint() is True
    payload = writer.payloads[-1]
    assert payload["format"] == "jarvis-hep.sampler-pickle"
    assert payload["method"] == "adaptivebridson"
    assert "persistence" not in payload
    point = payload["sampler_state"]["points"][0]
    assert point["uuid"] == "p1"
    assert point["f"] == 0.04
    assert payload["sampler_state"]["pending_uuids"] == ["p1"]
    assert payload["sampler_state"]["last_densify_new"] == 4
    assert writer.waits == 1

    restored = AdaptiveBridsonSampler()
    restored.set_config(sampler.config)
    restored.restore_checkpoint(payload["sampler_state"])
    assert restored._points[0].f == 0.04
    assert restored._points[0].uuid == "p1"
    assert restored._pending_uuids == {"p1"}
    assert restored._radius == 0.05
    assert restored._last_densify_new == 4
    assert restored._last_densify_attempted is True
    assert restored._accepted_index == 1


def test_adaptive_bridson_resume_reads_runtime_pickle_not_hdf5(tmp_path, monkeypatch):
    import jarvishep2.runtime._resume_service as resume
    from jarvishep2.Sampling.adaptive_bridson import AdaptiveBridsonSampler, LevelSetPoint
    from jarvishep2.runtime.worker_config import build_worker_config
    from jarvishep2.sampling.sample_checkpoint import uses_completion_feedback, uses_sample_checkpoint

    def forbidden(*_args, **_kwargs):
        pytest.fail("AdaptiveBridson resume inspected the HDF5 database")

    monkeypatch.setattr(resume, "read_persisted_sample_index_state", forbidden)
    monkeypatch.setattr(resume, "read_persisted_outcome_counts", forbidden)
    cfg = _adaptive_config(tmp_path)
    assert uses_sample_checkpoint(cfg) is True
    assert uses_completion_feedback(cfg) is False
    blueprint = build_worker_config(
        cfg, task_result_dir=str(tmp_path),
        extra={"publish_feedback": False, "completion_feedback": False},
    )
    assert blueprint["publish_feedback"] is True
    assert blueprint["completion_feedback"] is False

    core = Jarvis2Core(cfg)
    assert core.checkpoint_file().endswith("AdaptiveBridson/runtime.pickle")
    sampler = AdaptiveBridsonSampler()
    sampler.set_config(core.config)
    sampler._points = [LevelSetPoint(u=[0.2, 0.3], x={"x": 0.2, "y": 0.3}, f=0.11, uuid="kept")]
    sampler._uuid_to_index = {"kept": 0}
    sampler._accepted_index = 1
    database = tmp_path / "DATABASE" / "samples.hdf5"
    database.parent.mkdir()
    database.write_bytes(b"not an hdf5 file")
    atomic_pickle(core.checkpoint_file(), {
        "format": "jarvis-hep.sampler-pickle",
        "version": 1,
        "method": "adaptivebridson",
        "config_id": sampler._config_id,
        "integrity": {},
        "sampler_state": sampler.checkpoint_state(),
    })
    core.prepare_resume(resume=True)
    core._load_persisted_database_state(str(database.parent))
    assert database.read_bytes() == b"not an hdf5 file"
    assert core._resume_policy == "resume"
    assert core._persisted_uuids == set()
    assert core._persisted_index_prefix == 0
    restored = AdaptiveBridsonSampler()
    restored.set_config(core.config)
    assert core._import_sampler_checkpoint_state(restored) is True
    assert restored._points[0].uuid == "kept"
    assert restored._points[0].f == 0.11
