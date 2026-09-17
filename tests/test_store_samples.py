"""Storage opt-out preserves physics and durable results without SAMPLE artifacts."""

import json
import os
from pathlib import Path
from unittest import mock

import pytest

from jarvishep2.command_parser import CommandParser
from jarvishep2.core import Jarvis2Core
from jarvishep2.database import SimpleHDF5Writer, StreamingHDF5Writer
from jarvishep2.io.archiver import ArchiveProcessor
from jarvishep2.io.database import read_persisted_sample_index_state
from jarvishep2.io_portal import apply_hep_io_save
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.runtime.worker_config import build_worker_config
from jarvishep2.runtime_config import (
    get_archiver_config, get_runtime_block, get_sample_directory_config,
)
from jarvishep2.sample import Sample, materialize_failure_artifacts
from jarvishep2.sample_logger import NullSampleLogger
from jarvishep2.task_config import load_task_yaml
from jarvishep2.task_validation import validate_task_config
from jarvishep2.worker import Worker
from test_task_validation import _minimal_dynesty_config


def config(store=False):
    return {"EnvReqs": {"V2": {
        "store_samples": store,
        "sample_directory": {"enabled": True, "pack": True},
        "archiver": {"pack_buckets": True},
        "worker": {"sample_artifacts": "always"},
    }}}


def test_default_and_override_policy(tmp_path):
    assert get_runtime_block({})["store_samples"] is True
    cfg = config()
    cfg["Runtime"] = {"store_samples": True}
    blueprint = build_worker_config(cfg, task_result_dir=str(tmp_path))
    assert blueprint["sample_config"]["store_samples"] is False
    assert not get_sample_directory_config(cfg)["enabled"]
    assert not get_sample_directory_config(cfg)["pack"]
    assert not get_archiver_config(cfg)["pack_buckets"]
    assert get_archiver_config(config(True))["pack_buckets"]


@pytest.mark.parametrize("value", [False, True, "false", 0, None])
def test_strict_boolean_validation(value):
    cfg = _minimal_dynesty_config()
    cfg["EnvReqs"]["V2"]["store_samples"] = value
    if isinstance(value, bool):
        assert validate_task_config(cfg).ok
    else:
        assert not validate_task_config(cfg).ok
        with pytest.raises(ValueError, match="store_samples"):
            get_runtime_block(cfg)


def test_yaml_default_inheritance_and_task_override(tmp_path):
    (tmp_path / "jarvis.project.yaml").write_text("project: lightweight\n")
    (tmp_path / "deps").mkdir()
    defaults = tmp_path / "deps" / "environment_default.yaml"
    defaults.write_text("EnvReqs:\n  V2:\n    store_samples: false\n")
    task = tmp_path / "task.yaml"
    text = ("Scan:\n  name: light\nSampling:\n  Method: Random\n"
            "EnvReqs:\n  Check_default_dependencies:\n    required: true\n"
            "    default_yaml_path: '&J/deps/environment_default.yaml'\n")
    task.write_text(text)
    assert not get_runtime_block(load_task_yaml(str(task)))["store_samples"]
    task.write_text(text + "  V2:\n    store_samples: true\n")
    assert get_runtime_block(load_task_yaml(str(task)))["store_samples"]


def test_lazy_sdir_is_shared_silent_and_disposable(tmp_path):
    sample = Sample.from_params({"x": 1})
    sample.set_config({"store_samples": False, "task_result_dir": str(tmp_path),
                       "sample_artifacts": "never"})
    assert sample.save_dir is None
    logger = sample.info["logger"]
    assert isinstance(logger, NullSampleLogger)
    with mock.patch("jarvishep2.sample_logger._forward_sample_event") as forward:
        logger.error("ignored")
        logger.bind(module="child").error("also ignored")
        assert logger.event_count == 0
        forward.assert_not_called()
    parser = CommandParser(project_root=str(tmp_path))
    try:
        path = parser.resolve_sample("@Sdir/input", sample_info=sample.info, stage="execution")
        assert path == os.path.join(sample.save_dir, "input")
        assert parser.resolve_sample("@Sdir/input", sample_info=sample.info, stage="execution") == path
        assert not (Path(sample.save_dir) / "Sample_running.log").exists()
        assert materialize_failure_artifacts(sample.info, error="failure") is None
        assert "save_dir" not in sample.to_info_dict()
        scratch = sample.save_dir
    finally:
        sample.close()
        sample.cleanup_transient_directory()
    assert not os.path.exists(scratch)
    assert not (tmp_path / "SAMPLE").exists()


@pytest.mark.parametrize("failure", [None, "compute", "submit"])
def test_worker_cleans_scratch_on_success_and_failures(tmp_path, failure):
    queue = make_fakeredis_queue()
    blueprint = build_worker_config(config(), task_result_dir=str(tmp_path))
    worker = Worker(0, {}, blueprint)
    worker._redis = queue
    scratch_dirs = []

    def execute(sample):
        scratch = sample.materialize()
        scratch_dirs.append(scratch)
        Path(scratch, "needed-input").write_text("working data")
        sample.merge_observables({"answer": 42})
        if failure == "compute":
            raise RuntimeError("deliberate failure")
        sample.set_status("Completed")

    with mock.patch.object(worker._get_executor(), "process", side_effect=execute):
        if failure == "submit":
            with mock.patch.object(queue, "submit_result", side_effect=RuntimeError("unavailable")):
                worker.process_task({"uuid": "s", "u_coords": [], "execution_plan": []})
        else:
            worker.process_task({"uuid": "s", "u_coords": [], "execution_plan": []})
            result = queue.pull_result(timeout=1)
            assert result["observables"]["answer"] == 42
            assert result["status"] == ("Failed" if failure else "Completed")
            assert "save_dir" not in result
            assert "bucket_id" not in result
    assert scratch_dirs and all(not os.path.exists(path) for path in scratch_dirs)
    assert not (tmp_path / "SAMPLE").exists()


def test_numeric_worker_never_allocates_directory(tmp_path):
    worker = Worker(0, {}, build_worker_config(config(), task_result_dir=str(tmp_path)))
    worker._redis = make_fakeredis_queue()
    with mock.patch("jarvishep2.sample.tempfile.mkdtemp") as mkdir:
        worker.process_task({"uuid": "numeric", "u_coords": [], "execution_plan": []})
        mkdir.assert_not_called()
    assert worker._redis.pull_result(timeout=1)["status"] == "Completed"
    assert not (tmp_path / "SAMPLE").exists()


def test_database_retains_failure_and_resume_prefix_without_sample_tree(tmp_path):
    db = tmp_path / "DATABASE" / "samples.hdf5"
    writer = StreamingHDF5Writer(str(db))
    processor = ArchiveProcessor.from_config(
        writer, sample_root=str(tmp_path / "SAMPLE"), archiver_config=get_archiver_config(config()),
    )
    try:
        for index, status in enumerate(["Completed", "Failed"]):
            processor.ingest({"uuid": f"u{index}", "sample_index": index,
                              "status": status, "observables": {"x": index}})
        processor.flush_batch(force=True)
        assert processor.persistence_state()["persisted_prefix"] == 2
    finally:
        writer.close()
    assert [r["status"] for r in SimpleHDF5Writer(str(db)).read_records()] == ["Completed", "Failed"]
    assert read_persisted_sample_index_state(str(db.parent))[0] == 2
    assert not (tmp_path / "SAMPLE").exists()


@pytest.mark.parametrize("save", [True, False])
def test_output_copy_survives_pack_reuse_only_until_sample_finishes(tmp_path, save):
    sample = Sample.from_params({"x": 1})
    sample.set_config({"store_samples": False, "task_result_dir": str(tmp_path)})
    source = tmp_path / "calculator-output.txt"
    source.write_text("first result")
    # Existing retained results must never be removed by the opt-out.
    retained = tmp_path / "SAMPLE" / "previous-run.txt"
    retained.parent.mkdir()
    retained.write_text("keep")
    try:
        scratch = sample.materialize()
        copied = apply_hep_io_save(source_path=str(source), sample_info=sample.info,
                                   module="Calculator", spec={"save": save}, direction="output")
        assert Path(copied).is_relative_to(Path(scratch).resolve())
        source.write_text("next sample reused calculator pack")
        assert Path(copied).read_text() == "first result"
        assert not list(Path(scratch).rglob("*.log"))
    finally:
        sample.close()
        sample.cleanup_transient_directory()
    assert not Path(copied).exists()
    assert retained.read_text() == "keep"


def test_real_calculator_check_keeps_physics_without_sample_tree(tmp_path):
    from test_core_run_distributed import CHECK_MODULES_YAML, FIXTURES, _stop_factory_workers
    from test_worker_calculator import _start_tcp_fakeredis, _normalize_database_records

    _stop_factory_workers()
    server, redis_config = _start_tcp_fakeredis()
    core = Jarvis2Core()
    try:
        core.load_task_yaml(CHECK_MODULES_YAML, check_modules=True)
        core.config["EnvReqs"].setdefault("V2", {})["store_samples"] = False
        core.config["task_result_dir"] = str(tmp_path)
        core.config["Runtime"]["redis"] = redis_config
        core.runtime = core.config["Runtime"]
        core._ensure_managed_redis = mock.Mock()
        core._populate_info_from_config()
        count = core.run(check_modules=True, write_run_summary=False)
        assert count == 10
        rows = SimpleHDF5Writer(str(tmp_path / "DATABASE/test/samples.hdf5")).read_records()
        assert len(rows) == 10
        assert all(row["status"] == "Completed" and "z" in row and "LogL" in row for row in rows)
        expected = json.loads(Path(FIXTURES, "expected_calculator_records.json").read_text())
        assert _normalize_database_records(rows) == _normalize_database_records(expected)
        assert not (tmp_path / "SAMPLE").exists()
        assert not list(tmp_path.rglob("Sample_running.log"))
    finally:
        _stop_factory_workers()
        server.shutdown()
        server.server_close()
