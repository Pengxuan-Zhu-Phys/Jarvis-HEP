"""A fully durable accepted prefix resumes through rejection-only candidates."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys

import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry

from jarvishep2.database import SimpleHDF5Writer, read_persisted_sample_index_state
from jarvishep2.sampling.plugin_sampler import PluginSampler, preflight_plugin
from jarvishep2.task_config import load_task_yaml

from test_plugin_resume_fault_matrix import (
    ACCEPTED_CANDIDATES, DRIVER, PHYSICS, POINTS, SOURCE, _card, _json_lines, _semantic, save_fixture,
)

TAIL_DRIVER = '''
from jarvishep2.runtime.core import Jarvis2Core
_run_distributed_scan = Jarvis2Core.run_distributed_scan
def observe_new_submissions(self, *args, **kwargs):
    count = _run_distributed_scan(self, *args, **kwargs)
    Path(sys.argv[3] + ".submissions").write_text(str(count))
    return count
Jarvis2Core.run_distributed_scan = observe_new_submissions
''' + DRIVER


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required")
@pytest.mark.parametrize("database_only", [False, True], ids=["checkpoint", "database-only"])
def test_all_accepted_durable_resume_evaluates_only_rejected_tail(tmp_path, monkeypatch, database_only):
    root = tmp_path.resolve()
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(root / "registry"))
    monkeypatch.setenv("JARVIS_SCAN_REF_REGISTRY", str(root / "scan-refs.json"))
    monkeypatch.setenv("JARVIS_OPERAS_CURVE_CACHE_ROOT", str(root / "curve-cache"))
    (root / ".jarvis-project.json").write_text("{}", encoding="utf-8")
    (root / "matrix_physics.py").write_text(PHYSICS, encoding="utf-8")
    source = root / "FaultMatrix.py"
    source.write_text(SOURCE, encoding="utf-8")
    SamplerRegistry().load(source)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    workers, batch = 2, 5
    tag = hashlib.sha256(str(root).encode()).hexdigest()[:12]
    reference_scan, resumed_scan = f"tail_{tag}_reference", f"tail_{tag}_resumed"

    def run(scan, resume):
        task = root / f"{scan}.yaml"
        task.write_text(yaml.safe_dump(_card(scan, workers, batch, port)), encoding="utf-8")
        outcome = root / f"{scan}-outcome.json"
        phase_env = dict(
            env, JHS_MATRIX_PHYSICS_LOG=str(root / f"{scan}-physics.jsonl"),
            JHS_MATRIX_DRAW_LOG=str(root / f"{scan}-draw.jsonl"),
        )
        process = subprocess.run(
            [sys.executable, "-c", TAIL_DRIVER, str(task), "resume" if resume else "fresh", str(outcome)],
            cwd=root, env=phase_env, text=True, capture_output=True, timeout=60,
        )
        (root / f"{scan}.log").write_text(process.stdout + process.stderr, encoding="utf-8")
        with socket.socket() as probe:
            assert probe.connect_ex(("127.0.0.1", port)) != 0
        result = json.loads(outcome.read_text(encoding="utf-8"))
        result["exit_code"] = process.returncode
        observed = root / f"{scan}-outcome.json.submissions"
        result["actual_new_submitted"] = int(observed.read_text()) if observed.exists() else None
        return result

    reference = run(reference_scan, False)
    assert reference.get("outcome", {}).get("status") == "partial_failure", reference
    reference_database = root / "outputs" / reference_scan / "DATABASE/samples.hdf5"
    reference_rows = SimpleHDF5Writer(str(reference_database)).read_records()
    assert len(reference_rows) == len(ACCEPTED_CANDIDATES) == 22
    assert sum(row["status"] == "Failed" for row in reference_rows) == 5
    assert sorted(row["candidate"] for row in _json_lines(root / f"{reference_scan}-physics.jsonl")) == (
        ACCEPTED_CANDIDATES
    )

    task = root / f"{resumed_scan}.yaml"
    task.write_text(yaml.safe_dump(_card(resumed_scan, workers, batch, port)), encoding="utf-8")
    config = load_task_yaml(str(task))
    prepared = preflight_plugin(config, task_root=str(root), scan_name=resumed_scan)
    database = root / "outputs" / resumed_scan / "DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records(reference_rows)
    assert read_persisted_sample_index_state(str(database.parent)) == (22, set(), set())
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    by_index = {int(row["sample_index"]): row for row in reference_rows}
    for index in range(22):
        point = sampler.propose_next()
        sampler.accept_completion({
            "uuid": point.uuid, "sample_index": index, "terminal_status": by_index[index]["status"],
        })
    assert sampler._index == 26
    assert sampler._accepted_index == 22
    checkpoint = prepared.source_copy.parent.parent / "runtime.pickle"
    if not database_only:
        save_fixture(checkpoint, sampler)
    assert checkpoint.exists() is not database_only
    binding = database.parent / "sampler_identity.json"
    binding_before = binding.read_bytes()
    rows_before = _semantic(SimpleHDF5Writer(str(database)).read_records())
    resumed = run(resumed_scan, True)
    rows_after = SimpleHDF5Writer(str(database)).read_records()
    calls = _json_lines(root / f"{resumed_scan}-physics.jsonl")
    blocks = [row["block"] for row in _json_lines(root / f"{resumed_scan}-draw.jsonl")]
    assert calls == []
    assert len(rows_after) == 22 and _semantic(rows_after) == rows_before
    assert binding.read_bytes() == binding_before
    if database_only:
        assert resumed["exit_code"] == 1 and "checkpoint missing" in resumed["exception"]
        assert blocks == []
        return
    assert resumed.get("outcome", {}).get("status") == "partial_failure", resumed
    assert reference["exit_code"] == resumed["exit_code"] == 1
    assert resumed["actual_new_submitted"] == 0
    assert resumed["outcome"]["completed"] == 17
    assert resumed["outcome"]["failed"] == 5
    assert resumed["sampler_state"]["index"] == POINTS == 37
    assert resumed["sampler_state"]["accepted_index"] == 22
    assert blocks == [3, 4]
