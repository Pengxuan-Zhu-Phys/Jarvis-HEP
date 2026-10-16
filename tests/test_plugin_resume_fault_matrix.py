"""Joint recovery faults: sparse durable rows, failed points and a rejected tail.

The sparse DATABASE is an explicit crash-state fixture, not a simulated SIGINT.
Both its reference scan and its recovery use real Redis and spawn Workers.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from jarvis_hep_sampler import SamplerRegistry
from jarvishep2.database import SimpleHDF5Writer, read_persisted_sample_index_state
from jarvishep2.sampling.plugin_sampler import PluginSampler, preflight_plugin
from jarvishep2.sampling.runtime_checkpoint import (
    build_payload, build_run_spec, save_checkpoint,
)
from jarvishep2.task_config import load_task_yaml

POINTS = 37
ACCEPTED_CANDIDATES = list(range(4, 26))
RETAINED_INDICES = [0, 1, 3, 6, 8]

SOURCE = '''
import json
import os
import numpy as np
PLUGIN = {"api_version": 1, "method": "FaultMatrix", "version": "1",
          "profile": "independent", "block_size": 8}
def draw_block(context, rng, plan):
    first = context.block_index * context.block_size
    path = os.environ.get("JHS_MATRIX_DRAW_LOG")
    if path:
        with open(path, "a") as handle:
            handle.write(json.dumps({"block": context.block_index}) + "\\n")
    values = np.minimum((first + np.arange(context.block_size) + .5)
                        / context.point_number, 1)
    return np.repeat(values[:, None], context.dimension, axis=1)
'''

PHYSICS = '''
import json
import os
def evaluate(x, y):
    candidate = int(round(float(x) * 37 - .5))
    fd = os.open(os.environ["JHS_MATRIX_PHYSICS_LOG"],
                 os.O_APPEND | os.O_CREAT | os.O_WRONLY, 0o600)
    try:
        os.write(fd, (json.dumps({"candidate": candidate, "pid": os.getpid()})
                      + "\\n").encode())
    finally:
        os.close(fd)
    if candidate % 5 == 0:
        raise RuntimeError("intentional terminal physics failure")
    return {"score": float(x) + float(y)}
'''

DRIVER = '''
import json
import os
import sys
import traceback
from pathlib import Path
from jarvishep2.runtime.core import Jarvis2Core
import jarvishep2
import jarvis_hep_sampler
core = Jarvis2Core()
record = {"python": sys.executable, "pid": os.getpid(),
          "hep_file": jarvishep2.__file__, "sdk_file": jarvis_hep_sampler.__file__}
code = 1
try:
    core.load_task_yaml(sys.argv[1], resume=sys.argv[2] == "resume")
    core.config["skip_draw_flowchart"] = True
    core.config["skip_render_plots"] = True
    outcome = core.run(resume=sys.argv[2] == "resume")
    record["outcome"] = outcome.to_dict()
    code = outcome.exit_code
except Exception as exc:
    record["exception"] = f"{type(exc).__name__}: {exc}"
    traceback.print_exc()
finally:
    if core.sampler is not None:
        record["sampler_state"] = core.sampler.export_runtime_state()
    Path(sys.argv[3]).write_text(json.dumps(record, indent=2) + "\\n")
raise SystemExit(code)
'''


def _semantic(rows):
    return sorted([
        [int(row["sample_index"]), row["uuid"], float(row["x"]), float(row["y"]),
         row["status"], row.get("score")]
        for row in rows
    ])


def _json_lines(path):
    return [json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []


def _card(scan, workers, batch, port):
    return {
        "Scan": {"name": scan},
        "EnvReqs": {"V2": {
            "workers": workers, "batch_size": batch, "store_samples": False,
            "redis": {"host": "127.0.0.1", "port": port, "db": 0},
        }},
        "Sampling": {
            "Method": "FaultMatrix", "Bounds": {"point_number": POINTS, "seed": 7},
            "Variables": [{"name": name, "distribution": {
                "type": "Flat", "parameters": {"min": 0.0, "max": 1.0},
            }} for name in ("x", "y")],
            "selection": "(x > 0.1) & (x < 0.7)",
        },
        "Operas": {"Modules": [{
            "name": "Physics", "operator": "matrix_physics.evaluate",
            "input": [{"name": name, "expression": name} for name in ("x", "y")],
            "output": [{"name": "score", "entry": "score"}],
        }]},
    }


@pytest.fixture
def prepared_unit(tmp_path, monkeypatch):
    root = tmp_path.resolve()
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(root / "registry"))
    monkeypatch.setenv("JARVIS_OPERAS_CURVE_CACHE_ROOT", str(root / "curve-cache"))
    monkeypatch.setenv("JARVIS_OPERAS_PERSIST_FILE", str(root / "operas/legacy.json"))
    monkeypatch.setenv("JARVIS_OPERAS_SOURCES_FILE", str(root / "operas/sources.json"))
    monkeypatch.setenv("JARVIS_OPERAS_OVERRIDES_FILE", str(root / "operas/overrides.json"))
    monkeypatch.setenv("JARVIS_OPERAS_USER_SOURCE_DIR", str(root / "operas/user_sources"))
    source = root / "FaultMatrix.py"
    source.write_text(SOURCE)
    SamplerRegistry().load(source)
    config = _card("unit", 1, 3, 6379)
    config.pop("Operas")
    config.update(task_root=str(root), project_root=str(root), scan_name="unit",
                  task_result_dir=str(root / "outputs/unit"))
    prepared = preflight_plugin(config, task_root=str(root), scan_name="unit")
    return config, prepared


def test_database_reader_retains_pending_terminal_indices(prepared_unit):
    from jarvishep2.runtime.core import Jarvis2Core

    config, _ = prepared_unit
    database = Path(config["task_result_dir"]) / "DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records([
        {"sample_index": index, "uuid": str(index),
         "status": "Failed" if index in {1, 6} else "Completed"}
        for index in RETAINED_INDICES
    ])
    core = Jarvis2Core(config)
    core._load_persisted_database_state(str(database.parent))
    assert core._persisted_index_prefix == 2
    assert getattr(core, "_persisted_pending_indices", set()) == {3, 6, 8}
    assert core._persisted_records_count == 5
    assert core._persisted_failed_count == 2


@pytest.mark.parametrize("workers,batch", [(1, 3), (2, 5)])
def test_submission_skips_only_durable_indices_and_counts_pushes(prepared_unit, workers, batch):
    from jarvishep2.redis_queue import make_fakeredis_queue

    config, prepared = prepared_unit
    config["EnvReqs"]["V2"].update(workers=workers, batch_size=batch)
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    sampler.assert_checkpoint_attribute_contract()
    sampler.set_resume_mode(True)
    sampler.set_persisted_prefix(2)
    sampler.advance_to_persisted_prefix(2)
    # Feed the DATABASE-derived out-of-order facts. The old implementation
    # ignores these facts and re-submits three durable terminal rows.
    sampler._persisted_pending_indices = {3, 6, 8}
    queue = make_fakeredis_queue()
    sampler.set_redis(queue)
    sampler._max_inflight = 64
    try:
        pushed = sampler.run_distributed()
        assert pushed == 17
        submitted = [queue.pull_task(timeout=1) for _ in range(pushed)]
        assert [task["sample_index"] for task in submitted] == [
            index for index in range(22) if index not in RETAINED_INDICES
        ]
        assert sampler._index == 37 and sampler._accepted_index == 22
        assert not {"pending_indices", "persisted_pending_indices", "_persisted_pending_indices"}.intersection(
            sampler.export_runtime_state()
        )
    finally:
        queue.close()


def test_whole_durable_batch_still_captures_checkpoint(prepared_unit):
    from jarvishep2.redis_queue import make_fakeredis_queue

    config, prepared = prepared_unit
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    sampler.set_resume_mode(True)
    sampler.set_persisted_prefix(2)
    sampler.advance_to_persisted_prefix(2)
    sampler.record_persisted_proposals()
    sampler.set_persisted_pending_indices({3, 4, 5})
    queue = make_fakeredis_queue()
    sampler.set_redis(queue)
    try:
        # The prefix hole is submitted normally; only the following whole
        # batch is already durable, as a DATABASE reader can actually report.
        assert sampler.flush_batch([sampler.propose_next()]) == 1
        assert queue.pull_task(timeout=1)["sample_index"] == 2
        batch = [sampler.propose_next() for _ in range(3)]
        sampler.request_checkpoint_capture()
        assert sampler.flush_batch(batch) == 0
        assert not sampler._checkpoint_requested
        assert sampler._checkpoint_snapshots[-1][0] == 6
        assert queue.get_queue_lengths()["task_queue_length"] == 0
        assert sampler._accepted_index == 6
        assert sampler.checkpoint_runtime_state(safe=False)["accepted_index"] == 2
        sampler.set_persisted_prefix(6)
        assert sampler.checkpoint_runtime_state(safe=False)["accepted_index"] == 6
    finally:
        queue.close()


def test_pending_indices_are_copied_and_replaced(prepared_unit):
    config, prepared = prepared_unit
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    pending = {3, 6, 8}
    sampler.set_persisted_pending_indices(pending)
    pending.clear()
    assert sampler._persisted_pending_indices == {3, 6, 8}
    baseline_state = json.dumps(sampler.export_runtime_state(), sort_keys=True)
    sampler.set_persisted_pending_indices(set(range(10000)))
    assert json.dumps(sampler.export_runtime_state(), sort_keys=True) == baseline_state
    sampler.assert_checkpoint_attribute_contract()
    sampler.set_persisted_pending_indices(set())
    assert sampler._persisted_pending_indices == set()


def test_builtin_random_holes_keep_single_slot_backpressure(prepared_unit, monkeypatch):
    from jarvishep2.sampling.randoms import RandomS

    config, _ = prepared_unit
    config["Sampling"]["Method"] = "Random"
    config["Sampling"]["Bounds"]["point_number"] = 8
    config["Sampling"].pop("selection")
    sampler = RandomS()
    sampler.set_config(config)
    sampler.assert_checkpoint_attribute_contract()
    sampler.set_resume_mode(True)
    sampler.set_persisted_prefix(2)
    sampler.advance_to_persisted_prefix(2)
    sampler.set_persisted_pending_indices({3, 5})
    submitted = []
    monkeypatch.setattr(sampler, "_submit", submitted.append)
    monkeypatch.setattr(sampler, "_finished_count", lambda: len(submitted))
    sampler._max_inflight = 1
    original_wait = sampler._wait_for_inflight_room
    def checked_wait(*args, **kwargs):
        pushed = kwargs.get("pushed", args[0] if args else None)
        assert pushed == len(submitted)
        return original_wait(*args, **kwargs)
    monkeypatch.setattr(sampler, "_wait_for_inflight_room", checked_wait)
    assert sampler.run_distributed() == 4
    assert [sample.sample_index for sample in submitted] == [2, 4, 6, 7]
    assert sampler._accepted_index == sampler._index == 8


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required")
@pytest.mark.parametrize("workers,batch", [(1, 3), (2, 5)], ids=["w1-b3", "w2-b5"])
@pytest.mark.parametrize("database_only", [False, True], ids=["checkpoint", "database-only"])
def test_sparse_failed_batch_and_rejected_tail(tmp_path, monkeypatch, workers, batch, database_only):
    root = tmp_path.resolve()
    scan_token = hashlib.sha256(str(root).encode()).hexdigest()[:12]
    scan_names = {name: f"matrix_{scan_token}_{name}" for name in ("reference", "fault")}
    monkeypatch.setenv("JARVIS_SCAN_REF_REGISTRY", str(root / "scan-refs.json"))
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(root / "registry"))
    monkeypatch.setenv("JARVIS_OPERAS_CURVE_CACHE_ROOT", str(root / "curve-cache"))
    monkeypatch.setenv("JARVIS_OPERAS_PERSIST_FILE", str(root / "operas/legacy.json"))
    monkeypatch.setenv("JARVIS_OPERAS_SOURCES_FILE", str(root / "operas/sources.json"))
    monkeypatch.setenv("JARVIS_OPERAS_OVERRIDES_FILE", str(root / "operas/overrides.json"))
    monkeypatch.setenv("JARVIS_OPERAS_USER_SOURCE_DIR", str(root / "operas/user_sources"))
    (root / ".jarvis-project.json").write_text("{}")
    (root / "matrix_physics.py").write_text(PHYSICS)
    source = root / "FaultMatrix.py"
    source.write_text(SOURCE)
    SamplerRegistry().load(source)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    env.pop("PYTHONPATH", None)
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]

    def run(scan, resume):
        card = root / f"{scan}.yaml"
        card.write_text(yaml.safe_dump(_card(scan_names[scan], workers, batch, port)))
        output = root / f"{scan}-outcome.json"
        phase_env = dict(env, JHS_MATRIX_PHYSICS_LOG=str(root / f"{scan}-physics.jsonl"),
                         JHS_MATRIX_DRAW_LOG=str(root / f"{scan}-draw.jsonl"))
        process = subprocess.run(
            [sys.executable, "-c", DRIVER, str(card), "resume" if resume else "fresh", str(output)],
            cwd=root, env=phase_env, text=True, capture_output=True, timeout=60,
        )
        (root / f"{scan}.log").write_text(process.stdout + process.stderr)
        record = json.loads(output.read_text())
        record["exit_code"] = process.returncode
        return record

    reference = run("reference", False)
    assert reference.get("outcome", {}).get("status") == "partial_failure", reference
    reference_rows = SimpleHDF5Writer(str(root / f"outputs/{scan_names['reference']}/DATABASE/samples.hdf5")).read_records()
    assert len(reference_rows) == len(ACCEPTED_CANDIDATES) == 22
    by_index = {int(row["sample_index"]): row for row in reference_rows}
    assert sum(row["status"] == "Failed" for row in reference_rows) == 5
    # The first runtime batch mixes success with terminal failure; failure is
    # not a hole. Keep failed rows at accepted indices 1 and 6 in the fixture.
    assert {by_index[i]["status"] for i in range(batch)} == {"Completed", "Failed"}
    assert by_index[1]["status"] == by_index[6]["status"] == "Failed"
    task = root / "fault.yaml"
    task.write_text(yaml.safe_dump(_card(scan_names["fault"], workers, batch, port)))
    config = load_task_yaml(str(task))
    prepared = preflight_plugin(config, task_root=str(root), scan_name=scan_names["fault"])
    database = root / f"outputs/{scan_names['fault']}/DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records([by_index[i] for i in RETAINED_INDICES])
    prefix, pending, legacy = read_persisted_sample_index_state(str(database.parent))
    assert (prefix, pending, legacy) == (2, {3, 6, 8}, set())
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    for _ in range(prefix):
        assert sampler.propose_next() is not None
    checkpoint_state = sampler.export_runtime_state()
    assert checkpoint_state["accepted_index"] <= prefix
    checkpoint = prepared.source_copy.parent.parent / "state.pkl"
    if not database_only:
        run_spec = build_run_spec(config=config, scan_name=scan_names["fault"], task_root=str(root),
                                  task_result_dir=str(database.parent.parent),
                                  sampler_name="FaultMatrix", worker_parallel=workers)
        save_checkpoint(str(checkpoint), build_payload(
            run_spec=run_spec, sampler_state=checkpoint_state,
            persistence={"persisted_prefix": prefix}, reason="joint_fault_fixture",
            safe_barrier_confirmed=True,
        ))
    assert checkpoint.exists() is not database_only
    binding = database.parent / "sampler_identity.json"
    bound_bytes = binding.read_bytes()
    resumed = run("fault", True)
    resumed_rows = SimpleHDF5Writer(str(database)).read_records()
    reference_calls = _json_lines(root / "reference-physics.jsonl")
    resumed_calls = _json_lines(root / "fault-physics.jsonl")
    evidence = {
        "workers": workers, "batch_size": batch, "database_only": database_only,
        "scan_names": scan_names,
        "fixture": {"kind": "injected sparse durable crash state", "prefix": prefix,
                    "pending": sorted(pending), "retained_indices": RETAINED_INDICES,
                    "retained_failed_indices": [1, 6], "checkpoint": checkpoint_state if not database_only else None},
        "reference": reference, "resumed": resumed,
        "reference_records": _semantic(reference_rows), "resumed_records": _semantic(resumed_rows),
        "reference_physics_calls": reference_calls, "resumed_physics_calls": resumed_calls,
        "resumed_draw_blocks": _json_lines(root / "fault-draw.jsonl"),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    }
    (root / "fault-matrix-evidence.json").write_text(json.dumps(evidence, indent=2) + "\n")
    assert resumed.get("outcome", {}).get("status") == "partial_failure", resumed
    assert reference["exit_code"] == resumed["exit_code"] == 1
    assert binding.read_bytes() == bound_bytes
    assert _semantic(resumed_rows) == _semantic(reference_rows)
    assert len({row["uuid"] for row in resumed_rows}) == 22
    assert [row[0] for row in _semantic(resumed_rows)] == list(range(22))
    for accepted, candidate in enumerate(ACCEPTED_CANDIDATES):
        assert by_index[accepted]["x"] == by_index[accepted]["y"] == (candidate + .5) / POINTS
    assert sorted(row["candidate"] for row in reference_calls) == ACCEPTED_CANDIDATES
    expected_missing = [candidate for index, candidate in enumerate(ACCEPTED_CANDIDATES)
                        if index not in RETAINED_INDICES]
    assert sorted(row["candidate"] for row in resumed_calls) == expected_missing
    # Rejected prefix and last 11 candidates never become Worker tasks.
    assert not any(row["candidate"] >= 26 for row in reference_calls + resumed_calls)
    assert resumed["sampler_state"]["index"] == POINTS
    assert resumed["sampler_state"]["accepted_index"] == 22
    assert resumed["outcome"]["submitted"] == len(expected_missing) == 17
    assert resumed["outcome"]["completed"] == 17
    assert resumed["outcome"]["failed"] == 5
    assert resumed["outcome"]["archived"] == 22
    assert [row["block"] for row in evidence["resumed_draw_blocks"]][-1] == 4
