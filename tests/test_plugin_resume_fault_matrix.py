"""Joint recovery faults: a sparse archive, failed points and a rejected tail.

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
from jarvishep2.sampling.sample_checkpoint import FORMAT, VERSION, atomic_pickle, config_identity
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
def checkpoint(plan):
    return {}
def resume(setup, saved):
    return None
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
    if getattr(core.sampler, "supports_checkpoint", False):
        record["sampler_state"] = {"index": core.sampler._index,
                                   "accepted_index": core.sampler._accepted_index}
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


def test_resume_counts_come_from_checkpoint_and_ignore_sparse_archive(prepared_unit, monkeypatch):
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.runtime import _resume_service

    config, _ = prepared_unit
    core = Jarvis2Core(config)
    core._resume_checkpoint_payload = {"sampler_state": {"completed_count": 1, "failed_count": 1}}
    monkeypatch.setattr(_resume_service, "read_persisted_sample_index_state",
                        lambda path: pytest.fail("checkpoint resume scanned DATABASE"))
    core._load_persisted_database_state(config["task_result_dir"])
    assert core._persisted_index_prefix == 0
    assert core._persisted_pending_indices == set()
    assert core._persisted_records_count == 2
    assert core._persisted_failed_count == 1


def save_fixture(path, sampler):
    state = sampler.checkpoint_state()
    atomic_pickle(str(path), {
        "format": FORMAT, "version": VERSION, "method": sampler.method.lower(),
        "config_id": config_identity(sampler.config),
        "integrity": sampler._preflight.checkpoint_integrity, "sampler_state": state,
    })


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
    for index in range(prefix + 2):
        point = sampler.propose_next()
        if index < prefix:
            sampler.accept_completion({
                "uuid": point.uuid, "sample_index": index,
                "terminal_status": by_index[index]["status"],
            })
    checkpoint = prepared.source_copy.parent.parent / "runtime.pickle"
    if not database_only:
        save_fixture(checkpoint, sampler)
    assert checkpoint.exists() is not database_only
    binding = database.parent / "sampler_identity.json"
    bound_bytes = binding.read_bytes()
    resumed = run("fault", True)
    resumed_rows = SimpleHDF5Writer(str(database)).read_records()
    reference_calls = _json_lines(root / "reference-physics.jsonl")
    resumed_calls = _json_lines(root / "fault-physics.jsonl")
    assert binding.read_bytes() == bound_bytes
    assert sorted(row["candidate"] for row in reference_calls) == ACCEPTED_CANDIDATES
    if database_only:
        assert resumed["exit_code"] == 1 and "checkpoint missing" in resumed["exception"]
        assert resumed_calls == []
        assert _semantic(resumed_rows) == _semantic([by_index[i] for i in RETAINED_INDICES])
        return
    assert resumed.get("outcome", {}).get("status") == "partial_failure", resumed
    assert reference["exit_code"] == resumed["exit_code"] == 1
    # Only the two consumed completions are skipped. Durable archive rows
    # beyond the snapshot may repeat, including an unfinished failed row.
    expected_missing = ACCEPTED_CANDIDATES[prefix:]
    assert sorted(row["candidate"] for row in resumed_calls) == expected_missing
    assert _semantic(resumed_rows) == _semantic(
        [by_index[i] for i in RETAINED_INDICES] + [by_index[i] for i in range(prefix, 22)]
    )
    assert resumed["sampler_state"]["index"] == POINTS
    assert resumed["sampler_state"]["accepted_index"] == 22
    assert resumed["outcome"]["submitted"] == len(ACCEPTED_CANDIDATES)
    assert resumed["outcome"]["archived"] == len(expected_missing)
    assert resumed["outcome"]["completed"] == 17
    assert resumed["outcome"]["failed"] == 5
    assert [row["block"] for row in _json_lines(root / "fault-draw.jsonl")][-1] == 4
