"""Checkpoint authoritative resume protocol, including real runtime interruption."""
from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor
from jarvishep2.sampling.runtime_checkpoint import load_checkpoint
from jarvishep2.task_config import load_task_yaml


DRIVER = r'''
import json, sys, time
from pathlib import Path
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.distributor import Distributor
from jarvishep2.sampling.checkpointed_sampler import CheckpointHeartbeat
from jarvishep2.sampling.runtime_checkpoint import load_checkpoint

card, output, phase, timing = sys.argv[1:]
core = Jarvis2Core()
core.load_task_yaml(card)
core.config.update(skip_draw_flowchart=True, skip_render_plots=True)
if phase == "resume" and core.config["Sampling"]["Method"] == "Bridson":
    from jarvishep2.sampling.bridson import Bridson
    def no_generation(self, **kwargs):
        raise AssertionError("Bridson regenerated points during resume")
    Bridson.initialize = no_generation
if phase == "interrupt":
    klass = type(Distributor.set_method(core.config["Sampling"]["Method"]))
    original = klass.propose_next
    original_heartbeat = CheckpointHeartbeat.__init__
    def fast_heartbeat(self, *, interval_sec, save_callback):
        original_heartbeat(self, interval_sec=1.0, save_callback=save_callback)
    CheckpointHeartbeat.__init__ = fast_heartbeat
    paused = False
    def propose(self):
        global paused
        accepted = int(getattr(self, "_accepted_index", 0))
        if timing == "after" and accepted >= 2 and not paused:
            paused = True
            time.sleep(1.3)
        if accepted >= (5 if timing == "after" else 2):
            if timing == "after":
                deadline = time.monotonic() + 6
                while time.monotonic() < deadline:
                    state = load_checkpoint(core.checkpoint_file())
                    if int(state["sampler_state"].get("checkpoint_index", 0)) > 0:
                        break
                    self.request_checkpoint_capture()
                    self.capture_checkpoint_at_batch_boundary()
                    time.sleep(0.1)
                else:
                    raise AssertionError("heartbeat did not advance durable checkpoint")
            raise KeyboardInterrupt
        return original(self)
    klass.propose_next = propose
outcome = core.run(resume=phase == "resume", write_run_summary=False)
Path(output).write_text(json.dumps(outcome.to_dict()))
raise SystemExit(outcome.exit_code)
'''


def _card(tmp_path: Path, method: str, *, count: int = 16) -> tuple[Path, dict]:
    card = yaml.safe_load((Path(__file__).parent / "parity_project/random_opera.yaml").read_text())
    card["Scan"]["name"] = f"protocol_{method}_{tmp_path.name}"
    card["Sampling"]["Method"] = method
    card["Sampling"]["Bounds"] = {
        "Grid": {"seed": 11},
        "Random": {"point_number": count, "seed": 11},
        "Bridson": {"radius": 0.22, "max_attempt": 20, "seed": 11},
        "CSV": {},
    }[method]
    for variable in card["Sampling"]["Variables"]:
        variable["distribution"]["parameters"]["num"] = 4
    if method == "CSV":
        csv = tmp_path / "points.csv"
        csv.write_text("uuid,x,y\n" + "".join(f"row-{i},{(i + 1) / (count + 1):.6f},0.3\n" for i in range(count)))
        card["Sampling"]["Bounds"] = {"path": str(csv)}
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    card["EnvReqs"]["V2"]["redis"] = {"host": "127.0.0.1", "port": port, "db": 0}
    card["EnvReqs"]["V2"]["batch_size"] = 1
    path = tmp_path / "task.yaml"
    path.write_text(yaml.safe_dump(card))
    (tmp_path / ".jarvis-project.json").write_text("{}")
    return path, card


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("method", ["Grid", "Random", "Bridson", "CSV"])
@pytest.mark.parametrize("timing", ["before", "after"])
def test_fresh_interrupt_resume_replaces_old_rows(tmp_path, monkeypatch, method, timing):
    source, card = _card(tmp_path, method)
    config = load_task_yaml(str(source))
    sampler = Distributor.set_method(method)
    sampler.set_config(config)
    stale = []
    while (sample := sampler.propose_next()) is not None:
        stale.append({"uuid": sample.uuid, "sample_index": sample.sample_index,
                      "status": "Completed", "x": 0., "y": 0., "z": -999.})
    assert len(stale) > 5
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records(stale + [
        {"uuid": "old-tail", "sample_index": len(stale), "status": "Completed", "z": -999.}
    ])
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    for phase in ("interrupt", "resume"):
        output = tmp_path / f"{phase}.json"
        result = subprocess.run([sys.executable, "-c", DRIVER, str(source), str(output), phase, timing],
                                cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
        (tmp_path / f"{phase}.log").write_text(result.stdout + result.stderr)
        assert output.exists(), result.stdout + result.stderr
        if phase == "interrupt":
            assert result.returncode == 130, result.stdout + result.stderr
            checkpoint = tmp_path / f"checkpoints/{card['Scan']['name']}/{method}/state.pkl"
            payload = load_checkpoint(str(checkpoint))
            assert payload["persistence"]["replacing_existing"] is True
            assert (int(payload["sampler_state"].get("checkpoint_index", 0)) > 0) == (timing == "after")
        else:
            assert result.returncode == 0, result.stdout + result.stderr
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert len(rows) == len(stale)
    assert {int(row["sample_index"]) for row in rows} == set(range(len(stale)))
    assert all(float(row["z"]) != -999 for row in rows)


def test_bridson_million_points_resume_uses_mmap_without_generation(tmp_path, monkeypatch):
    from jarvishep2.sampling.bridson import Bridson
    from jarvishep2.sampling.runtime_checkpoint import array_fingerprint

    sampler = Bridson()
    sampler.vars = [object(), object()]
    points = np.lib.format.open_memmap(tmp_path / "points.npy", mode="w+", dtype="float64", shape=(1_000_000, 2))
    points[:, 0] = np.arange(len(points)) / len(points)
    points[:, 1] = 0.5
    digest = array_fingerprint(points)
    del points
    checkpoint = tmp_path / "state.pkl"
    monkeypatch.setattr(sampler, "initialize", lambda **kwargs: pytest.fail("Bridson regenerated points"))
    sampler.prepare_resume_points(str(checkpoint), {"points_sha256": digest})
    assert isinstance(sampler._P, np.memmap)
    assert len(sampler._P) == 1_000_000


def test_bridson_point_file_missing_damaged_or_mismatched(tmp_path):
    from jarvishep2.sampling.bridson import Bridson
    from jarvishep2.sampling.runtime_checkpoint import array_fingerprint

    sampler = Bridson()
    sampler.vars = [object(), object()]
    checkpoint = tmp_path / "state.pkl"
    with pytest.raises(ValueError, match="missing"):
        sampler.prepare_resume_points(str(checkpoint), {"points_sha256": "x"})
    path = tmp_path / "points.npy"
    path.write_bytes(b"damaged")
    with pytest.raises(ValueError, match="damaged"):
        sampler.prepare_resume_points(str(checkpoint), {"points_sha256": "x"})
    points = np.arange(20, dtype=np.float64).reshape(10, 2)
    np.save(path, points)
    with pytest.raises(ValueError, match="hash mismatch"):
        sampler.prepare_resume_points(str(checkpoint), {"points_sha256": "x"})
    sampler.prepare_resume_points(str(checkpoint), {"points_sha256": array_fingerprint(points)})
    assert isinstance(sampler._P, np.memmap)


def test_bridson_points_file_written_once_across_cursor_snapshots(tmp_path):
    from jarvishep2.sampling.bridson import Bridson

    config = {"Sampling": {"Method": "Bridson", "Bounds": {"seed": 11, "radius": 0.4, "max_attempt": 20},
              "Variables": [{"name": name, "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1, "length": 1}}}
                            for name in ("x", "y")]}}
    sampler = Bridson(); sampler.set_config(config)
    checkpoint = tmp_path / "state.pkl"
    sampler.prepare_fresh_points(str(checkpoint))
    path = tmp_path / "points.npy"
    stamp = path.stat().st_mtime_ns
    digest = sampler.export_runtime_state()["points_sha256"]
    for _ in range(3):
        assert sampler.export_runtime_state()["points_sha256"] == digest
    assert path.stat().st_mtime_ns == stamp
    assert len(list(tmp_path.glob("points.*.tmp"))) == 0


@pytest.mark.parametrize("damage", ["missing", "malformed", "hash"])
def test_bridson_incomplete_checkpoint_falls_back_to_fresh(tmp_path, damage):
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.sampling.bridson import Bridson
    from jarvishep2.sampling.runtime_checkpoint import build_payload, build_run_spec, save_checkpoint

    config = {"task_root": str(tmp_path), "scan_name": "bridson-damage",
              "Sampling": {"Method": "Bridson", "Bounds": {"seed": 11, "radius": 0.4, "max_attempt": 20},
                  "Variables": [{"name": name, "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1, "length": 1}}}
                                for name in ("x", "y")]},
              "EnvReqs": {"V2": {"checkpoint": {"enabled": True}}}}
    core = Jarvis2Core(config)
    original = Bridson(); original.set_config(core.config)
    original.prepare_fresh_points(core.checkpoint_file())
    spec = build_run_spec(config=core.config, scan_name="bridson-damage", task_root=str(tmp_path),
                          task_result_dir=str(tmp_path), sampler_name="Bridson")
    save_checkpoint(core.checkpoint_file(), build_payload(run_spec=spec, sampler_state=original.export_runtime_state()))
    points = Path(core.checkpoint_file()).with_name("points.npy")
    if damage == "missing":
        points.unlink()
    elif damage == "malformed":
        points.write_bytes(b"broken")
    else:
        with points.open("r+b") as handle:
            handle.seek(-8, os.SEEK_END)
            handle.write(b"tampered")
    core.prepare_resume(resume=True)
    resumed = Bridson(); resumed.set_config(core.config)
    core.set_sampler(resumed)
    assert core._resume_policy == "fresh"


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
def test_check_mode_writes_no_checkpoint_or_points(tmp_path):
    source, card = _card(tmp_path, "Bridson")
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    code = DRIVER.replace('outcome = core.run(resume=phase == "resume", write_run_summary=False)',
                          'outcome = core.run(check_modules=True, write_run_summary=False)')
    counts = []
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/test/samples.hdf5"
    for attempt in (1, 2):
        output = tmp_path / f"check-{attempt}.json"
        result = subprocess.run([sys.executable, "-c", code, str(source), str(output), "check", "before"],
                                cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
        counts.append(len(SimpleHDF5Writer(str(database)).read_records()))
    assert counts[1] > counts[0]
    checkpoint_dir = tmp_path / f"checkpoints/{card['Scan']['name']}/Bridson"
    assert not (checkpoint_dir / "state.pkl").exists()
    assert not (checkpoint_dir / "points.npy").exists()


def test_missing_or_disabled_checkpoint_starts_fresh(tmp_path, caplog):
    from jarvishep2.runtime.core import Jarvis2Core
    config = {"task_root": str(tmp_path), "scan_name": "no-state",
              "Sampling": {"Method": "Random", "Bounds": {"point_number": 10, "seed": 1}},
              "EnvReqs": {"V2": {"checkpoint": {"enabled": True}}}}
    core = Jarvis2Core(config)
    core.prepare_resume(resume=True)
    assert core._resume_policy == "fresh"
    assert "without state.pkl" in caplog.text
    config["EnvReqs"]["V2"]["checkpoint"]["enabled"] = False
    core = Jarvis2Core(config)
    Path(core.checkpoint_file()).parent.mkdir(parents=True, exist_ok=True)
    Path(core.checkpoint_file()).write_bytes(b"ignored")
    core.prepare_resume(resume=True)
    assert core._resume_policy == "fresh"
    assert "checkpoint.enabled: false" in caplog.text


def test_dynesty_nested_engine_resume_remains_outside_protocol(tmp_path):
    from jarvishep2.runtime.core import Jarvis2Core

    core = Jarvis2Core({"task_root": str(tmp_path), "scan_name": "nested",
                        "Sampling": {"Method": "Dynesty", "Variables": []}})
    core.prepare_resume(resume=True)
    assert core._resume_policy == "resume"


@pytest.mark.parametrize("change,expect_resume", [
    ({"workers": 4}, True),
    ({"seed": 2}, False),
    ({"point_number": 11}, False),
    ({"selection": "x > 0.2"}, False),
])
def test_sampling_fingerprint_refuses_card_drift(tmp_path, change, expect_resume):
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.sampling.runtime_checkpoint import build_payload, build_run_spec, save_checkpoint

    config = {"task_root": str(tmp_path), "scan_name": "card",
              "Sampling": {"Method": "Random", "Bounds": {"point_number": 10, "seed": 1}, "Variables": []},
              "EnvReqs": {"V2": {"workers": 1, "checkpoint": {"enabled": True}}}}
    core = Jarvis2Core(config)
    spec = build_run_spec(config=core.config, scan_name="card", task_root=str(tmp_path),
                          task_result_dir=str(tmp_path), sampler_name="Random")
    save_checkpoint(core.checkpoint_file(), build_payload(run_spec=spec, sampler_state={}))
    changed = json.loads(json.dumps(config))
    if "workers" in change:
        changed["EnvReqs"]["V2"].update(change)
    elif "selection" in change:
        changed["Sampling"].update(change)
    else:
        changed["Sampling"]["Bounds"].update(change)
    resumed = Jarvis2Core(changed)
    resumed.prepare_resume(resume=True)
    assert (resumed._resume_policy == "resume") is expect_resume


@pytest.mark.parametrize("field", ["seed", "radius", "max_attempt", "length"])
def test_bridson_bounds_change_rejected_before_loading_points(tmp_path, field):
    from jarvishep2.runtime.core import Jarvis2Core
    from jarvishep2.sampling.runtime_checkpoint import build_payload, build_run_spec, save_checkpoint

    variables = [{"name": name, "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1, "length": 1}}}
                 for name in ("x", "y")]
    config = {"task_root": str(tmp_path), "scan_name": "bounds",
              "Sampling": {"Method": "Bridson", "Bounds": {"seed": 1, "radius": 0.4, "max_attempt": 20},
                           "Variables": variables},
              "EnvReqs": {"V2": {"checkpoint": {"enabled": True}}}}
    core = Jarvis2Core(config)
    spec = build_run_spec(config=core.config, scan_name="bounds", task_root=str(tmp_path),
                          task_result_dir=str(tmp_path), sampler_name="Bridson")
    save_checkpoint(core.checkpoint_file(), build_payload(run_spec=spec, sampler_state={}))
    changed = json.loads(json.dumps(config))
    if field == "length":
        changed["Sampling"]["Variables"][0]["distribution"]["parameters"]["length"] = 2
    else:
        changed["Sampling"]["Bounds"][field] = 2 if field != "radius" else 0.2
    resumed = Jarvis2Core(changed)
    resumed.prepare_resume(resume=True)
    assert resumed._resume_policy == "fresh"


def test_grid_and_csv_point_fingerprints_reject_changed_sources(tmp_path):
    from jarvishep2.sampling.grid import Grid
    from jarvishep2.sampling.csv_sampler import CSVSampler

    config = {"Sampling": {"Method": "Grid", "Bounds": {"seed": 1},
              "Variables": [{"name": "x", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1, "num": 2}}}]}}
    grid = Grid(); grid.set_config(config)
    state = grid.export_runtime_state()
    config["Sampling"]["Variables"][0]["distribution"]["parameters"]["num"] = 3
    changed = Grid(); changed.set_config(config)
    with pytest.raises(ValueError, match="Grid point set changed"):
        changed.import_runtime_state(state)
    csv = tmp_path / "points.csv"
    csv.write_text("uuid,x\na,1\n")
    csv_config = {"Sampling": {"Method": "CSV", "Bounds": {"path": str(csv)}}}
    original = CSVSampler(); original.set_config(csv_config)
    state = original.export_runtime_state()
    csv.write_text("uuid,x\na,2\n")
    changed_csv = CSVSampler(); changed_csv.set_config(csv_config)
    with pytest.raises(ValueError, match="CSV source file changed"):
        changed_csv.import_runtime_state(state)


def test_indexed_checkpoint_size_does_not_scale_with_submissions(tmp_path):
    from jarvishep2.sampling.randoms import RandomS
    from jarvishep2.sampling.runtime_checkpoint import build_payload, build_run_spec, save_checkpoint

    sizes = []
    for point_number in (10, 100_000):
        sampler = RandomS()
        config = {"Sampling": {"Method": "Random", "Bounds": {"point_number": point_number, "seed": 1},
                  "Variables": [{"name": "x", "distribution": {"type": "Flat", "parameters": {"min": 0, "max": 1}}}]}}
        sampler.set_config(config)
        sampler.initialize()
        sampler.record_submitted_batch([f"id-{i}" for i in range(point_number)])
        assert sampler.submitted_uuids == frozenset()
        spec = build_run_spec(config=config, scan_name="size", task_root=str(tmp_path),
                              task_result_dir=str(tmp_path), sampler_name="Random")
        path = tmp_path / f"state-{point_number}.pkl"
        save_checkpoint(str(path), build_payload(run_spec=spec, sampler_state=sampler.export_runtime_state()))
        sizes.append(path.stat().st_size)
    assert abs(sizes[1] - sizes[0]) < 100


def test_final_trim_removes_stale_sealed_shard_rows(tmp_path):
    from jarvishep2.io.database import trim_scan_records

    database = tmp_path / "DATABASE"
    sealed = database / "samples.0001.hdf5"
    live = database / "samples.hdf5"
    SimpleHDF5Writer(str(sealed)).add_records([
        {"uuid": "old-0", "sample_index": 0, "z": -999},
        {"uuid": "old-tail", "sample_index": 3, "z": -999},
    ])
    SimpleHDF5Writer(str(live)).add_records([
        {"uuid": "new-0", "sample_index": 0, "z": 1},
        {"uuid": "new-1", "sample_index": 1, "z": 2},
    ])
    assert trim_scan_records(str(database), 2) == 2
    assert not sealed.exists()
    rows = SimpleHDF5Writer(str(live)).read_records()
    assert [(row["sample_index"], row["z"]) for row in rows] == [(0, 1), (1, 2)]


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("sharded", [False, True], ids=["live", "rolling-shards"])
def test_short_scan_removes_tail_of_completed_long_scan(tmp_path, sharded):
    source, card = _card(tmp_path, "Grid")
    if sharded:
        card["EnvReqs"]["V2"]["archiver"] = {"max_hdf5_bytes": 1, "batch_size": 1}
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    for size, expected in ((4, 16), (2, 4)):
        for variable in card["Sampling"]["Variables"]:
            variable["distribution"]["parameters"]["num"] = size
        source.write_text(yaml.safe_dump(card))
        output = tmp_path / f"outcome-{size}.json"
        result = subprocess.run([sys.executable, "-c", DRIVER, str(source), str(output), "fresh", "before"],
                                cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
        assert result.returncode == 0, result.stdout + result.stderr
        from jarvishep2.io.database import discover_samples_hdf5_shards
        paths = discover_samples_hdf5_shards(str(database.parent))
        rows = [row for path in paths for row in SimpleHDF5Writer(path).read_records()]
        if sharded:
            assert len(paths) > 1
        assert len(rows) == expected
        assert {int(row["sample_index"]) for row in rows} == set(range(expected))


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("disabled", [False, True], ids=["missing-state", "disabled-checkpoint"])
def test_resume_without_usable_checkpoint_restarts_from_zero(tmp_path, disabled):
    source, card = _card(tmp_path, "Bridson")
    card["EnvReqs"]["V2"]["checkpoint"] = {"enabled": not disabled}
    source.write_text(yaml.safe_dump(card))
    config = load_task_yaml(str(source))
    sampler = Distributor.set_method("Bridson")
    sampler.set_config(config)
    stale = []
    while (sample := sampler.propose_next()) is not None:
        stale.append({"uuid": sample.uuid, "sample_index": sample.sample_index,
                      "status": "Completed", "x": 0., "y": 0., "z": -999.})
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records(stale)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    output = tmp_path / "resume-without-state.json"
    code = DRIVER.replace('if phase == "resume" and core.config["Sampling"]["Method"] == "Bridson":', 'if False:')
    result = subprocess.run([sys.executable, "-c", code, str(source), str(output), "resume", "before"],
                            cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "starting fresh" in result.stdout + result.stderr
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert len(rows) == len(stale)
    assert all(float(row["z"]) != -999 for row in rows)
    checkpoint_dir = tmp_path / f"checkpoints/{card['Scan']['name']}/Bridson"
    assert (checkpoint_dir / "state.pkl").exists() is not disabled
    assert (checkpoint_dir / "points.npy").exists() is not disabled
