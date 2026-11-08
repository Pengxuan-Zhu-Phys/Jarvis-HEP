"""Checkpoint authoritative resume protocol, including real runtime interruption."""
from __future__ import annotations

import json
import os
import shutil
from collections import Counter
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
import json, sys
from pathlib import Path
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.distributor import Distributor

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
    def propose(self):
        if self._accepted_index >= 5:
            if timing == "after":
                self.save_checkpoint(force=True)
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
def test_fresh_interrupt_resume_recalculates_only_snapshot_unfinished_and_new_points(tmp_path, method, timing):
    from jarvishep2.sampling.sample_checkpoint import read_pickle

    source, card = _card(tmp_path, method)
    sampler = Distributor.set_method(method)
    sampler.set_config(load_task_yaml(str(source)))
    stale = []
    while (sample := sampler.propose_next()) is not None:
        stale.append({"uuid": sample.uuid, "sample_index": sample.sample_index,
                      "status": "Completed", "x": 0., "y": 0., "z": -999.})
    total = len(stale)
    assert total > 5
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records(stale)
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
            checkpoint = tmp_path / f"checkpoints/{card['Scan']['name']}/{method}/runtime.pickle"
            state = read_pickle(str(checkpoint))["sampler_state"]
            assert (state["next_sample_index"] > 0) == (timing == "after")
            redo = set(state["unfinished"]) | set(range(state["next_sample_index"], total))
            rows_at_interrupt = len(SimpleHDF5Writer(str(database)).read_records())
        else:
            assert result.returncode == 0, result.stdout + result.stderr
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert rows[:total] == stale
    assert {int(row["sample_index"]) for row in rows[total:]} == set(range(total))
    assert {int(row["sample_index"]) for row in rows[rows_at_interrupt:]} == redo
    assert all(float(row["z"]) != -999 for row in rows[total:])
    assert not checkpoint.with_name("state.pkl").exists()
    assert checkpoint.with_name("all.pickle").exists() == (method in {"Grid", "Bridson"})


NESTED_DRIVER = r"""
import json, os, sys, time
from pathlib import Path
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.dynesty_sampler import DynestySampler
from jarvishep2.sampling.sample_checkpoint import read_pickle

card, output, phase = sys.argv[1:]
core = Jarvis2Core()
core.load_task_yaml(card)
core.config.update(skip_draw_flowchart=True, skip_render_plots=True)
original_config = DynestySampler.set_config
def fast_engine_saves(self, config_info):
    original_config(self, config_info)
    self._checkpoint_every_sec = 0.2
DynestySampler.set_config = fast_engine_saves
if phase == "interrupt":
    original_save = DynestySampler.save_checkpoint
    landed = None
    def save(self, *, force=False):
        global landed
        if landed is not None and not force and self._accepted_index > landed:
            Path(output).with_suffix(".tail").write_text(str(self._accepted_index))
            raise KeyboardInterrupt
        done = original_save(self, force=force)
        if landed is not None or not done or self._accepted_index <= 0:
            return done
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline:
            path = core.checkpoint_file()
            if os.path.isfile(path):
                state = read_pickle(path)["sampler_state"]
                cursor = int(state.get("checkpoint_index", 0))
                if cursor > 0:
                    landed = cursor
                    return done
            time.sleep(0.05)
        raise AssertionError("engine snapshot was never landed")
    DynestySampler.save_checkpoint = save
outcome = core.run(resume=phase == "resume", write_run_summary=False)
Path(output).write_text(json.dumps(outcome.to_dict()))
raise SystemExit(outcome.exit_code)
"""


MCMC_DRIVER = r"""
import json, sys
from pathlib import Path
from jarvishep2.runtime import _resume_service
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.mcmc_sampler import MCMCBaseSampler

card, output, phase = sys.argv[1:]
core = Jarvis2Core()
core.load_task_yaml(card)
core.config.update(skip_draw_flowchart=True, skip_render_plots=True)
def forbidden(*args, **kwargs):
    raise AssertionError("MCMC checkpoint attempted to read HDF5 progress")
_resume_service.read_persisted_sample_index_state = forbidden
_resume_service.read_persisted_outcome_counts = forbidden
if phase == "interrupt":
    original = MCMCBaseSampler.save_checkpoint
    def save(self, *, force=False):
        if not force and self._generation >= 2:
            assert original(self, force=True)
            assert self._pickle_writer.wait()
            raise KeyboardInterrupt
        return original(self, force=force)
    MCMCBaseSampler.save_checkpoint = save
outcome = core.run(resume=phase == "resume", write_run_summary=False)
Path(output).write_text(json.dumps(outcome.to_dict()))
raise SystemExit(outcome.exit_code)
"""


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("method", ["ToyMCMC", "PTEnsemble"])
def test_mcmc_runtime_interrupt_and_resume(tmp_path, method):
    from jarvishep2.sampling.sample_checkpoint import read_pickle

    source, card = _card(tmp_path, "Random")
    card["Sampling"]["Method"] = method
    card["Sampling"]["Bounds"] = {"seed": 11, "num_chains": 4, "num_iters": 7, "proposal_scale": 0.1}
    if method == "PTEnsemble":
        card["Sampling"]["Bounds"]["temperature_ladder"] = [1.0, 2.0, 4.0, 8.0]
    source.write_text(yaml.safe_dump(card))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    checkpoint = tmp_path / f"checkpoints/{card['Scan']['name']}/{method}/runtime.pickle"
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    for phase, exit_code in (("interrupt", 130), ("resume", 0)):
        output = tmp_path / f"{phase}.json"
        result = subprocess.run([sys.executable, "-c", MCMC_DRIVER, str(source), str(output), phase],
                                cwd=tmp_path, env=env, text=True, capture_output=True, timeout=90)
        assert result.returncode == exit_code, result.stdout + result.stderr
        state = read_pickle(str(checkpoint))["sampler_state"]
        assert not checkpoint.with_name("state.pkl").exists()
        if phase == "interrupt":
            assert 2 <= state["progress"]["generation"] < 7
            assert bool(state["pending"]) == (method == "ToyMCMC")
            pending = set(state["pending"])
            rows_at_interrupt = SimpleHDF5Writer(str(database)).read_records()
        else:
            assert state["progress"]["finished"]
            assert [c["iter"] for c in state["chains"]] == [7] * 4
            assert not state["pending"]
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert rows[:len(rows_at_interrupt)] == rows_at_interrupt
    assert pending <= {row["uuid"] for row in rows[len(rows_at_interrupt):]}


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("method", ["Dynesty", "MultiNest"])
def test_nested_interrupt_resume_covers_every_index_once_past_the_redo_window(tmp_path, method):
    from jarvishep2.sampling.Source.Dynesty.py.dynesty.utils import restore_sampler
    from jarvishep2.sampling.sample_checkpoint import read_pickle

    source, card = _card(tmp_path, "Random")
    card["Sampling"]["Method"] = method
    run_nested = ({"maxiter": 150, "maxbatch": 1, "print_progress": False}
                  if method == "Dynesty" else {"maxiter": 150, "print_progress": False})
    card["Sampling"]["Bounds"] = {"nlive": 8, "seed": 11, "dlogz": 0.5, "run_nested": run_nested}
    source.write_text(yaml.safe_dump(card))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"),
               JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"))
    checkpoint = tmp_path / f"checkpoints/{card['Scan']['name']}/{method}/runtime.pickle"
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE"
    for phase in ("interrupt", "resume"):
        output = tmp_path / f"{phase}.json"
        result = subprocess.run([sys.executable, "-c", NESTED_DRIVER, str(source), str(output), phase],
                                cwd=tmp_path, env=env, text=True, capture_output=True, timeout=180)
        (tmp_path / f"{phase}.log").write_text(result.stdout + result.stderr)
        assert output.exists(), result.stdout + result.stderr
        if phase == "interrupt":
            assert result.returncode == 130, result.stdout + result.stderr
            assert not checkpoint.with_name("state.pkl").exists()
            payload = read_pickle(str(checkpoint))
            assert payload["format"] == "jarvis-hep.sampler-pickle"
            cursor = int(payload["sampler_state"]["checkpoint_index"])
            tail = int(output.with_suffix(".tail").read_text())
            # The interrupt keeps the landed snapshot, not the live engine.
            assert 0 < cursor < tail
            engine = restore_sampler(str(checkpoint.with_name("nested_engine.pkl")))
            assert int(engine.ncall) > 0
        else:
            assert result.returncode == 0, result.stdout + result.stderr
    rows = []
    for shard in sorted(database.glob("*.hdf5")):
        rows.extend(SimpleHDF5Writer(str(shard)).read_records())
    counts = Counter(int(row["sample_index"]) for row in rows)
    final = read_pickle(str(checkpoint))["sampler_state"]
    assert final["finished"] is True
    final_cursor = int(final["checkpoint_index"])
    assert final_cursor > tail
    # Resume continues from the landed cursor. Indices before it are not
    # redone. Indices at or after the interrupt tail are new. The open
    # interval is the only place a row can appear twice.
    assert set(range(cursor, final_cursor)).issubset(set(counts))
    assert {index for index, n in counts.items() if n > 1} <= set(range(cursor, tail))
    assert max(counts.values()) <= 2
    assert all(counts[index] == 1 for index in range(tail, final_cursor))


@pytest.mark.parametrize("method", ["Dynesty", "MultiNest"])
def test_nested_resume_without_runtime_pickle_refuses(tmp_path, method):
    from jarvishep2.runtime.core import Jarvis2Core

    core = Jarvis2Core({"task_root": str(tmp_path), "scan_name": "nested",
                        "Sampling": {"Method": method, "Variables": []}})
    # A leftover engine copy is never a resume source.
    engine_copy = Path(core.checkpoint_file()).with_name("nested_engine.pkl")
    engine_copy.parent.mkdir(parents=True)
    engine_copy.write_bytes(b"engine")
    with pytest.raises(FileNotFoundError, match="runtime.pickle"):
        core.prepare_resume(resume=True)
    assert engine_copy.read_bytes() == b"engine"
    assert core._resume_policy == "fresh"
    assert core._resume_checkpoint_payload is None
