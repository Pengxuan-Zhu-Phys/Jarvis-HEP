"""Units through persisted JO sources, real Redis/spawn, DATABASE and restore."""

import csv
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import sys

import pytest
import yaml

from jarvishep2.operas_units import DATABASE_KEY
from jarvishep2.io.database import SimpleHDF5Writer

FIXTURE = Path(__file__).parent / "fixtures/unit_ops.py"


class _RuntimeEnv(dict):
    def __repr__(self):
        return "<inherited runtime environment>"


DRIVER = """
import json,sys
from pathlib import Path
from jarvis_operas import persist_user_ops
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.operas_units import UnitBridgeError
card,fixture,output,phase=sys.argv[1:]
persist_user_ops(fixture,namespace="hep_unit_fixture")
core=Jarvis2Core()
core.load_task_yaml(card)
core.config.update(skip_draw_flowchart=True,skip_render_plots=True)
if phase=="interrupt":
    from jarvishep2.distributor import Distributor
    Random=type(Distributor.set_method(core.config["Sampling"]["Method"]))
    original=Random.propose_next
    def propose(self):
        if self._accepted_index>=2:
            raise KeyboardInterrupt
        return original(self)
    Random.propose_next=propose
try:
    result=core.run(resume=phase=="resume",check_modules=phase=="check",write_run_summary=False)
    Path(output).write_text(json.dumps(result.to_dict()))
    raise SystemExit(result.exit_code)
except UnitBridgeError as exc:
    Path(output).write_text(json.dumps(dict(code=exc.code,error=str(exc))))
    raise SystemExit(2)
"""


def _setup(tmp_path, *, store=False, mode="call"):
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    cfg = {
        "Scan": {"name": "unit_scan"},
        "EnvReqs": {
            "Check_default_dependencies": {"required": False},
            "V2": {
                "workers": 2,
                "batch_size": 1,
                "store_samples": store,
                "redis": {"host": "127.0.0.1", "port": port, "db": 0},
            },
        },
        "Sampling": {
            "Method": "Random",
            "Bounds": {"point_number": 4, "seed": 31},
            "Variables": [
                {
                    "name": "mass",
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.199, "max": 0.201},
                    },
                }
            ],
        },
        "Operas": {
            "Modules": [
                {
                    "name": "Limit",
                    "operator": "hep_unit_fixture.limit",
                    "call_mode": mode,
                    "input": [{"name": "mass_GeV", "entry": "mass", "unit": "TeV"}],
                    "output": [
                        {"name": "rate_fb", "entry": "value", "unit": "fb"},
                        {"name": "record", "entry": ""},
                    ],
                }
            ]
        },
    }
    path = tmp_path / "task.yaml"
    path.write_text(yaml.safe_dump(cfg))
    (tmp_path / ".jarvis-project.json").write_text("{}")
    env = _RuntimeEnv(
        os.environ,
        PYTHONDONTWRITEBYTECODE="1",
        JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "refs.json"),
        JARVIS_OPERAS_SOURCES_FILE=str(tmp_path / "jo/sources.json"),
        JARVIS_OPERAS_PERSIST_FILE=str(tmp_path / "jo/ops.json"),
        JARVIS_OPERAS_OVERRIDES_FILE=str(tmp_path / "jo/overrides.json"),
        JARVIS_OPERAS_USER_SOURCE_DIR=str(tmp_path / "jo/sources"),
        JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curves"),
        JARVIS_SAMPLERS_HOME=str(tmp_path / "samplers"),
    )
    return path, cfg, env


def _run(tmp_path, path, env, phase):
    output = tmp_path / (phase + ".json")
    process = subprocess.run(
        [
            sys.executable,
            "-c",
            DRIVER,
            str(path),
            str(FIXTURE.resolve()),
            str(output),
            phase,
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        capture_output=True,
        timeout=100,
    )
    (tmp_path / (phase + ".log")).write_text(process.stdout + process.stderr)
    assert output.exists(), (process.stdout + process.stderr)[-7000:]
    return process, json.loads(output.read_text())


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("store,mode", [(False, "call"), (True, "acall")])
def test_real_units_yaml_spawn_archive_and_resume(tmp_path, store, mode):
    path, cfg, env = _setup(tmp_path, store=store, mode=mode)
    process, result = _run(tmp_path, path, env, "run")
    assert process.returncode == 0, result
    directory = tmp_path / "outputs/unit_scan/DATABASE"
    rows = SimpleHDF5Writer(str(directory / "samples.hdf5")).read_records()
    assert len(rows) == 4
    for row in rows:
        assert row["rate_fb"] == pytest.approx(38.332)
        assert row["record"]["mass_received"] == pytest.approx(row["mass"] * 1000)
        assert row["record"]["unit"] == "pb"
        assert row[DATABASE_KEY]["observable_units"] == {"/rate_fb": "fb"}
        assert len(row[DATABASE_KEY]["calls"]) == 1
    manifest = json.loads((directory / "unit_execution.json").read_text())
    assert all(
        r[DATABASE_KEY]["execution_identity"] == manifest["execution_identity"]
        for r in rows
    )
    exported = list(csv.DictReader((directory / "samples.csv").open()))
    assert json.loads(exported[0][DATABASE_KEY]) == rows[0][DATABASE_KEY]
    if not store:
        assert not (tmp_path / "outputs/unit_scan/SAMPLE").exists()
    process, result = _run(tmp_path, path, env, "resume")
    assert process.returncode == 0, result
    assert len(SimpleHDF5Writer(str(directory / "samples.hdf5")).read_records()) == 4
    before = (directory / "samples.hdf5").read_bytes()
    cfg["Operas"]["Modules"][0]["output"][0]["unit"] = "pb"
    path.write_text(yaml.safe_dump(cfg))
    process, result = _run(tmp_path, path, env, "resume")
    assert process.returncode == 2 and result["code"] == "JV2-UNIT-004"
    assert (directory / "samples.hdf5").read_bytes() == before
    # check is an independent diagnostic batch, even after the formal unit
    # policy changed. It must leave the formal manifest untouched.
    process, result = _run(tmp_path, path, env, "check")
    assert process.returncode == 0, result
    assert json.loads((directory / "unit_execution.json").read_text()) == manifest
    tests = SimpleHDF5Writer(str(directory / "test/samples.hdf5")).read_records()
    assert tests and tests[0]["rate_fb"] == pytest.approx(0.038332)


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis required")
@pytest.mark.parametrize("method", ["Random", "Grid", "CSV", "Bridson"])
def test_real_unit_interrupt_resume(tmp_path, method):
    path, cfg, env = _setup(tmp_path)
    cfg["Sampling"]["Method"] = method
    if method == "Grid":
        cfg["Sampling"]["Bounds"] = {"seed": 31}
        cfg["Sampling"]["Variables"][0]["distribution"]["parameters"]["num"] = 4
    elif method == "CSV":
        source = tmp_path / "points.csv"
        source.write_text("uuid,mass\n" + "".join(f"point{i},0.2\n" for i in range(4)))
        cfg["Sampling"]["Bounds"] = {"path": str(source)}
    elif method == "Bridson":
        cfg["Sampling"]["Bounds"] = {"seed": 31, "radius": 0.2, "max_attempt": 20}
        cfg["Sampling"]["Variables"][0]["distribution"]["parameters"]["length"] = 1.0
        cfg["Sampling"]["Variables"].append(
            {
                "name": "aux",
                "distribution": {
                    "type": "Flat",
                    "parameters": {"min": 0.0, "max": 1.0, "length": 1.0},
                },
            }
        )
    path.write_text(yaml.safe_dump(cfg))
    process, result = _run(tmp_path, path, env, "interrupt")
    assert process.returncode == 130, result
    process, result = _run(tmp_path, path, env, "resume")
    assert process.returncode == 0, result
    rows = SimpleHDF5Writer(
        str(tmp_path / "outputs/unit_scan/DATABASE/samples.hdf5")
    ).read_records()
    expected = 4
    if method == "Bridson":
        import numpy as np

        expected = len(np.load(tmp_path / "checkpoints/unit_scan/Bridson/points.npy"))
    assert len(rows) == expected and len({r["sample_index"] for r in rows}) == expected
    assert all(r["rate_fb"] == pytest.approx(38.332) for r in rows)
