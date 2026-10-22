"""Native load/run/feedback/resume against private Redis and spawn Workers."""

from __future__ import annotations

import csv
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import time
from pathlib import Path

import pytest
import yaml

from jarvishep2.database import SimpleHDF5Writer

EXAMPLES = Path(__file__).resolve().parents[2] / "Jarvis-HEP-Sampler" / "examples"
pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(not shutil.which("redis-server"), reason="Redis required"),
]
ENTRY = [
    sys.executable,
    "-c",
    "from jarvishep2.cli.client import main; raise SystemExit(main())",
]


def setup(tmp_path, *, budget=61):
    env = dict(
        os.environ,
        JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"),
        JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "opera-cache"),
        JARVIS_OPERAS_CURVE_CACHE_PERSIST_DIR=str(tmp_path / "opera-cache/persist"),
        JARVIS_OPERAS_CURVE_CACHE_SOURCES_DIR=str(tmp_path / "opera-cache/sources"),
        JARVIS_OPERAS_CURVE_CACHE_OVERRIDES_DIR=str(tmp_path / "opera-cache/overrides"),
    )
    card = yaml.safe_load((EXAMPLES / "SliceMCMC.eggbox.task.yaml").read_text())
    card["Scan"]["name"] = "native-slice"
    card["Sampling"]["Bounds"].update(max_evaluations=budget, burn_in=0, steps=1000)
    card["EnvReqs"]["V2"]["workers"] = 2
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    card["EnvReqs"]["V2"]["redis"] = {"host": "127.0.0.1", "port": port, "db": 0}
    return env, card


def cli(tmp_path, env, label, *args, expected=0):
    result = subprocess.run(
        ENTRY + list(args),
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=90,
    )
    (tmp_path / f"{label}.log").write_text(result.stdout + result.stderr)
    assert result.returncode == expected, (result.stdout + result.stderr)[-6000:]
    return result


def rows(tmp_path):
    database = tmp_path / "outputs/native-slice/DATABASE"
    data = SimpleHDF5Writer(str(database / "samples.hdf5")).read_records()
    chain = list(csv.DictReader((database / "chain.csv").open()))
    return data, chain


@pytest.mark.parametrize("case", ["opera", "zero_density", "selection", "failure"])
def test_native_terminal_protocol_and_exact_budget(tmp_path, case):
    env, card = setup(tmp_path, budget=13)
    if case == "zero_density":
        card["Sampling"]["LogLikelihood"][0]["expression"] = "-oo"
    elif case == "selection":
        card["Sampling"]["selection"] = "x < 0"
    elif case == "failure":
        card["Sampling"]["LogLikelihood"][0]["expression"] = "nan"
    task = tmp_path / "task.yaml"
    task.write_text(yaml.safe_dump(card))
    cli(tmp_path, env, "load", "sampler", "load", str(EXAMPLES / "SliceMCMC.py"))
    cli(tmp_path, env, "validate", "validate", str(task))
    cli(
        tmp_path,
        env,
        "run",
        "run",
        str(task),
        "--skip-draw-flowchart",
        expected=1 if case == "failure" else 0,
    )
    data, chain = rows(tmp_path)
    assert len(data) == (0 if case == "selection" else 4 if case == "failure" else 13)
    assert len({row["uuid"] for row in data}) == len(data)
    if case in {"zero_density", "selection", "failure"}:
        assert chain == []
    if case == "opera":
        assert chain
        for row in data:
            assert float(row["LogL"]) == pytest.approx(
                (
                    __import__("math").sin(float(row["x"]))
                    * __import__("math").cos(float(row["y"]))
                    + 2
                )
                ** 5
            )
    checkpoint = tmp_path / "checkpoints/native-slice/SliceMCMC/state.pkl"
    assert checkpoint.exists()


def test_interrupted_generation_replays_after_remove_and_changes_worker_count(tmp_path):
    env, card = setup(tmp_path, budget=601)
    task = tmp_path / "task.yaml"
    task.write_text(yaml.safe_dump(card))
    cli(tmp_path, env, "load", "sampler", "load", str(EXAMPLES / "SliceMCMC.py"))
    with (tmp_path / "interrupted.log").open("w") as log:
        process = subprocess.Popen(
            ENTRY + ["run", str(task), "--skip-draw-flowchart"],
            cwd=tmp_path,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            chain = tmp_path / "outputs/native-slice/DATABASE/chain.csv"
            deadline = time.monotonic() + 45
            while time.monotonic() < deadline:
                assert process.poll() is None, (
                    tmp_path / "interrupted.log"
                ).read_text()[-6000:]
                if chain.exists() and len(chain.read_text().splitlines()) > 8:
                    break
                time.sleep(0.02)
            else:
                pytest.fail("no chain progress")
            process.send_signal(signal.SIGINT)
            assert process.wait(timeout=45) == 130
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait(timeout=45)
    cli(tmp_path, env, "remove", "sampler", "remove", "SliceMCMC", "--yes")
    card["EnvReqs"]["V2"]["workers"] = 3
    task.write_text(yaml.safe_dump(card))
    cli(tmp_path, env, "resume", "run", str(task), "--resume", "--skip-draw-flowchart")
    data, chain = rows(tmp_path)
    assert len(data) == len({r["uuid"] for r in data}) == 601
    assert len({(r["chain"], r["step"]) for r in chain}) == len(chain)
    # Compare the complete native posterior trajectory with the same SDK state machine.
    from jarvis_hep_sampler import GenerationalSession, load_plugin_file
    import math

    s = GenerationalSession.from_bounds(
        load_plugin_file(EXAMPLES / "SliceMCMC.py"),
        card["Sampling"]["Bounds"],
        variable_names=["x", "y"],
    )
    s.prepare()
    expected = []
    while (batch := s.ask()) is not None:
        feedback = [
            {
                "id": p["id"],
                "status": "ok",
                "logL": (
                    math.sin(p["u"][0] * 10 * math.pi)
                    * math.cos(p["u"][1] * 10 * math.pi)
                    + 2
                )
                ** 5,
            }
            for p in batch
        ]
        expected.extend(s.tell(feedback))
    assert len(chain) == len(expected)
    for actual, event in zip(chain, expected, strict=True):
        assert [float(actual["u_x"]), float(actual["u_y"])] == event["u"]
        assert float(actual["logL"]) == pytest.approx(event["logL"])
    summary = json.loads(
        (tmp_path / "outputs/native-slice/DATABASE/sampler_summary.json").read_text()
    )
    assert summary["proposals"] == 601 and summary["retained_draws"] == len(chain)
