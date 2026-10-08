"""Real Redis/spawn interruption and DATABASE-only recovery for indexed plug-ins."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import signal
import socket
import subprocess
import sys
import time

import pytest
import yaml

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(
        shutil.which("redis-server") is None, reason="real Redis is required"
    ),
]

SOURCE = """
import os
import time
import numpy as np
if os.environ.get("JHS_TEST_MARKER"):
    with open(os.environ["JHS_TEST_MARKER"], "a") as handle:
        handle.write(f"{os.getpid()}\\n")
PLUGIN = {"api_version": 1, "method": "RecoverIndexed", "version": "1.0",
          "profile": "independent", "block_size": 8}
def prepare(setup):
    return None
def draw_block(context, rng, plan):
    time.sleep(.12)
    i = context.block_index * context.block_size + np.arange(context.block_size)
    # Unused rows in the last block still obey the output range contract.
    u = np.minimum((i + .5) / context.point_number, 1)
    return np.repeat(u[:, None], context.dimension, axis=1)
"""


@pytest.mark.parametrize("database_only", [False, True])
def test_real_interruption_reload_remove_and_resume(tmp_path, database_only):
    entry = [
        sys.executable,
        "-c",
        "from jarvishep2.cli.client import main; raise SystemExit(main())",
    ]
    env = dict(
        os.environ,
        JARVIS_SAMPLERS_HOME=str(tmp_path / "registry"),
        PYTHONDONTWRITEBYTECODE="1",
        HOME=str(tmp_path / "home"),
        JHS_TEST_MARKER=str(tmp_path / "imports.txt"),
    )
    env.pop("PYTHONPATH", None)
    source = tmp_path / "RecoverIndexed.py"
    source.write_text(SOURCE)
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    card = {
        "Scan": {"name": "recover"},
        "EnvReqs": {
            "V2": {
                "workers": 2,
                "batch_size": 4,
                "redis": {"host": "127.0.0.1", "port": port, "db": 0},
            }
        },
        "Sampling": {
            "Method": "RecoverIndexed",
            "Variables": [
                {
                    "name": key,
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                }
                for key in ("x", "y")
            ],
            "Bounds": {"seed": 7, "point_number": 211},
            "selection": "x > 0.2",
            "LogLikelihood": [{"name": "LogL_Z", "expression": "z"}],
        },
        "Operas": {
            "Modules": [
                {
                    "name": "Eggbox",
                    "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                    "input": [{"name": key, "expression": key} for key in ("x", "y")],
                    "output": [{"name": "z", "entry": "z"}],
                }
            ]
        },
    }
    task = tmp_path / "task.yaml"
    task.write_text(yaml.safe_dump(card))

    def cli(label, *args):
        result = subprocess.run(
            entry + list(args),
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=60,
        )
        (tmp_path / f"{label}.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 0, (result.stdout + result.stderr)[-5000:]

    cli("load", "sampler", "load", str(source))
    database = tmp_path / "outputs/recover/DATABASE/samples.hdf5"
    with (tmp_path / "interrupt.log").open("w") as log:
        process = subprocess.Popen(
            entry + ["run", str(task), "--skip-draw-flowchart"],
            cwd=tmp_path,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        try:
            deadline = time.monotonic() + 45
            while time.monotonic() < deadline:
                assert process.poll() is None, (tmp_path / "interrupt.log").read_text()[
                    -5000:
                ]
                if database.exists():
                    try:
                        if len(SimpleHDF5Writer(str(database)).read_records()) >= 8:
                            break
                    except OSError:
                        pass
                time.sleep(0.1)
            else:
                pytest.fail("no durable prefix before interrupt")
            process.send_signal(signal.SIGINT)
            assert process.wait(timeout=40) == 130
        finally:
            if process.poll() is None:
                process.send_signal(signal.SIGTERM)
                process.wait(timeout=40)
    checkpoint = tmp_path / "checkpoints/recover/RecoverIndexed/state.pkl"
    assert checkpoint.exists()
    if database_only:
        checkpoint.unlink()
    source.write_text(SOURCE.replace('"version": "1.0"', '"version": "2.0"'))
    cli("reload", "sampler", "load", str(source))
    cli("remove", "sampler", "remove", "RecoverIndexed", "--yes")

    def refuse(label, message):
        protected = list(database.parent.iterdir()) + (
            [checkpoint] if checkpoint.exists() else []
        )
        before = {
            p: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in protected
            if p.is_file()
        }
        imports = (tmp_path / "imports.txt").read_bytes()
        result = subprocess.run(
            entry + ["run", str(task), "--resume", "--skip-draw-flowchart"],
            cwd=tmp_path,
            env=env,
            text=True,
            capture_output=True,
            timeout=30,
        )
        (tmp_path / f"{label}.log").write_text(result.stdout + result.stderr)
        assert result.returncode == 1 and message in result.stdout + result.stderr
        assert {p: hashlib.sha256(p.read_bytes()).hexdigest() for p in before} == before
        assert (tmp_path / "imports.txt").read_bytes() == imports
        with socket.socket() as sock:
            assert sock.connect_ex(("127.0.0.1", port)) != 0

    original_card = task.read_bytes()
    card["Sampling"]["Bounds"]["seed"] = 8
    task.write_text(yaml.safe_dump(card))
    refuse("changed-seed", "science settings differ")
    task.write_bytes(original_card)
    binding = json.loads((database.parent / "sampler_identity.json").read_text())
    bound_source = tmp_path / binding["source_file"]
    original_source = bound_source.read_bytes()
    bound_source.write_text("raise AssertionError('must not execute damaged source')")
    refuse("damaged-source", "SHA-256 mismatch")
    bound_source.write_bytes(original_source)
    if not database_only:
        original_checkpoint = checkpoint.read_bytes()
        checkpoint.write_bytes(b"this is not a pickle")
        refuse("damaged-checkpoint", "checkpoint rejected")
        checkpoint.write_bytes(original_checkpoint)
    cli("resume", "run", str(task), "--resume", "--skip-draw-flowchart")
    assert len((tmp_path / "imports.txt").read_text().splitlines()) == 4
    records = SimpleHDF5Writer(str(database)).read_records()
    candidates = [i for i in range(211) if (i + 0.5) / 211 > 0.2]
    ordered = sorted(records, key=lambda r: int(r["sample_index"]))
    assert len(ordered) == len(candidates) == 169
    binding = json.loads((database.parent / "sampler_identity.json").read_text())
    prefix = f"plugin:RecoverIndexed:{binding['hard_identity'][:16]}"
    for accepted, (record, candidate) in enumerate(zip(ordered, candidates)):
        assert float(record["x"]) == float(record["y"]) == (candidate + 0.5) / 211
        assert int(record["sample_index"]) == accepted
        assert record["uuid"] == deterministic_sampler_uuid(
            prefix=prefix, seed=7, sample_index=accepted
        )
    assert len({record["uuid"] for record in records}) == len(records)
    assert binding["metadata"]["version"] == "1.0"
