"""Verify an installed complete HEP wheel outside its source checkout."""

from pathlib import Path
import importlib.metadata as metadata
import json
import os
import socket
import subprocess
import sys
import tarfile
import tempfile

import yaml
import jarvishep2
import jarvis_hep_sampler
from jarvishep2.project_packager import create_project_package
from jarvishep2.project_samplers import INDEX, inspect_packaged_samplers

expected = sys.argv[1]
assert metadata.version("Jarvis-HEP") == jarvishep2.__version__ == expected
for module in (jarvishep2, jarvis_hep_sampler):
    assert Path(module.__file__).resolve().is_relative_to(Path(sys.prefix).resolve()), (
        module.__file__
    )
with tempfile.TemporaryDirectory() as tmp:
    base = Path(tmp)
    root = base / "project"
    (root / "bin").mkdir(parents=True)
    (root / ".jarvis-project.json").write_text("{}")
    (base / "home").mkdir()
    env = dict(
        os.environ, HOME=str(base / "home"), JARVIS_SAMPLERS_HOME=str(base / "registry")
    )
    # Project packaging reads the same isolated registry as the CLI subprocesses.
    os.environ["JARVIS_SAMPLERS_HOME"] = env["JARVIS_SAMPLERS_HOME"]
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    card = {
        "Scan": {"name": "release-smoke"},
        "EnvReqs": {
            "V2": {
                "workers": 2,
                "batch_size": 4,
                "redis": {"host": "127.0.0.1", "port": port, "db": 0},
            }
        },
        "Sampling": {
            "Method": "ReleaseSmoke",
            "Variables": [
                {
                    "name": n,
                    "distribution": {
                        "type": "Flat",
                        "parameters": {"min": 0.0, "max": 1.0},
                    },
                }
                for n in ("x", "y")
            ],
            "Bounds": {"seed": 7, "point_number": 16},
            "LogLikelihood": [{"name": "LogL_Z", "expression": "z"}],
        },
        "Operas": {
            "Modules": [
                {
                    "name": "Eggbox",
                    "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
                    "input": [{"name": n, "expression": n} for n in ("x", "y")],
                    "output": [{"name": "z", "entry": "z"}],
                }
            ]
        },
    }
    task = root / "bin/smoke.yaml"
    task.write_text(yaml.safe_dump(card))
    cli = [str(Path(sys.prefix) / "bin/Jarvis")]

    def run(args, cwd=root):
        result = subprocess.run(
            cli + args, cwd=cwd, env=env, text=True, capture_output=True, timeout=90
        )
        if result.returncode:
            raise SystemExit(
                f"Command failed: {args}\n" + result.stdout + result.stderr
            )

    for args in (
        ["sampler", "init", "ReleaseSmoke"],
        ["sampler", "check", "ReleaseSmoke.py", "--card", str(task)],
        ["sampler", "load", "ReleaseSmoke.py"],
        ["validate", str(task)],
        ["man", "sampler.ReleaseSmoke"],
        ["run", str(task), "--skip-draw-flowchart"],
    ):
        run(args)
    summary = json.loads((root / "outputs/release-smoke/run_summary.json").read_text())
    assert (
        summary["total_points_finished"] == 16 and summary["total_points_failed"] == 0
    ), summary
    assert (root / "outputs/release-smoke/DATABASE/sampler_identity.json").is_file()

    archive_path = create_project_package(str(root), "share").archive_path
    unpacked = base / "unpacked"
    with tarfile.open(archive_path) as archive:
        archive.extractall(unpacked, filter="data")
    fetched = unpacked / "project"
    asset = fetched / json.loads((fetched / INDEX).read_text())["entries"][0]["file"]
    # A conflict is visible; loading the packaged source is an explicit user step.
    source = root / "ReleaseSmoke.py"
    source.write_text(
        source.read_text().replace('"version": "0.1.0"', '"version": "9.0.0"')
    )
    run(["sampler", "load", str(source)])
    assert "different source" in inspect_packaged_samplers(str(fetched))[0]
    run(["sampler", "remove", "ReleaseSmoke", "--yes"])
    assert "not loaded" in inspect_packaged_samplers(str(fetched))[0]
    run(["sampler", "load", str(asset)])
    card["Scan"]["name"] = "fetched-smoke"
    (fetched / "bin/smoke.yaml").write_text(yaml.safe_dump(card))
    run(["validate", str(fetched / "bin/smoke.yaml")], cwd=fetched)
    run(["run", str(fetched / "bin/smoke.yaml"), "--skip-draw-flowchart"], cwd=fetched)
    summary = json.loads(
        (fetched / "outputs/fetched-smoke/run_summary.json").read_text()
    )
    assert (
        summary["total_points_finished"] == 16 and summary["total_points_failed"] == 0
    ), summary
print(
    "Complete HEP/SDK installed CLI, project reload and real Redis/spawn smoke passed."
)
