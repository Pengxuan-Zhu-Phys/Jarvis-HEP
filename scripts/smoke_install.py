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
import time

import yaml
import jarvishep2
import jarvis_hep_sampler
import jarvis_operas
from jarvishep2.project_packager import create_project_package
from jarvishep2.project_samplers import INDEX, inspect_packaged_samplers
from jarvishep2.sampler_sdk import require_sampler_sdk

expected = sys.argv[1]
assert metadata.version("Jarvis-HEP") == jarvishep2.__version__ == expected
assert require_sampler_sdk() == jarvis_hep_sampler.__version__ == metadata.version("Jarvis-HEP-Sampler")
for module in (jarvishep2, jarvis_hep_sampler, jarvis_operas):
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
    for key, directory in (("JARVIS_SCAN_REF_REGISTRY", "scan-refs.json"),
                           ("JARVIS_OPERAS_CURVE_CACHE_ROOT", "operas-cache"),
                           ("MPLCONFIGDIR", "matplotlib")):
        env[key] = str(base / directory)
        os.environ[key] = env[key]
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

    def command(args):
        return subprocess.run(cli + args, cwd=root, env=env, text=True,
                              capture_output=True, timeout=90)

    def run(args, cwd=root):
        result = subprocess.run(cli + args, cwd=cwd, env=env, text=True,
                                capture_output=True, timeout=90)
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
    ):
        run(args)

    bad = json.loads(json.dumps(card))
    bad["Sampling"]["LogLikelihood"][0]["expression"] = "nosuch.f(z)"
    bad_path = root / "bin/bad.yaml"
    bad_path.write_text(yaml.safe_dump(bad))
    result = command(["validate", str(bad_path), "--strict", "--json"])
    assert result.returncode != 0
    payload = json.loads(result.stdout)
    assert not payload["ok"] and any(i["code"] == "JV2-EXP-001" for i in payload["issues"])
    result = command(["run", str(bad_path), "--skip-draw-flowchart"])
    assert result.returncode != 0 and "nosuch.f(z)" in result.stderr + result.stdout
    result = command(["plot", str(task)])
    assert result.returncode == 2 and "gen-plot-yaml" in result.stderr
    for args in (["--monitor", "validate", str(task)], ["validate", str(task), "--monitor"]):
        result = command(args)
        assert result.returncode == 2 and "conflicting CLI intents" in result.stderr
    result = command(["operas", "init", "-h"])
    assert result.returncode == 0 and "Force recompilation" in result.stdout
    assert "<rich." not in result.stdout

    # Leave enough work for two real monitor commands to observe the live scan.
    card["Operas"]["Modules"][0]["operator"] = "jarvishep2.testing.resume.slow_sum"
    card["Operas"]["Modules"][0]["input"].append({"name": "delay", "expression": "1.0"})
    task.write_text(yaml.safe_dump(card))
    from jarvishep2.process_cleanup import list_active_scans
    from jarvishep2.monitor.scans import worker_inventory
    with (base / "live-smoke.log").open("w") as log:
        live = subprocess.Popen(cli + ["run", str(task), "--skip-draw-flowchart"],
                                cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT)
        try:
            deadline = time.monotonic() + 60
            scan = None
            while live.poll() is None and time.monotonic() < deadline:
                scan = next((s for s in list_active_scans() if s.name == "release-smoke"), None)
                if scan and len(worker_inventory(scan)) == 2:
                    break
                time.sleep(0.1)
            assert scan and len(worker_inventory(scan)) == 2, "No two-worker inventory"
            for selector in (scan.reference, scan.name):
                result = command(["monitor", "--once", "--json", selector])
                assert result.returncode == 0, result.stderr + result.stdout
                view = json.loads(result.stdout)
                assert view["ok"] and view["factory"]["workers_alive"] == 2
                assert {w["pid"] for w in view["workers"]} == {w["pid"] for w in worker_inventory(scan)}
                assert abs(time.time() - view["timestamp"]) < 10
            assert live.wait(timeout=90) == 0, (base / "live-smoke.log").read_text()
        finally:
            if live.poll() is None:
                live.terminate()
                try:
                    live.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    live.kill()
                    live.wait(timeout=5)
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
    card["Operas"]["Modules"][0]["operator"] = "jarvishep2.testing.eggbox.eggbox2d_numpy"
    card["Operas"]["Modules"][0]["input"].pop()
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
