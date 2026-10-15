"""Native external-sampler workflows across independent CLI processes and real Redis."""

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

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid

pytestmark = [
    pytest.mark.slow,
    pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required"),
]

# Redirect the library's legacy cache constant before entering the actual CLI;
# leave HOME and Python's import paths unchanged for normal installed discovery.
CLI_BOOTSTRAP = """
import os
import jarvishep2.official_project_library as library
library.USER_CATALOG_CACHE_PATH = os.environ["JHS_TEST_CATALOG_CACHE"]
from jarvishep2.cli.client import main
raise SystemExit(main())
"""

SOURCE = """
import os
from pathlib import Path
import numpy as np
with Path(os.environ["JHS_TEST_MARKER"]).open("a", encoding="utf-8") as handle:
    handle.write(str(os.getpid()) + "\\n")
PLUGIN = {"api_version": 1, "method": "MySampler", "version": "VERSION",
          "profile": "independent", "description": "CLI lifecycle test sampler"}
VALUE = COORDINATE

def draw(context, rng, plan):
    return [VALUE] * context.dimension

def status(plan, stats):
    return {"ready": np.bool_(True), "coordinate": np.float32(VALUE)}
"""


class SamplerCLI:
    def __init__(self, root: Path):
        self.root = root
        root.mkdir(parents=True, exist_ok=True)
        self.logs = root / "cli-logs"
        self.logs.mkdir()
        self.registry = root / "registry"
        self.marker = root / "plugin-imports.txt"
        catalog = root / "catalog.json"
        catalog.write_text(json.dumps({"schema_version": 1, "projects": []}), encoding="utf-8")
        self.env = dict(
            os.environ,
            JARVIS_SAMPLERS_HOME=str(self.registry),
            JARVIS_SCAN_REF_REGISTRY=str(root / "scan-refs.json"),
            JARVIS_REDIS_NOTICE_STATE=str(root / "redis-notice"),
            JHS_TEST_CATALOG_CACHE=str(root / "cache/catalog.json"),
            JHS_TEST_MARKER=str(self.marker),
            JARVIS_OFFICIAL_LIBRARY_INDEX_URL=catalog.as_uri(),
            JARVIS_OPERAS_PERSIST_FILE=str(root / "operas/legacy.json"),
            JARVIS_OPERAS_SOURCES_FILE=str(root / "operas/sources.json"),
            JARVIS_OPERAS_OVERRIDES_FILE=str(root / "operas/overrides.json"),
            JARVIS_OPERAS_USER_SOURCE_DIR=str(root / "operas/user_sources"),
            JARVIS_OPERAS_CURVE_CACHE_ROOT=str(root / "cache/curves"),
            XDG_CACHE_HOME=str(root / "cache"),
            MPLCONFIGDIR=str(root / "cache/matplotlib"),
            PYTHONDONTWRITEBYTECODE="1",
        )
        self.env.pop("PYTHONPATH", None)

    def __call__(self, label, *args, expected=0, cwd=None):
        command = [sys.executable, "-c", CLI_BOOTSTRAP, *map(str, args)]
        result = subprocess.run(
            command, cwd=cwd or self.root, env=self.env, text=True,
            capture_output=True, timeout=90,
        )
        log = self.logs / f"{label}.log"
        rendered = result.stdout + result.stderr
        log.write_text(rendered, encoding="utf-8")
        assert result.returncode == expected, f"{log}\n{rendered[-8000:]}"
        return rendered


def _source(path, version, coordinate):
    payload = SOURCE.replace("VERSION", version).replace("COORDINATE", str(coordinate))
    path.write_text(payload, encoding="utf-8")
    return hashlib.sha256(payload.encode()).hexdigest()


def _free_port():
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def _project(root):
    (root / "bin").mkdir(parents=True)
    (root / ".jarvis-project.json").write_text("{}", encoding="utf-8")
    return root


def _card(project, scan, port):
    card = {
        "Scan": {"name": scan},
        "EnvReqs": {"V2": {
            "workers": 2, "batch_size": 2,
            "redis": {"host": "127.0.0.1", "port": port, "db": 0},
        }},
        "Sampling": {
            "Method": "MySampler",
            "Variables": [
                {"name": name, "distribution": {
                    "type": "Flat", "parameters": {"min": 0.0, "max": 1.0},
                }} for name in ("x", "y")
            ],
            "Bounds": {"seed": 7, "point_number": 6},
            "LogLikelihood": [{"name": "LogL_Z", "expression": "z"}],
        },
        "Operas": {"Modules": [{
            "name": "Eggbox", "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
            "input": [{"name": name, "expression": name} for name in ("x", "y")],
            "output": [{"name": "z", "entry": "z"}],
        }]},
    }
    path = project / "bin" / f"{scan}.yaml"
    path.write_text(yaml.safe_dump(card), encoding="utf-8")
    return path


def _scan_records(project, scan, version, digest, coordinate):
    database = project / "outputs" / scan / "DATABASE"
    binding = json.loads((database / "sampler_identity.json").read_text(encoding="utf-8"))
    assert binding["metadata"]["version"] == version
    assert binding["source_sha256"] == digest
    assert hashlib.sha256((project / binding["source_file"]).read_bytes()).hexdigest() == digest
    rows = sorted(
        SimpleHDF5Writer(str(database / "samples.hdf5")).read_records(),
        key=lambda row: int(row["sample_index"]),
    )
    assert len(rows) == 6
    prefix = f"plugin:MySampler:{binding['hard_identity'][:16]}"
    result = []
    for index, row in enumerate(rows):
        assert float(row["x"]) == float(row["y"]) == coordinate
        assert int(row["sample_index"]) == index
        assert row["uuid"] == deterministic_sampler_uuid(prefix=prefix, seed=7, sample_index=index)
        result.append((row["uuid"], float(row["x"]), float(row["y"])))
    assert len({item[0] for item in result}) == 6
    return result


def test_native_cli_load_update_remove_and_bound_resume(tmp_path):
    cli = SamplerCLI(tmp_path)
    project = _project(tmp_path / "native-project")
    source = tmp_path / "MySampler.py"
    digest_v1 = _source(source, "1.0", .25)
    port = _free_port()
    card_a = _card(project, "scan_a", port)

    # No sampler check is required to register or use a valid external file.
    cli("load-v1", "sampler", "load", source)
    imported = cli.marker.read_bytes()
    assert "MySampler" in cli("list", "sampler", "list")
    assert "1.0" in cli("info", "sampler", "info", "MySampler")
    assert "MySampler" in cli("man", "man", "sampler.MySampler")
    cli("validate", "validate", card_a)
    assert cli.marker.read_bytes() == imported

    cli("check-first", "check", card_a)
    cli("check-second", "check", card_a)
    cli("run-a", "run", card_a, "--skip-draw-flowchart")
    records_a = _scan_records(project, "scan_a", "1.0", digest_v1, .25)

    digest_v2 = _source(source, "2.0", .75)
    assert "source changed" in cli("changed-working-source", "sampler", "list")
    card_unloaded = _card(project, "before_reload", port)
    cli("run-without-reload", "run", card_unloaded, "--skip-draw-flowchart")
    _scan_records(project, "before_reload", "1.0", digest_v1, .25)

    cli("reload-v2", "sampler", "load", source)
    card_b = _card(project, "scan_b", port)
    cli("run-b", "run", card_b, "--skip-draw-flowchart")
    records_b = _scan_records(project, "scan_b", "2.0", digest_v2, .75)

    bad_source = "raise RuntimeError('bad reload must not replace v2')\n"
    source.write_text(bad_source, encoding="utf-8")
    cli("failed-reload", "sampler", "load", source, expected=2)
    assert "2.0" in cli("info-after-failure", "sampler", "info", "MySampler")
    card_kept = _card(project, "after_failed_reload", port)
    cli("run-after-failure", "run", card_kept, "--skip-draw-flowchart")
    _scan_records(project, "after_failed_reload", "2.0", digest_v2, .75)

    cli("remove", "sampler", "remove", "MySampler", "--yes")
    assert source.read_text(encoding="utf-8") == bad_source
    fresh_card = _card(project, "after_remove", port)
    cli("validate-after-remove", "validate", fresh_card, expected=2)
    for label, card, version, digest, coordinate, original in (
        ("a", card_a, "1.0", digest_v1, .25, records_a),
        ("b", card_b, "2.0", digest_v2, .75, records_b),
    ):
        cli(f"resume-{label}", "run", card, "--resume", "--skip-draw-flowchart")
        assert _scan_records(project, f"scan_{label}", version, digest, coordinate) == original


def test_cli_project_share_fetch_and_explicit_load(tmp_path):
    sender = SamplerCLI(tmp_path / "sender")
    project = _project(sender.root / "original")
    source = sender.root / "MySampler.py"
    digest = _source(source, "1.0", .25)
    _card(project, "shared_scan", _free_port())
    sender("load", "sampler", "load", source)
    imported = sender.marker.read_bytes()
    sender("pack", "project", "pack", project, "--share")
    assert sender.marker.read_bytes() == imported
    archives = list(sender.root.glob("*.tar.gz"))
    assert len(archives) == 1

    receiver = SamplerCLI(tmp_path / "receiver")
    catalog = receiver.root / "catalog.json"
    catalog.write_text(json.dumps({"schema_version": 1, "projects": [{
        "name": "SamplerDemo", "access": "public", "archive_url": archives[0].as_uri(),
        "archive_root": "original", "entrypoint": "bin/shared_scan.yaml",
    }]}), encoding="utf-8")
    output = receiver("fetch", "project", "fetch", "SamplerDemo")
    assert "Jarvis sampler load" in output
    assert not receiver.marker.exists()
    assert not (receiver.registry / "sources.json").exists()
    fetched = receiver.root / "SamplerDemo"
    index = json.loads((fetched / ".jarvis-pack/samplers.json").read_text(encoding="utf-8"))
    assert len(index["entries"]) == 1
    entry = index["entries"][0]
    assert entry["source_sha256"] == digest
    payload = fetched / entry["file"]
    assert hashlib.sha256(payload.read_bytes()).hexdigest() == digest
    card = fetched / "bin/shared_scan.yaml"
    receiver("validate-before-load", "validate", card, expected=2)
    assert not receiver.marker.exists()

    receiver("load-packaged-source", "sampler", "load", payload)
    imported = receiver.marker.read_bytes()
    receiver("validate-after-load", "validate", card)
    assert receiver.marker.read_bytes() == imported
    receiver("check", "check", card)
    receiver("run", "run", card, "--skip-draw-flowchart")
    _scan_records(fetched, "shared_scan", "1.0", digest, .25)
