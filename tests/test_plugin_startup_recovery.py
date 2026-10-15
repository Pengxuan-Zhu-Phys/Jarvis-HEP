"""A published plug-in binding survives startup failure and supports native resume."""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path
from unittest import mock

import psutil
import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.errors import PluginConfigError, PluginResumeMismatch

from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.plugin_identity import identity_path
from jarvishep2.sampling.stateless_batch import deterministic_sampler_uuid

SOURCE = '''
import os
from pathlib import Path
marker = os.environ.get("JHS_STARTUP_IMPORTS")
if marker:
    with Path(marker).open("a", encoding="utf-8") as handle:
        handle.write(str(os.getpid()) + "\\n")
PLUGIN = {"api_version": 1, "method": "StartupBound", "version": "1.0",
          "profile": "independent"}
def prepare(setup):
    if os.environ.get("JHS_STARTUP_FAIL_PREPARE"):
        raise RuntimeError("injected prepare failure before source binding")
def draw(context, rng, plan):
    return [.25] * context.dimension
'''

CLI_BOOTSTRAP = '''
import json
import os
from pathlib import Path
import psutil
import jarvishep2.official_project_library as library
library.USER_CATALOG_CACHE_PATH = os.environ["JHS_STARTUP_CATALOG_CACHE"]
if os.environ.get("JHS_STARTUP_FAIL_AFTER_WORKERS"):
    from jarvishep2.runtime._runtime_supervisor import _RuntimeSupervisor
    def fail_archiver(self, *args, **kwargs):
        core = self._core
        managed = core._managed_redis
        workers = [{"pid": worker.pid, "created": psutil.Process(worker.pid).create_time(),
                    "alive": worker.is_alive()} for worker in core.factory.workers]
        redis = {"pid": managed.process.pid,
                 "created": psutil.Process(managed.process.pid).create_time(),
                 "alive": managed.process.poll() is None,
                 "started_by_us": managed.started_by_us, "port": managed.port}
        Path(os.environ["JHS_STARTUP_RESOURCES"]).write_text(
            json.dumps({"workers": workers, "redis": redis}), encoding="utf-8")
        raise RuntimeError("injected archiver startup failure after Worker spawn")
    _RuntimeSupervisor.init_archiver = fail_archiver
from jarvishep2.cli.client import main
raise SystemExit(main())
'''


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def card(root, port=6380):
    scan = "startup_" + hashlib.sha256(str(root.resolve()).encode()).hexdigest()[:12]
    return {
        "Scan": {"name": scan},
        "task_root": str(root),
        "task_result_dir": str(root / "outputs" / scan),
        "scan_name": scan,
        "Sampling": {
            "Method": "StartupBound", "Bounds": {"seed": 7, "point_number": 8},
            "Variables": [{"name": name, "distribution": {
                "type": "Flat", "parameters": {"min": 0.0, "max": 1.0},
            }} for name in ("x", "y")],
        },
        "EnvReqs": {"V2": {
            "workers": 2, "batch_size": 2,
            "redis": {"host": "127.0.0.1", "port": port, "db": 0},
        }},
        "Operas": {"Modules": [{
            "name": "Eggbox", "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
            "input": [{"name": name, "expression": name} for name in ("x", "y")],
            "output": [{"name": "z", "entry": "z"}],
        }]},
    }


@pytest.fixture
def registered(tmp_path, monkeypatch):
    registry_home = tmp_path / "registry"
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(registry_home))
    monkeypatch.setenv("JARVIS_SCAN_REF_REGISTRY", str(tmp_path / "scan-refs.json"))
    monkeypatch.setenv("JHS_STARTUP_IMPORTS", str(tmp_path / "imports.txt"))
    # Fast tests do not turn the pytest process into an OS-visible running scan.
    monkeypatch.setattr("jarvishep2.proc_title.set_process_title", lambda _: None)
    source = tmp_path / "StartupBound.py"
    source.write_text(SOURCE, encoding="utf-8")
    SamplerRegistry(registry_home).load(source)
    return card(tmp_path), source


def forbid_resources(core):
    core.init_redis = mock.Mock(side_effect=AssertionError("Redis must not start"))
    core.init_factory = mock.Mock(side_effect=AssertionError("Workers must not start"))
    core.init_archiver = mock.Mock(side_effect=AssertionError("Archiver must not start"))


def assert_no_science(config):
    database = identity_path(config).parent
    assert not list(database.glob("*.hdf5"))
    assert not list(database.glob("*.h5"))
    assert not list(database.glob("*.csv"))


def published_failure(config):
    core = Jarvis2Core(copy.deepcopy(config))
    core.check_environment_requirements = mock.Mock()
    core.prepare_libraries = mock.Mock()
    core.init_command_parser = mock.Mock()
    forbid_resources(core)
    core.init_redis.side_effect = RuntimeError("injected Redis startup failure after binding")
    with pytest.raises(RuntimeError, match="injected Redis startup failure"):
        core.run()
    core.init_redis.assert_called_once()
    core.init_factory.assert_not_called()
    core.init_archiver.assert_not_called()
    assert core.redis is core.factory is core.archiver is core._managed_redis is None
    assert core._shutdown_done
    binding = json.loads(identity_path(config).read_text(encoding="utf-8"))
    bound_source = Path(config["task_root"]) / binding["source_file"]
    assert binding["source_sha256"] == digest(bound_source)
    assert_no_science(config)
    return binding, bound_source


def test_prepare_failure_publishes_no_binding_or_runtime_resources(registered, monkeypatch):
    config, _ = registered
    monkeypatch.setenv("JHS_STARTUP_FAIL_PREPARE", "1")
    core = Jarvis2Core(copy.deepcopy(config))
    forbid_resources(core)
    with pytest.raises(PluginConfigError, match="injected prepare failure"):
        core.run()
    assert not identity_path(config).exists()
    assert not list(Path(config["task_root"]).glob("checkpoints/**/*.py"))
    assert_no_science(config)
    core.init_redis.assert_not_called()
    core.init_factory.assert_not_called()
    core.init_archiver.assert_not_called()


def test_startup_failure_binding_requires_explicit_resume_and_keeps_its_source(registered):
    config, source = registered
    binding, bound_source = published_failure(config)
    before = {path: path.read_bytes() for path in (identity_path(config), bound_source)}
    retry = Jarvis2Core(copy.deepcopy(config))
    forbid_resources(retry)
    with pytest.raises(PluginResumeMismatch, match="use --resume"):
        retry.run()
    retry.init_redis.assert_not_called()
    source.write_text(SOURCE.replace('"1.0"', '"2.0"').replace("[.25]", "[.75]"))
    registry = SamplerRegistry()
    registry.load(source)
    registry.remove("StartupBound")
    resumed = Jarvis2Core(copy.deepcopy(config))
    forbid_resources(resumed)
    resumed.bootstrap_distributed_runtime = mock.Mock(
        side_effect=RuntimeError("verified binding reached resource bootstrap")
    )
    with pytest.raises(RuntimeError, match="verified binding reached resource bootstrap"):
        resumed.run(resume=True)
    assert resumed._plugin_preflight.plugin.metadata.version == "1.0"
    assert resumed._plugin_preflight.plugin.source_sha256 == binding["source_sha256"]
    assert resumed._plugin_preflight.resumed_from_copy
    assert {path: path.read_bytes() for path in before} == before
    resumed.init_redis.assert_not_called()
    assert_no_science(config)


@pytest.mark.parametrize("corruption", ["science", "source"])
def test_startup_retry_rejects_drift_before_resource_initialization(registered, corruption):
    config, _ = registered
    _, bound_source = published_failure(config)
    changed = copy.deepcopy(config)
    if corruption == "science":
        changed["Sampling"]["Bounds"]["seed"] = 8
        message = "science settings differ"
    else:
        bound_source.write_text("raise AssertionError('damaged source must not execute')\n")
        message = "SHA-256 mismatch"
    before = {path: path.read_bytes() for path in (identity_path(config), bound_source)}
    imports = (Path(config["task_root"]) / "imports.txt").read_bytes()
    core = Jarvis2Core(changed)
    forbid_resources(core)
    with pytest.raises(PluginResumeMismatch, match=message):
        core.run(resume=True)
    core.init_redis.assert_not_called()
    core.init_factory.assert_not_called()
    core.init_archiver.assert_not_called()
    assert {path: path.read_bytes() for path in before} == before
    assert (Path(config["task_root"]) / "imports.txt").read_bytes() == imports
    assert_no_science(config)


def same_process_is_alive(record):
    try:
        process = psutil.Process(record["pid"])
        return process.create_time() == record["created"] and process.is_running()
    except psutil.NoSuchProcess:
        return False


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required")
def test_cli_worker_startup_failure_cleans_resources_then_bound_resume_succeeds(tmp_path):
    root = tmp_path.resolve()
    (root / ".jarvis-project.json").write_text("{}", encoding="utf-8")
    (root / "bin").mkdir()
    source = root / "StartupBound.py"
    source.write_text(SOURCE, encoding="utf-8")
    catalog = root / "catalog.json"
    catalog.write_text('{"schema_version": 1, "projects": []}', encoding="utf-8")
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    config = card(root, port)
    yaml_card = {key: value for key, value in config.items()
                 if key not in {"task_root", "task_result_dir", "scan_name"}}
    task = root / "bin/task.yaml"
    task.write_text(yaml.safe_dump(yaml_card), encoding="utf-8")
    imports_path = root / "imports.txt"
    resources_path = root / "startup-resources.json"
    env = dict(
        os.environ, JARVIS_SAMPLERS_HOME=str(root / "registry"),
        JARVIS_SCAN_REF_REGISTRY=str(root / "scan-refs.json"),
        JARVIS_REDIS_NOTICE_STATE=str(root / "redis-notice"),
        JARVIS_OFFICIAL_LIBRARY_INDEX_URL=catalog.as_uri(),
        JHS_STARTUP_CATALOG_CACHE=str(root / "cache/catalog.json"),
        JHS_STARTUP_IMPORTS=str(imports_path), JHS_STARTUP_RESOURCES=str(resources_path),
        JARVIS_OPERAS_PERSIST_FILE=str(root / "operas/legacy.json"),
        JARVIS_OPERAS_SOURCES_FILE=str(root / "operas/sources.json"),
        JARVIS_OPERAS_OVERRIDES_FILE=str(root / "operas/overrides.json"),
        JARVIS_OPERAS_USER_SOURCE_DIR=str(root / "operas/user_sources"),
        JARVIS_OPERAS_CURVE_CACHE_ROOT=str(root / "cache/curves"),
        XDG_CACHE_HOME=str(root / "cache"), MPLCONFIGDIR=str(root / "cache/matplotlib"),
        PYTHONDONTWRITEBYTECODE="1",
    )
    env.pop("PYTHONPATH", None)
    exits = {}

    def cli(label, *args, expected=0, fail_startup=False):
        call_env = dict(env)
        if fail_startup:
            call_env["JHS_STARTUP_FAIL_AFTER_WORKERS"] = "1"
        process = subprocess.run(
            [sys.executable, "-c", CLI_BOOTSTRAP, *map(str, args)],
            cwd=root, env=call_env, text=True, capture_output=True, timeout=90,
        )
        text = process.stdout + process.stderr
        (root / f"{label}.log").write_text(text, encoding="utf-8")
        exits[label] = process.returncode
        assert process.returncode == expected, text[-8000:]
        return text

    def port_is_closed():
        with socket.socket() as sock:
            return sock.connect_ex(("127.0.0.1", port)) != 0

    cli("load-v1", "sampler", "load", source)
    failed = cli("failed-startup", "run", task, "--skip-draw-flowchart",
                 expected=1, fail_startup=True)
    assert "injected archiver startup failure after Worker spawn" in failed
    resources = json.loads(resources_path.read_text(encoding="utf-8"))
    assert resources["redis"]["alive"] and resources["redis"]["started_by_us"]
    assert len(resources["workers"]) == 2 and all(row["alive"] for row in resources["workers"])
    assert all(not same_process_is_alive(row) for row in resources["workers"])
    assert not same_process_is_alive(resources["redis"])
    assert port_is_closed()
    assert_no_science(config)
    binding_path = identity_path(config)
    binding = json.loads(binding_path.read_text(encoding="utf-8"))
    bound_source = root / binding["source_file"]
    before = {path: path.read_bytes() for path in (binding_path, bound_source)}
    assert digest(bound_source) == binding["source_sha256"] == digest(source)
    fresh = cli("fresh-retry-refused", "run", task, "--skip-draw-flowchart", expected=1)
    assert "use --resume" in fresh and "existing scan" in fresh
    assert port_is_closed()
    source.write_text(SOURCE.replace('"1.0"', '"2.0"').replace("[.25]", "[.75]"))
    cli("reload-v2", "sampler", "load", source)
    cli("remove", "sampler", "remove", "StartupBound", "--yes")
    assert {path: path.read_bytes() for path in before} == before

    def refuse(label, message):
        imports_before = imports_path.read_bytes()
        text = cli(label, "run", task, "--resume", "--skip-draw-flowchart", expected=1)
        assert message in text
        assert imports_path.read_bytes() == imports_before
        assert port_is_closed()
        assert_no_science(config)

    yaml_card["Sampling"]["Bounds"]["seed"] = 8
    task.write_text(yaml.safe_dump(yaml_card), encoding="utf-8")
    refuse("science-drift-refused", "science settings differ")
    yaml_card["Sampling"]["Bounds"]["seed"] = 7
    task.write_text(yaml.safe_dump(yaml_card), encoding="utf-8")
    bound_source.write_text("raise AssertionError('damaged source must not execute')\n")
    refuse("source-corruption-refused", "SHA-256 mismatch")
    bound_source.write_bytes(before[bound_source])
    cli("resume-v1", "run", task, "--resume", "--skip-draw-flowchart")
    assert port_is_closed()
    assert {path: path.read_bytes() for path in before} == before
    rows = sorted(SimpleHDF5Writer(str(binding_path.parent / "samples.hdf5")).read_records(),
                  key=lambda row: int(row["sample_index"]))
    assert len(rows) == 8 and len({row["uuid"] for row in rows}) == 8
    prefix = f"plugin:StartupBound:{binding['hard_identity'][:16]}"
    for index, row in enumerate(rows):
        assert float(row["x"]) == float(row["y"]) == .25
        assert int(row["sample_index"]) == index
        assert row["uuid"] == deterministic_sampler_uuid(prefix=prefix, seed=7, sample_index=index)
    assert source.exists()
    evidence = {
        "scenario": "published binding, actual Redis/Worker startup failure, native resume",
        "exit_codes": exits, "resources_at_failure": resources,
        "workers_cleaned": True, "redis_cleaned": True,
        "scientific_rows_after_failure": 0, "scientific_rows_after_resume": len(rows),
        "bound_version": binding["metadata"]["version"],
        "binding_sha256": digest(binding_path), "source_sha256": digest(bound_source),
        "uuid_and_coordinate_reference": "constant .25, indexed UUID seed 7, 8 candidates",
        "logs_sha256": {path.name: digest(path) for path in root.glob("*.log")},
    }
    (root / "verification.json").write_text(json.dumps(evidence, indent=2) + "\n")
