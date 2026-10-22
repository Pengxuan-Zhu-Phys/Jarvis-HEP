"""HEP-side boundaries: fresh sampling, resume, and SDK CLI forwarding."""
from __future__ import annotations

import copy
import json
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path
from unittest import mock

import pytest
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.cli import main as sdk_sampler_main

from jarvishep2.cli.sampler import dispatch_sampler
from jarvishep2.database import SimpleHDF5Writer
from jarvishep2.distributor import Distributor
from jarvishep2.redis_queue import make_fakeredis_queue
from jarvishep2.runtime._resume_service import _ResumeService
from jarvishep2.runtime._runtime_supervisor import _RuntimeSupervisor


PLUGIN = '''
PLUGIN = {"api_version": 1, "method": "ReviewDemo", "version": "1",
          "profile": "independent"}
def draw(context, rng, plan):
    return [0.25] * context.dimension
'''


@pytest.fixture
def registry(tmp_path, monkeypatch):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    return SamplerRegistry()


def plugin_file(tmp_path, method="ReviewDemo"):
    path = tmp_path / f"{method}.py"
    path.write_text(PLUGIN.replace('"ReviewDemo"', repr(method)))
    return path


@pytest.mark.parametrize("method", ["Grid", "GRID", "Random", "RANDOM"])
def test_hep_load_rejects_reserved_names_without_registry_changes(tmp_path, registry, capsys, method):
    source = plugin_file(tmp_path, method)
    before = registry.entries()
    assert dispatch_sampler(["load", str(source), "--json"]) == 2
    payload = json.loads(capsys.readouterr().out)
    assert payload["ok"] is False
    assert registry.entries() == before


def test_hep_load_option_separator_cannot_bypass_reserved_names(tmp_path, registry, capsys):
    source = plugin_file(tmp_path, "Grid")
    assert dispatch_sampler(["load", "--json", "--", str(source)]) == 2
    assert json.loads(capsys.readouterr().out)["ok"] is False
    assert registry.entries() == {}


@pytest.mark.parametrize("operation", ["info", "check", "list"])
@pytest.mark.parametrize("as_json", [False, True])
def test_hep_registry_errors_do_not_escape(tmp_path, registry, capsys, operation, as_json):
    registry.load(plugin_file(tmp_path))
    (registry.home / "sources.json").write_text("{broken")
    args = [operation] + ([] if operation == "list" else ["ReviewDemo"])
    assert dispatch_sampler(args + (["--json"] if as_json else [])) == 2
    captured = capsys.readouterr()
    assert "Traceback" not in captured.out + captured.err
    if as_json:
        assert json.loads(captured.out)["ok"] is False


def test_hep_builtin_info_rejects_a_damaged_registry(registry, capsys):
    registry.home.mkdir(parents=True)
    registry.store_path.write_text("{broken", encoding="utf-8")
    assert dispatch_sampler(["info", "Grid", "--json"]) == 2
    payload = json.loads(capsys.readouterr().out)
    assert payload["command"] == "sampler info"
    assert payload["error"]["type"] == "RegistryError"


def test_hep_sampler_sdk_import_failure_uses_json_fallback(capsys):
    with mock.patch.dict(sys.modules, {"jarvis_hep_sampler.cli": None}):
        assert dispatch_sampler(["list", "--json"]) == 2
    payload = json.loads(capsys.readouterr().out)
    assert payload["command"] == "sampler list"
    assert payload["error"]["type"] == "ModuleNotFoundError"


def test_hep_sampler_sdk_version_error_uses_shared_text_panel(capsys):
    from jarvishep2.sampler_sdk import SamplerSDKError

    with mock.patch(
        "jarvishep2.sampler_sdk.require_sampler_sdk",
        side_effect=SamplerSDKError("SDK version is unavailable"),
    ):
        assert dispatch_sampler(["list"]) == 2
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "Jarvis sampler" in captured.err
    assert "SDK version is unavailable" in captured.err
    assert "[Jarvis]" not in captured.err


@pytest.mark.parametrize("dimension", ["0", "-1"])
def test_hep_check_rejects_nonpositive_dimension(tmp_path, registry, capsys, dimension):
    source = plugin_file(tmp_path)
    with mock.patch("jarvis_hep_sampler.testing.check_plugin") as checker:
        assert dispatch_sampler(["check", str(source), "--dimension", dimension, "--json"]) == 2
        checker.assert_not_called()
    assert json.loads(capsys.readouterr().out)["ok"] is False


def test_hep_check_rejects_negative_seed_before_loading(tmp_path, registry, capsys):
    source = plugin_file(tmp_path)
    with mock.patch("jarvis_hep_sampler.testing.check_plugin") as checker:
        assert dispatch_sampler(["check", str(source), "--seed", "-1", "--json"]) == 2
        checker.assert_not_called()
    payload = json.loads(capsys.readouterr().out)
    assert payload["command"] == "sampler check"
    assert payload["error"]["type"] == "ValueError"


@pytest.mark.parametrize("command", ["load", "check"])
def test_sdk_and_hep_sampler_forwarding_agree_on_reserved_file(
    tmp_path, registry, capsys, command,
):
    source = plugin_file(tmp_path, "Grid")
    args = [command, str(source), "--json"]
    sdk_code = sdk_sampler_main(args)
    sdk_payload = json.loads(capsys.readouterr().out)
    hep_code = dispatch_sampler(args)
    hep_payload = json.loads(capsys.readouterr().out)
    assert (sdk_code, hep_code) == (2, 2)
    assert sdk_payload["ok"] is hep_payload["ok"] is False
    assert sdk_payload["error"]["type"] == hep_payload["error"]["type"]
    assert hep_payload["command"] == f"sampler {command}"
    assert registry.entries() == {}


def test_hep_sampler_uses_sdk_help_with_list_info_extensions(capsys):
    assert sdk_sampler_main(["check", "-h"]) == 0
    sdk_help = capsys.readouterr().out
    assert dispatch_sampler(["check", "-h"]) == 0
    hep_help = capsys.readouterr().out
    assert hep_help.replace("Jarvis sampler", "jsampler") == sdk_help
    assert dispatch_sampler(["-h"]) == 0
    root_help = capsys.readouterr().out
    assert "List loaded samplers and methods shipped" in root_help
    assert "Show one loaded sampler or method shipped" in root_help
    assert root_help.count("with Jarvis-HEP") == 2


@pytest.mark.parametrize("method", ["Grid", "Random", "Bridson", "CSV"])
@pytest.mark.parametrize("policy", ["fresh", "resume", "rejected-checkpoint", "check"])
def test_durable_rows_filter_only_resume(tmp_path, monkeypatch, method, policy):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    csv = tmp_path / "points.csv"
    csv.write_text("uuid,x,y\nrow-0,0.1,0.2\nrow-1,0.3,0.4\nrow-2,0.5,0.6\n")
    config = {
        "Scan": {"name": "review"}, "task_root": str(tmp_path),
        "Runtime": {"workers": 1, "batch_size": 2, "max_inflight": 64,
                    "checkpoint": {"enabled": False}},
        "Sampling": {"Method": method,
            "Bounds": {"seed": 7, "point_number": 4, "radius": 0.4, "max_attempt": 5},
            "Variables": [{"name": n, "distribution": {"type": "Flat",
                "parameters": {"min": 0., "max": 1., "num": 2, "length": 1.}}}
                for n in ("x", "y")]},
    }
    if method == "CSV":
        config["Sampling"]["Bounds"] = {"path": str(csv)}
    if policy == "check":
        config["_check_modules_sample_layout"] = True
    sampler = Distributor.set_method(method)
    sampler.set_config(config)
    if hasattr(sampler, "initialize"):
        sampler.initialize()
    if hasattr(sampler, "_max_inflight"):
        sampler._max_inflight = 10000
    reference = Distributor.set_method(method)
    reference.set_config(config)
    if hasattr(reference, "initialize"):
        reference.initialize()
    indices = []
    while (sample := reference.propose_next()) is not None:
        indices.append(sample.sample_index)
    database = tmp_path / "DATABASE"
    database.mkdir()
    SimpleHDF5Writer(str(database / "samples.hdf5")).add_records([
        {"uuid": "old-0", "sample_index": 0, "status": "Completed", "z": -999},
        {"uuid": "old-2", "sample_index": 2, "status": "Completed", "z": -999},
        {"uuid": "row-2", "status": "Completed", "z": -999},
    ])
    core = mock.MagicMock()
    core.config = copy.deepcopy(config)
    core.info = {"scan_name": "review"}
    core._environment_summary = {}
    core.sampler = sampler
    core._resume_policy = "resume" if policy == "resume" else "fresh"
    core._persisted_index_prefix = core._persisted_records_count = 0
    core._persisted_failed_count = 0
    core._persisted_pending_indices = set()
    core._persisted_uuids = set()
    core._resolve_database_dir.return_value = str(database)
    service = _ResumeService(core)
    core._load_persisted_database_state.side_effect = service._load_persisted_database_state
    if policy == "rejected-checkpoint":
        core._resume_policy = "resume"
        checkpoint = tmp_path / "state.pkl"
        checkpoint.write_bytes(b"old")
        core.checkpoint_file.return_value = str(checkpoint)
        service._discard_resume_checkpoint(reason="changed card")
        assert not checkpoint.exists()
    queue = make_fakeredis_queue()
    sampler.set_redis(queue)
    supervisor = _RuntimeSupervisor(core)
    supervisor._publish_core_proc_board = mock.Mock()
    try:
        with mock.patch("jarvishep2.process_cleanup.ensure_scan_name_available"), \
             mock.patch("jarvishep2.runtime._runtime_supervisor.build_command_parser"), \
             mock.patch("jarvishep2.runtime._runtime_supervisor.prepare_hdf5_database_for_writer"):
            supervisor.bootstrap_distributed_runtime()
        pushed = sampler.run_distributed()
        actual = [queue.pull_task(timeout=1)["sample_index"] for _ in range(pushed)]
        expected = [i for i in indices if i not in {0, 2}] if policy == "resume" else indices
        assert actual == expected
    finally:
        sampler.shutdown_checkpointing()
        queue.close()


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required")
@pytest.mark.parametrize("method", ["Grid", "Random", "Bridson", "CSV"])
def test_fresh_runtime_replaces_old_results(tmp_path, monkeypatch, method):
    """Exercise real Redis, spawn Workers and HDF5 replacement, not only submission."""
    import yaml
    from jarvishep2.task_config import load_task_yaml
    from jarvishep2.testing.eggbox import eggbox2d_numpy

    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    card = yaml.safe_load((Path(__file__).parent / "parity_project/random_opera.yaml").read_text())
    card["Scan"]["name"] = f"fresh_{method}_{tmp_path.name}"
    card["Sampling"]["Method"] = method
    card["Sampling"]["Bounds"] = {
        "Grid": {"seed": 11},
        "Random": {"point_number": 6, "seed": 11},
        "Bridson": {"radius": 0.4, "max_attempt": 5, "seed": 11},
        "CSV": {},
    }[method]
    for variable in card["Sampling"]["Variables"]:
        variable["distribution"]["parameters"]["num"] = 2
    if method == "CSV":
        csv = tmp_path / "points.csv"
        csv.write_text("uuid,x,y\nrow-0,0.1,0.2\nrow-1,0.3,0.4\nrow-2,0.5,0.6\n")
        card["Sampling"]["Bounds"] = {"path": str(csv)}
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    card["EnvReqs"]["V2"]["redis"] = {"host": "127.0.0.1", "port": port, "db": 0}
    source = tmp_path / "task.yaml"
    source.write_text(yaml.safe_dump(card))
    (tmp_path / ".jarvis-project.json").write_text("{}")
    sampler = Distributor.set_method(method)
    sampler.set_config(load_task_yaml(str(source)))
    stale = []
    while (sample := sampler.propose_next()) is not None:
        stale.append({"uuid": sample.uuid, "sample_index": sample.sample_index,
                      "status": "Completed", "x": 0., "y": 0., "z": -999.})
    database = tmp_path / f"outputs/{card['Scan']['name']}/DATABASE/samples.hdf5"
    SimpleHDF5Writer(str(database)).add_records(stale)
    output = tmp_path / "outcome.json"
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1",
               JARVIS_SCAN_REF_REGISTRY=str(tmp_path / "scan-refs.json"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(tmp_path / "curve-cache"))
    driver = '''
import json, sys
from pathlib import Path
from jarvishep2.runtime.core import Jarvis2Core
core = Jarvis2Core()
core.load_task_yaml(sys.argv[1])
core.config.update(skip_draw_flowchart=True, skip_render_plots=True)
outcome = core.run(write_run_summary=False)
Path(sys.argv[2]).write_text(json.dumps(outcome.to_dict()))
raise SystemExit(outcome.exit_code)
'''
    result = subprocess.run([sys.executable, "-c", driver, str(source), str(output)],
                            cwd=tmp_path, env=env, text=True, capture_output=True, timeout=60)
    (tmp_path / "runtime.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads(output.read_text())["completed"] == len(stale)
    rows = SimpleHDF5Writer(str(database)).read_records()
    assert len(rows) == len(stale)
    assert {row["uuid"] for row in rows} == {row["uuid"] for row in stale}
    for row in rows:
        assert float(row["z"]) != -999.
        assert float(row["z"]) == pytest.approx(eggbox2d_numpy(float(row["x"]), float(row["y"]))["z"])


def test_old_sdk_host_contract_has_one_json_dependency_error(monkeypatch, capsys):
    import jarvis_hep_sampler
    from jarvishep2.cli.sampler import dispatch_sampler
    monkeypatch.setattr(jarvis_hep_sampler, "HOST_CONTRACT_VERSION", 0)
    assert dispatch_sampler(["list", "--json"]) == 2
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert payload["ok"] is False
    assert payload["error"]["type"] == "SamplerSDKError"
    assert "host contract 3" in payload["error"]["message"]
    assert "Traceback" not in captured.err
