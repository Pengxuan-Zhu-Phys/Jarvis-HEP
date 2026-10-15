"""Calculator smoke must never read or overwrite a formal plug-in scan."""

from __future__ import annotations

import copy
import hashlib
import os
import shutil
import socket
import subprocess
import sys
from pathlib import Path
from unittest import mock

import pytest
import yaml
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.errors import PluginConfigError, PluginResumeMismatch

from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampling.plugin_identity import identity_path
from jarvishep2.sampling.plugin_sampler import preflight_plugin
from jarvishep2.sampling.runtime_checkpoint import load_checkpoint


SOURCE = '''
PLUGIN = {"api_version": 1, "method": "CheckLayout", "version": "1",
          "profile": "independent"}
def prepare(setup):
    if setup.options.get("fail"):
        raise RuntimeError("deliberate smoke initialization failure")
def draw(context, rng, plan):
    return rng.random(context.dimension)
'''


@pytest.fixture(params=[False, True], ids=["default-output", "custom-output"])
def config(tmp_path, monkeypatch, request):
    root = tmp_path.resolve()
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(root / "registry"))
    source = root / "CheckLayout.py"
    source.write_text(SOURCE)
    SamplerRegistry().load(source)
    return {
        "Scan": {"name": "layout"},
        "scan_name": "layout",
        "task_root": str(root),
        "task_result_dir": str(root / ("custom-output" if request.param else "outputs/layout")),
        "Sampling": {
            "Method": "CheckLayout",
            "Bounds": {"point_number": 12, "seed": 7},
            "Variables": [{"name": "x", "distribution": {
                "type": "Flat", "parameters": {"min": 0.0, "max": 1.0},
            }}],
        },
        "EnvReqs": {"V2": {"workers": 1, "batch_size": 2}},
    }


def prepare(config, **kwargs):
    return preflight_plugin(config, task_root=config["task_root"],
                            scan_name=config["scan_name"], **kwargs)


def smoke_core(config):
    core = Jarvis2Core(copy.deepcopy(config))
    core._apply_check_modules_runtime_policy()
    return core


def checksums(files):
    return {path: hashlib.sha256(path.read_bytes()).hexdigest() for path in files}


def test_check_paths_do_not_create_formal_scan_state(config):
    core = smoke_core(config)
    prepared = prepare(core.config)
    result = Path(config["task_result_dir"])
    expected_checkpoint = (Path(config["task_root"]) / "checkpoints/layout/CheckLayout"
                           / "check/state.pkl")
    assert identity_path(core.config) == result / "DATABASE/test/sampler_identity.json"
    assert identity_path(core.config).is_file()
    assert not identity_path(config).exists()
    assert core._resolve_database_dir() == str(result / "DATABASE/test")
    assert Path(core.checkpoint_file()) == expected_checkpoint
    assert prepared.source_copy.parent == expected_checkpoint.parent / "plugin"
    assert not (expected_checkpoint.parent.parent / "state.pkl").exists()
    # A new formal run is still fresh after the calculator smoke preflight.
    assert prepare(config).source_copy.parent == expected_checkpoint.parent.parent / "plugin"


def test_repeated_checks_preserve_existing_formal_scan(config):
    formal = prepare(config)
    checkpoint = formal.source_copy.parent.parent / "state.pkl"
    checkpoint.write_bytes(b"formal checkpoint must stay untouched")
    database = identity_path(config).parent / "samples.hdf5"
    database.write_bytes(b"formal database must stay untouched")
    protected = [identity_path(config), formal.source_copy, checkpoint, database]
    before = checksums(protected)
    for seed in (8, 9):
        core = smoke_core(config)
        core.config["Sampling"]["Bounds"]["seed"] = seed
        prepared = prepare(core.config)
        # Smoke rows deliberately survive successive checks; they are not a
        # request to resume a scientific scan with its previous configuration.
        (identity_path(core.config).parent / "samples.hdf5").write_bytes(b"smoke rows")
        Path(core.checkpoint_file()).write_bytes(b"smoke checkpoint")
        assert prepared.identity != formal.identity
        assert checksums(protected) == before
    assert prepare(config, resume=True).identity == formal.identity
    changed = copy.deepcopy(config)
    changed["Sampling"]["Bounds"]["seed"] = 8
    with pytest.raises(PluginResumeMismatch, match="science settings differ"):
        prepare(changed, resume=True)
    assert checksums(protected) == before


def test_failed_check_cannot_rebind_formal_source(config):
    formal = prepare(config)
    before = checksums([identity_path(config), formal.source_copy])
    core = smoke_core(config)
    core.config["Sampling"]["Bounds"]["fail"] = True
    with pytest.raises(PluginConfigError, match="deliberate smoke initialization failure"):
        prepare(core.config)
    assert not identity_path(core.config).exists()
    assert checksums(list(before)) == before


def test_interrupt_saves_the_durable_plugin_cursor(config):
    core = Jarvis2Core(copy.deepcopy(config))
    core.preflight_sampler_plugin()
    sampler = core.init_sampler_from_config()
    core.start_runtime_checkpoint()
    try:
        for _ in range(2):
            assert sampler.propose_next() is not None
        sampler.request_checkpoint_capture()
        assert sampler.capture_checkpoint_at_batch_boundary()
        # The Archiver acknowledges this complete batch; the next two
        # proposals remain in flight when the interrupt arrives.
        core._persisted_index_prefix = 2
        sampler.set_persisted_prefix(2)
        for _ in range(2):
            assert sampler.propose_next() is not None
        path = core.save_runtime_checkpoint(reason="interrupt")
        saved = load_checkpoint(path)["sampler_state"]
        assert saved["accepted_index"] == 2
        assert saved["index"] == 2
        # The safe snapshot can be imported and reconciled under the same
        # strict source/identity/prefix checks used by a real resume.
        sampler.import_runtime_state(saved)
        assert sampler.advance_to_persisted_prefix(2) == 2
    finally:
        sampler.shutdown_checkpointing()


@pytest.mark.parametrize("reason", ["interrupt", "dynesty_engine_checkpoint", "dynesty_finished"])
def test_nonindexed_engine_keeps_live_checkpoint_behavior(config, reason):
    class Engine:
        uses_indexed_resume = False

        def export_runtime_state(self):
            return {"engine_state": "live"}

        def at_safe_barrier(self):
            return False

        def checkpoint_runtime_state(self, *, safe):
            return {"engine_state": "live" if safe else "old"}

    core = Jarvis2Core(copy.deepcopy(config))
    core.sampler = Engine()
    core.info["sampler_name"] = "FeedbackEngine"
    saved = load_checkpoint(core.save_runtime_checkpoint(reason=reason))
    assert saved["sampler_state"]["engine_state"] == "live"


@pytest.mark.parametrize("indexed", [True, False], ids=["indexed", "feedback"])
def test_bootstrap_keeps_indexed_cursor_after_reconciliation(config, indexed):
    from jarvishep2.runtime import _runtime_supervisor
    from jarvishep2.sampling.plugin_sampler import PluginSampler

    prepared = prepare(config)
    sampler = PluginSampler(prepared.spec)
    sampler.attach_preflight(prepared)
    sampler.set_config(config)
    sampler.uses_indexed_resume = indexed
    old_state = sampler.export_runtime_state()
    core = mock.MagicMock()
    core.config = config
    core.info = {"scan_name": config["scan_name"]}
    core._environment_summary = {}
    core._plugin_preflight = prepared
    core._resume_policy = "resume"
    core.sampler = sampler
    core._persisted_index_prefix = core._persisted_records_count = 2
    core._persisted_failed_count = 0
    core._persisted_uuids = set()
    core._resolve_database_dir.return_value = str(identity_path(config).parent)
    core.init_sampler_from_config.side_effect = lambda: sampler.import_runtime_state(old_state)
    core._import_sampler_checkpoint_state.side_effect = lambda target: target.import_runtime_state(old_state)
    supervisor = _runtime_supervisor._RuntimeSupervisor(core)
    supervisor._publish_core_proc_board = mock.Mock()
    with (
        mock.patch("jarvishep2.process_cleanup.ensure_scan_name_available"),
        mock.patch.object(_runtime_supervisor, "build_command_parser"),
        mock.patch.object(_runtime_supervisor, "prepare_hdf5_database_for_writer"),
    ):
        supervisor.bootstrap_distributed_runtime()
    if indexed:
        assert sampler._accepted_index == 2
        core._import_sampler_checkpoint_state.assert_not_called()
    else:
        core._import_sampler_checkpoint_state.assert_called_once_with(sampler)


@pytest.mark.slow
@pytest.mark.skipif(shutil.which("redis-server") is None, reason="real Redis is required")
@pytest.mark.parametrize("custom_output", [False, True], ids=["cli", "custom-output-api"])
def test_real_checks_then_run_preserve_science(tmp_path, custom_output):
    """Run the CLI lifecycle; custom result paths use the native Python API."""
    from jarvishep2.database import SimpleHDF5Writer

    root = tmp_path.resolve()
    (root / ".jarvis-project.json").write_text("{}")
    (root / "bin").mkdir()
    source = root / "CheckLayout.py"
    source.write_text(SOURCE)
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        port = sock.getsockname()[1]
    task = root / "bin/task.yaml"
    card = {
        "Scan": {"name": "check-layout"},
        "EnvReqs": {"V2": {
            "workers": 1, "batch_size": 2,
            "redis": {"host": "127.0.0.1", "port": port, "db": 0},
            "check_modules": {"n_samples": 2},
        }},
        "Sampling": {
            "Method": "CheckLayout", "Bounds": {"point_number": 12, "seed": 7},
            "Variables": [{"name": name, "distribution": {
                "type": "Flat", "parameters": {"min": 0.0, "max": 1.0},
            }} for name in ("x", "y")],
        },
        "Operas": {"Modules": [{
            "name": "Eggbox", "operator": "jarvishep2.testing.eggbox.eggbox2d_numpy",
            "input": [{"name": name, "expression": name} for name in ("x", "y")],
            "output": [{"name": "z", "entry": "z"}],
        }]},
    }
    task.write_text(yaml.safe_dump(card))
    output = root / ("custom-results" if custom_output else "outputs/check-layout")
    env = dict(os.environ, JARVIS_SAMPLERS_HOME=str(root / "registry"),
               JARVIS_OPERAS_CURVE_CACHE_ROOT=str(root / "curve-cache"),
               PYTHONDONTWRITEBYTECODE="1", JHS_CHECK_LAYOUT_OUTPUT=str(output))
    cli = [sys.executable, "-c",
           "from jarvishep2.cli.client import main; raise SystemExit(main())"]
    api = [sys.executable, "-c", '''
import os
import sys
from jarvishep2.runtime.core import Jarvis2Core
core = Jarvis2Core()
core.load_task_yaml(sys.argv[2])
core.config["task_result_dir"] = os.environ["JHS_CHECK_LAYOUT_OUTPUT"]
core.config["skip_draw_flowchart"] = True
core._populate_info_from_config()
outcome = core.check_modules(timeout=15) if sys.argv[1] == "check" else core.run()
raise SystemExit(outcome.exit_code)
''']

    def run(label, command, *args):
        process = subprocess.run(command + list(args), cwd=root, env=env,
                                 text=True, capture_output=True, timeout=60)
        text = process.stdout + process.stderr
        (root / f"{label}.log").write_text(text)
        assert process.returncode == 0, text[-5000:]
        with socket.socket() as sock:
            assert sock.connect_ex(("127.0.0.1", port)) != 0

    def execute(label, operation):
        if custom_output:
            run(label, api, operation, str(task))
        elif operation == "check":
            run(label, cli, operation, str(task), "--timeout", "15")
        else:
            run(label, cli, operation, str(task), "--skip-draw-flowchart")

    run("load", cli, "sampler", "load", str(source))
    for index in (1, 2):
        execute(f"check-{index}", "check")
        assert not (output / "DATABASE/sampler_identity.json").exists()
        assert (output / "DATABASE/test/sampler_identity.json").is_file()
    smoke_database = output / "DATABASE/test/samples.hdf5"
    assert len(SimpleHDF5Writer(str(smoke_database)).read_records()) == 4
    execute("run", "run")
    database = output / "DATABASE/samples.hdf5"
    assert len(SimpleHDF5Writer(str(database)).read_records()) == 12
    checkpoint_dir = root / "checkpoints/check-layout/CheckLayout"
    # Even an opaque formal checkpoint must never be opened by a calculator
    # smoke. It belongs only to the formal resume path.
    checkpoint = checkpoint_dir / "state.pkl"
    checkpoint.write_bytes(b"formal checkpoint must not be read or overwritten by check")
    protected = [p for p in (output / "DATABASE").iterdir() if p.is_file()]
    protected += list((checkpoint_dir / "plugin").glob("*.py")) + [checkpoint]
    before = checksums(protected)
    card["Sampling"]["Bounds"]["seed"] = 8
    task.write_text(yaml.safe_dump(card))
    execute("check-after-run", "check")
    assert checksums(protected) == before
    assert (checkpoint_dir / "check/plugin").is_dir()
