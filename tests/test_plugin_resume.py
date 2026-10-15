"""S2a: scan binding and refusal are independent of current user registration."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.errors import PluginConfigError, PluginResumeMismatch
from jarvis_hep_sampler.loader import load_plugin_file

from jarvishep2 import sampler_catalog
from jarvishep2.sampling.plugin_identity import (
    binding_record,
    identity_path,
    read_binding,
)
from jarvishep2.sampling.plugin_sampler import (
    PluginSampler,
    plugin_identity,
    preflight_plugin,
)

from test_sampler_plugins import HALF_CUBE, _config


@pytest.fixture
def setup(tmp_path, monkeypatch):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    source = tmp_path / "HalfCube.py"
    source.write_text(HALF_CUBE)
    SamplerRegistry().load(source)
    cfg = _config(tmp_path, "HalfCube", {"point_number": 16, "seed": 7, "upper": 0.8})
    return cfg, source


def prepare(cfg, **kwargs):
    return preflight_plugin(
        cfg, task_root=cfg["task_root"], scan_name=cfg["scan_name"], **kwargs
    )


def digest_tree(root):
    return {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in root.rglob("*")
        if p.is_file()
    }


@pytest.mark.parametrize(
    "field", ["seed", "point_number", "options", "Variables", "Mapper", "selection"]
)
def test_science_drift_refuses_without_mutating_scan(setup, field):
    cfg, _ = setup
    prepare(cfg)
    scan = Path(cfg["task_root"])
    before = digest_tree(scan)
    changed = copy.deepcopy(cfg)
    if field in {"seed", "point_number"}:
        changed["Sampling"]["Bounds"][field] += 1
    elif field == "options":
        changed["Sampling"]["Bounds"]["upper"] = 0.4
    elif field == "Variables":
        changed["Sampling"][field].reverse()
    elif field == "Mapper":
        changed["Sampling"][field] = [{"name": "mass", "expression": "2*x"}]
    else:
        changed["Sampling"][field] = ["x > .2"]
    with pytest.raises(PluginResumeMismatch, match="science settings differ"):
        prepare(changed, resume=True)
    assert digest_tree(scan) == before


def test_block_size_and_generation_mode_are_hard_identity(setup):
    cfg, source = setup
    plugin = load_plugin_file(source)
    kwargs = dict(config=cfg, seed=7, point_number=16, options={"upper": 0.8})
    from dataclasses import replace

    meta4 = replace(plugin.metadata, block_size=4)
    meta8 = replace(plugin.metadata, block_size=8)
    one = plugin_identity(metadata=meta4, generation_mode="draw_block", **kwargs)[0]
    assert (
        one
        != plugin_identity(metadata=meta8, generation_mode="draw_block", **kwargs)[0]
    )
    assert one != plugin_identity(metadata=meta4, generation_mode="draw", **kwargs)[0]


def test_failed_initialization_cannot_publish_a_binding(setup):
    cfg, _ = setup
    cfg["Sampling"]["Bounds"]["upper"] = 2
    with pytest.raises(PluginConfigError):
        prepare(cfg)
    assert not identity_path(cfg).exists()
    assert not list((Path(cfg["task_root"]) / "checkpoints").rglob("*.py"))


def test_reload_remove_and_newer_orphan_do_not_change_bound_source(setup):
    cfg, source = setup
    first = prepare(cfg)
    source.write_text(HALF_CUBE.replace('"version": "1.0"', '"version": "2.0"'))
    SamplerRegistry().load(source)
    orphan = first.source_copy.parent / "HalfCube-newer.py"
    orphan.write_text("raise RuntimeError('orphan must never be imported')")
    resumed = prepare(cfg, resume=True)
    assert resumed.plugin.source_sha256 == first.plugin.source_sha256
    SamplerRegistry().remove("HalfCube")
    with sampler_catalog.scan_plugin_context(binding_record(cfg)):
        assert sampler_catalog.is_stateless("HalfCube")
        assert sampler_catalog.plugin_record("HalfCube").version == "1.0"
        assert prepare(cfg, resume=True).identity == first.identity
    assert sampler_catalog.spec_for("HalfCube") is None


def test_corrupt_source_is_refused_before_import(setup):
    cfg, _ = setup
    first = prepare(cfg)
    first.source_copy.write_text("raise AssertionError('must not execute')")
    before = digest_tree(Path(cfg["task_root"]))
    with pytest.raises(PluginResumeMismatch, match="SHA-256 mismatch"):
        prepare(cfg, resume=True)
    assert digest_tree(Path(cfg["task_root"])) == before


def test_missing_copy_only_repairs_from_identical_snapshot(setup):
    cfg, source = setup
    first = prepare(cfg)
    first.source_copy.unlink()
    resumed = prepare(cfg, resume=True)
    assert resumed.source_copy.read_bytes() == first.plugin.source_bytes
    resumed.source_copy.unlink()
    source.write_text(HALF_CUBE + "\n# changed\n")
    SamplerRegistry().load(source)
    with pytest.raises(PluginResumeMismatch, match="SHA-256 mismatch"):
        prepare(cfg, resume=True)
    assert not resumed.source_copy.exists()


@pytest.mark.parametrize("kind", ["missing", "corrupt"])
def test_database_without_a_valid_binding_is_refused(setup, kind):
    cfg, _ = setup
    first = prepare(cfg)
    (identity_path(cfg).parent / "samples.hdf5").write_bytes(b"durable data")
    if kind == "missing":
        identity_path(cfg).unlink()
    else:
        identity_path(cfg).write_text('{"broken": true}')
    with pytest.raises(PluginResumeMismatch):
        prepare(cfg, resume=True)
    assert first.source_copy.exists()


@pytest.mark.parametrize(
    "key,value",
    [
        ("maxp", 8),
        ("seed", 0),
        ("index", 100),
        ("accepted_index", -1),
        ("selectionexp", ["x > .2"]),
        ("plugin", ["invalid"]),
    ],
)
def test_checkpoint_cannot_overwrite_prepared_configuration(setup, key, value):
    cfg, _ = setup
    preflight = prepare(cfg)
    sampler = PluginSampler(preflight.spec)
    sampler.attach_preflight(preflight)
    sampler.set_config(cfg)
    sampler.propose_next()
    before = sampler.export_runtime_state()
    state = copy.deepcopy(before)
    state[key] = value
    with pytest.raises(PluginResumeMismatch):
        sampler.import_runtime_state(state)
    assert sampler.export_runtime_state() == before


def test_unreachable_durable_prefix_is_a_failure(setup):
    cfg, _ = setup
    preflight = prepare(cfg)
    sampler = PluginSampler(preflight.spec)
    sampler.attach_preflight(preflight)
    sampler.set_config(cfg)
    with pytest.raises(PluginResumeMismatch, match="cannot reach"):
        sampler.advance_to_persisted_prefix(17)


def test_environment_warning_is_soft_but_identity_record_is_unchanged(setup, caplog):
    cfg, _ = setup
    first = prepare(cfg)
    data = read_binding(cfg)
    data["environment"]["python"] = "3.10.0"
    identity_path(cfg).write_text(json.dumps(data))
    before = identity_path(cfg).read_bytes()
    assert prepare(cfg, resume=True).identity == first.identity
    assert "environment changed" in caplog.text
    assert identity_path(cfg).read_bytes() == before
