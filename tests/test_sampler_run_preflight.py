"""The real Jarvis run/check entrypoints reject ambiguous sampler registries."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from unittest import mock

import pytest
from jarvis_hep_sampler.registry import RegistryEntry, SamplerRegistry

from jarvishep2.cli.client import main
from jarvishep2.cli.sampler import dispatch_sampler
from jarvishep2.logging import shutdown_jarvis_logging


def _legacy_grid_registration(registry: SamplerRegistry, tmp_path: Path, method: str = "Grid") -> None:
    source = tmp_path / f"{method}.py"
    source.write_text(
        'PLUGIN = {"api_version": 1, "method": ' + repr(method) + ', "version": "1", '
        '"profile": "independent"}\n'
        'def draw(context, rng, plan): return [0.25] * context.dimension\n'
    )
    payload = source.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    snapshot = f"user_sources/Grid-{digest[:12]}.py"
    target = registry.home / snapshot
    target.parent.mkdir(parents=True)
    target.write_bytes(payload)
    entry = RegistryEntry(
        method=method,
        source_path=str(source),
        snapshot=snapshot,
        source_sha256=digest,
        loaded_at="",
        sdk_version="legacy",
        plugin={"api_version": 1, "method": method, "version": "1", "profile": "independent"},
        options_schema=None,
        hooks=("draw",),
    )
    registry.store_path.write_text(json.dumps({"version": 1, "entries": {method: entry.to_dict()}}))


@pytest.mark.parametrize("command", ["run", "check"])
def test_jarvis_rejects_legacy_collision_before_redis(tmp_path, monkeypatch, capsys, command):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    registry = SamplerRegistry()
    _legacy_grid_registration(registry, tmp_path)
    card = tmp_path / "task.yaml"
    card.write_text(
        (Path(__file__).parent / "parity_project/random_opera.yaml").read_text()
        .replace("name: random-opera", "name: legacy-collision")
    )
    (tmp_path / ".jarvis-project.json").write_text("{}")

    try:
        with mock.patch("jarvishep2.runtime.core.Jarvis2Core.init_redis") as init_redis:
            assert main([command, str(card)]) == 2
            init_redis.assert_not_called()
    finally:
        shutdown_jarvis_logging()
    stderr = capsys.readouterr().err
    assert "Grid" in stderr
    assert "Jarvis sampler remove Grid --yes" in stderr


@pytest.mark.parametrize("method", ["Grid", "grid"])
def test_sampler_list_info_show_legacy_collision_and_remove_repairs_it(
    tmp_path, monkeypatch, capsys, method,
):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    _legacy_grid_registration(SamplerRegistry(), tmp_path, method)

    assert dispatch_sampler(["list", "--json"]) == 0
    rows = json.loads(capsys.readouterr().out)["samplers"]
    grid = next(row for row in rows if row["method"] == method)
    assert grid["status"] == "name-conflict"
    assert grid["repair_command"] == f"Jarvis sampler remove {method} --yes"

    assert dispatch_sampler(["info", method, "--json"]) == 0
    info = json.loads(capsys.readouterr().out)
    assert info["status"] == "name-conflict"
    assert info["repair_command"] == f"Jarvis sampler remove {method} --yes"

    assert dispatch_sampler(["list"]) == 0
    listing = capsys.readouterr().out
    assert "name-conflict" in listing
    assert f"Jarvis sampler remove {method} --yes" in listing
    assert dispatch_sampler(["info", method]) == 0
    detail = capsys.readouterr().out
    assert detail.count("STATUS") == 1
    assert "name-conflict" in detail
    assert f"Jarvis sampler remove {method} --yes" in detail

    assert dispatch_sampler(["remove", method, "--yes", "--json"]) == 0
    capsys.readouterr()
    assert dispatch_sampler(["info", method, "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["status"] == "built-in"
