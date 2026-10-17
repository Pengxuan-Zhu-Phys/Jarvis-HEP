"""N5: merged registry report and the runtime's early registry gate."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from unittest import mock

import pytest
from jarvis_hep_sampler import SamplerRegistry

from jarvishep2.distributor import Distributor
from jarvishep2.runtime.core import Jarvis2Core
from jarvishep2.sampler_catalog import (
    SOURCE_USER,
    SamplerRegistryPreflightError,
    builtin_plugin_specs,
    builtin_specs,
    registry_report,
)


@pytest.fixture
def registry(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> SamplerRegistry:
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "samplers"))
    return SamplerRegistry()


def _load(registry: SamplerRegistry, tmp_path: Path, method: str = "External") -> None:
    source = tmp_path / f"{method}.py"
    source.write_text(
        f'PLUGIN = {{"api_version": 1, "method": "{method}", '
        '"version": "1", "profile": "independent"}\n'
        'def draw(context, rng, plan):\n'
        '    return rng.random(context.dimension)\n',
        encoding="utf-8",
    )
    registry.load(source)


def _rename_registration(registry: SamplerRegistry, original: str, renamed: str) -> None:
    raw = json.loads(registry.store_path.read_text(encoding="utf-8"))
    entry = raw["entries"].pop(original)
    entry["plugin"]["method"] = renamed
    raw["entries"][renamed] = entry
    registry.store_path.write_text(json.dumps(raw), encoding="utf-8")


def _core(tmp_path: Path, method: str) -> Jarvis2Core:
    return Jarvis2Core({
        "Sampling": {"Method": method},
        "task_root": str(tmp_path),
        "scan_name": "registry-test",
    })


def test_report_keeps_old_conflicts_visible_and_distributor_keeps_builtin_priority(
    registry: SamplerRegistry, tmp_path: Path,
) -> None:
    _load(registry, tmp_path)
    _rename_registration(registry, "External", "grid")
    report = registry_report()
    rows = report.matching("Grid")
    assert len(rows) == 2
    assert rows[0].source == "builtin"
    assert rows[1].source == SOURCE_USER
    assert rows[1].status == "name-conflict"
    assert rows[1].repair_command == "Jarvis sampler remove grid --yes"
    assert type(Distributor.set_method("Grid")).__name__ == "Grid"
    with pytest.raises(SamplerRegistryPreflightError, match="remove grid --yes"):
        _core(tmp_path, "Random").preflight_sampler_plugin()
    assert registry.remove("grid") is not None
    _core(tmp_path, "Random").preflight_sampler_plugin()


def test_report_rechecks_snapshot_without_registry_change_and_blocks_selected_plugin(
    registry: SamplerRegistry, tmp_path: Path,
) -> None:
    _load(registry, tmp_path)
    with mock.patch("jarvis_hep_sampler.loader.load_plugin_file", side_effect=AssertionError("plugin imported")):
        assert registry_report().for_method("External", source=SOURCE_USER).status == "ok"
    snapshot = registry.snapshot_path(registry.get("External"))
    snapshot.write_text(snapshot.read_text(encoding="utf-8") + "\n# changed\n", encoding="utf-8")
    broken = registry_report().for_method("External", source=SOURCE_USER)
    assert broken is not None and broken.status == "snapshot-broken"
    assert "SHA-256" in broken.problems[0]
    with mock.patch("jarvishep2.sampling.plugin_sampler.preflight_plugin") as plugin_preflight:
        with pytest.raises(SamplerRegistryPreflightError, match="snapshot is unusable"):
            _core(tmp_path, "External").preflight_sampler_plugin()
    plugin_preflight.assert_not_called()


def test_preflight_lists_all_conflicts_before_redis(
    registry: SamplerRegistry, tmp_path: Path,
) -> None:
    for old, collision in (("External", "Grid"), ("Another", "CSV")):
        _load(registry, tmp_path, old)
        _rename_registration(registry, old, collision)
    core = _core(tmp_path, "Random")
    with mock.patch.object(core, "init_logger"), \
         mock.patch.object(core, "check_environment_requirements"), \
         mock.patch.object(core, "init_redis") as init_redis:
        with pytest.raises(SamplerRegistryPreflightError) as caught:
            core.bootstrap_distributed_runtime()
    assert "remove Grid --yes" in str(caught.value)
    assert "remove CSV --yes" in str(caught.value)
    init_redis.assert_not_called()


def test_unrelated_broken_snapshot_warns_but_builtin_can_start(
    registry: SamplerRegistry, tmp_path: Path, caplog: pytest.LogCaptureFixture,
) -> None:
    _load(registry, tmp_path)
    registry.snapshot_path(registry.get("External")).unlink()
    with caplog.at_level(logging.WARNING):
        _core(tmp_path, "Random").preflight_sampler_plugin()
    assert "Unrelated sampler External" in caplog.text


def test_unreadable_registry_warns_for_builtin_and_blocks_plugin(
    registry: SamplerRegistry, tmp_path: Path, caplog: pytest.LogCaptureFixture,
) -> None:
    registry.home.mkdir(parents=True)
    registry.store_path.write_text("{broken", encoding="utf-8")
    report = registry_report()
    assert report.error is not None
    assert report.entries[0].status == "built-in"
    with caplog.at_level(logging.WARNING):
        _core(tmp_path, "Random").preflight_sampler_plugin()
    assert "Sampler registry cannot be read" in caplog.text
    with pytest.raises(SamplerRegistryPreflightError, match="Cannot read sampler registry"):
        _core(tmp_path, "External").preflight_sampler_plugin()


def test_sdk_reserved_names_cover_hep_builtins() -> None:
    from jarvis_hep_sampler.metadata import RESERVED_METHOD_NAMES

    expected = {spec.name.casefold() for spec in (*builtin_specs(), *builtin_plugin_specs())}
    actual = {name.casefold() for name in RESERVED_METHOD_NAMES}
    assert expected <= actual


@pytest.mark.parametrize(("registered", "query"), [("grid", "Grid"), ("Grid", "grid"), ("grid", "GRID")])
def test_builtin_info_reports_case_insensitive_legacy_collision(
    registry, tmp_path, capsys, registered, query,
):
    from jarvishep2.cli.sampler import dispatch_sampler
    _load(registry, tmp_path)
    _rename_registration(registry, "External", registered)
    assert dispatch_sampler(["info", query, "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["method"] == "Grid"
    assert payload["conflicts"][0]["method"] == registered
    assert payload["conflicts"][0]["repair_command"] == f"Jarvis sampler remove {registered} --yes"
    assert dispatch_sampler(["info", query]) == 0
    assert f"remove {registered} --yes" in capsys.readouterr().out


def test_hep_list_and_info_are_rendered_by_the_sdk(registry, tmp_path, capsys):
    from jarvis_hep_sampler import cli
    from jarvishep2.cli.sampler import dispatch_sampler
    with mock.patch.object(cli, "_render_list_section", wraps=cli._render_list_section) as listing:
        assert dispatch_sampler(["list"]) == 0
        assert listing.call_count == 2
    output = capsys.readouterr().out
    assert "No samplers loaded." in output and "Jarvis-HEP" in output
    with mock.patch.object(cli, "_render_info", wraps=cli._render_info) as info:
        assert dispatch_sampler(["info", "Grid"]) == 0
        info.assert_called_once()
