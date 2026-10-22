"""Jarvis-HEP requires the SDK release that externalizes LatinHypercube (1.1.2)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from jarvishep2.sampler_sdk import SDK_SPECIFIER, SamplerSDKError, require_sampler_sdk


def test_installed_sdk_meets_minimum() -> None:
    import jarvis_hep_sampler
    from importlib.metadata import version

    installed = require_sampler_sdk()
    assert installed == version("Jarvis-HEP-Sampler") == jarvis_hep_sampler.__version__
    assert SDK_SPECIFIER == ">=1.1.2"


@pytest.mark.parametrize("installed", ["0.1.0", "1.0.0", "1", "1.0.2", "1.0.3", "1.1.0", "1.1.1", "1.1.2rc1"])
def test_sdk_requirement_rejects_older_releases(monkeypatch: pytest.MonkeyPatch, installed: str) -> None:
    monkeypatch.setattr("jarvishep2.sampler_sdk.version", lambda name: installed)
    with pytest.raises(SamplerSDKError, match=SDK_SPECIFIER):
        require_sampler_sdk()


@pytest.mark.parametrize("installed", ["1.1.2", "1.1.3", "1.2.0", "2.0.0"])
def test_sdk_requirement_accepts_minimum_and_newer(
    monkeypatch: pytest.MonkeyPatch, installed: str
) -> None:
    monkeypatch.setattr("jarvishep2.sampler_sdk.version", lambda name: installed)
    assert require_sampler_sdk() == installed


def test_user_registry_reports_an_unsupported_sdk(monkeypatch: pytest.MonkeyPatch) -> None:
    from jarvishep2.sampler_catalog import user_registry_error

    monkeypatch.setattr("jarvishep2.sampler_sdk.version", lambda name: "0.1.0")
    problem = user_registry_error()
    assert problem is not None
    assert "0.1.0" in problem and SDK_SPECIFIER in problem


def test_default_environments_require_sampler_minimum() -> None:
    root = Path(__file__).resolve().parents[1] / "jarvishep2"
    paths = (
        root / "card" / "environment_default.yaml",
        root / "project_template" / "deps" / "environment_default.yaml",
    )
    for path in paths:
        document = yaml.safe_load(path.read_text(encoding="utf-8"))
        dependencies = document["EnvReqs"]["Python"]["Dependencies"]
        match = next(item for item in dependencies if item["name"] == "Jarvis-HEP-Sampler")
        assert match["required"] is True
        assert match["version"] == SDK_SPECIFIER
        packaging = next(item for item in dependencies if item["name"] == "packaging")
        assert packaging["required"] is True
        assert packaging["version"] == ">=23.0"


def test_published_sdk_without_new_host_contract_gets_actionable_error(monkeypatch):
    import jarvis_hep_sampler
    monkeypatch.setattr(jarvis_hep_sampler, "HOST_CONTRACT_VERSION", 0)
    with pytest.raises(SamplerSDKError, match="host contract 2"):
        require_sampler_sdk()
