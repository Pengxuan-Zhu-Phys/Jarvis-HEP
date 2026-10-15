"""Jarvis-HEP requires a Jarvis-HEP-Sampler release newer than 1."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from jarvishep2.sampler_sdk import SDK_SPECIFIER, SamplerSDKError, require_sampler_sdk


def test_installed_sdk_is_newer_than_v1() -> None:
    import jarvis_hep_sampler
    from importlib.metadata import version

    installed = require_sampler_sdk()
    assert installed == version("Jarvis-HEP-Sampler") == jarvis_hep_sampler.__version__
    assert SDK_SPECIFIER == ">1"


@pytest.mark.parametrize("installed", ["0.1.0", "1.0.0", "1"])
def test_sdk_requirement_rejects_v1_and_older(monkeypatch: pytest.MonkeyPatch, installed: str) -> None:
    monkeypatch.setattr("jarvishep2.sampler_sdk.version", lambda name: installed)
    with pytest.raises(SamplerSDKError, match=SDK_SPECIFIER):
        require_sampler_sdk()


@pytest.mark.parametrize("installed", ["1.0.1", "1.2.0", "2.0.0"])
def test_sdk_requirement_accepts_anything_newer_than_v1(
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


def test_default_environments_require_sampler_newer_than_v1() -> None:
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
        assert match["version"] == ">1"
