"""Local version inventory, missing packages, and terminal output contracts."""
from __future__ import annotations

import json
import sys
from importlib import metadata

import pytest

from jarvishep2.cli.client import dispatch_version, main
from jarvishep2.versioning import ecosystem_versions, get_runtime_version


def test_ecosystem_missing_package_uses_null(monkeypatch):
    def installed(name):
        if name == "Jarvis-Operas":
            raise metadata.PackageNotFoundError(name)
        return "9.8.7"
    monkeypatch.setattr(metadata, "version", installed)
    rows = ecosystem_versions()
    assert [r["name"] for r in rows] == [
        "Jarvis-HEP", "Jarvis-HEP-Sampler", "Jarvis-Operas", "Jarvis-PLOT", "Jarvis-HEP-Portal",
    ]
    assert rows[0]["version"] == get_runtime_version()
    assert next(r for r in rows if r["name"] == "Jarvis-Operas")["version"] is None


@pytest.mark.parametrize("is_tty,no_color", [(False, ""), (True, "1")])
def test_root_report_keeps_banner_and_aligns_ecosystem(monkeypatch, capsys, is_tty, no_color):
    monkeypatch.setenv("NO_COLOR", no_color)
    monkeypatch.setattr(sys.stdout, "isatty", lambda: is_tty)
    monkeypatch.setattr("jarvishep2.versioning.ecosystem_versions", lambda: [
        {"name": "Jarvis-HEP", "version": get_runtime_version()},
        {"name": "Jarvis-Operas", "version": None},
    ])
    assert dispatch_version() == 0
    output = capsys.readouterr().out
    assert "Author:" in output and "Version:" in output
    assert output.index("Author:") < output.index("Package")
    assert "Jarvis-HEP" in output and "not installed" in output
    assert "\x1b[" not in output


@pytest.mark.parametrize("tokens", [["--version", "--json"], ["--json", "-v"]])
def test_root_json_versions(tokens, capsys):
    assert main(tokens) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["name"] == "Jarvis-HEP"
    assert payload["version"] == get_runtime_version()
    assert [r["name"] for r in payload["ecosystem"]] == [r["name"] for r in ecosystem_versions()]
