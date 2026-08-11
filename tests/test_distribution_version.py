from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from unittest import mock

from jarvishep2.project_packager import _jarvis_version
from jarvishep2.versioning import get_runtime_version


def _installed_version(name: str) -> str:
    versions = {
        "Jarvis-HEP": "1.7.5",
        "jarvishep2": "2.0.0",
    }
    try:
        return versions[name]
    except KeyError as exc:
        raise PackageNotFoundError(name) from exc


def test_runtime_version_ignores_installed_v1_distribution() -> None:
    with mock.patch("importlib.metadata.version", side_effect=_installed_version):
        assert get_runtime_version() == "2.0.0"


def test_project_packager_version_ignores_installed_v1_distribution() -> None:
    with mock.patch(
        "jarvishep2.project_packager.importlib_metadata.version",
        side_effect=_installed_version,
    ):
        assert _jarvis_version() == "2.0.0"


def test_runtime_version_prefers_jarvis_hep_v2_distribution() -> None:
    def installed_version(name: str) -> str:
        return {"Jarvis-HEP": "2.1.0", "jarvishep2": "2.0.0"}[name]

    with mock.patch("importlib.metadata.version", side_effect=installed_version):
        assert get_runtime_version() == "2.1.0"
