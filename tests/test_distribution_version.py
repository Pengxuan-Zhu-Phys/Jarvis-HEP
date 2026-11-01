from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from unittest import mock

from jarvishep2 import __version__ as package_version
from jarvishep2.project_packager import _jarvis_version
from jarvishep2.versioning import get_runtime_version


def _installed_version(name: str) -> str:
    versions = {
        "Jarvis-HEP": "1.7.5",
        "jarvishep2": "0.0.0",
    }
    try:
        return versions[name]
    except KeyError as exc:
        raise PackageNotFoundError(name) from exc


def test_package_constant_wins_over_stale_metadata() -> None:
    assert get_runtime_version() == package_version
    with mock.patch("importlib.metadata.version", side_effect=_installed_version):
        assert get_runtime_version() == package_version
        assert _jarvis_version() == package_version


def test_without_constant_runtime_version_uses_installed_v2(subtests) -> None:
    cases = (
        # An installed V1 Jarvis-HEP is ignored.
        ({"Jarvis-HEP": "1.7.5", "jarvishep2": "2.0.0"}, "2.0.0"),
        ({"Jarvis-HEP": "2.1.0", "jarvishep2": "2.0.0"}, "2.1.0"),
    )
    for installed, expected in cases:
        with (
            subtests.test(installed=installed),
            mock.patch(
                "jarvishep2.versioning._read_package_version_constant",
                return_value=None,
            ),
            mock.patch("importlib.metadata.version", side_effect=installed.__getitem__),
        ):
            assert get_runtime_version() == expected

