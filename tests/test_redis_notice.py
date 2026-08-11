"""First-run Redis installation guidance for the Jarvis CLI."""

from __future__ import annotations

import io
import tempfile
from pathlib import Path
from unittest import mock

from jarvishep2.redis_notice import (
    emit_first_run_redis_notice,
    find_redis_server,
    redis_install_command,
)


def test_find_redis_server_uses_redis_server_on_path() -> None:
    with mock.patch("jarvishep2.redis_notice.shutil.which", return_value="/usr/bin/redis-server") as which:
        assert find_redis_server() == "/usr/bin/redis-server"
    which.assert_called_once_with("redis-server")


def test_find_redis_server_accepts_redis6_and_valkey_binary_names() -> None:
    with mock.patch(
        "jarvishep2.redis_notice.shutil.which",
        side_effect=[None, "/usr/bin/redis6-server"],
    ) as which:
        assert find_redis_server() == "/usr/bin/redis6-server"
    assert which.call_args_list == [
        mock.call("redis-server"),
        mock.call("redis6-server"),
    ]


def test_install_commands_cover_macos_and_debian() -> None:
    assert redis_install_command(system="Darwin") == "brew install redis"
    assert redis_install_command(
        system="Linux", distro={"ID": "ubuntu"}, which=lambda _: None
    ) == "sudo apt install -y redis-server"


def test_install_commands_cover_common_linux_families() -> None:
    cases = {
        "fedora": "sudo dnf install -y valkey",
        "rhel": "sudo dnf install -y redis",
        "amzn": "sudo dnf install -y valkey",
        "arch": "sudo pacman -S valkey",
        "manjaro": "sudo pacman -S valkey",
        "opensuse-tumbleweed": "sudo zypper install -y redis",
        "alpine": "sudo apk add redis",
        "gentoo": "sudo emerge --ask dev-db/redis",
        "void": "sudo xbps-install -S redis",
        "nixos": "nix profile install nixpkgs#redis",
        "solus": "sudo eopkg install redis",
        "mageia": "sudo urpmi redis",
    }
    for distro_id, expected in cases.items():
        assert redis_install_command(
            system="Linux", distro={"ID": distro_id}, which=lambda _: "/usr/bin/dnf"
        ) == expected
    assert redis_install_command(
        system="Linux",
        distro={"ID": "amzn", "VERSION_ID": "2"},
        which=lambda _: None,
    ) == "sudo amazon-linux-extras install redis6 -y"


def test_install_commands_cover_bsd_hosts() -> None:
    assert redis_install_command(system="FreeBSD") == "sudo pkg install redis"
    assert redis_install_command(system="OpenBSD") == "doas pkg_add redis"
    assert redis_install_command(
        system="Linux", distro={"ID": "debian"}, which=lambda _: None
    ) == "sudo apt install -y redis-server"


def test_linux_fallback_uses_available_package_manager() -> None:
    assert redis_install_command(
        system="Linux",
        distro={"ID": "unknown", "ID_LIKE": ""},
        which=lambda name: "/usr/bin/dnf" if name == "dnf" else None,
    ) == "sudo dnf install -y redis"
    assert redis_install_command(
        system="Linux",
        distro={"ID": "unknown", "ID_LIKE": ""},
        which=lambda name: "/usr/bin/apk" if name == "apk" else None,
    ) == "sudo apk add redis"


def test_notice_is_emitted_once_and_keeps_normal_output_stream_clean() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        state_path = Path(tmp) / "redis-notice.marker"
        output = io.StringIO()
        with (
            mock.patch.dict("os.environ", {"JARVIS_REDIS_NOTICE_STATE": str(state_path)}),
            mock.patch("jarvishep2.redis_notice.find_redis_server", return_value=None),
            mock.patch("jarvishep2.redis_notice.platform.system", return_value="Darwin"),
        ):
            assert emit_first_run_redis_notice(stream=output)
            assert not emit_first_run_redis_notice(stream=output)

        text = output.getvalue()
        assert "Redis server was not found" in text
        assert "brew install redis" in text
        assert text.count("Redis server was not found") == 1
        assert state_path.is_file()


def test_cli_emits_notice_after_the_normal_command_returns() -> None:
    events: list[str] = []

    def display() -> int:
        events.append("display")
        return 0

    def notice() -> bool:
        events.append("notice")
        return True

    from jarvishep2.client import main

    with (
        mock.patch("jarvishep2.client.dispatch_version", side_effect=display),
        mock.patch("jarvishep2.redis_notice.emit_first_run_redis_notice", side_effect=notice),
    ):
        assert main(["--version"]) == 0
    assert events == ["display", "notice"]
