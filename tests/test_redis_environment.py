"""Redis onboarding reports and CLI queries stay useful without running a service."""

from __future__ import annotations

import argparse
import io
import json
import subprocess
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

import pytest

from jarvishep2.queue import redis_notice


def _report(**kwargs):
    options = {
        "system": "Linux",
        "distro": {"ID": "ubuntu", "PRETTY_NAME": "Ubuntu"},
        "which": lambda _name: None,
        "probe_version": False,
    }
    options.update(kwargs)
    return redis_notice.build_redis_environment_report(**options)


def _invoke(*argv):
    from jarvishep2.cli.client import main

    output, errors = io.StringIO(), io.StringIO()
    with (
        redirect_stdout(output),
        redirect_stderr(errors),
        mock.patch("jarvishep2.redis_notice.emit_first_run_redis_notice"),
    ):
        status = main(list(argv))
    return status, output.getvalue(), errors.getvalue()


def test_missing_report_distinguishes_executable_from_service_status():
    with mock.patch("socket.create_connection") as connect:
        report = _report()
    assert report["kind"] == "redis-environment"
    assert report["ok"] is True
    assert report["system"]["distro_id"] == "ubuntu"
    assert report["server"]["found"] is False
    assert report["server"]["path"] is None
    assert report["install"]["command"] == "sudo apt install -y redis-server"
    assert report["start"]["required_for_jarvis"] is False
    assert report["start"]["command"]
    assert report["runtime_supported"] is True
    connect.assert_not_called()


def test_report_accepts_every_server_executable_used_by_runtime(subtests):
    for executable in ("redis-server", "redis6-server", "valkey-server"):
        with subtests.test(executable=executable):
            binary = "/opt/example/" + executable
            report = _report(which=lambda name, e=executable, b=binary: b if name == e else None)
            assert report["server"]["found"] is True
            assert report["server"]["path"] == binary
            assert report["server"]["executable"] == executable


def test_report_installation_guidance_follows_host_system(subtests):
    apt = ["sudo apt update", "sudo apt install -y redis-server"]
    cases = (
        # (system, os-release, package manager on PATH, install commands)
        ("Darwin", {}, None, ["brew install redis"]),
        ("Linux", {"ID": "debian"}, None, apt),
        ("Linux", {"ID": "fedora"}, None, ["sudo dnf install -y valkey"]),
        ("Linux", {"ID": "arch"}, None, ["sudo pacman -S valkey"]),
        ("Linux", {"ID": "alpine", "VERSION_ID": "3.19"}, None, ["sudo apk update", "sudo apk add redis"]),
        ("Linux", {"ID": "opensuse-tumbleweed"}, None, ["sudo zypper install -y redis"]),
        ("Linux", {"ID": "amzn", "VERSION_ID": "2023"}, None, ["sudo dnf install -y valkey"]),
        ("Linux", {"ID": "amzn", "VERSION_ID": "2"}, None, ["sudo amazon-linux-extras install redis6 -y"]),
        # The exact identity takes precedence over Fedora ancestry.
        ("Linux", {"ID": "rocky", "ID_LIKE": "rhel centos fedora"}, "dnf", ["sudo dnf install -y redis"]),
        ("Linux", {"ID": "rhel", "VERSION_ID": "9"}, "dnf", ["sudo dnf install -y redis"]),
        ("Linux", {"ID": "rhel", "VERSION_ID": "10"}, "dnf", ["sudo dnf install -y valkey"]),
        # Derivatives use their distribution family.
        ("Linux", {"ID": "custom", "ID_LIKE": "debian"}, None, apt),
        ("Linux", {"ID": "endeavouros", "ID_LIKE": "arch"}, None, ["sudo pacman -S valkey"]),
        # Modern apt releases use the existing compatible distro package.
        ("Linux", {"ID": "ubuntu", "VERSION_ID": "24.04"}, None, apt),
        ("Linux", {"ID": "debian", "VERSION_ID": "12"}, None, apt),
        # The fallback uses the package manager actually found.
        ("Linux", {}, "apt-get", ["sudo apt-get update", "sudo apt-get install -y redis-server"]),
    )
    for system, distro, tool, commands in cases:
        with subtests.test(system=system, distro=distro):
            report = _report(
                system=system,
                distro=distro,
                which=lambda name, tool=tool: f"/usr/bin/{tool}" if name == tool else None,
            )
            assert report["install"]["commands"] == commands
            assert report["install"]["command"] == commands[-1]



def test_rhel_eight_selects_updated_redis_six_stream_instead_of_default_redis_five(subtests):
    for distro_id in ("rhel", "rocky", "almalinux"):
        with subtests.test(distro_id=distro_id):
            report = _report(
                distro={"ID": distro_id, "VERSION_ID": "8.10", "ID_LIKE": "rhel centos fedora"},
                which=lambda name: "/usr/bin/dnf" if name == "dnf" else None,
            )
            assert report["install"]["commands"] == ["sudo dnf module install -y redis:6"]
            notes = " ".join(report["install"]["notes"])
            assert "6.2.7" in notes
            assert "AppStream" in notes
            assert "redis:5" in notes
            assert "Jarvis redis check" in notes
            assert all(
                url.startswith("https://docs.redhat.com/") for url in report["install"]["documentation"]
            )




def test_legacy_apt_distributions_get_compatible_redis_from_official_repository(subtests):
    cases = (
        ({"ID": "ubuntu", "VERSION_ID": "20.04"}, "focal"),
        ({"ID": "ubuntu", "VERSION_ID": "22.04"}, "jammy"),
        ({"ID": "debian", "VERSION_ID": "11"}, "bullseye"),
    )
    for metadata, codename in cases:
        with subtests.test(codename=codename):
            with mock.patch.object(redis_notice.subprocess, "run") as run:
                report = _report(distro=metadata)
            commands = report["install"]["commands"]
            text = "\n".join(commands)
            assert report["install"]["command"] == "sudo apt-get install -y redis"
            assert commands[-1] == "sudo apt-get install -y redis"
            assert "https://packages.redis.io/gpg" in text
            assert "signed-by=/usr/share/keyrings/redis-archive-keyring.gpg" in text
            assert f"packages.redis.io/deb {codename} main" in text
            assert commands.count("sudo apt-get update") == 2
            assert "repository" in " ".join(report["install"]["notes"]).lower()
            run.assert_not_called()



def test_hosts_without_a_verified_recipe_get_guidance_but_no_commands(subtests):
    cases = (
        # Unverified legacy apt recipes would advise a known incompatible package.
        ("Linux", {"ID": "ubuntu", "VERSION_ID": "18.04"}, None),
        ("Linux", {"ID": "debian", "VERSION_ID": "10"}, None),
        ("UnknownOS", {}, "package manager"),
        ("Linux", {"ID": "custom"}, "package manager"),
    )
    for system, distro, fragment in cases:
        with subtests.test(system=system, distro=distro):
            report = _report(system=system, distro=distro)
            assert report["install"]["commands"] == []
            assert report["install"]["notes"]
            assert report["install"]["documentation"]
            if fragment:
                assert fragment in report["install"]["command"]





def test_native_windows_guides_users_to_wsl_for_both_jarvis_and_redis():
    report = _report(system="Windows", distro={})
    assert report["runtime_supported"] is False
    text = json.dumps(report).lower()
    assert "wsl" in text
    assert "winget install --id redis.redis" not in text
    assert "jarvis" in " ".join(report["notes"]).lower()


def test_native_windows_manual_start_does_not_reuse_a_windows_executable_path():
    report = _report(
        system="Windows", distro={},
        which=lambda name: r"C:\Program Files\Redis\redis-server.exe" if name == "redis-server" else None,
    )
    assert report["server"]["found"] is True
    assert report["start"]["command"].startswith("redis-server ")
    assert "Program Files" not in report["start"]["command"]
    assert "WSL" in " ".join(report["start"]["notes"])


def test_wsl_is_detected_as_supported_linux_environment(subtests):
    cases = (
        ("environment", {"WSL_DISTRO_NAME": "Ubuntu"}, "6.8.0"),
        ("kernel", {}, "5.15.0-microsoft-standard-WSL2"),
    )
    for source, environment, release in cases:
        with subtests.test(source=source):
            with (
                mock.patch.dict("os.environ", environment, clear=True),
                mock.patch.object(redis_notice.platform, "release", return_value=release),
            ):
                report = _report()
            assert report["system"]["wsl"] is True
            assert report["runtime_supported"] is True
            assert report["install"]["command"] == "sudo apt install -y redis-server"


def test_version_probe_is_bounded_and_only_runs_the_version_command(subtests):
    cases = (
        ("redis-server", "Redis server v=7.2.8 sha=00000000:0 malloc=libc bits=64 build=test", "", "7.2.8"),
        ("valkey-server", "Valkey server v=8.1.2 sha=00000000:0 malloc=libc bits=64 build=test", "", "8.1.2"),
        # Redis 6.2 is the minimum (BLMOVE); some builds print the version on stderr.
        ("redis-server", "", "Redis server v=6.2.0 sha=test", "6.2.0"),
        ("redis-server", "Redis server v=6.0.20 sha=test", "", "6.0.20"),
    )
    for executable, stdout, stderr, version in cases:
        with subtests.test(version=version):
            binary = "/opt/example binaries/" + executable
            completed = subprocess.CompletedProcess([binary, "--version"], 0, stdout=stdout, stderr=stderr)
            with (
                mock.patch.object(redis_notice.subprocess, "run", return_value=completed) as run,
                mock.patch("socket.create_connection") as connect,
            ):
                report = _report(
                    which=lambda name, e=executable, b=binary: b if name == e else None,
                    probe_version=True,
                )
            usable = version != "6.0.20"
            assert report["server"]["version"] == version
            assert report["server"]["version_status"] == "ok"
            assert report["server"]["supports_blmove"] is usable
            assert report["server"]["usable"] is usable
            if not usable:
                assert "6.2" in json.dumps(report)
                assert "BLMOVE" in json.dumps(report)
            run.assert_called_once()
            assert run.call_args.args[0] == [binary, "--version"]
            assert 0 < run.call_args.kwargs["timeout"] <= 5
            assert not run.call_args.kwargs.get("shell", False)
            connect.assert_not_called()


def test_no_version_process_runs_without_probe_or_without_a_server(subtests):
    installed = lambda name: "/usr/bin/redis-server" if name == "redis-server" else None  # noqa: E731
    cases = (
        ("found, no probe", installed, False, {"found": True, "version": None}),
        ("missing, probe", lambda _name: None, True, {"version_status": "not_found", "usable": False}),
    )
    for label, which, probe, expected in cases:
        with subtests.test(label):
            with mock.patch.object(redis_notice.subprocess, "run") as run:
                report = _report(which=which, probe_version=probe)
            assert {key: report["server"][key] for key in expected} == expected
            run.assert_not_called()



def test_unusable_version_probe_reports_a_found_but_unverified_executable(subtests):
    failures = (
        subprocess.TimeoutExpired(["redis-server", "--version"], timeout=2),
        PermissionError("redis-server is not executable"),
        subprocess.CompletedProcess(["redis-server", "--version"], 1, stdout="", stderr="broken executable"),
        subprocess.CompletedProcess(["redis-server", "--version"], 0, stdout="unrecognized output", stderr=""),
    )
    for failure in failures:
        with subtests.test(failure=repr(failure)):
            kwargs = {"side_effect": failure} if isinstance(failure, Exception) else {"return_value": failure}
            with mock.patch.object(redis_notice.subprocess, "run", **kwargs):
                report = _report(
                    which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None,
                    probe_version=True,
                )
            assert report["ok"] is True
            assert report["server"]["found"] is True
            assert report["server"]["version_error"]
            assert report["server"]["version_status"] not in {"compatible", "ok", "supported"}
            assert report["server"]["usable"] is False




def test_manual_start_command_quotes_binary_paths_and_uses_local_optional_port():
    report = _report(which=lambda name: "/opt/example binaries/redis-server" if name == "redis-server" else None)
    command = report["start"]["command"]
    assert command.startswith("'/opt/example binaries/redis-server' ")
    assert "--bind 127.0.0.1" in command
    assert "--port 6380" in command
    assert report["start"]["required_for_jarvis"] is False


def test_redis_argv_routes_to_its_dispatcher_unchanged(subtests):
    from jarvishep2.cli.client import normalize_argv

    # Bare `redis` and `redis info` are the same overview alias.
    cases = [
        (("redis",), "info", False),
        (("redis", "--json"), "info", True),
        (("redis", "info"), "info", False),
        (("redis", "info", "--json"), "info", True),
    ]
    # --json is kept before, after or on both sides of the child command.
    for command in ("info", "check", "install", "start"):
        cases += [
            (("redis", "--json", command), command, True),
            (("redis", command, "--json"), command, True),
            (("redis", "--json", command, "--json"), command, True),
        ]
    for argv, topic, as_json in cases:
        with subtests.test(argv=argv):
            assert normalize_argv(list(argv)) == list(argv)
            with mock.patch("jarvishep2.cli.redis.dispatch_redis", return_value=0) as dispatch:
                status, _output, _errors = _invoke(*argv)
            assert status == 0
            dispatch.assert_called_once_with(topic=topic, as_json=as_json)


def test_explicit_redis_report_claims_first_run_notice_without_duplicate_guidance():
    from jarvishep2.cli.client import main

    output, errors = io.StringIO(), io.StringIO()
    with (
        redirect_stdout(output),
        redirect_stderr(errors),
        mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=_report()),
        mock.patch("jarvishep2.redis_notice.emit_first_run_redis_notice") as notice,
    ):
        status = main(["redis", "--json"])
    assert status == 0
    assert json.loads(output.getvalue())["server"]["found"] is False
    assert errors.getvalue() == ""
    notice.assert_called_once_with(show_hint=False)


def test_explicit_redis_report_writes_marker_and_does_not_repeat_next_notice(tmp_path):
    from jarvishep2.cli.client import main

    marker = tmp_path / "onboarding.marker"
    output, errors = io.StringIO(), io.StringIO()
    with (
        mock.patch.dict("os.environ", {"JARVIS_REDIS_NOTICE_STATE": str(marker)}),
        mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=_report()),
        mock.patch.object(redis_notice, "find_redis_server", return_value=None),
        redirect_stdout(output), redirect_stderr(errors),
    ):
        status = main(["redis", "--json"])
        assert marker.is_file()
        assert redis_notice.emit_first_run_redis_notice(stream=errors) is False
    assert status == 0
    assert json.loads(output.getvalue())["server"]["found"] is False
    assert errors.getvalue() == ""


def test_redis_help_and_agent_catalog_expose_read_only_environment_queries():
    from jarvishep2.cli.agent import command_catalog
    from jarvishep2.cli.help import build_parser

    parser = build_parser()
    assert "redis" in parser.format_help().lower()
    choices = next(action.choices for action in parser._actions if isinstance(action, argparse._SubParsersAction))
    help_text = choices["redis"].format_help()
    for topic in ("info", "check", "install", "start", "--json"):
        assert topic in help_text
    entry = next(row for row in command_catalog() if row["command"] == "redis")
    assert entry["effect"] == "query"
    assert "json" in entry["output_formats"]


def test_redis_setup_actions_are_real_nested_subcommands_with_standard_group_help():
    from jarvishep2.cli.help import build_parser

    root = build_parser()
    root_commands = next(action.choices for action in root._actions if isinstance(action, argparse._SubParsersAction))
    redis = root_commands["redis"]
    nested_actions = [action for action in redis._actions if isinstance(action, argparse._SubParsersAction)]
    assert len(nested_actions) == 1
    assert set(nested_actions[0].choices) == {"info", "check", "install", "start"}
    assert not any(action.dest == "topic" and action.choices for action in redis._actions)
    text = redis.format_help()
    assert "Usage: Jarvis redis [OPTIONS] COMMAND [ARGS]..." in text
    assert "╭─ Commands" in text
    assert "╭─ Scan workflow" not in text
    assert "╭─ Arguments" not in text
    for command in ("info", "check", "install", "start"):
        assert any(line.lstrip().startswith("│ " + command + " ") for line in text.splitlines())


def test_each_redis_subcommand_has_focused_standard_help_without_querying_environment(subtests):
    from jarvishep2.cli.client import main

    for command in ("info", "check", "install", "start"):
        with subtests.test(command=command):
            output, errors = io.StringIO(), io.StringIO()
            with (
                redirect_stdout(output), redirect_stderr(errors),
                mock.patch("jarvishep2.redis_notice.emit_first_run_redis_notice"),
                mock.patch("jarvishep2.cli.redis.build_redis_environment_report") as build,
                pytest.raises(SystemExit) as stopped,
            ):
                main(["redis", command, "--help"])
            assert stopped.value.code == 0
            text = output.getvalue()
            assert f"Usage: Jarvis redis {command} [OPTIONS]" in text
            assert "╭─ Options" in text
            assert "--json" in text
            assert "--help" in text
            assert "╭─ Commands" not in text
            assert "╭─ Arguments" not in text
            assert errors.getvalue() == ""
            build.assert_not_called()




def test_agent_catalog_describes_each_redis_child_as_an_inspectable_query():
    from jarvishep2.cli.agent import command_catalog

    entries = {row["command"]: row for row in command_catalog()}
    for command in ("info", "check", "install", "start"):
        entry = entries["redis " + command]
        assert entry["argv_prefix"] == ["Jarvis", "redis", command]
        assert entry["help_argv"] == ["Jarvis", "redis", command, "--help"]
        assert entry["effect"] == "query"
        assert entry["arguments_complete"] is True
        assert "json" in entry["output_formats"]
        assert any("--json" in argument["flags"] for argument in entry["arguments"])


def test_bounded_agent_search_discovers_the_focused_redis_child_help(subtests):
    for command in ("info", "check", "install", "start"):
        with subtests.test(command=command):
            with mock.patch("jarvishep2.cli.redis.build_redis_environment_report") as build:
                status, output, _errors = _invoke(
                    "agent", "search", "redis " + command, "--limit", "3", "--json"
                )
            assert status == 0
            result = json.loads(output)
            assert len(result["results"]) <= 3
            matches = [row for row in result["results"] if row["name"] == "redis " + command]
            assert matches
            assert matches[0]["next_argv"] == ["Jarvis", "redis", command, "--help"]
            assert matches[0]["effect"] == "query"
            build.assert_not_called()


def test_missing_redis_is_one_json_report_and_only_check_fails(subtests):
    for topic, expected in (("info", 0), ("install", 0), ("start", 0), ("check", 1)):
        with subtests.test(topic=topic):
            with (
                mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=_report()),
                mock.patch.object(redis_notice.subprocess, "run") as run,
                mock.patch("socket.create_connection") as connect,
            ):
                status, output, _errors = _invoke("redis", topic, "--json")
            assert status == expected
            report = json.loads(output)
            assert report["kind"] == "redis-environment"
            assert report["server"]["found"] is False
            run.assert_not_called()
            connect.assert_not_called()



def test_check_uses_executable_compatibility_without_needing_a_running_service(subtests):
    for usable in (True, False):
        with subtests.test(usable=usable):
            report = _report(which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None)
            report["server"].update(
                usable=usable, supports_blmove=usable, version="7.2.8" if usable else "6.0.20"
            )
            with mock.patch(
                "jarvishep2.cli.redis.build_redis_environment_report", return_value=report
            ) as build:
                status, output, _errors = _invoke("redis", "check", "--json")
            assert status == (0 if usable else 1)
            assert json.loads(output)["ok"] is usable
            build.assert_called_once_with(probe_version=True)


def test_cli_internal_report_error_is_single_json_with_failure_status():
    with mock.patch("jarvishep2.cli.redis.build_redis_environment_report", side_effect=OSError("unreadable platform")):
        status, output, _errors = _invoke("redis", "info", "--json")
    assert status == 1
    report = json.loads(output)
    assert report["ok"] is False
    assert "unreadable platform" in report["error"]["message"]


def test_unknown_system_installation_help_keeps_manual_instruction_visible(subtests):
    for system in ("Linux", "UnknownOS"):
        with subtests.test(system=system):
            report = _report(system=system, distro={"ID": "custom"})
            with mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=report):
                status, output, _errors = _invoke("redis", "install")
            assert status == 0
            assert report["install"]["command"] in output
            assert "package manager" in output


def test_install_and_start_only_print_instructions(subtests):
    for topic in ("install", "start"):
        with subtests.test(topic=topic):
            with (
                mock.patch.object(redis_notice, "find_redis_server", return_value="/usr/bin/redis-server"),
                mock.patch.object(redis_notice.subprocess, "run") as run,
                mock.patch("subprocess.Popen") as popen,
                mock.patch("socket.create_connection") as connect,
            ):
                status, output, _errors = _invoke("redis", topic, "--json")
            assert status == 0
            assert json.loads(output)["start"]["required_for_jarvis"] is False
            run.assert_not_called()
            popen.assert_not_called()
            connect.assert_not_called()


def test_redis_machine_readable_usage_errors_are_one_json_object(subtests):
    tails = (
        ("unknown", "--json"),
        ("info", "--unknown", "--json"),
        ("check", "--plot", "--json"),
        ("--json=bad",),
        ("install", "unexpected", "--json"),
        ("--json", "check", "--unknown"),
        ("--json", "invalid-command"),
    )
    for tail in tails:
        with subtests.test(tail=tail):
            status, output, _errors = _invoke("redis", *tail)
            assert status == 2
            assert json.loads(output)["ok"] is False
