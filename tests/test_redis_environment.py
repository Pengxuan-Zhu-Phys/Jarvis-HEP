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


@pytest.mark.parametrize("executable", ["redis-server", "redis6-server", "valkey-server"])
def test_report_accepts_every_server_executable_used_by_runtime(executable):
    binary = "/opt/example/" + executable
    report = _report(which=lambda name: binary if name == executable else None)
    assert report["server"]["found"] is True
    assert report["server"]["path"] == binary
    assert report["server"]["executable"] == executable


@pytest.mark.parametrize(
    ("system", "distro", "expected"),
    [
        ("Darwin", {}, "brew install redis"),
        ("Linux", {"ID": "debian"}, "sudo apt install -y redis-server"),
        ("Linux", {"ID": "fedora"}, "sudo dnf install -y valkey"),
        ("Linux", {"ID": "arch"}, "sudo pacman -S valkey"),
        ("Linux", {"ID": "alpine", "VERSION_ID": "3.19"}, "sudo apk add redis"),
        ("Linux", {"ID": "opensuse-tumbleweed"}, "sudo zypper install -y redis"),
        ("Linux", {"ID": "amzn", "VERSION_ID": "2023"}, "sudo dnf install -y valkey"),
        ("Linux", {"ID": "amzn", "VERSION_ID": "2"}, "sudo amazon-linux-extras install redis6 -y"),
    ],
)
def test_report_installation_guidance_follows_host_system(system, distro, expected):
    report = _report(system=system, distro=distro)
    assert report["install"]["command"] == expected
    assert expected in report["install"]["commands"]


def test_linux_exact_identity_takes_precedence_over_fedora_ancestry():
    report = _report(
        distro={"ID": "rocky", "ID_LIKE": "rhel centos fedora"},
        which=lambda name: "/usr/bin/dnf" if name == "dnf" else None,
    )
    assert report["install"]["command"] == "sudo dnf install -y redis"


@pytest.mark.parametrize("distro_id", ["rhel", "rocky", "almalinux"])
def test_rhel_eight_selects_updated_redis_six_stream_instead_of_default_redis_five(distro_id):
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
    assert all(url.startswith("https://docs.redhat.com/") for url in report["install"]["documentation"])


@pytest.mark.parametrize(("version", "package"), [("9", "redis"), ("10", "valkey")])
def test_modern_rhel_releases_select_their_current_package(version, package):
    report = _report(
        distro={"ID": "rhel", "VERSION_ID": version},
        which=lambda name: "/usr/bin/dnf" if name == "dnf" else None,
    )
    assert report["install"]["command"] == f"sudo dnf install -y {package}"


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        ({"ID": "custom", "ID_LIKE": "debian"}, "sudo apt install -y redis-server"),
        ({"ID": "endeavouros", "ID_LIKE": "arch"}, "sudo pacman -S valkey"),
    ],
)
def test_linux_derivatives_use_distribution_family(metadata, expected):
    assert _report(distro=metadata)["install"]["command"] == expected


@pytest.mark.parametrize(
    ("metadata", "codename"),
    [
        ({"ID": "ubuntu", "VERSION_ID": "20.04"}, "focal"),
        ({"ID": "ubuntu", "VERSION_ID": "22.04"}, "jammy"),
        ({"ID": "debian", "VERSION_ID": "11"}, "bullseye"),
    ],
)
def test_legacy_apt_distributions_get_compatible_redis_from_official_repository(metadata, codename):
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


@pytest.mark.parametrize(
    "metadata",
    [
        {"ID": "ubuntu", "VERSION_ID": "24.04"},
        {"ID": "debian", "VERSION_ID": "12"},
    ],
)
def test_modern_apt_distributions_use_the_existing_compatible_distro_package(metadata):
    report = _report(distro=metadata)
    assert report["install"]["command"] == "sudo apt install -y redis-server"
    assert report["install"]["commands"] == ["sudo apt update", "sudo apt install -y redis-server"]


@pytest.mark.parametrize(
    "metadata",
    [
        {"ID": "ubuntu", "VERSION_ID": "18.04"},
        {"ID": "debian", "VERSION_ID": "10"},
    ],
)
def test_unverified_legacy_apt_recipes_do_not_advise_a_known_incompatible_package(metadata):
    report = _report(distro=metadata)
    assert report["install"]["commands"] == []
    assert report["install"]["documentation"]
    assert report["install"]["notes"]


def test_unknown_system_does_not_invent_an_executable_package_manager_command():
    report = _report(system="UnknownOS", distro={})
    assert not report["install"]["commands"]
    assert report["install"]["notes"]


def test_unknown_linux_without_package_manager_preserves_honest_install_guidance():
    report = _report(distro={"ID": "custom"})
    assert report["install"]["commands"] == []
    assert "package manager" in report["install"]["command"]
    assert report["server"]["found"] is False


def test_apt_get_fallback_uses_the_package_manager_actually_found():
    report = _report(distro={}, which=lambda name: "/usr/bin/apt-get" if name == "apt-get" else None)
    assert report["install"]["command"] == "sudo apt-get install -y redis-server"
    assert "sudo apt-get update" in report["install"]["commands"]


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


@pytest.mark.parametrize("source", ["environment", "kernel"])
def test_wsl_is_detected_as_supported_linux_environment(source):
    environment = {"WSL_DISTRO_NAME": "Ubuntu"} if source == "environment" else {}
    release = "5.15.0-microsoft-standard-WSL2" if source == "kernel" else "6.8.0"
    with mock.patch.dict("os.environ", environment, clear=True), mock.patch.object(redis_notice.platform, "release", return_value=release):
        report = _report()
    assert report["system"]["wsl"] is True
    assert report["runtime_supported"] is True
    assert report["install"]["command"] == "sudo apt install -y redis-server"


@pytest.mark.parametrize(
    ("executable", "version_text", "version"),
    [
        ("redis-server", "Redis server v=7.2.8 sha=00000000:0 malloc=libc bits=64 build=test", "7.2.8"),
        ("valkey-server", "Valkey server v=8.1.2 sha=00000000:0 malloc=libc bits=64 build=test", "8.1.2"),
    ],
)
def test_version_probe_is_bounded_and_only_runs_the_version_command(executable, version_text, version):
    binary = "/opt/example binaries/" + executable
    completed = subprocess.CompletedProcess([binary, "--version"], 0, stdout=version_text, stderr="")
    with (
        mock.patch.object(redis_notice.subprocess, "run", return_value=completed) as run,
        mock.patch("socket.create_connection") as connect,
    ):
        report = _report(
            which=lambda name: binary if name == executable else None,
            probe_version=True,
        )
    assert report["server"]["version"] == version
    assert report["server"]["version_status"] == "ok"
    assert report["server"]["supports_blmove"] is True
    assert report["server"]["usable"] is True
    run.assert_called_once()
    assert run.call_args.args[0] == [binary, "--version"]
    assert 0 < run.call_args.kwargs["timeout"] <= 5
    assert not run.call_args.kwargs.get("shell", False)
    connect.assert_not_called()


def test_guidance_without_probe_does_not_execute_an_installed_program():
    with mock.patch.object(redis_notice.subprocess, "run") as run:
        report = _report(which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None)
    assert report["server"]["found"] is True
    assert report["server"]["version"] is None
    run.assert_not_called()


def test_server_missing_skips_version_process_even_when_probe_requested():
    with mock.patch.object(redis_notice.subprocess, "run") as run:
        report = _report(probe_version=True)
    assert report["server"]["version_status"] == "not_found"
    assert report["server"]["usable"] is False
    run.assert_not_called()


@pytest.mark.parametrize(
    "failure",
    [
        subprocess.TimeoutExpired(["redis-server", "--version"], timeout=2),
        PermissionError("redis-server is not executable"),
        subprocess.CompletedProcess(["redis-server", "--version"], 1, stdout="", stderr="broken executable"),
        subprocess.CompletedProcess(["redis-server", "--version"], 0, stdout="unrecognized output", stderr=""),
    ],
)
def test_unusable_version_probe_reports_a_found_but_unverified_executable(failure):
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


def test_old_redis_version_explains_the_blmove_requirement():
    old = subprocess.CompletedProcess([], 0, stdout="Redis server v=6.0.20 sha=test", stderr="")
    with mock.patch.object(redis_notice.subprocess, "run", return_value=old):
        report = _report(
            which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None,
            probe_version=True,
        )
    assert report["server"]["version"] == "6.0.20"
    assert report["server"]["supports_blmove"] is False
    assert report["server"]["usable"] is False
    assert "6.2" in json.dumps(report)
    assert "BLMOVE" in json.dumps(report)


def test_redis_minimum_version_is_accepted_and_stderr_version_output_is_supported():
    completed = subprocess.CompletedProcess([], 0, stdout="", stderr="Redis server v=6.2.0 sha=test")
    with mock.patch.object(redis_notice.subprocess, "run", return_value=completed):
        report = _report(
            which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None,
            probe_version=True,
        )
    assert report["server"]["version"] == "6.2.0"
    assert report["server"]["usable"] is True


def test_manual_start_command_quotes_binary_paths_and_uses_local_optional_port():
    report = _report(which=lambda name: "/opt/example binaries/redis-server" if name == "redis-server" else None)
    command = report["start"]["command"]
    assert command.startswith("'/opt/example binaries/redis-server' ")
    assert "--bind 127.0.0.1" in command
    assert "--port 6380" in command
    assert report["start"]["required_for_jarvis"] is False


def test_redis_command_routes_to_its_dispatcher_and_normalizes_unchanged():
    from jarvishep2.cli.client import normalize_argv

    assert normalize_argv(["redis", "install", "--json"]) == ["redis", "install", "--json"]
    with mock.patch("jarvishep2.cli.redis.dispatch_redis", return_value=0) as dispatch:
        status, _output, _errors = _invoke("redis", "install", "--json")
    assert status == 0
    dispatch.assert_called_once_with(topic="install", as_json=True)


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


@pytest.mark.parametrize("command", ["info", "check", "install", "start"])
def test_each_redis_subcommand_has_focused_standard_help_without_querying_environment(command):
    from jarvishep2.cli.client import main

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


@pytest.mark.parametrize("command", ["info", "check", "install", "start"])
@pytest.mark.parametrize("placement", ["group", "child", "both"])
def test_redis_json_flag_is_preserved_before_or_after_the_child_command(command, placement):
    argv = ["redis"]
    if placement in {"group", "both"}:
        argv.append("--json")
    argv.append(command)
    if placement in {"child", "both"}:
        argv.append("--json")
    with mock.patch("jarvishep2.cli.redis.dispatch_redis", return_value=0) as dispatch:
        status, _output, _errors = _invoke(*argv)
    assert status == 0
    dispatch.assert_called_once_with(topic=command, as_json=True)


@pytest.mark.parametrize(
    ("argv", "as_json"),
    [
        (("redis",), False),
        (("redis", "--json"), True),
        (("redis", "info"), False),
        (("redis", "info", "--json"), True),
    ],
)
def test_bare_redis_and_explicit_info_remain_the_same_overview_alias(argv, as_json):
    with mock.patch("jarvishep2.cli.redis.dispatch_redis", return_value=0) as dispatch:
        status, _output, _errors = _invoke(*argv)
    assert status == 0
    dispatch.assert_called_once_with(topic="info", as_json=as_json)


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


@pytest.mark.parametrize("command", ["info", "check", "install", "start"])
def test_bounded_agent_search_discovers_the_focused_redis_child_help(command):
    with mock.patch("jarvishep2.cli.redis.build_redis_environment_report") as build:
        status, output, _errors = _invoke("agent", "search", "redis " + command, "--limit", "3", "--json")
    assert status == 0
    result = json.loads(output)
    assert len(result["results"]) <= 3
    matches = [row for row in result["results"] if row["name"] == "redis " + command]
    assert matches
    assert matches[0]["next_argv"] == ["Jarvis", "redis", command, "--help"]
    assert matches[0]["effect"] == "query"
    build.assert_not_called()


@pytest.mark.parametrize("topic", ["info", "install", "start"])
def test_missing_redis_guidance_is_query_success_and_single_json(topic):
    with mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=_report()):
        status, output, _errors = _invoke("redis", topic, "--json")
    assert status == 0
    report = json.loads(output)
    assert report["kind"] == "redis-environment"
    assert report["server"]["found"] is False


def test_check_missing_redis_reports_failure_without_installing_or_starting():
    with (
        mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=_report()),
        mock.patch.object(redis_notice.subprocess, "run") as run,
        mock.patch("socket.create_connection") as connect,
    ):
        status, output, _errors = _invoke("redis", "check", "--json")
    assert status == 1
    assert json.loads(output)["server"]["found"] is False
    run.assert_not_called()
    connect.assert_not_called()


@pytest.mark.parametrize("usable", [True, False])
def test_check_uses_executable_compatibility_without_needing_a_running_service(usable):
    report = _report(which=lambda name: "/usr/bin/redis-server" if name == "redis-server" else None)
    report["server"].update(usable=usable, supports_blmove=usable, version="7.2.8" if usable else "6.0.20")
    with mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=report) as build:
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


@pytest.mark.parametrize("system", ["Linux", "UnknownOS"])
def test_unknown_system_installation_help_keeps_manual_instruction_visible(system):
    report = _report(system=system, distro={"ID": "custom"})
    with mock.patch("jarvishep2.cli.redis.build_redis_environment_report", return_value=report):
        status, output, _errors = _invoke("redis", "install")
    assert status == 0
    assert report["install"]["command"] in output
    assert "package manager" in output


@pytest.mark.parametrize("topic", ["install", "start"])
def test_install_and_start_only_print_instructions(topic):
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


@pytest.mark.parametrize(
    "tail",
    [
        ("unknown", "--json"),
        ("info", "--unknown", "--json"),
        ("check", "--plot", "--json"),
        ("--json=bad",),
        ("install", "unexpected", "--json"),
        ("--json", "check", "--unknown"),
        ("--json", "invalid-command"),
    ],
)
def test_redis_machine_readable_usage_errors_are_one_json_object(tail):
    status, output, _errors = _invoke("redis", *tail)
    assert status == 2
    assert json.loads(output)["ok"] is False
