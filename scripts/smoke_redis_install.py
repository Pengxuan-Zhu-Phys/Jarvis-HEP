"""Verify Redis onboarding from an already installed complete HEP wheel.

Run with the release venv's Python and ``-I`` outside the source checkout.
Every CLI invocation also uses isolated Python with a temporary home/caches.
Child guards permit only a server's ``--version`` process and reject endpoint
connections; no package installation, service startup, or scan is performed.
"""

from __future__ import annotations

import importlib.metadata as metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import jarvishep2


_GUARDED_CLI = r"""
import json, pathlib, runpy, socket, subprocess, sys
cli, arguments, allow_version = sys.argv[1], json.loads(sys.argv[2]), sys.argv[3] == '1'
original_popen = subprocess.Popen
def guarded_popen(argv, *args, **kwargs):
    names = {'redis-server', 'redis6-server', 'valkey-server'}
    if (allow_version and not kwargs.get('shell') and isinstance(argv, (list, tuple))
            and len(argv) == 2 and argv[1] == '--version'
            and pathlib.Path(argv[0]).name in names):
        return original_popen(argv, *args, **kwargs)
    raise AssertionError('Redis query attempted an unexpected process: ' + repr(argv))
def reject_connection(*args, **kwargs):
    raise AssertionError('Redis query attempted a network connection')
subprocess.Popen = guarded_popen
socket.create_connection = reject_connection
socket.socket.connect = reject_connection
sys.argv = [cli, *arguments]
runpy.run_path(cli, run_name='__main__')
"""


def main(expected: str) -> None:
    prefix = Path(sys.prefix).resolve()
    package_path = Path(jarvishep2.__file__).resolve()
    assert metadata.version("Jarvis-HEP") == jarvishep2.__version__ == expected
    assert package_path.is_relative_to(prefix), (
        f"Smoke must import the installed wheel, got {package_path} outside {prefix}"
    )
    cli = prefix / "bin" / "Jarvis"
    assert cli.is_file(), f"Installed CLI is missing: {cli}"

    with tempfile.TemporaryDirectory(prefix="jarvis-redis-install-smoke-") as temporary:
        base = Path(temporary)
        (base / "home").mkdir()
        (base / "empty-bin").mkdir()
        env = dict(os.environ)
        for name in ("PYTHONPATH", "PYTHONHOME", "JARVIS_HEP_TASK_ROOT", "JHEP_TASK_ROOT"):
            env.pop(name, None)
        env.update({
            "HOME": str(base / "home"),
            "JARVIS_REDIS_NOTICE_STATE": str(base / "redis-notice.marker"),
            "JARVIS_SCAN_REF_REGISTRY": str(base / "scan-refs.json"),
            "JARVIS_OPERAS_CURVE_CACHE_ROOT": str(base / "operas-cache"),
            "MPLCONFIGDIR": str(base / "matplotlib"),
            "XDG_CONFIG_HOME": str(base / "config"),
            "XDG_CACHE_HOME": str(base / "cache"),
            "NO_COLOR": "1",
        })
        checks = 0

        def command(argv: list[str], *, expected_exit: int | None = 0,
                    allow_version: bool = False, environment=None):
            nonlocal checks
            result = subprocess.run(
                [sys.executable, "-I", "-c", _GUARDED_CLI, str(cli),
                 json.dumps(argv), "1" if allow_version else "0"],
                cwd=base, env=environment if environment is not None else env,
                text=True, capture_output=True, timeout=45,
            )
            if expected_exit is not None:
                assert result.returncode == expected_exit, (
                    argv, result.returncode, result.stdout, result.stderr
                )
            checks += 1
            return result

        def query(argv: list[str], **kwargs):
            result = command(argv, **kwargs)
            try:
                value = json.loads(result.stdout)
            except ValueError:
                raise AssertionError(
                    f"Expected one JSON object from {argv}: {result.stdout}\n{result.stderr}"
                ) from None
            assert isinstance(value, dict), (argv, value)
            return value, result

        group = command(["redis", "--help"]).stdout
        assert "Usage: Jarvis redis [OPTIONS] COMMAND [ARGS]..." in group
        assert "╭─ Commands" in group and "╭─ Arguments" not in group
        for topic in ("info", "check", "install", "start"):
            focused = command(["redis", topic, "--help"]).stdout
            assert f"Usage: Jarvis redis {topic} [OPTIONS]" in focused
            assert "--json" in focused and "╭─ Commands" not in focused
            assert focused != group
            entry, _ = query([
                "agent", "--section", "commands", "--command", "redis " + topic, "--json"
            ])
            catalog = entry["command_catalog"]
            assert len(catalog) == 1 and catalog[0]["command"] == "redis " + topic
            assert catalog[0]["effect"] == "query"
            assert catalog[0]["help_argv"] == ["Jarvis", "redis", topic, "--help"]
            found, _ = query(["agent", "search", "redis " + topic, "--limit", "3", "--json"])
            matches = [row for row in found["results"] if row["name"] == "redis " + topic]
            assert matches and matches[0]["next_argv"] == ["Jarvis", "redis", topic, "--help"]

            for argv in (["redis", topic, "--json"], ["redis", "--json", topic],
                         ["redis", "--json", topic, "--json"]):
                report, result = query(
                    argv, expected_exit=None if topic == "check" else 0,
                    allow_version=topic in {"info", "check"},
                )
                assert report["kind"] == "redis-environment" and report["topic"] == topic
                assert report["start"]["required_for_jarvis"] is False
                assert result.stderr == "", (argv, result.stderr)
                if topic == "check":
                    usable = bool(report["server"]["usable"]) and report["runtime_supported"] is not False
                    assert report["ok"] is usable and result.returncode == (0 if usable else 1)
                else:
                    assert report["ok"] is True
                if topic in {"install", "start"} and report["server"]["found"]:
                    assert report["server"]["version_status"] == "not_checked"

        implicit, _ = query(["redis", "--json"], allow_version=True)
        explicit, _ = query(["redis", "info", "--json"], allow_version=True)
        assert implicit == explicit
        assert command(["redis"], allow_version=True).stdout == command(
            ["redis", "info"], allow_version=True
        ).stdout

        for argv in (["redis", "unknown", "--json"], ["redis", "check", "--unknown", "--json"],
                     ["redis", "--json=bad"], ["redis", "install", "unexpected", "--json"]):
            error, _ = query(argv, expected_exit=2)
            assert error["ok"] is False

        missing_env = dict(env)
        missing_env["PATH"] = str(base / "empty-bin")
        missing_marker = base / "missing-redis.marker"
        missing_env["JARVIS_REDIS_NOTICE_STATE"] = str(missing_marker)
        first = command(["--version"], environment=missing_env)
        assert expected in first.stdout
        assert "Redis server was not found on PATH" in first.stderr
        assert "Jarvis redis install" in first.stderr and "Jarvis redis start" in first.stderr
        assert missing_marker.is_file()
        assert "Redis server was not found" not in command(["--version"], environment=missing_env).stderr
        for topic in ("info", "check", "install", "start"):
            report, result = query(
                ["redis", topic, "--json"], expected_exit=1 if topic == "check" else 0,
                environment=missing_env,
            )
            assert report["server"]["found"] is False and report["server"]["version_status"] == "not_found"
            assert report["ok"] is (topic != "check")
            assert result.stderr == ""

        assert not any((base / name).exists() for name in ("outputs", "logs", "checkpoints"))
        print(json.dumps({
            "ok": True, "installed_version": expected, "package_path": str(package_path),
            "checks": checks, "process_policy": "server --version only",
            "scan_execution": False, "package_installation": False, "service_start": False,
        }))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python -I smoke_redis_install.py VERSION")
    main(sys.argv[1])
