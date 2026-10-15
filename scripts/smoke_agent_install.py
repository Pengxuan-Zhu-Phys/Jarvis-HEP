"""Verify agent authoring interfaces from an already installed complete wheel.

Run with the release venv's Python and ``-I`` outside the source checkout. All
fixtures and dependency caches are temporary; this script does not run scans,
install packages or send mail.
"""

from __future__ import annotations

import importlib.metadata as metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import yaml

import jarvishep2


def main(expected: str) -> None:
    prefix = Path(sys.prefix).resolve()
    package_path = Path(jarvishep2.__file__).resolve()
    assert metadata.version("Jarvis-HEP") == jarvishep2.__version__ == expected
    assert package_path.is_relative_to(prefix), (
        f"Smoke must import the installed wheel, got {package_path} outside {prefix}"
    )
    cli = prefix / "bin" / "Jarvis"
    assert cli.is_file(), f"Installed CLI is missing: {cli}"

    with tempfile.TemporaryDirectory(prefix="jarvis-agent-install-smoke-") as temporary:
        base = Path(temporary)
        (base / "home").mkdir()
        env = dict(os.environ)
        for name in ("PYTHONPATH", "PYTHONHOME", "JARVIS_HEP_TASK_ROOT", "JHEP_TASK_ROOT"):
            env.pop(name, None)
        env.update({
            "HOME": str(base / "home"),
            "JARVIS_SCAN_REF_REGISTRY": str(base / "scan-refs.json"),
            "JARVIS_OPERAS_CURVE_CACHE_ROOT": str(base / "operas-cache"),
            "MPLCONFIGDIR": str(base / "matplotlib"),
            "XDG_CONFIG_HOME": str(base / "config"),
            "XDG_CACHE_HOME": str(base / "cache"),
        })
        checks = []

        def command(argv: list[str], *, expected_exit: int = 0):
            result = subprocess.run(
                [str(cli), *argv], cwd=base, env=env, text=True,
                capture_output=True, timeout=90,
            )
            assert result.returncode == expected_exit, (
                argv, result.returncode, result.stdout, result.stderr
            )
            checks.append(argv)
            return result

        def query(argv: list[str], *, expected_exit: int = 0):
            result = command(argv, expected_exit=expected_exit)
            try:
                return json.loads(result.stdout)
            except ValueError:
                raise AssertionError(
                    f"Expected one JSON value from {argv}: {result.stdout}\n{result.stderr}"
                ) from None

        def follow_manual(argv: list[str]):
            assert argv[:2] == ["Jarvis", "man"], argv
            return query(argv[1:])

        full = query(["agent", "--json"])
        assert full["installed_version"] == expected and full["contract_version"] == 1
        brief = query(["agent", "--brief", "--json"])
        assert brief["rules"] == full["rules"] and brief["workflows"] == full["workflows"]
        assert brief == query(["man", "agent-readme", "--brief", "--json"])
        run_entry = query([
            "agent", "--section", "commands", "--command", "run", "--json"
        ])["command_catalog"][0]
        assert run_entry["effect"] == "execute"

        search = query(["agent", "search", "checkpoint heartbeat", "--limit", "3", "--json"])
        assert search["results"]
        for argv in {tuple(row["next_argv"]) for row in search["results"]}:
            assert follow_manual(list(argv))["diagnostics"] == []

        bundle = query(["agent", "schema", "--json"])
        assert bundle["root_schema_id"] in bundle["schemas"]
        assert len(bundle["schemas"]) == len(bundle["manifest"]["schema_files"])
        index = query(["agent", "schema", "--list", "--json"])
        assert set(index["schema_ids"]) == set(bundle["schemas"])
        csv_manual = query(["man", "sampler.CSV", "--json"])
        csv_bounds = query(csv_manual["bounds_schema_argv"][1:])
        assert "path" in csv_bounds["required"] and "uuid_column" in csv_bounds["properties"]

        fragment = query(["agent", "template", "dynesty", "--json"])
        assert fragment["artifact"]["kind"] == "sampling_fragment"
        composition = next(row for row in fragment["prerequisites"] if row["kind"] == "composition")
        assert composition["operation"]["target_path"] == "$.Sampling"
        assert any("<ComposedTask.yaml>" in row.get("argv", []) for row in fragment["prerequisites"])
        template = query(["agent", "template", "random-operas", "--json"])
        operator = next(row for row in template["prerequisites"] if row["kind"] == "operator")
        assert operator["reference_kind"] == "python_import_path" and not operator.get("query_argv")
        assert template["body"] == command(["agent", "template", "random-operas"]).stdout

        command(["project", "create", "AgentInstallSmoke"])
        project = base / "AgentInstallSmoke"
        card = project / "bin" / "authored.yaml"
        document = yaml.safe_load(template["body"])
        document["Scan"]["name"] = "installed-agent-smoke"
        document["EnvReqs"]["V2"] = {"workers": 3, "checkpoint": {"heartbeat": 60}}
        card.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        inspected = query(["agent", "inspect", str(card), "--section", "runtime", "--json"])
        assert inspected["runtime"]["workers"] == 3
        assert inspected["runtime"]["checkpoint"]["heartbeat"] == 60
        assert inspected["metadata"]["execution_performed"] is False
        sources = query(["agent", "inspect", str(card), "--section", "sources", "--json"])
        assert sources["sources"]["runtime_fields"]["$.Runtime.workers"]["layer"] == "task"
        assert sources["sources"]["runtime_fields"]["$.Runtime.batch_size"]["layer"] == "default_environment"
        assert query(["validate", str(card), "--json"])["ok"]

        duplicate = project / "bin" / "duplicate.yaml"
        duplicate.write_text(
            card.read_text(encoding="utf-8").replace("Scan:\n", "Scan:\n  name: first_name\n", 1),
            encoding="utf-8",
        )
        bad = query(["validate", str(duplicate), "--json"], expected_exit=2)
        issue = next(row for row in bad["issues"] if row["code"] == "JV2-YAML-003")
        assert issue["path"] == "$.Scan.name" and issue["location"] and issue["related_locations"]
        assert follow_manual(issue["manual_argv"])["diagnostics"] == []
        rejected = query([
            "agent", "inspect", str(duplicate), "--section", "validation", "--json"
        ], expected_exit=2)
        assert rejected["issues"][0]["code"] == "JV2-YAML-003"
        unknown = query(["man", "calculator.execution.JSON.typo", "--json"], expected_exit=2)
        assert follow_manual(unknown["parent_manual_argv"])["diagnostics"] == []

        draft = query([
            "agent", "request", "--title", "Installed agent smoke",
            "--description", "Verify literal request evidence fields",
            "--yaml", "n/a: documentation request", "--commands", "Jarvis agent --brief --json",
            "--diagnostics", "exit 0", "--actual", "Current behavior",
            "--expected", "Requested behavior", "--constraints", "Preserve scientific assumptions",
            "--installation-method", "Complete release wheel",
            "--json",
        ])["data"]
        assert draft["draft_ready"] and draft["missing_fields"] == []
        assert draft["automatic_send"] is False and draft["installed_version"] == expected
        assert draft["recipient"] == "Jarvis.HEP.AI@mail.grokbot.com"
        assert "Resquest@Jarvis-HEP" in draft["subject"]
        assert query(["agent", "request", "--send", "--json"], expected_exit=2)["ok"] is False
        assert not any((project / name).exists() for name in ("outputs", "logs", "images", "checkpoints"))

        print(json.dumps({
            "ok": True, "installed_version": expected, "package_path": str(package_path),
            "checks": len(checks), "scan_execution": False, "mail_sent": False,
        }))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python -I smoke_agent_install.py VERSION")
    main(sys.argv[1])
