"""Agents can discover and close an authoring loop without starting a scan."""

import argparse
import io
import json
from contextlib import redirect_stdout
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from jarvishep2.cli.agent import agent_contract
from jarvishep2.cli.client import main, normalize_argv
from jarvishep2.cli.help import build_parser
from jarvishep2.cli.man import resolve_man_request
from jarvishep2.project_scaffold import create_project_scaffold
from jarvishep2.versioning import get_runtime_version


def invoke(*argv):
    output = io.StringIO()
    with redirect_stdout(output), patch("jarvishep2.redis_notice.emit_first_run_redis_notice"):
        status = main(list(argv))
    return status, output.getvalue()


def test_bootstrap_alias_is_same_versioned_contract_and_discoverable():
    status, text = invoke("agent", "--json")
    alias_status, alias_text = invoke("man", "agent-readme", "--json")
    contract = json.loads(text)
    assert status == alias_status == 0
    assert contract == json.loads(alias_text)
    assert contract["contract_version"] == 1
    assert contract["installed_version"] == get_runtime_version()
    assert "Jarvis agent --json" in build_parser().format_help()
    center = resolve_man_request([])
    assert any(row["man"] == "Jarvis man agent-readme" for row in center["list_rows"])
    assert normalize_argv(["agent", "--json"]) == ["agent", "--json"]


def test_catalog_covers_actual_cli_and_describes_execution_boundaries():
    contract = agent_contract()
    entries = {entry["command"]: entry for entry in contract["command_catalog"]}
    actual = next(a for a in build_parser()._actions if isinstance(a, argparse._SubParsersAction)).choices
    assert set(actual) == {name for name in entries if " " not in name}
    assert entries["check"]["effect"] == entries["run"]["effect"] == "execute"
    assert entries["sampler check"]["effect"] == "execute"
    assert entries["sampler load"]["effect"] == "execute_write"
    assert entries["validate"]["effect"] == "validate"
    assert entries["kill"]["effect"] == "terminate"
    for command in ("man", "plot", "portal", "operas"):
        assert entries[command]["arguments_complete"] is False
    assert "caches" in contract["output_contract"]["query_side_effects"]
    assert "already authorized" in " ".join(contract["rules"])


def test_human_contract_includes_workflow_boundaries_and_request_protocol():
    status, text = invoke("man", "agent-readme")
    assert status == 0
    for evidence in ("Jarvis validate", "manual_argv", "Bare Jarvis TASK.yaml", "Schema zones",
                     "Sampling fragments", "Command effects", "Resquest@Jarvis-HEP",
                     "Jarvis.HEP.AI@mail.grokbot.com", "No automatic sending"):
        assert evidence in text


@pytest.mark.parametrize("argv", [
    ("agent", "--json", "--unknown"),
    ("agent", "--json=bad"),
    ("man", "--json", "--type"),
    ("man", "--json=bad"),
    ("man", "--json", "--direction", "wrong"),
    ("validate", "--json"),
    ("validate", "missing.yaml", "--json", "--unknown"),
    ("validate", "missing.yaml", "--json=bad"),
    ("validate", "missing.yaml", "--json", "--plot"),
    ("monitor", "--json", "--unknown"),
    ("agent", "template", "--json", "missing"),
    ("agent", "schema", "--json", "--id", "https://example.invalid/schema"),
    ("agent", "request", "--json", "--send"),
])
def test_machine_readable_usage_failures_reach_agent_from_top_level(argv):
    status, text = invoke(*argv)
    assert status == 2
    payload = json.loads(text)
    assert payload["ok"] is False


def test_man_internal_failure_is_single_json(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("catalog unavailable")
    monkeypatch.setattr("jarvishep2.cli.man.resolve_man_request", fail)
    status, text = invoke("man", "yaml", "--json")
    assert status == 1
    assert json.loads(text)["diagnostics"][0]["message"].endswith("catalog unavailable")


def test_cli_authoring_diagnosis_repair_loop_without_runtime(tmp_path):
    project = Path(create_project_scaffold("AgentScan", cwd=str(tmp_path)))
    before = sorted(str(p.relative_to(project)) for p in project.rglob("*"))
    with patch("jarvishep2.core.Jarvis2Core") as core, patch("jarvishep2.redis_queue.RedisQueue") as redis:
        status, text = invoke("agent", "template", "random-operas", "--json")
        assert status == 0
        template = json.loads(text)
        assert template["artifact"]["kind"] == "complete_task_card"
        # Inspection/export did not create any task or runtime files.
        assert before == sorted(str(p.relative_to(project)) for p in project.rglob("*"))
        task = project / "bin" / "agent.yaml"
        document = yaml.safe_load(template["body"])
        bounds = document["Sampling"]["Bounds"]
        original_budget = bounds.pop("point_number")
        task.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        status, text = invoke("validate", str(task), "--json")
        report = json.loads(text)
        assert status == 2 and report["ok"] is False
        issue = next(i for i in report["issues"] if i["code"] == "JV2-MTH-020")
        status, text = invoke(*issue["manual_argv"][1:])
        assert status == 0
        manual = json.loads(text)
        assert manual["context"]["method"] == "Random"
        bounds["point_number"] = original_budget
        task.write_text(yaml.safe_dump(document, sort_keys=False), encoding="utf-8")
        status, text = invoke("validate", str(task), "--json")
        repaired = json.loads(text)
        assert status == 0 and repaired["ok"] is True
        assert repaired["metadata"]["execution_performed"] is False
        assert repaired["metadata"]["physics_correctness_verified"] is False
        core.assert_not_called()
        redis.assert_not_called()
    assert not any((project / name).exists() for name in ("outputs", "logs", "checkpoints"))


def test_fragment_identity_is_consistent_in_man_and_template():
    status, text = invoke("agent", "template", "dynesty", "--json")
    assert status == 0
    assert json.loads(text)["artifact"]["kind"] == "sampling_fragment"
    assert resolve_man_request(["example.dynesty"])["artifact"]["kind"] == "sampling_fragment"
    row = next(r for r in resolve_man_request(["example"])["list_rows"] if r["name"] == "dynesty")
    assert row["kind"] == "sampling_fragment"


def test_registered_plugin_schema_is_exposed_without_importing_plugin(tmp_path, monkeypatch):
    from types import SimpleNamespace
    schema = {"type": "object", "properties": {"step": {"type": "number", "minimum": 0, "maximum": 1}}, "additionalProperties": False}
    record = SimpleNamespace(options_schema=schema, description="Test metadata", profile="independent",
                             version="1.0", source="user", source_path=str(tmp_path / "unimportable.py"))
    monkeypatch.setattr("jarvishep2.sampler_catalog.plugin_record", lambda name: record)
    from jarvishep2.cli.man import _plugin_method_page, _print_page
    page = _plugin_method_page("MetadataOnly")
    assert page["options_schema"] == schema
    out = io.StringIO()
    with redirect_stdout(out):
        _print_page(page, as_json=True)
    assert json.loads(out.getvalue())["options_schema"]["properties"]["step"]["maximum"] == 1
