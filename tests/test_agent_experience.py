"""Regression scenarios discovered while using the installed CLI as an agent."""

import io
import json
from contextlib import redirect_stdout
from unittest.mock import patch

import pytest
import yaml

from jarvishep2.cli.client import main
from jarvishep2.task_config import TaskCardLoadError, load_task_yaml
from jarvishep2.cardload.yaml_loader import UniqueKeySafeLoader


def invoke(*argv):
    out = io.StringIO()
    with redirect_stdout(out), patch("jarvishep2.redis_notice.emit_first_run_redis_notice"):
        try:
            code = main(list(argv))
        except SystemExit as exc:
            code = exc.code
    return code, out.getvalue()


def test_short_bootstrap_keeps_all_rules_and_workflows_with_far_less_output():
    status, text = invoke("agent", "--brief", "--json")
    full_status, full_text = invoke("agent", "--json")
    alias_status, alias = invoke("man", "agent-readme", "--brief", "--json")
    short, full = json.loads(text), json.loads(full_text)
    assert status == full_status == alias_status == 0
    assert short == json.loads(alias)
    assert len(text.encode()) < len(full_text.encode()) / 3
    assert short["rules"] == full["rules"]
    assert short["workflows"] == full["workflows"]
    assert short["output_contract"] == full["output_contract"]
    assert "body_markdown" not in short and "command_catalog" not in short


def test_exact_command_section_is_a_small_usable_parser_description():
    status, text = invoke("agent", "--section", "commands", "--command", "run", "--json")
    data = json.loads(text)
    assert status == 0
    assert len(data["command_catalog"]) == 1
    entry = data["command_catalog"][0]
    assert entry["effect"] == "execute"
    assert any("--resume" in arg["flags"] for arg in entry["arguments"])
    assert invoke(*entry["help_argv"][1:])[0] == 0


@pytest.mark.parametrize("argv", [
    ("agent", "--brief", "--section", "commands", "--json"),
    ("agent", "--section", "rules", "--command", "run", "--json"),
    ("agent", "--command", "not-a-command", "--json"),
    ("agent", "search", "--json"),
    ("agent", "search", "checkpoint", "--limit", "0", "--json"),
])
def test_new_discovery_usage_errors_are_structured(argv):
    status, text = invoke(*argv)
    assert status == 2
    assert json.loads(text)["ok"] is False


def test_shared_json_position_before_helper_is_supported():
    status, text = invoke("agent", "--json", "request", "--title", "Question")
    assert status == 0
    assert json.loads(text)["data"]["title"] == "Question"


def test_search_can_answer_checkpoint_question_through_real_parent_manual():
    status, text = invoke("agent", "search", "checkpoint heartbeat", "--limit", "3", "--json")
    assert status == 0
    result = json.loads(text)
    assert result["results"]
    assert result["execution_performed"] is False
    found = False
    for entry in result["results"]:
        status, page_text = invoke(*entry["next_argv"][1:])
        assert status == 0
        page = json.loads(page_text)
        if page["path"] == "$.EnvReqs.V2.checkpoint":
            heartbeat = next(row for row in page["keys"] if row["name"] == "heartbeat")
            assert heartbeat["unit"] == "seconds"
            assert heartbeat["minimum"] == heartbeat["default"] == 30
            found = True
    assert found


def test_searching_execution_points_to_help_and_never_executes():
    with patch("jarvishep2.core.Jarvis2Core") as core:
        status, text = invoke("agent", "search", "run", "--json")
        assert status == 0
        result = json.loads(text)
        entry = next(row for row in result["results"] if row["kind"] == "command" and row["name"] == "run")
        assert entry["next_argv"] == ["Jarvis", "run", "--help"]
        core.assert_not_called()
    status, text = invoke("agent", "search", "unfindabletokenxyz", "--json")
    assert status == 0
    assert json.loads(text)["results"] == []


@pytest.mark.parametrize("topic", [
    "sampler.CSV.typo", "sampler.Dynesty.run_nested",
    "yaml.EnvReqs.V2.checkpoint.heartbeat", "tokens.no_such_token",
    "example.random.extra", "operas.unknown",
    "calculator.module.path", "calculator.execution.input.JSON.extra",
    "yaml.Sampling.Bounds.unknown", "yaml.Sampling.Variables.typo",
    "yaml.Calculators.Modules.typo",
    "yaml.Calculators.Modules.execution.typo",
    "yaml.Calculators.Modules.execution.typo.more",
    "calculator.execution.JSON.typo",
])
def test_unknown_or_leaf_topics_do_not_silently_return_parent_success(topic):
    status, text = invoke("man", topic, "--json")
    assert status == 2
    failure = json.loads(text)
    assert failure["ok"] is False
    assert failure["parent_manual_argv"]
    parent_status, parent_text = invoke(*failure["parent_manual_argv"][1:])
    assert parent_status == 0
    assert json.loads(parent_text)["diagnostics"] == []


def test_token_example_can_be_copied_as_yaml():
    status, text = invoke("man", "tokens", "--json")
    assert status == 0
    example = json.loads(text)["examples"][0]
    assert yaml.safe_load(example)["path"] == "&J/calculators/runtime"


@pytest.mark.parametrize("topic,flags", [
    ("calculator.execution.input", ["--direction", "output"]),
    ("calculator.execution.input.JSON", ["--type", "CSV"]),
    ("calculator.execution.JSON", ["--type", "CSV"]),
    ("yaml.Calculators.Modules.execution.output", ["--direction", "input"]),
    ("yaml.Calculators.Modules.execution.input.JSON", ["--type", "CSV"]),
])
def test_conflicting_manual_selectors_do_not_silently_choose_a_different_surface(topic, flags):
    status, text = invoke("man", topic, *flags, "--json")
    assert status == 2
    failure = json.loads(text)
    assert "Conflicting" in failure["summary"]
    assert invoke(*failure["parent_manual_argv"][1:])[0] == 0


@pytest.mark.parametrize("document,path", [
    ("Scan:\n  name: first\n  name: second\n", "$.Scan.name"),
    ("Calculators:\n  Modules:\n    - name: First\n      name: Second\n", "$.Calculators.Modules[0].name"),
    ("EnvReqs:\n  V2:\n    workers: 2\n    'workers': 4\n", "$.EnvReqs.V2.workers"),
])
def test_duplicate_explicit_keys_are_rejected_with_both_locations(tmp_path, document, path):
    card = tmp_path / "task.yaml"
    card.write_text(document, encoding="utf-8")
    with pytest.raises(TaskCardLoadError) as caught:
        load_task_yaml(str(card))
    error = caught.value
    assert error.code == "JV2-YAML-003"
    assert error.path == path
    assert error.location["line"] > error.related_locations[0]["line"]
    with patch("jarvishep2.core.Jarvis2Core") as core, patch("jarvishep2.redis_queue.RedisQueue") as redis:
        status, text = invoke("validate", str(card), "--json")
        assert status == 2
        item = json.loads(text)["issues"][0]
        assert item["path"] == path and item["code"] == "JV2-YAML-003"
        assert item["location"] == error.location
        core.assert_not_called()
        redis.assert_not_called()


@pytest.mark.parametrize("document", [
    "base: &base {workers: 2, batch_size: 256}\nEnvReqs: {V2: {<<: *base, workers: 4}}\n",
    "a: &a {x: 1}\nb: &b {x: 2, y: 3}\nsettings: {<<: [*a, *b], y: 4}\n",
    "a: &a {x: 1}\nb: *a\n",
    "a: &a {x: 1}\nsettings: {<<: *a, '<<': literal_data}\n",
])
def test_legal_anchor_merge_overrides_keep_standard_safe_loader_semantics(document):
    assert yaml.load(document, Loader=UniqueKeySafeLoader) == yaml.safe_load(document)


def test_duplicate_in_selected_default_yaml_is_diagnosed_before_merge(tmp_path):
    defaults = tmp_path / "environment.yaml"
    defaults.write_text("EnvReqs:\n  V2:\n    workers: 2\n    workers: 4\n", encoding="utf-8")
    card = tmp_path / "task.yaml"
    card.write_text("Scan: {name: test}\nEnvReqs:\n  Check_default_dependencies:\n    required: true\n    default_yaml_path: environment.yaml\n", encoding="utf-8")
    status, text = invoke("validate", str(card), "--json")
    assert status == 2
    item = json.loads(text)["issues"][0]
    assert item["location"]["file"] == str(defaults)
    assert item["path"] == "$.EnvReqs.V2.workers"
