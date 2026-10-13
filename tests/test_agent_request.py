"""Feedback drafts must carry evidence and remain explicit about sending boundaries."""

from __future__ import annotations

import json

import pytest

from jarvishep2.cli import agent_request


def _complete_request_values():
    return {
        "title": "新的似然接口",
        "description": "希望支持新接口\n保留现有模型。",
        "yaml": "Scan: {name: example}\nSampling:\n  Method: Random\n  Bounds: {point_number: 1}\n",
        "commands": "Jarvis validate example.yaml --json\nJarvis man yaml --json",
        "diagnostics": '{"ok": false, "issues": [{"code": "JV2-SCH-001"}]}\nExit: 2',
        "actual": "字段被拒绝\n没有执行扫描。",
        "expected": "支持所述接口\n返回明确结果。",
        "constraints": "不改动物理模型\n保持现有配置兼容。",
        "installation_method": "pip install Jarvis-HEP==2.0.15\nExisting dependency versions attached separately.",
    }


def test_feature_request_contract_preserves_plot_style_feedback_boundary():
    contract = agent_request.feature_request_contract()
    assert contract["recipient"] == "Jarvis.HEP.AI@mail.grokbot.com"
    assert contract["subject_marker"] == "Resquest@Jarvis-HEP"
    assert contract["automatic_send"] is False
    assert "current conversation" in contract["maintainer_session"]
    assert "user to submit" in contract["other_user_session"]
    assert "chat" in contract["unknown_session"]
    assert "paths or machine names" in contract["identity_rule"]
    assert "explicit user authorization" in contract["send_rule"]
    assert any("installed Jarvis-HEP version" in detail for detail in contract["required_details"])
    assert "credentials" in contract["sharing_rule"]
    assert "inline text" in contract["input_rule"]
    assert "never executed" in contract["input_rule"]
    assert "n/a" in contract["input_rule"]


def test_default_draft_includes_installed_environment_and_all_reproducer_sections(monkeypatch):
    monkeypatch.setattr("jarvishep2.versioning.get_runtime_version", lambda: "9.8.7-test")
    monkeypatch.setattr(agent_request.platform, "system", lambda: "TestOS")
    monkeypatch.setattr(agent_request.platform, "release", lambda: "42")
    monkeypatch.setattr(agent_request.platform, "python_version", lambda: "3.12.99")
    data = agent_request.build_request()
    assert data["installed_version"] == "9.8.7-test"
    assert data["environment"]["os"] == "TestOS 42"
    assert data["environment"]["python"] == "3.12.99"
    assert data["automatic_send"] is False
    assert data["draft_ready"] is False
    assert data["missing_fields"] == [
        "title", "description", "yaml", "commands", "diagnostics", "actual",
        "expected", "constraints", "installation_method",
    ]
    body = data["body"]
    for text in (
        "To: Jarvis.HEP.AI@mail.grokbot.com",
        "Subject: Resquest@Jarvis-HEP:",
        "Jarvis-HEP version: 9.8.7-test",
        "TestOS 42",
        "3.12.99",
        "Installation method:",
        "Goal / requested behavior:",
        "Minimal YAML:",
        "anonymized sample:",
        "Exact commands and diagnostics:",
        "Current result:",
        "Expected behavior:",
        "Constraints:",
        "remove credentials",
    ):
        assert text in body


def test_request_title_cannot_create_extra_headers_and_description_remains_multiline():
    data = agent_request.build_request("New sampler\nBcc: another@example.org", "目标：保留模型\n第二行")
    assert data["subject"] == "Resquest@Jarvis-HEP: New sampler Bcc: another@example.org"
    assert "\nBcc:" not in data["body"]
    assert "目标：保留模型\n第二行" in data["body"]


def test_contract_returns_independent_data():
    first = agent_request.feature_request_contract()
    first["required_details"].clear()
    assert agent_request.feature_request_contract()["required_details"]


def test_filled_evidence_is_in_its_own_section_with_multiline_unicode_unchanged():
    values = _complete_request_values()
    data = agent_request.build_request(**values)
    assert data["missing_fields"] == []
    assert data["draft_ready"] is True
    assert data["automatic_send"] is False
    assert "unverified" in data["draft_ready_meaning"]
    assert "email was sent" in data["draft_ready_meaning"]
    sections = {
        "description": "Goal / requested behavior:",
        "yaml": "Minimal YAML:",
        "commands": "Commands:",
        "diagnostics": "Diagnostics / exit codes:",
        "actual": "Current result:",
        "expected": "Expected behavior:",
        "constraints": "Constraints:",
    }
    for field, heading in sections.items():
        assert data[field] == values[field]
        assert f"{heading}\n{values[field]}\n" in data["body"]
    assert data["environment"]["installation_method"] == values["installation_method"]
    assert f"Installation method: {values['installation_method']}\n" in data["body"]
    assert "[Paste the smallest task card" not in data["body"]


def test_readiness_tracks_omitted_and_blank_evidence_without_guessing_applicability():
    values = _complete_request_values()
    values["yaml"] = " \n\t"
    values["commands"] = ""
    values["installation_method"] = " "
    data = agent_request.build_request(**values)
    assert data["draft_ready"] is False
    assert data["missing_fields"] == ["yaml", "commands", "installation_method"]
    assert "[Paste the smallest task card" in data["body"]
    assert "[Paste the exact CLI commands.]" in data["body"]


def test_explicit_not_applicable_text_is_retained_without_being_inferred():
    values = _complete_request_values()
    values["yaml"] = "n/a: this request concerns discoverability before a card exists."
    values["diagnostics"] = "n/a: missing documentation; no failing command exists."
    data = agent_request.build_request(**values)
    assert data["draft_ready"] is True
    assert data["yaml"] == values["yaml"]
    assert data["diagnostics"] == values["diagnostics"]
    assert values["yaml"] in data["body"]


def test_request_human_output_is_complete_copyable_draft(capsys):
    assert agent_request.dispatch_request([]) == 0
    captured = capsys.readouterr()
    assert captured.out.startswith("To: Jarvis.HEP.AI@mail.grokbot.com\nSubject: Resquest@Jarvis-HEP:")
    assert "Minimal YAML:" in captured.out
    assert captured.err == ""


def test_request_json_output_preserves_user_text_and_has_no_stderr(capsys):
    assert agent_request.dispatch_request(["--title", "采样器支持", "--description", "需要新能力\n约束说明", "--json"]) == 0
    captured = capsys.readouterr()
    response = json.loads(captured.out)
    assert response["ok"] is True
    assert response["kind"] == "agent.request"
    assert response["data"]["title"] == "采样器支持"
    assert response["data"]["description"] == "需要新能力\n约束说明"
    assert response["data"]["automatic_send"] is False
    assert response["data"]["draft_ready"] is False
    assert captured.err == ""


def test_cli_can_populate_every_required_section_without_postprocessing_the_draft(capsys):
    values = _complete_request_values()
    argv = ["--json"]
    for name, text in values.items():
        argv.extend([f"--{name.replace('_', '-')}", text])
    assert agent_request.dispatch_request(argv) == 0
    captured = capsys.readouterr()
    response = json.loads(captured.out)
    assert response["ok"] is True
    assert response["data"]["draft_ready"] is True
    assert response["data"]["missing_fields"] == []
    for name, text in values.items():
        assert response["data"][name] == text
    assert captured.err == ""


def test_path_and_command_inputs_are_literal_data_without_file_read_or_execution(tmp_path, capsys):
    source = tmp_path / "private.yaml"
    source.write_text("PRIVATE CONTENT THAT MUST NOT BE INCLUDED", encoding="utf-8")
    executed = tmp_path / "command-executed.txt"
    values = _complete_request_values()
    values["yaml"] = str(source)
    values["commands"] = f'python3 -c "from pathlib import Path; Path({str(executed)!r}).write_text(\'executed\')"'
    argv = ["--json"]
    for name, text in values.items():
        argv.extend([f"--{name.replace('_', '-')}", text])
    assert agent_request.dispatch_request(argv) == 0
    response = json.loads(capsys.readouterr().out)
    assert response["data"]["yaml"] == str(source)
    assert str(source) in response["data"]["body"]
    assert "PRIVATE CONTENT THAT MUST NOT BE INCLUDED" not in response["data"]["body"]
    assert response["data"]["commands"] == values["commands"]
    assert not executed.exists()
    assert list(tmp_path.iterdir()) == [source]


def test_repeated_title_keeps_last_value_on_one_header_line(capsys):
    assert agent_request.dispatch_request([
        "--title", "first", "--title", "last\nBcc: another@example.org", "--json",
    ]) == 0
    data = json.loads(capsys.readouterr().out)["data"]
    assert data["title"] == "last Bcc: another@example.org"
    assert "\nBcc:" not in data["body"]
    assert data["automatic_send"] is False


@pytest.mark.parametrize("argv", [
    ["--json", "--title"], ["--json", "--send"], ["--json=bad"],
    ["--json", "--yaml"], ["--json", "--diagnostics"], ["--json", "--installation-method"],
])
def test_request_json_usage_errors_are_single_objects_with_nonzero_status(argv, capsys):
    assert agent_request.dispatch_request(argv) == 2
    captured = capsys.readouterr()
    response = json.loads(captured.out)
    assert response["ok"] is False
    assert response["error"]["code"] == "usage_error"
    assert response["error"]["message"]
    assert captured.err == ""


def test_request_internal_error_is_json_and_does_not_claim_draft_success(monkeypatch, capsys):
    def fail(**kwargs):
        raise RuntimeError("environment unavailable")

    monkeypatch.setattr(agent_request, "build_request", fail)
    assert agent_request.dispatch_request(["--json"]) == 1
    captured = capsys.readouterr()
    response = json.loads(captured.out)
    assert response["ok"] is False
    assert response["error"] == {"code": "request_error", "message": "environment unavailable"}
    assert captured.err == ""


def test_request_json_help_keeps_stdout_parseable(capsys):
    assert agent_request.dispatch_request(["--help", "--json"]) == 0
    captured = capsys.readouterr()
    response = json.loads(captured.out)
    assert response["ok"] is True
    assert "No email is sent" in response["data"]["usage"]
    assert "Every input is inline text" in response["data"]["usage"]
    assert "--yaml TEXT" in response["data"]["usage"]
    assert "missing_fields" in response["data"]["usage"]
    assert response["data"]["feature_request"]["automatic_send"] is False
    assert captured.err == ""
