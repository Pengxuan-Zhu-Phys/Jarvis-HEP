"""Generate feature-request drafts without sending mail or modifying files."""

from __future__ import annotations

import argparse
import json
import platform
import sys
from collections.abc import Sequence
from typing import Any

__all__ = ["feature_request_contract", "build_request", "dispatch_request"]

_RECIPIENT = "Jarvis.HEP.AI@mail.grokbot.com"
# Preserve the spelling of the request marker used by the Jarvis-PLOT contract.
_SUBJECT_MARKER = "Resquest@Jarvis-HEP"
_FIELD_PLACEHOLDERS = {
    "title": "[Short title for the missing capability]",
    "description": "[Explain the goal and requested behavior.]",
    "yaml": "[Paste the smallest task card that demonstrates the request.]",
    "commands": "[Paste the exact CLI commands.]",
    "diagnostics": "[Paste the exit codes and relevant CLI diagnostics.]",
    "actual": "[Describe what the installed package actually does.]",
    "expected": "[Describe the desired result and why it matters.]",
    "constraints": (
        "[State compatibility requirements, scientific assumptions, execution limits,\n"
        "and which existing behavior or data must be preserved.]"
    ),
    "installation_method": "[pip / uv / source / other; include relevant package versions]",
}
_DRAFT_READY_MEANING = (
    "All required text sections were explicitly filled. Their contents remain unverified; "
    "this does not establish feature support, a valid reproducer, or that email was sent."
)
_REQUEST_TEMPLATE = """To: {recipient}
Subject: {subject}

Jarvis-HEP version: {installed_version}
Environment:
  OS: {os}
  Python: {python}
  Installation method: {installation_method}

Request title: {title}
Goal / requested behavior:
{description}

Minimal YAML:
{yaml}

Relevant data fields or anonymized sample:
[Optional: include only the fields or small sample needed to reproduce the behavior.]

Exact commands and diagnostics:
Commands:
{commands}

Diagnostics / exit codes:
{diagnostics}

Current result:
{actual}

Expected behavior:
{expected}

Constraints:
{constraints}

Before sharing, remove credentials, sensitive paths, and private data.
"""


def feature_request_contract() -> dict[str, Any]:
    """Return the installed CLI's feedback protocol and copyable draft template."""
    return {
        "recipient": _RECIPIENT,
        "subject_marker": _SUBJECT_MARKER,
        "automatic_send": False,
        "maintainer_session": "Report directly in the current conversation.",
        "other_user_session": "Provide a draft for the user to submit by email.",
        "unknown_session": "Report in chat and provide email submission instructions.",
        "identity_rule": "Do not infer maintainer identity from local paths or machine names.",
        "required_details": [
            "installed Jarvis-HEP version from Jarvis --version",
            "OS, Python version, installation method, and relevant dependency versions",
            "goal and required behavior",
            "minimal YAML and relevant data fields or anonymized sample",
            "exact commands, exit codes, and diagnostics",
            "current result and expected behavior",
            "scientific assumptions, execution limits, and compatibility constraints",
        ],
        "evidence_rule": (
            "Check the installed man, schema, and relevant dependency CLI before reporting "
            "a missing capability. Distinguish unsupported behavior, invalid YAML, and "
            "environment or dependency failures. Do not present a workaround as native support."
        ),
        "sharing_rule": "Remove credentials, sensitive paths, and private data before sharing.",
        "send_rule": (
            "This command only outputs a draft. It does not send email. An agent must not "
            "send the draft through another tool without explicit user authorization."
        ),
        "draft_command": ["Jarvis", "agent", "request", "--json"],
        "required_fields": list(_FIELD_PLACEHOLDERS),
        "input_fields": [
            {"name": name, "flag": f"--{name.replace('_', '-')}", "format": "inline_text"}
            for name in _FIELD_PLACEHOLDERS
        ],
        "input_rule": (
            "All input parameters are inline text. Paths are not read as files, YAML is not "
            "loaded, and command or diagnostic text is never executed. If a section is not "
            "applicable, explicitly enter n/a with an explanation; applicability is not guessed."
        ),
        "draft_ready_meaning": _DRAFT_READY_MEANING,
        "template": _REQUEST_TEMPLATE,
    }


def build_request(
    title: str = "",
    description: str = "",
    *,
    yaml: str = "",
    commands: str = "",
    diagnostics: str = "",
    actual: str = "",
    expected: str = "",
    constraints: str = "",
    installation_method: str = "",
) -> dict[str, Any]:
    """Populate a draft from inline text; never load input files or run commands."""
    from jarvishep2.versioning import get_runtime_version
    installed_version = get_runtime_version()

    # A title occupies one header line. Keep all other supplied text verbatim.
    raw_fields = {
        "title": title, "description": description, "yaml": yaml, "commands": commands,
        "diagnostics": diagnostics, "actual": actual, "expected": expected,
        "constraints": constraints, "installation_method": installation_method,
    }
    supplied = {name: str(value) if value is not None else "" for name, value in raw_fields.items()}
    missing_fields = [name for name, text in supplied.items() if not text.strip()]
    rendered = {
        name: text if text.strip() else _FIELD_PLACEHOLDERS[name]
        for name, text in supplied.items()
    }
    request_title = " ".join(rendered["title"].split())
    request_description = rendered["description"]
    subject = f"{_SUBJECT_MARKER}: {request_title}"
    environment = {
        "os": f"{platform.system()} {platform.release()}".strip(),
        "python": platform.python_version(),
        "installation_method": rendered["installation_method"],
    }
    draft = _REQUEST_TEMPLATE.format(
        recipient=_RECIPIENT,
        subject=subject,
        installed_version=installed_version,
        os=environment["os"],
        python=environment["python"],
        title=request_title,
        description=request_description,
        **{name: text for name, text in rendered.items() if name not in {"title", "description"}},
    )
    return {
        "installed_version": installed_version,
        "recipient": _RECIPIENT,
        "subject_marker": _SUBJECT_MARKER,
        "subject": subject,
        "automatic_send": False,
        "title": request_title,
        "description": request_description,
        **{name: text for name, text in supplied.items() if name not in {"title", "description"}},
        "environment": environment,
        "missing_fields": missing_fields,
        "draft_ready": not missing_fields,
        "draft_ready_meaning": _DRAFT_READY_MEANING,
        "body": draft,
        "feature_request": feature_request_contract(),
    }


class _RequestUsageError(ValueError):
    pass


class _RequestParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise _RequestUsageError(message)


def _build_parser() -> argparse.ArgumentParser:
    parser = _RequestParser(
        prog="Jarvis agent request",
        description=(
            "Print a feature-request email draft. No email is sent and no file is written. "
            "Every input is inline text; no input file is read and no command is executed."
        ),
        epilog=(
            "Omitted or blank required sections remain placeholders and are listed in missing_fields. "
            "Explicitly enter n/a with a reason where appropriate. draft_ready means the text "
            "sections are filled; it does not verify their contents or send the draft."
        ),
        add_help=False,
        allow_abbrev=False,
    )
    parser.add_argument("--title", default="", metavar="TEXT", help="Short title as inline text.")
    parser.add_argument("--description", default="", metavar="TEXT", help="Goal and requested behavior as inline text.")
    fields = {
        "yaml": "Minimal YAML as inline text; this is not a file path option.",
        "commands": "Exact CLI commands as inline text; commands are not executed.",
        "diagnostics": "Diagnostics and exit codes as inline text.",
        "actual": "Current result as inline text.",
        "expected": "Desired behavior as inline text.",
        "constraints": "Scientific, execution, and compatibility constraints as inline text.",
        "installation-method": "Installation method and relevant package versions as inline text.",
    }
    for flag, help_text in fields.items():
        parser.add_argument(f"--{flag}", default="", metavar="TEXT", help=help_text)
    parser.add_argument("--json", action="store_true", help="Print one JSON object, including failures.")
    parser.add_argument("-h", "--help", action="store_true", help="Show this command's usage.")
    return parser


def _emit_error(message: str, *, code: str, json_output: bool) -> None:
    if json_output:
        print(json.dumps(
            {"ok": False, "kind": "agent.request", "contract_version": 1, "error": {"code": code, "message": message}},
            ensure_ascii=False,
        ))
    else:
        print(f"Jarvis agent request: {message}", file=sys.stderr)


def dispatch_request(argv: Sequence[str] | None = None) -> int:
    """Print a copyable draft or a machine-readable response (0 / 1 / 2)."""
    raw_argv = list(argv) if argv is not None else []
    json_output = any(arg == "--json" or arg.startswith("--json=") for arg in raw_argv)
    try:
        parser = _build_parser()
        args = parser.parse_args(raw_argv)
        if args.help:
            if args.json:
                print(json.dumps(
                    {
                        "ok": True,
                        "kind": "agent.request.help",
                        "contract_version": 1,
                        "data": {"usage": parser.format_help(), "feature_request": feature_request_contract()},
                    },
                    ensure_ascii=False,
                ))
            else:
                print(parser.format_help(), end="")
            return 0
        data = build_request(
            title=args.title, description=args.description, yaml=args.yaml,
            commands=args.commands, diagnostics=args.diagnostics, actual=args.actual,
            expected=args.expected, constraints=args.constraints,
            installation_method=args.installation_method,
        )
        if args.json:
            print(json.dumps({"ok": True, "kind": "agent.request", "contract_version": 1, "data": data}, ensure_ascii=False))
        else:
            print(data["body"], end="")
        return 0
    except _RequestUsageError as exc:
        _emit_error(str(exc), code="usage_error", json_output=json_output)
        return 2
    except Exception as exc:
        _emit_error(str(exc), code="request_error", json_output=json_output)
        return 1
