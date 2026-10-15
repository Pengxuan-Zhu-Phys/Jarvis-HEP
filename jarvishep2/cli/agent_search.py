"""Bounded keyword discovery over the installed manual and CLI help metadata."""

from __future__ import annotations

import argparse
from collections import deque
import json
import shlex
import sys
from typing import Any

from jarvishep2.versioning import get_runtime_version


def _manual_args(command: str) -> list[str] | None:
    tokens = shlex.split(command)
    if tokens[:2] != ["Jarvis", "man"] or "agent-readme" in tokens:
        return None
    return [token for token in tokens[2:] if token != "--json"]


def search_index() -> tuple[list[dict[str, Any]], list[str]]:
    from jarvishep2.cli.agent import command_catalog
    from jarvishep2.cli.man import build_man_parser, resolve_man_request

    entries: list[dict[str, Any]] = []
    notes: list[str] = []
    pending = deque([[], ["yaml"], ["sampler"], ["calculator"], ["operas"], ["tokens"], ["example"], ["sampler.mcmc-runtime"]])
    visited: set[tuple[str, ...]] = set()
    emitted: set[tuple[str, str, str]] = set()
    parser = build_man_parser()
    while pending and len(visited) < 160:
        tokens = pending.popleft()
        if tuple(tokens) in visited:
            continue
        visited.add(tuple(tokens))
        try:
            args = parser.parse_args(tokens)
            page = resolve_man_request(args.tokens, io_type=args.io_type, direction=args.direction, full=args.full)
        except (KeyError, ValueError):
            continue
        argv = ["Jarvis", "man", *tokens, "--json"]
        path, context = str(page.get("path") or "$"), dict(page.get("context") or {})
        rows = [{"name": tokens[0] if tokens else "manual", "description": page.get("summary", ""), "kind": "topic"}]
        rows += [{**row, "kind": "field"} for row in page.get("keys", [])]
        rows += [{**row, "kind": "listed_item"} for row in page.get("list_rows", [])]
        rows += [{**row, "name": row["method"], "kind": "sampler", "man": f"Jarvis man sampler.{row['method']}"} for row in page.get("methods", [])]
        rows += [{**row, "name": row.get("type", ""), "kind": "file_type", "man": row.get("command")} for row in page.get("type_rows", [])]
        for row in rows:
            name = str(row.get("name") or "")
            if name == "agent-readme":
                continue
            target = _manual_args(str(row.get("man") or ""))
            if target is not None:
                pending.append(target)
            kind = row["kind"]
            entry_path = path + "." + name if kind == "field" else path
            identity = (kind, entry_path if kind == "field" else name, json.dumps(context, sort_keys=True))
            if identity in emitted:
                continue
            emitted.add(identity)
            description = str(row.get("description") or row.get("summary") or "")
            entries.append({
                "name": name, "kind": kind, "path": entry_path, "context": context,
                "description": description[:600], "status": page.get("status"),
                "next_argv": ["Jarvis", "man", *target, "--json"] if target is not None else argv,
                "source": "installed_manual",
            })
        for command in page.get("see_also", []):
            target = _manual_args(str(command))
            if target is not None:
                pending.append(target)
    if pending:
        notes.append("The installed manual traversal reached 160 pages; inspect additional delegated help directly.")
    for row in command_catalog():
        entries.append({"name": row["command"], "kind": "command", "path": None,
                        "context": {}, "description": row["effects_detail"],
                        "effect": row["effect"], "next_argv": row["help_argv"],
                        "source": "installed_cli_help"})
    return entries, notes


def search(query: str, limit: int = 8) -> dict[str, Any]:
    terms = list(dict.fromkeys(query.casefold().split()))
    entries, notes = search_index()
    matches = []
    for entry in entries:
        name = entry["name"].casefold()
        path = str(entry["path"] or "").casefold()
        description = entry["description"].casefold()
        context = json.dumps(entry["context"]).casefold()
        matched = [term for term in terms if term in name or term in path or term in description or term in context]
        if not matched:
            continue
        score = sum(8 * (term == name) + 4 * (term in name) + 2 * (term in path) + (term in description) + (term in context) for term in matched)
        matches.append({**entry, "score": score, "matched_terms": matched})
    matches.sort(key=lambda entry: (-len(entry["matched_terms"]), -entry["score"], entry["name"], json.dumps(entry["context"], sort_keys=True)))
    return {"ok": True, "kind": "agent.search", "contract_version": 1,
            "installed_version": get_runtime_version(), "query": query,
            "matching_entries": len(matches), "results": matches[:limit],
            "notes": notes,
            "scope": "Keyword matching over installed man/help metadata. Results are evidence pointers, not a claim that a requested behavior is supported. No match is not proof of unsupported behavior.",
            "execution_performed": False}


class _UsageError(ValueError):
    pass


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        raise _UsageError(message)


def dispatch_search(argv: list[str] | None = None) -> int:
    parser = _Parser(prog="Jarvis agent search", add_help=False, allow_abbrev=False,
                     description="Find installed manual topics/fields and command help by English keywords or literal YAML names; next_argv only queries documentation.")
    parser.add_argument("query", nargs="?", help="Quote multiple keywords, e.g. 'checkpoint heartbeat'")
    parser.add_argument("--limit", type=int, default=8, help="Maximum results, 1..50 (default 8)")
    parser.add_argument("--json", action="store_true", help="Machine-readable bounded results")
    parser.add_argument("-h", "--help", action="store_true")
    raw = list(argv or [])
    as_json = any(token == "--json" or token.startswith("--json=") for token in raw)
    try:
        args = parser.parse_args(raw)
        if args.help:
            print(json.dumps({"contract_version": 1, "help": parser.format_help()}) if args.json else parser.format_help())
            return 0
        if not args.query or not args.query.strip():
            raise _UsageError("A non-empty keyword query is required.")
        if not 1 <= args.limit <= 50:
            raise _UsageError("--limit must be between 1 and 50.")
        result = search(args.query, args.limit)
        if args.json:
            print(json.dumps(result, indent=2, sort_keys=True))
        else:
            print(result["scope"])
            for row in result["results"]:
                print(f"{row['name']} ({row['kind']}): {row['description']}\n  {shlex.join(row['next_argv'])}")
            if not result["results"]:
                print("No keyword matches; inspect the owning CLI help before concluding unsupported behavior.")
        return 0
    except Exception as exc:
        status = 2 if isinstance(exc, _UsageError) else 1
        message = str(exc)
        if as_json:
            print(json.dumps({"ok": False, "kind": "agent.search", "contract_version": 1,
                              "error": {"code": "usage_error" if status == 2 else "search_error", "message": message}}))
        else:
            print(f"Jarvis agent search: {message}", file=sys.stderr)
        return status
