#!/usr/bin/env python3
from __future__ import annotations

import os
import sys
import json

from jarvishep.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold
from jarvishep.versioning import render_logo_with_version


def _render_version_banner() -> str:
    logo_path = os.path.join(os.path.dirname(__file__), "card", "logo")
    banner = render_logo_with_version(logo_path)
    links_text = _render_document_links()
    if links_text:
        banner = f"{banner}\n\n{links_text}"
    return banner


def _normalize_public_value(value: object) -> str:
    text = str(value).strip() if value is not None else ""
    if not text:
        return ""
    if text.upper() in {"TODO", "N/A", "NA", "NONE", "NULL", "UNKNOWN", "TBD", "TBA"}:
        return ""
    return text


def _render_resources_block(documents: dict, paper: dict) -> str:
    candidates = [
        ("Online docs", documents.get("online_docs")),
        ("User manual", documents.get("user_manual")),
        ("Tutorials", documents.get("tutorials")),
        ("API reference", documents.get("api_reference")),
        ("Homepage", documents.get("homepage")),
        ("Paper", paper.get("url")),
    ]
    public_links = [
        (label, _normalize_public_value(value))
        for label, value in candidates
    ]
    public_links = [(label, value) for label, value in public_links if value]
    if not public_links:
        return ""

    lines = ["Resources:"]
    lines.extend(f"\t{label}:\t{value}" for label, value in public_links)
    return "\n".join(lines)


def _normalize_reference_entry(entry: object) -> dict:
    if not isinstance(entry, dict):
        return {}

    return {
        "caption": _normalize_public_value(entry.get("caption") or entry.get("summary")),
        "title": _normalize_public_value(entry.get("title")),
        "author": _normalize_public_value(entry.get("author") or entry.get("authors")),
        "arxiv": _normalize_public_value(entry.get("arxiv")),
        "doi": _normalize_public_value(entry.get("doi") or entry.get("doi_link")),
    }


def _render_references_block(payload: dict, paper: dict) -> str:
    references = payload.get("references", {})
    if not isinstance(references, dict):
        references = {}

    section_specs = [
        ("Jarvis-HEP", references.get("jarvis_hep")),
        ("Samplers", references.get("builtin_scanners")),
    ]
    rendered_sections = []
    for section_name, entries in section_specs:
        if not isinstance(entries, list):
            entries = []
        normalized = [_normalize_reference_entry(entry) for entry in entries]
        normalized = [entry for entry in normalized if entry.get("title")]
        if not normalized:
            continue

        block_lines = [f"\t{section_name}:"]
        for idx, entry in enumerate(normalized, start=1):
            title_line = f"\t\t[{idx}] "
            if entry["caption"]:
                title_line += entry["caption"]
            else:
                title_line += entry["title"]
            block_lines.append(title_line)
            block_lines.append(f"\t\t\tTitle:\t{entry['title']}")
            if entry["author"]:
                block_lines.append(f"\t\t\tAuthor:\t{entry['author']}")
            if entry["arxiv"]:
                block_lines.append(f"\t\t\tarXiv:\t{entry['arxiv']}")
            if entry["doi"]:
                block_lines.append(f"\t\t\tDOI:\t{entry['doi']}")
        rendered_sections.append("\n".join(block_lines))

    # Backward-compatibility with legacy "paper" schema.
    if not rendered_sections:
        legacy_title = _normalize_public_value(paper.get("title"))
        legacy_url = _normalize_public_value(paper.get("url"))
        if legacy_title or legacy_url:
            legacy_block = ["\tJarvis-HEP:"]
            legacy_block.append(f"\t\t[1] {legacy_title or 'Jarvis-HEP Paper'}")
            if legacy_url:
                legacy_block.append(f"\t\t\tLink:\t{legacy_url}")
            rendered_sections.append("\n".join(legacy_block))

    if not rendered_sections:
        return ""

    return "References:\n" + "\n".join(rendered_sections)


def _render_document_links() -> str:
    links_path = os.path.join(os.path.dirname(__file__), "card", "document_links.json")
    if not os.path.exists(links_path):
        return ""

    try:
        with open(links_path, "r", encoding="utf-8") as f:
            payload = json.load(f)
    except Exception:
        return ""

    if not isinstance(payload, dict):
        return ""

    documents = payload.get("documents", {})
    paper = payload.get("paper", {})
    if not isinstance(documents, dict):
        documents = {}
    if not isinstance(paper, dict):
        paper = {}

    sections = [
        _render_resources_block(documents, paper),
        _render_references_block(payload, paper),
    ]
    sections = [section for section in sections if section]
    return "\n\n".join(sections)


def _version_fast_path(argv) -> int | None:
    if "--version" not in argv and "-v" not in argv:
        return None

    print(_render_version_banner())
    return 0


def _mkproject_fast_path(argv) -> int | None:
    if "--mkproject" not in argv:
        return None

    conflicts = [
        flag for flag in (
            "--plot",
            "--convert",
            "--monitor",
            "--check-modules",
        )
        if flag in argv
    ]
    if conflicts:
        print(f"[Jarvis-HEP] --mkproject cannot be combined with: {' '.join(conflicts)}")
        return 2

    i = argv.index("--mkproject")
    if i + 1 >= len(argv) or argv[i + 1].startswith("-"):
        print("[Jarvis-HEP] Missing project name for --mkproject")
        return 2

    # Fast-path mode is intentionally strict: only `--mkproject <name>` is accepted.
    extra_tokens = argv[1:i] + argv[i + 2 :]
    if extra_tokens:
        unsupported_opts = [tok for tok in extra_tokens if tok.startswith("-")]
        if unsupported_opts:
            print(
                "[Jarvis-HEP] --mkproject does not accept option(s): "
                + " ".join(unsupported_opts)
            )
        else:
            print(
                "[Jarvis-HEP] --mkproject does not accept extra arguments: "
                + " ".join(extra_tokens)
            )
        return 2

    project_name = argv[i + 1]
    try:
        project_root = create_project_scaffold(project_name, cwd=os.getcwd())
    except ValueError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except FileExistsError as exc:
        print(f"[Jarvis-HEP] Project directory already exists: {exc}")
        return 1

    print(f"[Jarvis-HEP] Project scaffold created at: {project_root}")
    print(f"[Jarvis-HEP] Created folders: {', '.join(PROJECT_SUBDIRS)}")
    return 0


def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv

    version_code = _version_fast_path(argv)
    if version_code is not None:
        return version_code

    fast_path_code = _mkproject_fast_path(argv)
    if fast_path_code is not None:
        return fast_path_code

    from jarvishep.core import Core

    jc = Core()
    jc.initialization()
    if getattr(jc.args, "version", False):
        print(_render_version_banner())
    elif getattr(jc.args, "mkproject", None):
        jc.mkproject()
    elif jc.scan_mode:
        jc.run_sampling()
    elif jc.args.cvtDB:
        jc.convert()
    elif jc.args.plot:
        jc.plot()
    elif jc.args.monitor:
        jc.monitor()
    else:
        import argparse as _ap
        _ap.ArgumentParser(prog="Jarvis").print_help()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
