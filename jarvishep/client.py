#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys

from jarvishep.log_kv import format_two_column_log
from jarvishep.official_project_library import (
    OfficialLibraryError,
    OfficialProjectFetchError,
    OfficialProjectNotFoundError,
    fetch_official_project,
    get_official_project,
    list_official_projects,
)
from jarvishep.project_packager import (
    ProjectPackError,
    create_project_pack_manifest,
    create_project_package,
    create_project_package_from_manifest,
)
from jarvishep.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold
from jarvishep.versioning import render_logo_with_version

RESET = "\033[0m"
BOLD = "\033[1m"


_TOP_LEVEL_HELP_TEXT = """Usage:
  Jarvis [file] [options]
  Jarvis project <command> [arguments]
  Jarvis portal formats

Jarvis Program Help Center

Main entry points:
  file                  Run Jarvis with a YAML input file
  project               Manage Jarvis standalone projects
  portal formats        Show Jarvis-Portal supported calculator IO formats

General options:
  -h, --help            Show this help message and exit
  -d, --debug           Run Jarvis-HEP in debug mode
  -v, --version         Print version and runtime package information
  --refs                Print full reference information

Workflow options:
  --plot                Run plotting mode
  --convert             Convert sample.hdf5 into CSV format
  --monitor             Start a real-time resource monitor
  --resume              Resume from the latest checkpoint without prompting
  --check-modules       Run calculator/module checks
  --skip-library-installation
                        Skip library installation
  --skip-draw-flowchart
                        Skip JarvisPLOT flowchart.png rendering

Hint:
  Run `Jarvis project -h` to see project workflow commands.
  Run `Jarvis portal -h` to see Portal inspection commands.
"""

_PROJECT_HELP_TEXT = """Usage:
  Jarvis project <command> [arguments]

Manage Jarvis standalone projects.

Commands:
  create <name>      Create a new local project scaffold
  pack [path]        Pack a local project for sharing, reproduction, or full export
  browse             List verified projects in the official Jarvis library
  fetch <name>       Fetch an official project into a local directory
  info <name>        Show details for an official project

Scope:
  create and pack work on local projects.
  browse, fetch, and info work on the official Jarvis project library.

Examples:
  Jarvis project create MyProject
  Jarvis project browse
  Jarvis project fetch Example_Bridson
  Jarvis project info Example_Bridson
"""

_PROJECT_PACK_HELP_TEXT = """Usage:
  Jarvis project pack [path] [--share | --repro | --full] [--man]
  Jarvis project pack <pack_manifest.yaml>

Pack a local standalone project.

Modes:
  --share           Create a lightweight shareable bundle
  --repro           Create a reproducible bundle
  --full            Create a full bundle
  --man             Write a pack manifest only; do not create an archive

Notes:
  If no mode is specified, `--share` is used by default.
"""

_PROJECT_SUBCOMMAND_HELP = {
    "create": """Usage:
  Jarvis project create <name>

Create a new local project scaffold.
""",
    "browse": """Usage:
  Jarvis project browse

List verified projects in the official Jarvis library.
""",
    "fetch": """Usage:
  Jarvis project fetch <name>

Fetch an official project into a local directory.
""",
    "info": """Usage:
  Jarvis project info <name>

Show details for an official project.
""",
}

_PORTAL_HELP_TEXT = """Usage:
  Jarvis portal formats

Inspect Jarvis-Portal integration.

Commands:
  formats            Show calculator IO formats registered by Jarvis-Portal
"""

_HELP_FLAGS = {"-h", "--help"}
_PACK_MODE_FLAGS = {
    "--share": "share",
    "--repro": "repro",
    "--full": "full",
}
_PACK_MANIFEST_FLAG = "--man"


def _render_version_banner() -> str:
    logo_path = os.path.join(os.path.dirname(__file__), "card", "logo")
    banner = render_logo_with_version(logo_path)
    links_text = _render_document_links(
        include_resources=True,
        include_compact_reference=True,
        include_references=False,
    )
    if links_text:
        banner = f"{banner}\n\n{links_text}"
    return banner


def _render_refs_text() -> str:
    logo_path = os.path.join(os.path.dirname(__file__), "card", "logo")
    banner = render_logo_with_version(logo_path)
    links_text = _render_document_links(
        include_resources=False,
        include_compact_reference=False,
        include_references=True,
        include_builtin_scanners=True,
    )
    if links_text:
        return f"{banner}\n\n{links_text}"
    return banner


def _normalize_public_value(value: object) -> str:
    text = str(value).strip() if value is not None else ""
    if not text:
        return ""
    if text.upper() in {"TODO", "N/A", "NA", "NONE", "NULL", "UNKNOWN", "TBD", "TBA"}:
        return ""
    return text


def _render_resources_block(documents: dict, paper: dict, primary_reference: dict | None = None) -> str:
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
    if not public_links and not primary_reference:
        return ""

    reference_title = ""
    reference_fields: list[tuple[str, str]] = []
    if primary_reference:
        title = primary_reference.get("title", "")
        if title:
            reference_title = title
        if primary_reference.get("author"):
            reference_fields.append(("Author", primary_reference["author"]))
        if primary_reference.get("arxiv"):
            reference_fields.append(("arXiv", primary_reference["arxiv"]))
        if primary_reference.get("doi"):
            reference_fields.append(("DOI", primary_reference["doi"]))

    labels = [f"{label}:" for label, _value in public_links]
    if reference_title:
        labels.append("Reference:")
    label_width = max(len(label) for label in labels)
    lines = ["Resources:"]
    lines.extend(f"  {f'{label}:':<{label_width}}  {value}" for label, value in public_links)
    if reference_title:
        lines.append(f"  {'Reference:':<{label_width}}  {BOLD}{reference_title}{RESET}")
        metadata_indent = "  " + (" " * label_width) + "  "
        for label, value in reference_fields:
            lines.append(f"{metadata_indent}{label}: {value}")
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


def _render_references_block(payload: dict, paper: dict, *, include_builtin_scanners: bool = False) -> str:
    references = payload.get("references", {})
    if not isinstance(references, dict):
        references = {}

    section_specs = [
        ("Jarvis-HEP", references.get("jarvis_hep")),
    ]
    if include_builtin_scanners:
        section_specs.append(("Built-in Scanners", references.get("builtin_scanners")))
    rendered_sections = []
    for section_name, entries in section_specs:
        if not isinstance(entries, list):
            entries = []
        normalized = [_normalize_reference_entry(entry) for entry in entries]
        normalized = [entry for entry in normalized if entry.get("title")]
        if not normalized:
            continue

        block_lines = [f"  {section_name}:"]
        for idx, entry in enumerate(normalized, start=1):
            title_line = f"[{idx}] "
            if entry["caption"]:
                title_line += entry["caption"]
            else:
                title_line += entry["title"]
            block_lines.append(title_line)
            block_lines.append(f"  {'Title:':<7}  {entry['title']}")
            if entry["author"]:
                block_lines.append(f"  {'Author:':<7}  {entry['author']}")
            if entry["arxiv"]:
                block_lines.append(f"  {'arXiv:':<7}  {entry['arxiv']}")
            if entry["doi"]:
                block_lines.append(f"  {'DOI:':<7}  {entry['doi']}")
        rendered_sections.append("\n".join(block_lines))

    # Backward-compatibility with legacy "paper" schema.
    if not rendered_sections:
        legacy_title = _normalize_public_value(paper.get("title"))
        legacy_url = _normalize_public_value(paper.get("url"))
        if legacy_title or legacy_url:
            legacy_block = ["  Jarvis-HEP:"]
            legacy_block.append(f"[1] {legacy_title or 'Jarvis-HEP Paper'}")
            if legacy_url:
                legacy_block.append(f"  {'Link:':<7}  {legacy_url}")
            rendered_sections.append("\n".join(legacy_block))

    if not rendered_sections:
        return ""

    return "References:\n" + "\n".join(rendered_sections)


def _primary_jarvis_reference(payload: dict) -> dict:
    references = payload.get("references", {})
    if not isinstance(references, dict):
        return {}
    entries = references.get("jarvis_hep")
    if not isinstance(entries, list):
        return {}
    for entry in entries:
        normalized = _normalize_reference_entry(entry)
        if normalized.get("title"):
            return normalized
    return {}


def _render_document_links(
    *,
    include_resources: bool = True,
    include_compact_reference: bool = False,
    include_references: bool = True,
    include_builtin_scanners: bool = False,
) -> str:
    links_path = os.path.join(os.path.dirname(__file__), "card", "document_links.json")
    if not os.path.exists(links_path):
        return ""

    try:
        with open(links_path, "r", encoding="utf-8") as f1:
            payload = json.load(f1)
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

    sections = []
    if include_resources:
        primary_reference = _primary_jarvis_reference(payload) if include_compact_reference else None
        sections.append(_render_resources_block(documents, paper, primary_reference))
    if include_references:
        sections.append(_render_references_block(payload, paper, include_builtin_scanners=include_builtin_scanners))
    sections = [section for section in sections if section]
    return "\n\n".join(sections)


def _print_top_level_help() -> None:
    print(_TOP_LEVEL_HELP_TEXT, end="")


def _print_project_help() -> None:
    print(_PROJECT_HELP_TEXT, end="")


def _print_project_pack_help() -> None:
    print(_PROJECT_PACK_HELP_TEXT, end="")


def _print_project_subcommand_help(command: str) -> None:
    print(_PROJECT_SUBCOMMAND_HELP[command], end="")


def _print_portal_help() -> None:
    print(_PORTAL_HELP_TEXT, end="")


def _run_portal_formats() -> int:
    try:
        from jarvis_portal import available_formats
    except ImportError as exc:
        print(f"[Jarvis-HEP] Jarvis-Portal is not available: {exc}")
        return 1

    try:
        all_formats = available_formats()
        input_formats = available_formats("input")
        output_formats = available_formats("output")
    except Exception as exc:
        print(f"[Jarvis-HEP] Failed to query Jarvis-Portal IO registry: {exc}")
        return 1

    print(
        _render_portal_formats_table(
            [
                ("All", all_formats),
                ("Input", input_formats),
                ("Output", output_formats),
            ]
        )
    )
    return 0


def _format_names(values) -> str:
    names = [str(value) for value in values or []]
    return ", ".join(names) if names else "none"


def _render_portal_formats_table(rows) -> str:
    normalized = [(str(scope), list(formats or [])) for scope, formats in rows]
    format_names = []
    seen = set()
    for _, formats in normalized:
        for format_name in formats:
            text = str(format_name)
            key = text.casefold()
            if key in seen:
                continue
            seen.add(key)
            format_names.append(text)

    scopes = [scope for scope, _ in normalized]
    headers = ["Format", *scopes]
    table_rows = []
    supported_by_scope = {
        scope: {str(format_name).casefold() for format_name in formats}
        for scope, formats in normalized
    }
    for format_name in format_names:
        key = format_name.casefold()
        table_rows.append(
            [format_name, *["✓" if key in supported_by_scope[scope] else "" for scope in scopes]]
        )

    widths = [
        max(len(str(row[idx])) for row in [headers, *table_rows])
        for idx in range(len(headers))
    ]
    spacing = "  "
    rule = spacing.join("-" * width for width in widths)
    header_line = spacing.join(
        f"{header:<{width}}" for header, width in zip(headers, widths)
    )
    lines = [
        "Jarvis-Portal calculator IO formats",
        rule,
        header_line,
        rule,
    ]
    for row in table_rows:
        lines.append(spacing.join(f"{cell:<{width}}" for cell, width in zip(row, widths)))
    lines.append(rule)
    return "\n".join(lines)


def _human_bytes(value: int) -> str:
    units = ["B", "KB", "MB", "GB", "TB"]
    amount = float(max(0, int(value)))
    idx = 0
    while amount >= 1024.0 and idx < len(units) - 1:
        amount /= 1024.0
        idx += 1
    if idx == 0:
        return f"{int(amount)} {units[idx]}"
    return f"{amount:.2f} {units[idx]}"


def _run_project_create(project_name: str) -> int:
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


def _run_project_pack(project_path: str | None, profile: str) -> int:
    try:
        report = create_project_package(
            project_root=project_path,
            profile=profile,
        )
    except ProjectPackError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except Exception as exc:
        print(f"[Jarvis-HEP] Failed to package project: {exc}")
        return 1

    print(
        format_two_column_log(
            "Jarvis-HEP project package created",
            [
                ("Project root", report.project_root),
                ("Archive", report.archive_path),
                ("Mode", report.profile),
                ("Packed files", report.included_files),
                ("Skipped files", report.excluded_files),
                ("Payload size", _human_bytes(report.total_bytes)),
            ],
        )
    )
    return 0


def _run_project_pack_manifest(project_path: str | None, profile: str) -> int:
    try:
        report = create_project_pack_manifest(
            project_root=project_path,
            profile=profile,
        )
    except ProjectPackError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except Exception as exc:
        print(f"[Jarvis-HEP] Failed to write project pack manifest: {exc}")
        return 1

    print(
        format_two_column_log(
            "Jarvis-HEP project pack manifest created",
            [
                ("Project root", report.project_root),
                ("Manifest", report.manifest_path),
                ("Pack ID", report.pack_id),
                ("Mode", report.profile),
                ("Manifest files", report.included_files),
                ("Manifest excludes", report.excluded_files),
                ("Output", report.output),
            ],
        )
    )
    return 0


def _run_project_pack_from_manifest(manifest_path: str) -> int:
    try:
        report = create_project_package_from_manifest(manifest_path)
    except ProjectPackError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except Exception as exc:
        print(f"[Jarvis-HEP] Failed to package project from manifest: {exc}")
        return 1

    print(
        format_two_column_log(
            "Jarvis-HEP project package created from manifest",
            [
                ("Project root", report.project_root),
                ("Archive", report.archive_path),
                ("Mode", report.profile),
                ("Packed files", report.included_files),
                ("Skipped files", report.excluded_files),
                ("Payload size", _human_bytes(report.total_bytes)),
            ],
        )
    )
    return 0


def _looks_like_yaml_path(path: str | None) -> bool:
    if path is None:
        return False
    return path.lower().endswith((".yaml", ".yml"))


def _run_project_browse() -> int:
    try:
        projects = list_official_projects()
    except OfficialLibraryError as exc:
        print(f"[Jarvis-HEP] Failed to query the official Jarvis library: {exc}")
        return 1

    if not projects:
        print("[Jarvis-HEP] No verified projects are currently listed in the official Jarvis library.")
        return 0

    print("[Jarvis-HEP] Verified projects in the official Jarvis library:")
    for project in projects:
        name = project.get("name") or "<unnamed>"
        category = project.get("category") or "general"
        summary = project.get("summary") or "No summary available."
        print(f"[Jarvis-HEP] - {name} | {category} | {summary}")
    return 0


def _run_project_info(project_name: str) -> int:
    try:
        project = get_official_project(project_name)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except OfficialLibraryError as exc:
        print(f"[Jarvis-HEP] Failed to query the official Jarvis library: {exc}")
        return 1

    print(f"[Jarvis-HEP] Official project: {project['name']}")
    print(f"[Jarvis-HEP] Summary: {project.get('summary') or 'N/A'}")
    print(f"[Jarvis-HEP] Category: {project.get('category') or 'N/A'}")
    print(f"[Jarvis-HEP] Entrypoint: {project.get('entrypoint') or 'N/A'}")
    notes = project.get("compatibility_notes") or "None"
    print(f"[Jarvis-HEP] Compatibility notes: {notes}")
    return 0


def _run_project_fetch(project_name: str) -> int:
    try:
        report = fetch_official_project(project_name)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 2
    except OfficialProjectFetchError as exc:
        print(f"[Jarvis-HEP] {exc}")
        return 1
    except OfficialLibraryError as exc:
        print(f"[Jarvis-HEP] Failed to query the official Jarvis library: {exc}")
        return 1

    print(
        f"[Jarvis-HEP] Official project '{report.project_name}' fetched to: {report.target_dir}"
    )
    print(f"[Jarvis-HEP] Entrypoint: {report.entrypoint or 'N/A'}")
    return 0


def _parse_pack_arguments(tokens: list[str]) -> tuple[str | None, str, bool, bool] | int:
    path: str | None = None
    mode_flag: str | None = None
    manifest_only = False

    for tok in tokens:
        if tok in _HELP_FLAGS:
            _print_project_pack_help()
            return 0
        if tok == _PACK_MANIFEST_FLAG:
            manifest_only = True
            continue
        if tok in _PACK_MODE_FLAGS:
            if mode_flag is not None:
                print("[Jarvis-HEP] project pack modes are mutually exclusive: --share, --repro, --full")
                return 2
            mode_flag = tok
            continue
        if tok.startswith("-"):
            print(f"[Jarvis-HEP] Unsupported option for project pack: {tok}")
            return 2
        if path is None:
            path = tok
            continue
        print(f"[Jarvis-HEP] Unexpected argument for project pack: {tok}")
        return 2

    if manifest_only and _looks_like_yaml_path(path):
        print("[Jarvis-HEP] --man writes a new manifest from a project path, not from a manifest file")
        return 2

    is_manifest_input = _looks_like_yaml_path(path)
    if is_manifest_input and mode_flag is not None:
        print("[Jarvis-HEP] Manifest packing does not accept --share, --repro, or --full")
        return 2

    profile = _PACK_MODE_FLAGS.get(mode_flag, "share")
    return path, profile, manifest_only, is_manifest_input


def _handle_project_subcommand(command: str, args: list[str]) -> int:
    if command == "create":
        if len(args) == 1 and args[0] in _HELP_FLAGS:
            _print_project_subcommand_help("create")
            return 0
        if len(args) != 1 or args[0].startswith("-"):
            print("[Jarvis-HEP] Usage error: Jarvis project create <name>")
            return 2
        return _run_project_create(args[0])

    if command == "pack":
        parsed = _parse_pack_arguments(args)
        if isinstance(parsed, int):
            return parsed
        path, profile, manifest_only, is_manifest_input = parsed
        if manifest_only:
            return _run_project_pack_manifest(path, profile)
        if is_manifest_input:
            return _run_project_pack_from_manifest(str(path))
        return _run_project_pack(path, profile)

    if command == "browse":
        if len(args) == 1 and args[0] in _HELP_FLAGS:
            _print_project_subcommand_help("browse")
            return 0
        if args:
            print("[Jarvis-HEP] Usage error: Jarvis project browse")
            return 2
        return _run_project_browse()

    if command == "fetch":
        if len(args) == 1 and args[0] in _HELP_FLAGS:
            _print_project_subcommand_help("fetch")
            return 0
        if len(args) != 1 or args[0].startswith("-"):
            print("[Jarvis-HEP] Usage error: Jarvis project fetch <name>")
            return 2
        return _run_project_fetch(args[0])

    if command == "info":
        if len(args) == 1 and args[0] in _HELP_FLAGS:
            _print_project_subcommand_help("info")
            return 0
        if len(args) != 1 or args[0].startswith("-"):
            print("[Jarvis-HEP] Usage error: Jarvis project info <name>")
            return 2
        return _run_project_info(args[0])

    print(f"[Jarvis-HEP] Unknown project command: {command}")
    _print_project_help()
    return 2


def _project_fast_path(argv: list[str]) -> int | None:
    if len(argv) < 2 or argv[1] != "project":
        return None

    tokens = argv[2:]
    if not tokens:
        _print_project_help()
        return 0
    if tokens[0] in _HELP_FLAGS:
        if len(tokens) > 1:
            print("[Jarvis-HEP] Usage error: Jarvis project --help")
            return 2
        _print_project_help()
        return 0

    return _handle_project_subcommand(tokens[0], tokens[1:])


def _portal_fast_path(argv: list[str]) -> int | None:
    if len(argv) < 2 or argv[1] != "portal":
        return None

    tokens = argv[2:]
    if not tokens:
        _print_portal_help()
        return 0
    if tokens[0] in _HELP_FLAGS:
        if len(tokens) > 1:
            print("[Jarvis-HEP] Usage error: Jarvis portal --help")
            return 2
        _print_portal_help()
        return 0

    command = tokens[0]
    args = tokens[1:]
    if command == "formats":
        if len(args) == 1 and args[0] in _HELP_FLAGS:
            print("Usage:\n  Jarvis portal formats\n")
            return 0
        if args:
            print("[Jarvis-HEP] Usage error: Jarvis portal formats")
            return 2
        return _run_portal_formats()

    print(f"[Jarvis-HEP] Unknown portal command: {command}")
    _print_portal_help()
    return 2


def _top_level_help_fast_path(argv: list[str]) -> int | None:
    if len(argv) == 1:
        _print_top_level_help()
        return 0
    if len(argv) >= 2 and argv[1] in {"project", "portal"}:
        return None
    if any(token in _HELP_FLAGS for token in argv[1:]):
        _print_top_level_help()
        return 0
    return None


def _version_fast_path(argv: list[str]) -> int | None:
    if "--version" not in argv and "-v" not in argv:
        return None

    print(_render_version_banner())
    return 0


def _refs_fast_path(argv: list[str]) -> int | None:
    if "--refs" not in argv:
        return None

    print(_render_refs_text())
    return 0


def main(argv=None) -> int:
    argv = sys.argv if argv is None else list(argv)

    project_code = _project_fast_path(argv)
    if project_code is not None:
        return project_code

    portal_code = _portal_fast_path(argv)
    if portal_code is not None:
        return portal_code

    help_code = _top_level_help_fast_path(argv)
    if help_code is not None:
        return help_code

    refs_code = _refs_fast_path(argv)
    if refs_code is not None:
        return refs_code

    version_code = _version_fast_path(argv)
    if version_code is not None:
        return version_code

    from jarvishep.core import Core

    jc = Core()
    old_argv = sys.argv
    try:
        sys.argv = list(argv)
        jc.initialization()
    finally:
        sys.argv = old_argv

    if getattr(jc.args, "version", False):
        print(_render_version_banner())
    elif jc.scan_mode:
        jc.run_sampling()
    elif jc.args.cvtDB:
        jc.convert()
    elif jc.args.plot:
        jc.plot()
    elif jc.args.monitor:
        jc.monitor()
    else:
        _print_top_level_help()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
