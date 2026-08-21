#!/usr/bin/env python3
"""Jarvis project tools: create, pack, browse, fetch, info, encrypt (D25.9)."""

from __future__ import annotations

import os
import sys

from rich import box
from rich.console import Console, Group
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.cli.help import build_project_help_parser, _render_help_for_args
from jarvishep2.run_outcome import EXIT_OK, EXIT_RUN_FAILED, EXIT_USAGE

_PROJECT_COMMANDS = frozenset(
    {"create", "pack", "browse", "fetch", "info", "encrypt"}
)
_PACK_MODE_FLAGS = {
    "--share": "share",
    "--repro": "repro",
    "--full": "full",
}
_PACK_MANIFEST_FLAG = "--man"
_PACK_ENCRYPT_FLAGS = frozenset({"--encrypt", "--encrypted"})
_HELP_FLAGS = frozenset({"-h", "--help"})

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


def _print_kv_block(title: str, rows: list[tuple[str, object]]) -> None:
    print(title)
    width = max((len(str(key)) for key, _ in rows), default=8)
    for key, value in rows:
        print(f"  {str(key):<{width}}  {value}")


def _print_project_help() -> None:
    build_project_help_parser().print_help()


def _run_project_create(project_name: str) -> int:
    from jarvishep2.project_scaffold import PROJECT_SUBDIRS, create_project_scaffold

    try:
        project_root = create_project_scaffold(project_name, cwd=os.getcwd())
    except ValueError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except FileExistsError as exc:
        print(f"[Jarvis] Project directory already exists: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    print(f"[Jarvis] Project scaffold created at: {project_root}")
    print(f"[Jarvis] Created folders: {', '.join(PROJECT_SUBDIRS)}")
    print("[Jarvis] Try: Jarvis run bin/quickstart_bridson_operas.yaml")
    return EXIT_OK


def _run_project_pack(
    project_path: str | None,
    profile: str,
    *,
    encrypt: bool = False,
    key: str | None = None,
) -> int:
    from jarvishep2.project_crypto import ProjectCryptoError, encrypt_file, resolve_fetch_key
    from jarvishep2.project_packager import ProjectPackError, create_project_package

    if encrypt and not resolve_fetch_key(key):
        print(
            "[Jarvis] --encrypt requires --key KEY or JARVIS_PROJECT_FETCH_KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE

    try:
        report = create_project_package(project_root=project_path, profile=profile)
    except ProjectPackError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis] Failed to package project: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    archive_path = report.archive_path
    encrypted_path = None
    if encrypt:
        try:
            encrypted_path = archive_path + ".jenc"
            if archive_path.endswith(".tar.gz"):
                encrypted_path = archive_path[: -len(".tar.gz")] + ".tar.gz.jenc"
            encrypt_file(archive_path, encrypted_path, key=str(resolve_fetch_key(key)))
        except ProjectCryptoError as exc:
            print(f"[Jarvis] {exc}", file=sys.stderr)
            return EXIT_RUN_FAILED

    rows = [
        ("Project root", report.project_root),
        ("Archive", report.archive_path),
        ("Mode", report.profile),
        ("Packed files", report.included_files),
        ("Skipped files", report.excluded_files),
        ("Payload size", _human_bytes(report.total_bytes)),
    ]
    if encrypted_path:
        rows.append(("Encrypted", encrypted_path))
        rows.append(("Fetch with", f"Jarvis project fetch <name> --key <key>"))
    _print_kv_block("Jarvis project package created", rows)
    return EXIT_OK


def _run_project_encrypt(archive_path: str, *, key: str | None) -> int:
    from jarvishep2.project_crypto import ProjectCryptoError, encrypt_file, resolve_fetch_key

    resolved = resolve_fetch_key(key)
    if not resolved:
        print(
            "[Jarvis] encrypt requires --key KEY or JARVIS_PROJECT_FETCH_KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE
    path = os.path.abspath(os.path.expanduser(archive_path))
    if not os.path.isfile(path):
        print(f"[Jarvis] Archive not found: {path}", file=sys.stderr)
        return EXIT_USAGE
    if path.endswith(".jenc"):
        print("[Jarvis] File already looks encrypted (.jenc)", file=sys.stderr)
        return EXIT_USAGE
    out = path + ".jenc"
    if path.endswith(".tar.gz"):
        out = path[: -len(".tar.gz")] + ".tar.gz.jenc"
    try:
        encrypt_file(path, out, key=resolved)
    except ProjectCryptoError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    _print_kv_block(
        "Jarvis project archive encrypted",
        [
            ("Input", path),
            ("Encrypted", out),
            ("Scheme", "openssl-aes-256-cbc"),
            ("Users fetch with", "Jarvis project fetch NAME --key …"),
        ],
    )
    return EXIT_OK


def _run_project_pack_manifest(project_path: str | None, profile: str) -> int:
    from jarvishep2.project_packager import ProjectPackError, create_project_pack_manifest

    try:
        report = create_project_pack_manifest(project_root=project_path, profile=profile)
    except ProjectPackError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis] Failed to write project pack manifest: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    _print_kv_block(
        "Jarvis project pack manifest created",
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
    return EXIT_OK


def _run_project_pack_from_manifest(manifest_path: str) -> int:
    from jarvishep2.project_packager import (
        ProjectPackError,
        create_project_package_from_manifest,
    )

    try:
        report = create_project_package_from_manifest(manifest_path)
    except ProjectPackError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except Exception as exc:
        print(f"[Jarvis] Failed to package project from manifest: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    _print_kv_block(
        "Jarvis project package created from manifest",
        [
            ("Project root", report.project_root),
            ("Archive", report.archive_path),
            ("Mode", report.profile),
            ("Packed files", report.included_files),
            ("Skipped files", report.excluded_files),
            ("Payload size", _human_bytes(report.total_bytes)),
        ],
    )
    return EXIT_OK


def _run_project_browse() -> int:
    from jarvishep2.official_project_library import (
        PROJECT_FETCH_KEY_ENV,
        OfficialLibraryError,
        list_official_projects,
    )

    try:
        projects = list_official_projects()
    except OfficialLibraryError as exc:
        print(f"[Jarvis] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    if not projects:
        Console().print("[dim]No verified projects are listed in the official library.[/]")
        return EXIT_OK

    table = Table(box=box.SIMPLE_HEAVY, show_header=True, header_style="bold dim")
    table.add_column("NAME", style="bold #c8c8ff", no_wrap=True)
    table.add_column("ACCESS", no_wrap=True)
    table.add_column("KEY", no_wrap=True)
    table.add_column("CATEGORY", no_wrap=True)
    table.add_column("SUMMARY", ratio=1, overflow="fold")

    for project in projects:
        access = str(project.get("access") or "public")
        needs_key = bool(project.get("requires_key"))
        access_style = "bold yellow" if access == "restricted" else "bold green"
        key_style = "bold yellow" if needs_key else "bold green"
        table.add_row(
            Text(str(project.get("name") or ""), style="bold #c8c8ff"),
            Text(access, style=access_style),
            Text("required" if needs_key else "no", style=key_style),
            Text(str(project.get("category") or "-")),
            Text(str(project.get("summary") or "-"), style="dim"),
        )

    key_note = Text()
    key_note.append("Key: ", style="bold")
    key_note.append("required", style="bold yellow")
    key_note.append(" → fetch needs --key or ", style="dim")
    key_note.append(PROJECT_FETCH_KEY_ENV, style="bold cyan")
    key_note.append("; ", style="dim")
    key_note.append("no", style="bold green")
    key_note.append(" → open download.", style="dim")

    public_example = next(
        (project for project in projects if not bool(project.get("requires_key"))),
        None,
    )
    encrypted_example = next(
        (project for project in projects if bool(project.get("requires_key"))),
        None,
    )
    fetch_examples = Text()
    fetch_examples.append("Fetch examples:\n", style="bold")
    if public_example is not None:
        fetch_examples.append("  public:     ", style="dim")
        fetch_examples.append(
            f"Jarvis project fetch {public_example.get('name')}",
            style="bold cyan",
        )
    if encrypted_example is not None:
        if public_example is not None:
            fetch_examples.append("\n")
        fetch_examples.append("  encrypted:  ", style="dim")
        fetch_examples.append(
            f"Jarvis project fetch {encrypted_example.get('name')} --key 'YOUR_KEY'",
            style="bold cyan",
        )

    Console().print(
        Panel(
            Group(table, Text(""), key_note, Text(""), fetch_examples),
            title=(
                "[bold]Official library[/] [dim]·[/] "
                "[bold green]Jarvis-Examples GitHub JSON[/]"
            ),
            box=box.ROUNDED,
            border_style="dim",
            padding=(0, 1),
        )
    )
    return EXIT_OK


def _run_project_info(project_name: str) -> int:
    from jarvishep2.official_project_library import (
        OfficialLibraryError,
        OfficialProjectNotFoundError,
        get_official_project,
    )

    try:
        project = get_official_project(project_name)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except OfficialLibraryError as exc:
        print(f"[Jarvis] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    needs_key = bool(project.get("requires_key"))
    access = str(project.get("access") or "public")
    access_style = "bold yellow" if access == "restricted" else "bold green"
    key_style = "bold yellow" if needs_key else "bold green"

    details = Table.grid(padding=(0, 1))
    details.add_column(style="bold dim", no_wrap=True)
    details.add_column(ratio=1, overflow="fold")
    details.add_row("SUMMARY", Text(str(project.get("summary") or "N/A"), style="dim"))
    details.add_row("CATEGORY", Text(str(project.get("category") or "N/A")))
    details.add_row("ACCESS", Text(access, style=access_style))
    details.add_row(
        "KEY REQUIRED",
        Text("yes" if needs_key else "no", style=key_style),
    )
    if needs_key and project.get("encryption_hint"):
        details.add_row(
            "KEY HINT",
            Text(str(project.get("encryption_hint")), style="dim"),
        )
    if needs_key:
        details.add_row(
            "ENCRYPTION",
            Text(
                str(project.get("encryption_scheme") or "openssl-aes-256-cbc"),
                style="dim",
            ),
        )
    details.add_row("ENTRYPOINT", Text(str(project.get("entrypoint") or "N/A")))
    details.add_row(
        "COMPATIBILITY",
        Text(str(project.get("compatibility_notes") or "None"), style="dim"),
    )

    title = Text("Official project", style="bold")
    title.append(" · ", style="dim")
    title.append(str(project["name"]), style="bold #c8c8ff")
    Console().print(
        Panel(
            details,
            title=title,
            box=box.ROUNDED,
            border_style="dim",
            padding=(0, 1),
        )
    )
    return EXIT_OK


def _run_project_fetch(project_name: str, *, key: str | None = None) -> int:
    from jarvishep2.official_project_library import (
        OfficialLibraryError,
        OfficialProjectFetchError,
        OfficialProjectNotFoundError,
        fetch_official_project,
    )

    try:
        report = fetch_official_project(project_name, key=key)
    except OfficialProjectNotFoundError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_USAGE
    except OfficialProjectFetchError as exc:
        print(f"[Jarvis] {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED
    except OfficialLibraryError as exc:
        print(f"[Jarvis] Failed to query the official library: {exc}", file=sys.stderr)
        return EXIT_RUN_FAILED

    access = str(report.access or "public")
    access_style = "bold yellow" if access == "restricted" else "bold green"
    details = Table.grid(padding=(0, 1))
    details.add_column(style="bold dim", no_wrap=True)
    details.add_column(ratio=1, overflow="fold")
    details.add_row("SAVED TO", Text(report.target_dir, style="bold cyan"))
    details.add_row("ACCESS", Text(access, style=access_style))
    details.add_row(
        "ENTRYPOINT",
        Text(report.entrypoint or "N/A"),
    )
    if report.required_key:
        details.add_row("KEY", Text("used", style="bold yellow"))

    title = Text("Fetched", style="bold")
    title.append(" · ", style="dim")
    title.append(str(report.project_name), style="bold #c8c8ff")
    Console().print(
        Panel(
            details,
            title=title,
            box=box.ROUNDED,
            border_style="green",
            padding=(0, 1),
        )
    )
    return EXIT_OK


def _parse_fetch_arguments(tokens: list[str]) -> tuple[str, str | None] | int:
    """Parse ``fetch <name> [--key KEY]``."""
    if not tokens or (len(tokens) == 1 and tokens[0] in _HELP_FLAGS):
        print(
            "usage: Jarvis project fetch <name> [--key KEY]\n"
            "  restricted projects need --key or JARVIS_PROJECT_FETCH_KEY\n"
        )
        return 0 if tokens and tokens[0] in _HELP_FLAGS else EXIT_USAGE

    name: str | None = None
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis] Unsupported option for project fetch: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if name is None:
            name = tok
            i += 1
            continue
        print(f"[Jarvis] Unexpected argument for project fetch: {tok}", file=sys.stderr)
        return EXIT_USAGE

    if not name:
        print("[Jarvis] Usage: Jarvis project fetch <name> [--key KEY]", file=sys.stderr)
        return EXIT_USAGE
    return name, key


def _looks_like_yaml_path(path: str | None) -> bool:
    if path is None:
        return False
    return path.lower().endswith((".yaml", ".yml"))


def _parse_pack_arguments(
    tokens: list[str],
) -> tuple[str | None, str, bool, bool, bool, str | None] | int:
    path: str | None = None
    mode_flag: str | None = None
    manifest_only = False
    encrypt = False
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in _HELP_FLAGS:
            print(
                "usage: Jarvis project pack [path] [--share|--repro|--full] [--man]\n"
                "       Jarvis project pack [path] --encrypt --key KEY\n"
                "       Jarvis project pack <pack_manifest.yaml>\n"
            )
            return 0
        if tok == _PACK_MANIFEST_FLAG:
            manifest_only = True
            i += 1
            continue
        if tok in _PACK_ENCRYPT_FLAGS:
            encrypt = True
            i += 1
            continue
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok in _PACK_MODE_FLAGS:
            if mode_flag is not None:
                print(
                    "[Jarvis] project pack modes are mutually exclusive: "
                    "--share, --repro, --full",
                    file=sys.stderr,
                )
                return EXIT_USAGE
            mode_flag = tok
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis] Unsupported option for project pack: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if path is None:
            path = tok
            i += 1
            continue
        print(f"[Jarvis] Unexpected argument for project pack: {tok}", file=sys.stderr)
        return EXIT_USAGE

    if manifest_only and _looks_like_yaml_path(path):
        print(
            "[Jarvis] --man writes a new manifest from a project path, not from a manifest file",
            file=sys.stderr,
        )
        return EXIT_USAGE

    is_manifest_input = _looks_like_yaml_path(path)
    if is_manifest_input and mode_flag is not None:
        print(
            "[Jarvis] Manifest packing does not accept --share, --repro, or --full",
            file=sys.stderr,
        )
        return EXIT_USAGE
    if is_manifest_input and encrypt:
        print(
            "[Jarvis] --encrypt is for packing a project path; "
            "encrypt a finished archive with: Jarvis project encrypt FILE --key K",
            file=sys.stderr,
        )
        return EXIT_USAGE

    profile = _PACK_MODE_FLAGS.get(mode_flag or "", "share")
    return path, profile, manifest_only, is_manifest_input, encrypt, key


def _parse_encrypt_arguments(tokens: list[str]) -> tuple[str, str | None] | int:
    if not tokens or (len(tokens) == 1 and tokens[0] in _HELP_FLAGS):
        print("usage: Jarvis project encrypt <archive.tar.gz> --key KEY\n")
        return 0 if tokens and tokens[0] in _HELP_FLAGS else EXIT_USAGE
    path: str | None = None
    key: str | None = None
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if tok in {"--key", "-k"}:
            if i + 1 >= len(tokens):
                print("[Jarvis] --key requires a value", file=sys.stderr)
                return EXIT_USAGE
            key = tokens[i + 1]
            i += 2
            continue
        if tok.startswith("--key="):
            key = tok.split("=", 1)[1]
            i += 1
            continue
        if tok.startswith("-"):
            print(f"[Jarvis] Unsupported option for project encrypt: {tok}", file=sys.stderr)
            return EXIT_USAGE
        if path is None:
            path = tok
            i += 1
            continue
        print(f"[Jarvis] Unexpected argument: {tok}", file=sys.stderr)
        return EXIT_USAGE
    if not path:
        print(
            "[Jarvis] Usage: Jarvis project encrypt <archive.tar.gz> --key KEY",
            file=sys.stderr,
        )
        return EXIT_USAGE
    return path, key


def dispatch_project(project_argv: list[str] | None = None) -> int:
    """Handle ``Jarvis project create|pack|browse|fetch|info`` (D12.5)."""
    args = list(project_argv or [])
    if not args or args[0] in _HELP_FLAGS:
        _print_project_help()
        return EXIT_OK
    if any(arg in _HELP_FLAGS for arg in args):
        return _render_help_for_args(build_project_help_parser(), args)

    command = args[0]
    rest = args[1:]
    if command not in _PROJECT_COMMANDS:
        print(f"[Jarvis] Unknown project command: {command}", file=sys.stderr)
        _print_project_help()
        return EXIT_USAGE

    if command == "create":
        if len(rest) == 1 and rest[0] in _HELP_FLAGS:
            print("usage: Jarvis project create <name>\n")
            return EXIT_OK
        if len(rest) != 1 or rest[0].startswith("-"):
            print("[Jarvis] Usage: Jarvis project create <name>", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_create(rest[0])

    if command == "pack":
        parsed = _parse_pack_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        path, profile, manifest_only, is_manifest_input, encrypt, key = parsed
        if manifest_only:
            return _run_project_pack_manifest(path, profile)
        if is_manifest_input:
            return _run_project_pack_from_manifest(str(path))
        return _run_project_pack(path, profile, encrypt=encrypt, key=key)

    if command == "encrypt":
        parsed = _parse_encrypt_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        archive, key = parsed
        return _run_project_encrypt(archive, key=key)

    if command == "browse":
        if rest:
            if len(rest) == 1 and rest[0] in _HELP_FLAGS:
                print("usage: Jarvis project browse\n")
                return EXIT_OK
            print("[Jarvis] Usage: Jarvis project browse", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_browse()

    if command == "fetch":
        parsed = _parse_fetch_arguments(rest)
        if isinstance(parsed, int):
            return parsed
        name, key = parsed
        return _run_project_fetch(name, key=key)

    if command == "info":
        if len(rest) == 1 and rest[0] in _HELP_FLAGS:
            print("usage: Jarvis project info <name>\n")
            return EXIT_OK
        if len(rest) != 1 or rest[0].startswith("-"):
            print("[Jarvis] Usage: Jarvis project info <name>", file=sys.stderr)
            return EXIT_USAGE
        return _run_project_info(rest[0])

    return EXIT_USAGE

