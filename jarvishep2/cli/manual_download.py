"""User-facing PDF manual download command."""

from __future__ import annotations

import argparse
import json
import sys
import traceback

from rich.console import Console
from rich.panel import Panel
from rich.text import Text

from jarvishep2.cli.help import JarvisArgumentParser, add_manual_download_arguments
from jarvishep2.manual_download import ManualDownloadError, download_manual


class _Parser(JarvisArgumentParser):
    def error(self, message: str) -> None:
        raise ManualDownloadError(
            message, hint="See Jarvis download-manual --help.", exit_code=2,
        )


def dispatch_download_manual(argv: list[str]) -> int:
    as_json = any(arg == "--json" or arg.startswith("--json=") for arg in argv)
    try:
        parser = _Parser(
            prog="Jarvis download-manual",
            description="Download the PDF manual from a Jarvis-HEP GitHub Release.",
        )
        add_manual_download_arguments(parser)
        args = parser.parse_args(argv)
    except ManualDownloadError as exc:
        return _report_error(exc, as_json=as_json)
    return run_download_manual(args)


def run_download_manual(args: argparse.Namespace) -> int:
    """Execute a parsed download request, shared by both CLI dispatch paths."""
    as_json = args.json
    try:
        report = download_manual(
            output=args.output, release=args.release, force=args.force, timeout=args.timeout,
        )
    except KeyboardInterrupt:
        error = ManualDownloadError(
            "Manual download interrupted.", hint="Retry Jarvis download-manual.", exit_code=130,
        )
        return _report_error(error, as_json=as_json)
    except ManualDownloadError as exc:
        if args.debug:
            traceback.print_exc(file=sys.stderr)
        return _report_error(exc, as_json=as_json)
    if as_json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        print(f"Downloaded manual: {report['path']}")
        print(f"Release: {report['release']}\nSource: {report['source']}")
    return 0


def _report_error(error: ManualDownloadError, *, as_json: bool) -> int:
    if as_json:
        print(json.dumps({
            "ok": False, "command": "download-manual",
            "error": {"type": "usage_error" if error.exit_code == 2 else
                      "interrupted" if error.exit_code == 130 else "download_error",
                      "message": str(error), "hint": error.hint},
        }, indent=2, ensure_ascii=False))
    else:
        Console(stderr=True).print(Panel(
            Text(f"{error}\n\n{error.hint}"), title="Jarvis download-manual",
            border_style="red",
        ))
    return error.exit_code
