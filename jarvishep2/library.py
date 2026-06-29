#!/usr/bin/env python3
"""LibDeps helpers for calculator isolation (WP-D2.3)."""

from __future__ import annotations

import os
from collections.abc import Mapping
from typing import Any

from jarvishep2.command_parser import CommandParser, ResolvedExecutable


class LibraryManager:
    """Minimal LibDeps helper: symlink safe tools into Sample dirs."""

    def __init__(self, config: Mapping[str, Any] | None = None, *, task_root: str | None = None) -> None:
        self.config = dict(config or {})
        self.task_root = str(task_root or os.getcwd())

    @staticmethod
    def requires_shadow(clone_shadow: bool) -> bool:
        return bool(clone_shadow)

    def link_into_sample(self, source_path: str, sample_dir: str, link_name: str) -> str:
        """Symlink a concurrency-safe tool into a Sample directory (zero-copy)."""
        source = os.path.abspath(str(source_path))
        if not os.path.exists(source):
            raise FileNotFoundError(f"LibDeps source does not exist: {source}")
        sample_root = os.path.abspath(str(sample_dir))
        os.makedirs(sample_root, exist_ok=True)
        link_path = os.path.join(sample_root, str(link_name))
        if os.path.lexists(link_path):
            return link_path
        os.symlink(source, link_path)
        return link_path

    def resolve_registered(self, parser: CommandParser) -> dict[str, ResolvedExecutable]:
        """Return the Phase-1 registered executable map from a CommandParser."""
        return dict(parser.registered)


__all__ = ["LibraryManager"]