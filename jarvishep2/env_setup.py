#!/usr/bin/env python3
"""Capture and cache sourced shell environments for calculator subprocesses (WP-D3.2)."""

from __future__ import annotations

import hashlib
import os
import shlex
import subprocess
from collections.abc import Callable, Mapping, Sequence
from typing import Any

EnvRunner = Callable[..., subprocess.CompletedProcess[str]]


class EnvCapture:
    """Capture ``source``-based environments once per Worker process."""

    _cache: dict[str, dict[str, str]] = {}
    _runner: EnvRunner | None = None

    @classmethod
    def set_runner(cls, runner: EnvRunner | None) -> None:
        """Inject a subprocess runner (tests only)."""
        cls._runner = runner

    @classmethod
    def clear_cache(cls) -> None:
        cls._cache.clear()

    @classmethod
    def capture_from_source(
        cls,
        script: str,
        base_env: Mapping[str, str] | None = None,
    ) -> dict[str, str]:
        """Run ``source <script> && env`` and return the parsed environment."""
        script_path = os.path.abspath(str(script))
        if not os.path.isfile(script_path):
            raise FileNotFoundError(f"env_setup source script does not exist: {script_path}")

        cache_key = cls._cache_key(script_path, base_env)
        cached = cls._cache.get(cache_key)
        if cached is not None:
            return dict(cached)

        run_env = {str(k): str(v) for k, v in (base_env or os.environ).items()}
        command = f"source {shlex.quote(script_path)} && env"
        completed = cls._run_capture(["bash", "-c", command], env=run_env)
        if int(completed.returncode) != 0:
            stderr = (completed.stderr or "").strip()
            raise RuntimeError(
                f"env_setup failed for '{script_path}' (exit {completed.returncode}): {stderr}"
            )
        stdout = completed.stdout or ""
        if not stdout.strip():
            raise RuntimeError(f"env_setup produced empty environment for '{script_path}'")

        captured = cls._parse_env(stdout)
        cls._cache[cache_key] = dict(captured)
        return dict(captured)

    @classmethod
    def merged_env(
        cls,
        scripts: Sequence[str],
        base_env: Mapping[str, str] | None = None,
    ) -> dict[str, str]:
        """Merge captured environments from *scripts* in order over *base_env*."""
        merged = {str(k): str(v) for k, v in (base_env or os.environ).items()}
        for script in scripts:
            path = str(script).strip()
            if not path:
                continue
            captured = cls.capture_from_source(path, base_env=merged)
            merged.update(captured)
        return merged

    @staticmethod
    def _cache_key(script_path: str, base_env: Mapping[str, str] | None) -> str:
        if base_env is None:
            digest = "default"
        else:
            payload = repr(sorted((str(k), str(v)) for k, v in base_env.items())).encode()
            digest = hashlib.sha256(payload).hexdigest()[:16]
        return f"{os.path.abspath(script_path)}:{digest}"

    @classmethod
    def _run_capture(
        cls,
        command: list[str],
        *,
        env: Mapping[str, str],
    ) -> subprocess.CompletedProcess[str]:
        runner = cls._runner or subprocess.run
        return runner(
            command,
            capture_output=True,
            text=True,
            env=dict(env),
            check=False,
        )

    @staticmethod
    def _parse_env(text: str) -> dict[str, str]:
        parsed: dict[str, str] = {}
        for line in str(text).splitlines():
            if not line or line.startswith("#"):
                continue
            key, sep, value = line.partition("=")
            if not sep:
                continue
            parsed[key] = value
        return parsed


def resolve_env_setup_sources(
    env_setup: Sequence[Mapping[str, Any]] | None,
    *,
    command_parser: Any = None,
) -> list[str]:
    """Extract and optionally Phase-1-resolve ``env_setup`` source script paths."""
    scripts: list[str] = []
    for item in env_setup or []:
        if not isinstance(item, Mapping):
            continue
        source = str(item.get("source", "")).strip()
        if not source:
            continue
        if command_parser is not None and hasattr(command_parser, "resolve_static"):
            source = command_parser.resolve_static(source)
        scripts.append(source)
    return scripts


__all__ = ["EnvCapture", "resolve_env_setup_sources"]