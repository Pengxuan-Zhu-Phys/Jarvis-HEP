"""CLI results agree with the registry commit even when cleanup needs a retry."""

from __future__ import annotations

import os
import subprocess
import sys

import pytest
from jarvis_hep_sampler import SamplerRegistry

SOURCE = '''
PLUGIN = {"api_version": 1, "method": "ManagedDemo", "version": "1",
          "profile": "independent"}
def draw(context, rng, plan):
    return [0.25] * context.dimension
'''


@pytest.mark.parametrize("operation", ["load", "remove"])
def test_cleanup_warning_does_not_turn_a_committed_cli_operation_into_failure(
    tmp_path, monkeypatch, operation
):
    monkeypatch.setenv("JARVIS_SAMPLERS_HOME", str(tmp_path / "registry"))
    source = tmp_path / "ManagedDemo.py"
    source.write_text(SOURCE, encoding="utf-8")
    registry = SamplerRegistry()
    previous = registry.load(source).entry
    retired = registry.snapshot_path(previous)
    source.write_text(SOURCE.replace('"version": "1"', '"version": "2"'), encoding="utf-8")
    env = dict(
        os.environ,
        JARVIS_REDIS_NOTICE_STATE=str(tmp_path / "redis-notice"),
        PYTHONDONTWRITEBYTECODE="1",
    )
    # Inject only a filesystem cleanup failure in a new CLI interpreter. The
    # actual registry commit, logging setup, dispatcher and exit status run.
    bootstrap = '''
import sys
from pathlib import Path
from jarvis_hep_sampler import SamplerRegistry
retired = Path(sys.argv[1])
original = SamplerRegistry._delete_snapshot
def fail_cleanup(self, snapshot, entries):
    if self.home / snapshot == retired:
        raise PermissionError("injected cleanup failure")
    return original(self, snapshot, entries)
SamplerRegistry._delete_snapshot = fail_cleanup
from jarvishep2.cli.client import main
raise SystemExit(main(sys.argv[2:]))
'''
    args = ["load", str(source)] if operation == "load" else ["remove", "ManagedDemo", "--yes"]
    result = subprocess.run(
        [sys.executable, "-c", bootstrap, str(retired), "sampler", *args],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "already committed" in result.stderr
    assert "next successful sampler load/remove" in result.stderr
    current = registry.get("ManagedDemo")
    assert (current.version if current else None) == ("2" if operation == "load" else None)
    assert retired.exists() and source.exists()

    # A subsequent ordinary command completes the deferred cleanup.
    retry = subprocess.run(
        [sys.executable, "-c", "from jarvishep2.cli.client import main; "
         "raise SystemExit(main())", "sampler", "load", str(source)],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=30,
    )
    assert retry.returncode == 0, retry.stdout + retry.stderr
    assert not retired.exists()
    assert registry.load_snapshot("ManagedDemo").metadata.version == "2"
