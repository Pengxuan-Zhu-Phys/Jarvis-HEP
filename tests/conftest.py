"""Pytest hooks for the Jarvis-HEP V2 test suite."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Sequence


def drop_default_slow_markexpr(
    markexpr: str,
    invocation_args: Sequence[str],
    *,
    value_options: set[str] | None = None,
) -> str:
    """Keep the fast gate unless the user named a test file on the CLI.

    ``pyproject.toml`` sets ``addopts = -m not slow``. That filter still
    applies when a path is given, so ``pytest tests/test_mcmc_sampler.py``
    would collect nothing. Drop the *default* filter when the invocation
    names a ``.py`` file (or a node id) and the user did not pass ``-m``.
    ``pytest`` / ``pytest tests/`` stay on the fast gate.
    """
    expr = (markexpr or "").strip()
    if expr != "not slow":
        return expr
    named_paths: list[str] = []
    idx = 0
    args = [str(item) for item in invocation_args]
    value_options = _pytest_value_options() if value_options is None else value_options
    while idx < len(args):
        arg = args[idx]
        if arg in {"-m", "--markexpr"}:
            return expr
        if arg.startswith("-m") or arg.startswith("--markexpr="):
            return expr
        if arg == "--":
            named_paths.extend(args[idx + 1 :])
            break
        # Use pytest's registered options, including those added by plugins,
        # so option values (e.g. --deselect nodeid) never become test paths.
        if arg in value_options:
            idx += 2
            continue
        if not arg.startswith("-"):
            named_paths.append(arg)
        idx += 1
    for path in named_paths:
        node = path.split("::", 1)[0]
        if node.endswith(".py") or "::" in path:
            return ""
    return expr


def _pytest_value_options(parser=None) -> set[str]:
    from _pytest.config import get_config

    config = get_config() if parser is None else None
    try:
        if parser is None:
            parser = config._parser
        no_value = {"store_true", "store_false", "store_const", "count", "help", "version"}
        return {name for group in [*parser._groups, parser._anonymous]
                for option in group.options
                if option.attrs().get("nargs") != 0
                and option.attrs().get("action") not in no_value
                for name in option.names()}
    finally:
        if config is not None:
            config._ensure_unconfigure()


def _isolate_sampler_registry() -> None:
    """Never read the developer's own ~/.jarvis/samplers (Jarvis-HEP-Sampler registry).

    Spawned Workers inherit the variable. Tests that load plug-ins point
    ``JARVIS_SAMPLERS_HOME`` at their own temporary directory.
    """
    if os.environ.get("JARVIS_SAMPLERS_HOME"):
        return
    os.environ["JARVIS_SAMPLERS_HOME"] = tempfile.mkdtemp(prefix="jarvis-samplers-test-")


def pytest_configure(config) -> None:
    _isolate_sampler_registry()
    markexpr = getattr(config.option, "markexpr", "") or ""
    invocation = list(getattr(config.invocation_params, "args", ()) or ())
    value_options = _pytest_value_options(config._parser)
    config.option.markexpr = drop_default_slow_markexpr(
        markexpr, invocation, value_options=value_options,
    )


def pytest_sessionfinish(session, exitstatus: int) -> None:
    try:
        from test_distributed_acceptance import DistributedAcceptanceTests, _write_benchmark_report
    except ImportError:
        return
    if DistributedAcceptanceTests.metrics:
        _write_benchmark_report()
