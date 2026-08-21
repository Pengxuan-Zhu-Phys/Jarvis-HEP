"""Pytest hooks for the Jarvis-HEP V2 test suite."""

from __future__ import annotations

from collections.abc import Sequence


def drop_default_slow_markexpr(
    markexpr: str,
    invocation_args: Sequence[str],
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
    while idx < len(args):
        arg = args[idx]
        if arg in {"-m", "--markexpr"}:
            return expr
        if arg.startswith("-m") or arg.startswith("--markexpr="):
            return expr
        if arg == "--":
            named_paths.extend(args[idx + 1 :])
            break
        if not arg.startswith("-"):
            named_paths.append(arg)
        idx += 1
    for path in named_paths:
        node = path.split("::", 1)[0]
        if node.endswith(".py") or "::" in path:
            return ""
    return expr


def pytest_configure(config) -> None:
    markexpr = getattr(config.option, "markexpr", "") or ""
    invocation = list(getattr(config.invocation_params, "args", ()) or ())
    config.option.markexpr = drop_default_slow_markexpr(markexpr, invocation)


def pytest_sessionfinish(session, exitstatus: int) -> None:
    try:
        from test_distributed_acceptance import DistributedAcceptanceTests, _write_benchmark_report
    except ImportError:
        return
    if DistributedAcceptanceTests.metrics:
        _write_benchmark_report()