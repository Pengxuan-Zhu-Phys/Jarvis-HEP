"""Pytest hooks for the Jarvis-HEP V2 test suite."""

from __future__ import annotations


def pytest_sessionfinish(session, exitstatus: int) -> None:
    try:
        from test_distributed_acceptance import DistributedAcceptanceTests, _write_benchmark_report
    except ImportError:
        return
    if DistributedAcceptanceTests.metrics:
        _write_benchmark_report()