#!/usr/bin/env python3
"""Unit tests for jarvishep2.log_kv."""

from __future__ import annotations

from jarvishep2.log_kv import format_two_column_log


def _assert_max_width(block: str, max_width: int) -> None:
    for line in block.splitlines():
        assert len(line) <= max_width, f"line exceeds max width {max_width}: {line!r}"


def test_two_column_log_wraps_long_values_under_80_columns() -> None:
    long_value = "x" * 160
    rendered = format_two_column_log(
        "Subprocess scheduler configured",
        [("max_concurrency", 8), ("very_long_value", long_value)],
        max_width=80,
    )
    _assert_max_width(rendered, 80)
    assert "very_long_value" in rendered
    assert long_value[:40] in rendered


def test_two_column_log_wraps_long_keys_under_80_columns() -> None:
    long_key = "very_long_configuration_key_name_" * 3
    rendered = format_two_column_log(
        "Config dump",
        [(long_key, "ok")],
        max_width=80,
    )
    _assert_max_width(rendered, 80)
    assert "ok" in rendered


def test_two_column_log_wraps_title_under_80_columns() -> None:
    long_title = "Scheduler report " + ("#" * 120)
    rendered = format_two_column_log(
        long_title,
        [("status", "ready")],
        max_width=80,
    )
    _assert_max_width(rendered, 80)
    assert "status" in rendered