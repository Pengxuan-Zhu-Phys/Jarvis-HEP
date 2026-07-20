#!/usr/bin/env python3
"""Post-run sampler diagnostics under ``DATABASE/`` (D13.6).

Additive files next to ``samples.hdf5``:

* ``sampler_summary.json`` — method summary (accept rates, R-hat, logZ, …)
* ``chain_history.csv`` — MCMC/ensemble per-iteration accept/logl rows
  (design goal: ``chain_id``, ``step``, ``accepted``, ``weight``)

Nested methods already write ``dynesty_result.csv`` / ``multinest_result.csv``
with ``log_weight``; this module does not re-export nested dead-points.
"""

from __future__ import annotations

import csv
import json
import os
from collections.abc import Mapping, Sequence
from typing import Any


def _database_dir(task_result_dir: str) -> str:
    path = os.path.join(os.path.abspath(str(task_result_dir)), "DATABASE")
    os.makedirs(path, exist_ok=True)
    return path


def write_sampler_summary_json(
    task_result_dir: str,
    summary: Mapping[str, Any],
    *,
    filename: str = "sampler_summary.json",
) -> str:
    """Write ``DATABASE/sampler_summary.json`` (pretty JSON)."""
    out = os.path.join(_database_dir(task_result_dir), filename)
    with open(out, "w", encoding="utf-8") as handle:
        json.dump(dict(summary), handle, indent=2, sort_keys=False, default=str)
        handle.write("\n")
    return out


def write_chain_history_csv(
    task_result_dir: str,
    rows: Sequence[Mapping[str, Any]],
    *,
    filename: str = "chain_history.csv",
) -> str | None:
    """Write per-iteration MCMC history; returns path or None if *rows* empty."""
    if not rows:
        return None
    fieldnames = [
        "chain_id",
        "step",
        "accepted",
        "weight",
        "logl",
        "temperature",
        "state",
    ]
    # Preserve extra keys in stable order.
    extras: list[str] = []
    seen = set(fieldnames)
    for row in rows:
        for key in row.keys():
            name = str(key)
            if name not in seen:
                seen.add(name)
                extras.append(name)
    fields = fieldnames + extras
    out = os.path.join(_database_dir(task_result_dir), filename)
    with open(out, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            payload = {key: row.get(key, "") for key in fields}
            writer.writerow(payload)
    return out


def chain_history_rows_from_registry(registry: Any) -> list[dict[str, Any]]:
    """Flatten ChainRegistry histories into diagnostic CSV rows.

    ``weight`` is ``1`` when the Metropolis step was accepted, else ``0``
    (equal-weight chain view; nested importance weights live elsewhere).
    """
    rows: list[dict[str, Any]] = []
    if registry is None:
        return rows
    try:
        chains = list(registry.all())
    except Exception:
        return rows
    for chain in chains:
        cid = int(getattr(chain, "chain_id", -1))
        history = getattr(chain, "history", None)
        if history is None:
            continue
        try:
            events = history.all()
        except Exception:
            continue
        for ev in events:
            accepted = bool(getattr(ev, "accepted", False))
            rows.append(
                {
                    "chain_id": cid,
                    "step": int(getattr(ev, "iter", 0) or 0),
                    "accepted": int(accepted),
                    "weight": 1 if accepted else 0,
                    "logl": getattr(ev, "logl", ""),
                    "temperature": float(getattr(ev, "temperature", 1.0) or 1.0),
                    "state": str(getattr(ev, "state", "")),
                }
            )
    return rows


def export_mcmc_diagnostics(
    task_result_dir: str,
    *,
    summary: Mapping[str, Any] | None,
    registry: Any = None,
) -> dict[str, str]:
    """Write MCMC/ensemble diagnostics under DATABASE/. Returns written paths."""
    written: dict[str, str] = {}
    if summary:
        path = write_sampler_summary_json(task_result_dir, summary)
        written["sampler_summary"] = path
    rows = chain_history_rows_from_registry(registry)
    hist = write_chain_history_csv(task_result_dir, rows)
    if hist:
        written["chain_history"] = hist
    return written


def export_nested_diagnostics(
    task_result_dir: str,
    *,
    summary: Mapping[str, Any] | None,
) -> dict[str, str]:
    """Write nested-sampler summary JSON (CSV is method-specific elsewhere)."""
    written: dict[str, str] = {}
    if summary:
        path = write_sampler_summary_json(task_result_dir, summary)
        written["sampler_summary"] = path
    return written


# Column contract (documented for D13.6 DATABASE section)
DATABASE_CHAIN_DIAGNOSTIC_COLUMNS = (
    "chain_id",  # int — proposal / history row
    "step",  # int — Metropolis iteration index
    "stage",  # int — delayed-rejection stage (samples.hdf5 proposals)
    "accepted",  # 0|1 — chain_history.csv only (post-absorb)
    "weight",  # 0|1 chain; nested uses log_weight in *nest*_result.csv
)

DATABASE_NESTED_RESULT_COLUMNS = (
    "uuid",
    "log_weight",
    "log_Like",
    "log_PriorVolume",
    "log_Evidence",
    "log_Evidence_err",
    "samples_nlive",
    "ncall",
    "samples_it",
    "samples_id",
    "information",
)


__all__ = [
    "DATABASE_CHAIN_DIAGNOSTIC_COLUMNS",
    "DATABASE_NESTED_RESULT_COLUMNS",
    "chain_history_rows_from_registry",
    "export_mcmc_diagnostics",
    "export_nested_diagnostics",
    "write_chain_history_csv",
    "write_sampler_summary_json",
]
