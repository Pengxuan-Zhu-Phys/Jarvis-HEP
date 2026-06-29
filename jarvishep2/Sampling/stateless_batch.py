#!/usr/bin/env python3
"""Shared batch submission helpers for stateless V2 samplers."""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from typing import Any

from jarvishep2.sample import Sample


def deterministic_sampler_uuid(*, prefix: str, seed: int, sample_index: int) -> str:
    digest = hashlib.sha256(
        f"{prefix}:{int(seed)}:{int(sample_index)}".encode("utf-8")
    ).hexdigest()
    return (
        f"{digest[0:8]}-{digest[8:12]}-{digest[12:16]}-{digest[16:20]}-{digest[20:32]}"
    )


def flush_batch(sampler: Any, batch: list[Sample]) -> int:
    if not batch:
        return 0
    if len(batch) == 1:
        sampler._submit(batch[0])
    else:
        sampler._submit_group(batch)
    sampler._submitted_uuids.extend(sample.uuid for sample in batch)
    return len(batch)


def run_stateless_distributed(
    sampler: Any,
    *,
    propose_next: Callable[[], Sample | None],
) -> int:
    pushed = 0
    batch: list[Sample] = []
    batch_size = max(1, int(getattr(sampler, "_batch_size", 1) or 1))
    while True:
        sample = propose_next()
        if sample is None:
            break
        batch.append(sample)
        if len(batch) >= batch_size:
            pushed += flush_batch(sampler, batch)
            batch = []
    if batch:
        pushed += flush_batch(sampler, batch)
    return pushed


__all__ = ["deterministic_sampler_uuid", "flush_batch", "run_stateless_distributed"]