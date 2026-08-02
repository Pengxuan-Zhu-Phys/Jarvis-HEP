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
    """Delegate to ``FixedSetSampler.flush_batch`` when available."""
    method = getattr(sampler, "flush_batch", None)
    if callable(method):
        return int(method(batch))
    proposed = list(batch)
    should_submit = getattr(sampler, "should_submit_uuid", None)
    if callable(should_submit):
        batch = [sample for sample in batch if should_submit(sample.uuid)]
    if not batch:
        if proposed and callable(getattr(sampler, "record_persisted_proposals", None)):
            sampler.record_persisted_proposals()
        return 0
    if len(batch) == 1:
        sampler._submit(batch[0])
    else:
        sampler._submit_group(batch)
    record_batch = getattr(sampler, "record_submitted_batch", None)
    if callable(record_batch):
        record_batch([sample.uuid for sample in batch])
    capture = getattr(sampler, "capture_checkpoint_at_batch_boundary", None)
    if callable(capture):
        capture()
    return len(batch)


def run_stateless_distributed(
    sampler: Any,
    *,
    propose_next: Callable[[], Sample | None] | None = None,
) -> int:
    """Prefer ``sampler.run_distributed()``; fall back to propose loop."""
    if propose_next is None and hasattr(sampler, "run_distributed"):
        return int(sampler.run_distributed())
    if propose_next is None:
        raise TypeError("run_stateless_distributed requires propose_next or sampler.run_distributed")
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
