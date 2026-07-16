#!/usr/bin/env python3
"""Redis-backed evaluation pool for dynesty (D13.5).

Implements the same *shape* as V1 ``JarvisFactoryAsyncPool`` (``map`` / ``submit``)
but transports work over Redis tasks + ``hep:feedback`` barriers.
"""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable, Sequence
from typing import Any

import numpy as np

from jarvishep2.Sampling.stateless_batch import deterministic_sampler_uuid
from jarvishep2.logging import get_jarvis_logger
from jarvishep2.redis_queue import RedisQueue
from jarvishep2.sample import Sample


class RedisEvaluationPool:
    """Pool interface for dynesty: batch logL via Redis Workers.

    Parameters
    ----------
    redis :
        Connected :class:`RedisQueue`.
    build_sample :
        ``build_sample(u_or_v, uuid) -> Sample`` — create a task Sample.
    extract_logl :
        ``extract_logl(feedback_record) -> float``.
    batch_size :
        Max tasks pushed before waiting (pipeline size).
    timeout :
        Seconds to wait for a full generation barrier.
    seed :
        For deterministic uuids when the input vector has no trailing uuid.
    """

    def __init__(
        self,
        redis: RedisQueue,
        *,
        build_sample: Callable[[Any, str], Sample],
        extract_logl: Callable[[Mapping], float] | None = None,
        batch_size: int = 16,
        timeout: float = 3600.0,
        seed: int = 0,
        method: str = "Dynesty",
    ) -> None:
        from collections.abc import Mapping as _Mapping  # noqa: F401 — typing

        self.redis = redis
        self.build_sample = build_sample
        self.extract_logl = extract_logl or _default_extract_logl
        self.batch_size = max(1, int(batch_size))
        self.timeout = float(timeout)
        self.seed = int(seed)
        self.method = str(method)
        self.njobs = self.batch_size
        self._call_index = 0
        self._logger = get_jarvis_logger("sampler.dynesty.pool")

    @property
    def size(self) -> int:
        return self.njobs

    def map(self, func: Callable, iterable: Iterable) -> list[Any]:
        """Map *func* over *iterable*.

        Dynesty uses the same ``pool.map`` for **prior_transform** and
        **loglikelihood**. Prior transforms must run locally (they mint uuids);
        loglikelihoods go through Redis Workers + ``hep:feedback``.
        """
        items = list(iterable)
        if not items:
            return []
        # Dynesty reuses pool.map for:
        #   * prior_transform (unit cube → coords[+uuid])
        #   * loglikelihood / LogLikelihood.map
        #   * internal_sampler.sample (SamplerArgument)
        # Only logL batches go remote.
        from jarvishep2.Sampling.Source.Dynesty.py.dynesty.jarvis_uuid import (
            looks_uuid_augmented,
        )

        first = items[0]
        if type(first).__name__ in {"SamplerArgument"}:
            return [func(x) for x in items]
        remote = self.redis is not None and (
            self._is_loglikelihood_callable(func) or looks_uuid_augmented(first)
        )
        if remote:
            return self._redis_batch_logl(items)
        return [func(x) for x in items]

    @staticmethod
    def _is_loglikelihood_callable(func: Callable) -> bool:
        if type(func).__name__ == "LogLikelihood":
            return True
        name = str(getattr(func, "__name__", "") or "").lower()
        if "logl" in name or name in {"loglikelihood", "loglike"}:
            return True
        return False

    def _redis_batch_logl(self, items: list[Any]) -> list[Any]:
        pending: dict[str, int] = {}
        samples: list[Sample] = []
        for item in items:
            uuid, payload = _uuid_and_payload(item, seed=self.seed, index=self._call_index)
            self._call_index += 1
            sample = self.build_sample(payload, uuid)
            sample.uuid = uuid
            pending[uuid] = len(samples)
            samples.append(sample)

        for i in range(0, len(samples), self.batch_size):
            chunk = samples[i : i + self.batch_size]
            self.redis.push_many_tasks([s.to_task_dict() for s in chunk])

        results: list[Any | None] = [None] * len(samples)
        deadline = time.monotonic() + max(1.0, self.timeout)
        remaining = set(pending.keys())
        while remaining:
            if time.monotonic() > deadline:
                raise TimeoutError(
                    f"RedisEvaluationPool timed out with {len(remaining)} pending"
                )
            wait = max(1, min(5, int(deadline - time.monotonic())))
            record = self.redis.pull_feedback(timeout=wait)
            if record is None:
                continue
            uuid = str(record.get("uuid", ""))
            if uuid not in remaining:
                continue
            remaining.discard(uuid)
            idx = pending[uuid]
            results[idx] = self.extract_logl(record)

        # Dynesty's live-point init does ``[_.val for _ in mapper(logl_wrap, …)]``
        # so return objects with a ``.val`` attribute (same as LoglOutput).
        class _LoglVal:
            __slots__ = ("val", "blob")

            def __init__(self, v: float) -> None:
                self.val = float(v)
                self.blob = None

        out: list[Any] = []
        for i, val in enumerate(results):
            if val is None:
                raise RuntimeError(f"missing feedback for sample index {i}")
            # If *func* is a LogLikelihood wrapper, prefer its blob handling later;
            # live-init only needs .val.
            out.append(_LoglVal(float(val)))
        return out

    def submit(self, func: Callable, item: Any):
        """Single-item submit — returns a trivial future-like holder."""
        result = self.map(func, [item])[0]

        class _F:
            def result(self_nonlocal):
                return result

        return _F()

    def close(self) -> None:
        return None

    def join(self) -> None:
        return None


def _default_extract_logl(record: dict[str, Any]) -> float:
    status = str(record.get("status", "Completed"))
    if status == "Failed":
        return -1.0e300
    obs = dict(record.get("observables") or {})
    if "LogL" in obs:
        return float(obs["LogL"])
    terms = [
        float(v)
        for k, v in obs.items()
        if str(k).startswith("LogL") and k != "LogL"
    ]
    if terms:
        return float(sum(terms))
    return -1.0e300


def _uuid_and_payload(item: Any, *, seed: int, index: int) -> tuple[str, np.ndarray]:
    """Extract uuid (if trailing) and numeric payload for Sample.u_coords."""
    arr = np.asarray(item, dtype=object).reshape(-1)
    if arr.size >= 2:
        last = arr[-1]
        if isinstance(last, (str, bytes, np.str_)):
            uuid = str(last)
            payload = np.asarray(arr[:-1], dtype=np.float64)
            return uuid, payload
        try:
            float(last)
        except (TypeError, ValueError):
            uuid = str(last)
            payload = np.asarray(arr[:-1], dtype=np.float64)
            return uuid, payload
    payload = np.asarray(arr, dtype=np.float64)
    uuid = deterministic_sampler_uuid(prefix="dynesty", seed=seed, sample_index=index)
    return uuid, payload


# typing re-export for annotations above
from collections.abc import Mapping  # noqa: E402

__all__ = ["RedisEvaluationPool"]
