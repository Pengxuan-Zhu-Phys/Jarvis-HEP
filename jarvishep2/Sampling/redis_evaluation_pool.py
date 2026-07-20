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
        logger=None,
    ) -> None:
        self.redis = redis
        self.build_sample = build_sample
        self.extract_logl = extract_logl or _default_extract_logl
        self.batch_size = max(1, int(batch_size))
        self.timeout = float(timeout)
        self.seed = int(seed)
        self.method = str(method)
        self.njobs = self.batch_size
        self._call_index = 0
        if logger is not None:
            self._logger = logger
        else:
            self._logger = get_jarvis_logger(
                "sampler.dynesty.pool",
                module="Dynesty.Pool",
            )

    @property
    def size(self) -> int:
        return self.njobs

    def map(self, func: Callable, iterable: Iterable) -> list[Any]:
        """Map *func* over *iterable*.

        Dynesty uses the same ``pool.map`` for **prior_transform** and
        **loglikelihood**. Prior transforms must run locally (they mint uuids);
        loglikelihoods go through Redis Workers + ``hep:feedback``.

        D13.7: unknown call shapes raise instead of silently evaluating a
        candidate physics logL in the control process.
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

        is_logl = self._is_loglikelihood_callable(func)
        is_ptform = self._is_prior_transform_callable(func)
        uuid_aug = looks_uuid_augmented(first)
        wants_remote = is_logl or uuid_aug

        if wants_remote:
            if self.redis is None:
                raise RuntimeError(
                    "RedisEvaluationPool.map: loglikelihood batch requires a "
                    "connected Redis queue; refusing local physics evaluation "
                    f"(func={_callable_label(func)}, uuid_augmented={uuid_aug})"
                )
            return self._redis_batch_logl(items)

        if is_ptform:
            return [func(x) for x in items]

        raise RuntimeError(
            "RedisEvaluationPool.map: ambiguous dispatch — refusing silent "
            "local fallback that could bypass Workers/calculators/Operas. "
            f"func={_callable_label(func)}, first_type={type(first).__name__}, "
            f"uuid_augmented={uuid_aug}. Expected prior_transform, "
            "LogLikelihood, uuid-augmented logL vectors, or SamplerArgument."
        )

    @staticmethod
    def _is_loglikelihood_callable(func: Callable) -> bool:
        if type(func).__name__ == "LogLikelihood":
            return True
        # dynesty wraps user logL as _function_wrapper(name='loglikelihood')
        # and also exposes LogLikelihood.map as bound method.
        label = _callable_label(func).lower()
        if "logl" in label or label in {"loglikelihood", "loglike"}:
            return True
        return False

    @staticmethod
    def _is_prior_transform_callable(func: Callable) -> bool:
        label = _callable_label(func).lower()
        if "prior_transform" in label or label in {"ptform", "priortransform"}:
            return True
        if type(func).__name__ == "_function_wrapper":
            wname = str(getattr(func, "name", "") or "").lower()
            if "prior" in wname and "logl" not in wname:
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
                # Destructive BLPOP: log before discard so resume/uuid bugs
                # leave a greppable trail (D13.7b).
                self._logger.warning(
                    "dropping unmatched hep:feedback uuid=%s "
                    "(pending_batch=%d expected=%d logL=%s)",
                    uuid or "<empty>",
                    len(remaining),
                    len(pending),
                    record.get("logL", ""),
                )
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
    """Map flat feedback logL to a dynesty-safe finite (D13.8).

    Wire uses real −∞ for unusable points; dynesty paths historically expect a
    large negative finite, so convert non-finite logL to ``-1e300`` here only.
    """
    from jarvishep2.feedback_return import feedback_logl, is_unusable_logl

    logl = feedback_logl(record)
    if is_unusable_logl(logl):
        return -1.0e300
    return float(logl)


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


def _callable_label(func: Callable) -> str:
    """Stable short name for dispatch diagnostics."""
    if type(func).__name__ == "_function_wrapper":
        return str(getattr(func, "name", "") or type(func).__name__)
    name = getattr(func, "__name__", None)
    if name:
        return str(name)
    return type(func).__name__


# typing re-export for annotations above
from collections.abc import Mapping  # noqa: E402

__all__ = ["RedisEvaluationPool"]
