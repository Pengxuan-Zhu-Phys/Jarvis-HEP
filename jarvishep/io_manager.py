#!/usr/bin/env python3
from __future__ import annotations

import asyncio
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from typing import Any, Callable


class IOManager:
    """Internal calculator-workflow executor for blocking file work."""

    def __init__(self, max_workers: int = 8, *, thread_name_prefix: str = "jarvis-hep-io"):
        workers = max(1, int(max_workers))
        self._max_workers = workers
        self._executor = ThreadPoolExecutor(
            max_workers=workers,
            thread_name_prefix=thread_name_prefix,
        )
        self._closed = False

    @property
    def max_workers(self) -> int:
        return int(self._max_workers)

    async def run_blocking(self, fn: Callable[..., Any], *args, **kwargs) -> Any:
        if self._closed:
            raise RuntimeError("IOManager is already shut down.")
        loop = asyncio.get_running_loop()
        bound = partial(fn, *args, **kwargs)
        return await loop.run_in_executor(self._executor, bound)

    def shutdown(self, *, wait: bool = True, cancel_futures: bool = True) -> None:
        if self._closed:
            return
        self._closed = True
        self._executor.shutdown(wait=wait, cancel_futures=cancel_futures)
