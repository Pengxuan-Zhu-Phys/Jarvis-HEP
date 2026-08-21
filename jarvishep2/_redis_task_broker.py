"""Compatibility shim. Canonical module is ``jarvishep2.queue._redis_task_broker``."""

from __future__ import annotations

import sys

from jarvishep2.queue import _redis_task_broker as _impl

sys.modules[__name__] = _impl
