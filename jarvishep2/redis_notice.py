"""Compatibility shim. Canonical module is ``jarvishep2.queue.redis_notice``."""

from __future__ import annotations

import sys

from jarvishep2.queue import redis_notice as _impl

sys.modules[__name__] = _impl
