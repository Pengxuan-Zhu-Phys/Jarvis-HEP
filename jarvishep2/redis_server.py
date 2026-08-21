"""Compatibility shim. Canonical module is ``jarvishep2.queue.redis_server``."""

from __future__ import annotations

import sys

from jarvishep2.queue import redis_server as _impl

sys.modules[__name__] = _impl
