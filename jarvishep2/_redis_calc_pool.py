"""Compatibility shim. Canonical module is ``jarvishep2.queue._redis_calc_pool``."""

from __future__ import annotations

import sys

from jarvishep2.queue import _redis_calc_pool as _impl

sys.modules[__name__] = _impl
