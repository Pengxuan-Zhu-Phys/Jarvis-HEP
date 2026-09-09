"""Compatibility shim. Canonical module is ``jarvishep2.queue._redis_proc_board``."""

from __future__ import annotations

import sys

from jarvishep2.queue import _redis_proc_board as _impl

sys.modules[__name__] = _impl
