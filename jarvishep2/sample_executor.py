"""Compatibility shim. Canonical module is ``jarvishep2.runtime.sample_executor``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import sample_executor as _impl

sys.modules[__name__] = _impl
