"""Compatibility shim. Canonical module is ``jarvishep2.runtime._runtime_supervisor``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import _runtime_supervisor as _impl

sys.modules[__name__] = _impl
