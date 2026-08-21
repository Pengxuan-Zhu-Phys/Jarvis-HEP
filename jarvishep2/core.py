"""Compatibility shim. Canonical module is ``jarvishep2.runtime.core``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import core as _impl

sys.modules[__name__] = _impl
