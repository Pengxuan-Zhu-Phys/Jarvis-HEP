"""Compatibility shim. Canonical module is ``jarvishep2.runtime.worker``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import worker as _impl

sys.modules[__name__] = _impl
