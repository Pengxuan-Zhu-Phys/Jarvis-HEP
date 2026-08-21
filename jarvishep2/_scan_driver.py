"""Compatibility shim. Canonical module is ``jarvishep2.runtime._scan_driver``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import _scan_driver as _impl

sys.modules[__name__] = _impl
