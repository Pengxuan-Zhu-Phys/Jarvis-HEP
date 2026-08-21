"""Compatibility shim. Canonical module is ``jarvishep2.io.file_ops``."""

from __future__ import annotations

import sys

from jarvishep2.io import file_ops as _impl

sys.modules[__name__] = _impl
