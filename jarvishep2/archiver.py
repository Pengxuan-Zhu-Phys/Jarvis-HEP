"""Compatibility shim. Canonical module is ``jarvishep2.io.archiver``."""

from __future__ import annotations

import sys

from jarvishep2.io import archiver as _impl

sys.modules[__name__] = _impl
