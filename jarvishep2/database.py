"""Compatibility shim. Canonical module is ``jarvishep2.io.database``."""

from __future__ import annotations

import sys

from jarvishep2.io import database as _impl

sys.modules[__name__] = _impl
