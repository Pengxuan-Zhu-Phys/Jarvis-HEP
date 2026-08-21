"""Compatibility shim. Canonical module is ``jarvishep2.io.sample_bucket``."""

from __future__ import annotations

import sys

from jarvishep2.io import sample_bucket as _impl

sys.modules[__name__] = _impl
