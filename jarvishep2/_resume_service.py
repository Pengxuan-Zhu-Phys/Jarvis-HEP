"""Compatibility shim. Canonical module is ``jarvishep2.runtime._resume_service``."""

from __future__ import annotations

import sys

from jarvishep2.runtime import _resume_service as _impl

sys.modules[__name__] = _impl
