"""Compatibility shim. Canonical module is ``jarvishep2.cardload.task_schema``."""

from __future__ import annotations

import sys

from jarvishep2.cardload import task_schema as _impl

sys.modules[__name__] = _impl
