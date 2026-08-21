"""Compatibility shim. Canonical module is ``jarvishep2.cardload.task_config``."""

from __future__ import annotations

import sys

from jarvishep2.cardload import task_config as _impl

sys.modules[__name__] = _impl
