"""Compatibility shim. Canonical module is ``jarvishep2.cli.client``."""

from __future__ import annotations

import sys

from jarvishep2.cli import client as _impl

sys.modules[__name__] = _impl
