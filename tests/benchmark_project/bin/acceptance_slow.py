#!/usr/bin/env python3
"""Slow calculator fixture for D7.1 worker-scaling acceptance (WP-D7.1)."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

# Default sleep matches tests/fixtures/slow_calculators/slow_a.py (0.35 s).
SLEEP_SEC = float(sys.argv[2]) if len(sys.argv) > 2 else 0.35

time.sleep(max(0.05, SLEEP_SEC))

if len(sys.argv) < 2:
    raise SystemExit("usage: acceptance_slow.py <output.json> [sleep_sec]")

out = Path(sys.argv[1])
out.write_text(json.dumps({"z": 1.0, "acceptance": True}), encoding="utf-8")