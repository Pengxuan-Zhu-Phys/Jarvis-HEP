#!/usr/bin/env python3
import json
import sys
import time
from pathlib import Path

time.sleep(0.35)
if len(sys.argv) < 2:
    raise SystemExit("usage: slow_a.py <output.json>")
out = Path(sys.argv[1])
out.write_text(json.dumps({"a": 1}), encoding="utf-8")