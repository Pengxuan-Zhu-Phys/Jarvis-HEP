#!/usr/bin/env python3
import json
import sys
import time
from pathlib import Path

time.sleep(0.35)
if len(sys.argv) < 2:
    raise SystemExit("usage: slow_b.py <output.json>")
out = Path(sys.argv[1])
out.write_text(json.dumps({"b": 2}), encoding="utf-8")