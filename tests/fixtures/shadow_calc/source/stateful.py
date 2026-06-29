#!/usr/bin/env python3
import json
import sys
from pathlib import Path

if len(sys.argv) < 2:
    raise SystemExit("usage: stateful.py <output.json>")

marker = Path("isolation.marker")
count = int(marker.read_text(encoding="utf-8")) if marker.exists() else 0
count += 1
marker.write_text(str(count), encoding="utf-8")
Path(sys.argv[1]).write_text(json.dumps({"isolation_count": count}), encoding="utf-8")