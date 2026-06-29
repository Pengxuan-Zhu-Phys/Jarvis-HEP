#!/usr/bin/env python3
import json
import sys
from pathlib import Path

if len(sys.argv) < 2:
    raise SystemExit("usage: safe.py <output.json>")

Path(sys.argv[1]).write_text(json.dumps({"safe": True}), encoding="utf-8")