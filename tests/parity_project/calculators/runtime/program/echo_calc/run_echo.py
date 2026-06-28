#!/usr/bin/env python3
import json

payload = json.load(open("in.json", encoding="utf-8"))
json.dump(
    {"calc_z": float(payload["x"]) * 2.0},
    open("out.json", "w", encoding="utf-8"),
)