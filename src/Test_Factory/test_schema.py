#!/usr/bin/env python3 
from jsonschema import validate, ValidationError

schema = {
    "type": "object",
    "properties": {
        "type": {"const": "Dump"},
        "variables": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "expression": {"type": "string"},
                    "entry": {"type": "string"}
                },
                "required": ["name"],
                "additionalProperties": False
            }
        }
    },
    "required": ["type", "variables"]
}

data = {
    "type": "Dump",
    "variables": [
        {'name': 'xx', 'expression': 'x * Pi'}, 
        {'name': 'yy', 'expression': 'y * Pi'}, 
        {'name': 'cx', 'expression': '(x + y) * Pi', 'entry': 'this.is.an.entry'}]
}

try:
    validate(instance=data, schema=schema)
    print("Validation passed.")
except ValidationError as e:
    print(f"Validation failed: {e}")
