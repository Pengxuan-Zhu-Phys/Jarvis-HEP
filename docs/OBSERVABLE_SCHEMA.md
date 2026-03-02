# Observable Schema And CSV Flattening

Jarvis-HEP now stores observables with two layers:

1. Full-fidelity storage in HDF5 (`samples.N.hdf5`)
2. User-configurable CSV export via schema (`samples.schema.json`)

The schema is editable. You can change flatten rules and regenerate CSV without rerunning the scan.

---

## Files Produced

Given database root `.../DATABASE/samples.hdf5`, Jarvis-HEP will use:

- `.../DATABASE/samples.N.hdf5` for chunked raw records
- `.../DATABASE/samples.schema.json` for column metadata and flatten rules
- `.../DATABASE/samples.N.csv` for CSV exports

`N` is the HDF5 chunk index (`0`, `1`, ...).

---

## Schema Structure

Typical schema example:

```json
{
  "version": 1,
  "pathroot": "/path/to/DATABASE/samples",
  "created_at": "2026-03-01T10:00:00+00:00",
  "updated_at": "2026-03-01T10:30:00+00:00",
  "flatten_defaults": {
    "array": "json",
    "list": "json",
    "dict": "json",
    "mixed": "json"
  },
  "columns": {
    "uuid": {
      "kind": "str",
      "dtype": "str",
      "flatten": { "mode": "scalar" }
    },
    "LogL": {
      "kind": "float",
      "dtype": "float",
      "flatten": { "mode": "scalar" }
    },
    "arr": {
      "kind": "ndarray",
      "dtype": "float64",
      "shape": [2, 2],
      "ndim": 2,
      "flatten": { "mode": "json" }
    }
  }
}
```

Notes:

- `shape` or `ndim` can become `"variable"` if records are inconsistent.
- `kind` can become `"mixed"` if a column changes type across records.
- Unknown or missing fields are auto-normalized with fallback defaults.

---

## Flatten Modes

Each column uses `columns.<name>.flatten.mode`.

Supported modes:

- `scalar`
  - Intended for scalar fields.
  - If value is structured (list/dict), exporter falls back to JSON cell text.
- `json`
  - Structured value is written as JSON string in one CSV cell.
- `split`
  - Structured value is expanded into multiple columns.
  - Arrays/lists use index notation (`arr[0]`, `arr[1]`, `mat[0][1]`).
  - Dicts use dot notation (`meta.mass`, `meta.flags[0]`).
  - Supports optional column renaming via `flatten.name_map`.
- `drop`
  - Column is skipped in CSV output.

Default fallback logic:

- Missing/invalid mode -> fallback by kind (`ndarray/list/dict/mixed -> json`, others -> `scalar`)

---

## Recommended Workflow

1. Run scan normally.
2. Open `samples.schema.json`.
3. Modify `flatten.mode` for columns you want to reshape in CSV.
4. Regenerate CSV:

```bash
Jarvis /path/to/task.yaml --convert
```

You can repeat step 3 and step 4 without rerunning sampling.

---

## Example: Split One Array Column

Schema edit:

```json
{
  "columns": {
    "arr": {
      "flatten": { "mode": "split" }
    }
  }
}
```

Resulting CSV headers may include:

```text
uuid,LogL,arr[0],arr[1],...
```

Jarvis-HEP will auto-generate a default identity map after split export:

```json
{
  "columns": {
    "arr": {
      "flatten": {
        "mode": "split",
        "name_map": {
          "arr[0]": "arr[0]",
          "arr[1]": "arr[1]"
        }
      }
    }
  }
}
```

You can customize names by editing map values:

```json
{
  "columns": {
    "arr": {
      "flatten": {
        "mode": "split",
        "name_map": {
          "arr[0]": "arr_left",
          "arr[1]": "arr_right"
        }
      }
    }
  }
}
```

---

## Safety And Recovery

Jarvis-HEP applies safe fallbacks:

- Schema file missing:
  - A default schema is created automatically.
- Schema file broken (invalid JSON):
  - Conversion falls back to default schema and emits warning.
- Invalid flatten mode:
  - Mode is auto-corrected to a valid fallback and written back.

Warnings are emitted in:

- runtime logger (`hdf5writer`)
- Python warnings path (`utils.convert_hdf5_to_csv`)

---

## Compatibility

- Existing scalar workflows remain unchanged.
- Array/list/dict observables are preserved in HDF5.
- CSV behavior is configurable and reversible via schema edits.
