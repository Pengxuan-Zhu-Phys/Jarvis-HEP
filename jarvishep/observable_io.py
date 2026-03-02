from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from typing import Any

import numpy as np
import sympy

VALID_FLATTEN_MODES = {"json", "split", "drop", "scalar"}
SCHEMA_VERSION_CURRENT = 1
SCHEMA_VERSION_MIN_COMPATIBLE = 1


def _is_sympy_boolean(value: Any) -> bool:
    try:
        from sympy.logic.boolalg import Boolean

        return isinstance(value, Boolean)
    except Exception:
        return False


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def make_json_compatible(value: Any) -> Any:
    """Recursively convert values into JSON-serializable Python objects."""
    if isinstance(value, dict):
        return {str(k): make_json_compatible(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [make_json_compatible(v) for v in value]
    if isinstance(value, np.ndarray):
        return make_json_compatible(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.str_):
        return str(value)
    if isinstance(value, (sympy.Float, sympy.Integer)):
        return float(value)
    if _is_sympy_boolean(value):
        return bool(value)
    if isinstance(value, sympy.Basic):
        if value.is_real and value.is_number:
            try:
                return float(value)
            except Exception:
                return str(value)
        return str(value)
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def _infer_kind(value: Any) -> str:
    if isinstance(value, np.ndarray):
        return "ndarray"
    if isinstance(value, (list, tuple, set)):
        return "list"
    if isinstance(value, dict):
        return "dict"
    if isinstance(value, np.generic):
        return "numpy_scalar"
    if isinstance(value, (sympy.Float, sympy.Integer)) or _is_sympy_boolean(value):
        return "sympy_scalar"
    if value is None:
        return "null"
    if isinstance(value, bool):
        return "bool"
    if isinstance(value, int):
        return "int"
    if isinstance(value, float):
        return "float"
    if isinstance(value, str):
        return "str"
    return type(value).__name__


def _infer_list_shape(values: list[Any]) -> list[int]:
    shape: list[int] = []
    current: Any = values
    while isinstance(current, list):
        shape.append(len(current))
        if not current:
            break
        if all(isinstance(item, list) for item in current):
            first_len = len(current[0])
            if all(len(item) == first_len for item in current):
                current = current[0]
                continue
        break
    return shape


def _infer_list_element_kind(values: list[Any]) -> str | None:
    for item in values:
        if item is None:
            continue
        if isinstance(item, list):
            nested = _infer_list_element_kind(item)
            if nested:
                return f"list<{nested}>"
            return "list"
        return _infer_kind(item)
    return None


def infer_column_descriptor(value: Any) -> dict[str, Any]:
    kind = _infer_kind(value)
    desc: dict[str, Any] = {"kind": kind}

    if isinstance(value, np.ndarray):
        desc["shape"] = list(value.shape)
        desc["ndim"] = int(value.ndim)
        desc["dtype"] = str(value.dtype)
        return desc

    if isinstance(value, (list, tuple, set)):
        seq = list(value)
        desc["shape"] = _infer_list_shape(seq)
        desc["ndim"] = len(desc["shape"])
        element_kind = _infer_list_element_kind(seq)
        if element_kind:
            desc["element_kind"] = element_kind
        return desc

    if isinstance(value, dict):
        desc["keys_preview"] = [str(k) for k in list(value.keys())[:20]]
        return desc

    if isinstance(value, np.generic):
        desc["dtype"] = str(value.dtype)
        return desc

    if isinstance(value, (sympy.Float, sympy.Integer)):
        desc["dtype"] = type(value).__name__
        return desc

    if value is not None:
        desc["dtype"] = type(value).__name__
    return desc


def _default_mode_for_kind(kind: str) -> str:
    if kind in {"array", "ndarray", "list", "dict", "mixed"}:
        return "split"
    return "scalar"


def create_default_schema(pathroot: str) -> dict[str, Any]:
    now = utc_now_iso()
    return {
        "version": SCHEMA_VERSION_CURRENT,
        "pathroot": pathroot,
        "created_at": now,
        "updated_at": now,
        "flatten_defaults": {
            "array": "split",
            "list": "split",
            "dict": "split",
            "mixed": "split",
        },
        "columns": {},
    }


def _atomic_write_json(path: str, data: dict[str, Any]) -> None:
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
        f.flush()
    os.replace(tmp_path, path)


def load_schema(path: str, pathroot: str) -> dict[str, Any]:
    warnings: list[str] = []

    if not os.path.exists(path):
        warnings.append(f"Schema file not found; creating default schema -> {path}")
        schema = create_default_schema(pathroot)
        schema["_warnings"] = warnings
        return schema

    if os.path.exists(path):
        try:
            with open(path, "r") as f:
                loaded = json.load(f)
            if not isinstance(loaded, dict):
                raise ValueError("Schema root is not an object")
            if not isinstance(loaded.get("columns", {}), dict):
                raise ValueError("Schema 'columns' must be an object")

            normalized, changed, warn_msgs = sanitize_schema(loaded, pathroot=pathroot)
            warnings.extend(warn_msgs)
            if changed:
                warnings.append("Schema had invalid/missing fields; normalized with fallback defaults.")
            normalized["_warnings"] = warnings
            return normalized
        except Exception as exc:
            warnings.append(f"Failed to read schema ({exc}); fallback to default schema -> {path}")
            schema = create_default_schema(pathroot)
            schema["_warnings"] = warnings
            return schema

    schema = create_default_schema(pathroot)
    schema["_warnings"] = warnings
    return schema


def save_schema(path: str, schema: dict[str, Any]) -> None:
    schema["updated_at"] = utc_now_iso()
    to_save = dict(schema)
    to_save.pop("_warnings", None)
    _atomic_write_json(path, to_save)


def sanitize_schema(schema: dict[str, Any], pathroot: str) -> tuple[dict[str, Any], bool, list[str]]:
    """Normalize user-edited schema and collect warnings.

    Version policy:
    - Missing/invalid/legacy versions are migrated to ``SCHEMA_VERSION_CURRENT``.
    - Future versions are preserved (not downgraded) and normalized best-effort
      for known fields to keep backward-compatible behavior.
    """
    normalized = dict(schema)
    changed = False
    warnings: list[str] = []

    version, version_changed, version_warning = _normalize_schema_version(normalized.get("version"))
    normalized["version"] = version
    changed = changed or version_changed
    if version_warning:
        warnings.append(version_warning)

    if normalized.get("pathroot") != pathroot:
        normalized["pathroot"] = pathroot
        changed = True

    if "created_at" not in normalized:
        normalized["created_at"] = utc_now_iso()
        changed = True
    normalized.setdefault("updated_at", utc_now_iso())

    defaults = normalized.get("flatten_defaults")
    if not isinstance(defaults, dict):
        defaults = {}
        normalized["flatten_defaults"] = defaults
        changed = True
        warnings.append("schema.flatten_defaults is invalid; replaced with defaults.")

    for kind in ("array", "list", "dict", "mixed"):
        mode = str(defaults.get(kind, _default_mode_for_kind(kind))).strip().lower()
        if mode not in VALID_FLATTEN_MODES:
            fallback = _default_mode_for_kind(kind)
            defaults[kind] = fallback
            changed = True
            warnings.append(
                f"schema.flatten_defaults.{kind} has invalid mode '{mode}'; fallback to '{fallback}'."
            )
        else:
            if defaults.get(kind) != mode:
                defaults[kind] = mode
                changed = True

    columns = normalized.get("columns")
    if not isinstance(columns, dict):
        normalized["columns"] = {}
        changed = True
        warnings.append("schema.columns is invalid; reset to empty object.")
        return normalized, changed, warnings

    for name, meta in list(columns.items()):
        if not isinstance(meta, dict):
            columns[name] = {"kind": "unknown", "flatten": {"mode": "scalar"}}
            changed = True
            warnings.append(f"schema.columns.{name} is not an object; reset with fallback config.")
            continue

        kind = str(meta.get("kind", "unknown"))
        flatten_cfg = meta.get("flatten")
        if not isinstance(flatten_cfg, dict):
            flatten_cfg = {"mode": _default_mode_for_kind(kind)}
            meta["flatten"] = flatten_cfg
            changed = True
            warnings.append(f"schema.columns.{name}.flatten is invalid; fallback applied.")

        mode = str(flatten_cfg.get("mode", "")).strip().lower()
        if not mode:
            mode = _default_mode_for_kind(kind)
            flatten_cfg["mode"] = mode
            changed = True
            warnings.append(f"schema.columns.{name}.flatten.mode missing; fallback to '{mode}'.")
        elif mode not in VALID_FLATTEN_MODES:
            fallback = _default_mode_for_kind(kind)
            flatten_cfg["mode"] = fallback
            changed = True
            warnings.append(
                f"schema.columns.{name}.flatten.mode '{mode}' invalid; fallback to '{fallback}'."
            )
        elif flatten_cfg.get("mode") != mode:
            flatten_cfg["mode"] = mode
            changed = True

        name_map = flatten_cfg.get("name_map", None)
        if name_map is not None and not isinstance(name_map, dict):
            flatten_cfg["name_map"] = {}
            changed = True
            warnings.append(f"schema.columns.{name}.flatten.name_map is invalid; reset to empty object.")
        elif isinstance(name_map, dict):
            normalized_name_map: dict[str, str] = {}
            map_changed = False
            for kk, vv in name_map.items():
                nk = str(kk)
                nv = str(vv).strip()
                if not nv:
                    nv = nk
                    map_changed = True
                normalized_name_map[nk] = nv
                if nk != kk or nv != vv:
                    map_changed = True
            if map_changed:
                flatten_cfg["name_map"] = normalized_name_map
                changed = True

    return normalized, changed, warnings


def _normalize_schema_version(raw_version: Any) -> tuple[int, bool, str | None]:
    if raw_version is None:
        return (
            SCHEMA_VERSION_CURRENT,
            True,
            (
                f"schema.version missing; migrated to {SCHEMA_VERSION_CURRENT} "
                f"(legacy compatibility mode)."
            ),
        )

    parsed: int | None = None
    try:
        parsed = int(raw_version)
    except Exception:
        return (
            SCHEMA_VERSION_CURRENT,
            True,
            (
                f"schema.version '{raw_version}' invalid; migrated to "
                f"{SCHEMA_VERSION_CURRENT}."
            ),
        )

    changed = parsed != raw_version
    if parsed < SCHEMA_VERSION_MIN_COMPATIBLE:
        return (
            SCHEMA_VERSION_CURRENT,
            True,
            (
                f"schema.version {parsed} is below minimum compatible "
                f"{SCHEMA_VERSION_MIN_COMPATIBLE}; migrated to {SCHEMA_VERSION_CURRENT}."
            ),
        )

    if parsed > SCHEMA_VERSION_CURRENT:
        return (
            parsed,
            changed,
            (
                f"schema.version {parsed} is newer than supported {SCHEMA_VERSION_CURRENT}; "
                "applying best-effort compatibility normalization for known fields."
            ),
        )

    return parsed, changed, None


def _merge_descriptor(existing: dict[str, Any], current: dict[str, Any]) -> tuple[dict[str, Any], bool]:
    merged = dict(existing)
    changed = False

    old_kind = str(merged.get("kind", "unknown"))
    new_kind = str(current.get("kind", "unknown"))
    if old_kind != new_kind:
        kinds = set(merged.get("variants", []))
        kinds.update({old_kind, new_kind})
        merged["kind"] = "mixed"
        merged["variants"] = sorted(kinds)
        changed = True

    for key in ("dtype", "element_kind"):
        new_val = current.get(key)
        if new_val is None:
            continue
        old_val = merged.get(key)
        if old_val is None:
            merged[key] = new_val
            changed = True
        elif old_val != new_val and old_val != "mixed":
            merged[key] = "mixed"
            changed = True

    if "shape" in current:
        new_shape = current["shape"]
        old_shape = merged.get("shape")
        if old_shape is None:
            merged["shape"] = new_shape
            changed = True
        elif old_shape != new_shape and old_shape != "variable":
            merged["shape"] = "variable"
            changed = True

    if "ndim" in current:
        new_ndim = current["ndim"]
        old_ndim = merged.get("ndim")
        if old_ndim is None:
            merged["ndim"] = new_ndim
            changed = True
        elif old_ndim != new_ndim and old_ndim != "variable":
            merged["ndim"] = "variable"
            changed = True

    if "keys_preview" in current:
        merged_keys = set(merged.get("keys_preview", []))
        incoming_keys = [str(k) for k in current.get("keys_preview", [])]
        before = len(merged_keys)
        merged_keys.update(incoming_keys)
        if len(merged_keys) != before:
            merged["keys_preview"] = sorted(merged_keys)[:20]
            changed = True

    if "flatten" not in merged or not isinstance(merged.get("flatten"), dict):
        merged["flatten"] = {"mode": _default_mode_for_kind(str(merged.get("kind", "unknown")))}
        changed = True
    elif "mode" not in merged["flatten"]:
        merged["flatten"]["mode"] = _default_mode_for_kind(str(merged.get("kind", "unknown")))
        changed = True

    return merged, changed


def update_schema_with_record(schema: dict[str, Any], record: Any) -> bool:
    """Update schema using one observable record. Returns True if modified."""
    if not isinstance(record, dict):
        return False

    columns = schema.setdefault("columns", {})
    changed = False

    for raw_key, value in record.items():
        key = str(raw_key)
        current = infer_column_descriptor(value)
        existing = columns.get(key)
        if existing is None:
            current["flatten"] = {"mode": _default_mode_for_kind(current.get("kind", "unknown"))}
            columns[key] = current
            changed = True
            continue

        merged, merged_changed = _merge_descriptor(existing, current)
        if merged_changed:
            columns[key] = merged
            changed = True

    if changed:
        schema["updated_at"] = utc_now_iso()
    return changed


def _json_cell(value: Any) -> str:
    return json.dumps(make_json_compatible(value), ensure_ascii=False, separators=(",", ":"))


def _flatten_nested(prefix: str, value: Any, row: dict[str, Any]) -> None:
    if isinstance(value, dict):
        if not value:
            row[prefix] = "{}"
            return
        for key, item in value.items():
            _flatten_nested(f"{prefix}.{key}", item, row)
        return

    if isinstance(value, list):
        if not value:
            row[prefix] = "[]"
            return
        for idx, item in enumerate(value):
            _flatten_nested(f"{prefix}[{idx}]", item, row)
        return

    row[prefix] = value


def _resolve_flatten_mode(schema: dict[str, Any], key: str, value: Any) -> str:
    columns = schema.get("columns", {}) if isinstance(schema, dict) else {}
    col_meta = columns.get(key, {}) if isinstance(columns, dict) else {}

    flatten_cfg = col_meta.get("flatten", {}) if isinstance(col_meta, dict) else {}
    if isinstance(flatten_cfg, dict) and flatten_cfg.get("mode"):
        mode = str(flatten_cfg["mode"]).strip().lower()
        if mode in VALID_FLATTEN_MODES:
            return mode

    kind = str(col_meta.get("kind", _infer_kind(value)))
    if kind == "ndarray":
        kind = "array"
    defaults = schema.get("flatten_defaults", {}) if isinstance(schema, dict) else {}
    mode = str(defaults.get(kind, _default_mode_for_kind(kind))).strip().lower()
    if mode not in VALID_FLATTEN_MODES:
        return _default_mode_for_kind(kind)
    return mode


def _ensure_split_name_map(
    schema: dict[str, Any],
    column: str,
    generated_keys: list[str],
    value: Any,
    populate_name_map: bool,
) -> tuple[dict[str, str], bool]:
    if not isinstance(schema, dict):
        return {}, False

    changed = False
    columns = schema.setdefault("columns", {})
    if not isinstance(columns, dict):
        schema["columns"] = {}
        columns = schema["columns"]
        changed = True

    col_meta = columns.get(column)
    if not isinstance(col_meta, dict):
        if populate_name_map:
            col_meta = infer_column_descriptor(value)
            col_meta["flatten"] = {"mode": "split", "name_map": {}}
            columns[column] = col_meta
            changed = True
        else:
            col_meta = {}

    flatten_cfg = col_meta.get("flatten")
    if not isinstance(flatten_cfg, dict):
        flatten_cfg = {"mode": "split"}
        if populate_name_map:
            col_meta["flatten"] = flatten_cfg
            changed = True

    name_map = flatten_cfg.get("name_map")
    if not isinstance(name_map, dict):
        name_map = {}
        if populate_name_map:
            flatten_cfg["name_map"] = name_map
            changed = True

    if populate_name_map:
        for key in generated_keys:
            if key not in name_map:
                name_map[key] = key
                changed = True
            else:
                alias = str(name_map[key]).strip()
                if not alias:
                    name_map[key] = key
                    changed = True
                elif alias != name_map[key]:
                    name_map[key] = alias
                    changed = True

    normalized_map: dict[str, str] = {}
    for key in generated_keys:
        alias = str(name_map.get(key, key)).strip()
        if not alias:
            alias = key
        normalized_map[key] = alias
    return normalized_map, changed


def _flatten_record_for_csv_internal(
    record: dict[str, Any],
    schema: dict[str, Any],
    populate_name_map: bool,
) -> tuple[dict[str, Any], bool]:
    """Apply per-column flatten policies to a record for CSV export."""
    if not isinstance(record, dict):
        return {}, False

    row: dict[str, Any] = {}
    schema_changed = False
    normalized_record = make_json_compatible(record)

    for key, value in normalized_record.items():
        column = str(key)
        mode = _resolve_flatten_mode(schema, column, value)

        if mode == "drop":
            continue

        if mode == "split":
            split_items: dict[str, Any] = {}
            if isinstance(value, (dict, list)):
                _flatten_nested(column, value, split_items)
            else:
                split_items[column] = value

            generated_keys = list(split_items.keys())
            name_map, changed = _ensure_split_name_map(
                schema=schema,
                column=column,
                generated_keys=generated_keys,
                value=value,
                populate_name_map=populate_name_map,
            )
            schema_changed = schema_changed or changed
            for raw_key, raw_value in split_items.items():
                row[name_map.get(raw_key, raw_key)] = raw_value
            continue

        if isinstance(value, (dict, list)):
            row[column] = _json_cell(value)
        else:
            row[column] = value

    return row, schema_changed


def flatten_record_for_csv(record: dict[str, Any], schema: dict[str, Any] | None = None) -> dict[str, Any]:
    schema = schema or {}
    row, _ = _flatten_record_for_csv_internal(
        record=record,
        schema=schema,
        populate_name_map=False,
    )
    return row


def flatten_records_for_csv(
    records: list[Any],
    schema: dict[str, Any] | None = None,
    populate_name_map: bool = False,
) -> tuple[list[dict[str, Any]], bool]:
    schema = schema or {}
    rows: list[dict[str, Any]] = []
    changed = False

    for record in records:
        if not isinstance(record, dict):
            continue
        row, row_changed = _flatten_record_for_csv_internal(
            record=record,
            schema=schema,
            populate_name_map=populate_name_map,
        )
        rows.append(row)
        changed = changed or row_changed
    return rows, changed


def collect_csv_fieldnames(rows: list[dict[str, Any]]) -> list[str]:
    """Build ordered union of CSV column names from flattened rows."""
    seen: set[str] = set()
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    return fieldnames


def resolve_schema_path(pathroot_or_hdf5: str) -> str:
    """Resolve schema path from a root path or an indexed HDF5 file path."""
    normalized = str(pathroot_or_hdf5)
    if normalized.endswith(".hdf5.snap"):
        normalized = normalized[: -len(".snap")]

    root, ext = os.path.splitext(normalized)
    if ext.lower() == ".hdf5":
        root_path = root
    else:
        root_path = normalized

    root_path = re.sub(r"\.\d+$", "", root_path)
    return f"{root_path}.schema.json"
