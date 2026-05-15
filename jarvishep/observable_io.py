from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from typing import Any

import numpy as np
import sympy

SCHEMA_VERSION_CURRENT = 3
SCHEMA_VERSION_MIN_COMPATIBLE = 3
CSV_EXPORT_POLICY_DEFAULTS = {
    "scalar": "keep",
    "path_like_string": "keep",
    "small_list": "split",
    "small_flat_dict": "split",
    "large_list": "drop",
    "large_dict": "drop",
    "nested_dict": "drop",
    "mixed": "drop",
    "unsupported": "drop",
}
CSV_EXPORT_LIMIT_DEFAULTS = {
    "max_list_len": 20,
    "max_dict_keys": 20,
    "max_nested_depth": 1,
}
PATH_LIKE_FIELD_NAMES = {
    "gmfit_output",
    "spheno_input",
    "spheno_output",
    "rge_output",
}

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
        desc["length"] = len(seq)
        desc["shape"] = _infer_list_shape(seq)
        desc["ndim"] = len(desc["shape"])
        element_kind = _infer_list_element_kind(seq)
        if element_kind:
            desc["element_kind"] = element_kind
        return desc

    if isinstance(value, dict):
        desc["keys_preview"] = [str(k) for k in list(value.keys())[:20]]
        desc["num_keys"] = len(value)
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


def create_default_schema(pathroot: str) -> dict[str, Any]:
    now = utc_now_iso()
    return {
        "version": SCHEMA_VERSION_CURRENT,
        "pathroot": pathroot,
        "created_at": now,
        "updated_at": now,
        "csv_export_policy": dict(CSV_EXPORT_POLICY_DEFAULTS),
        "csv_export_limits": dict(CSV_EXPORT_LIMIT_DEFAULTS),
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
        msg = f"Schema file not found; creating default schema -> {path}"
        warnings.append(msg)
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
                msg = "Schema had invalid/missing fields; normalized with fallback defaults."
                warnings.append(msg)
            normalized["_warnings"] = warnings
            return normalized
        except Exception as exc:
            msg = f"Failed to read schema ({exc}); fallback to default schema -> {path}"
            warnings.append(msg)
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
    """Normalize schema metadata and simplify obsolete export blocks."""
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

    if "flatten_defaults" in normalized:
        normalized.pop("flatten_defaults", None)
        changed = True
        warnings.append("schema.flatten_defaults is obsolete; removed.")

    policy = normalized.get("csv_export_policy")
    if not isinstance(policy, dict):
        normalized["csv_export_policy"] = dict(CSV_EXPORT_POLICY_DEFAULTS)
        changed = True
        warnings.append("schema.csv_export_policy is invalid or missing; replaced with defaults.")
    else:
        normalized_policy = dict(CSV_EXPORT_POLICY_DEFAULTS)
        for key, value in policy.items():
            if key in CSV_EXPORT_POLICY_DEFAULTS and value in {"keep", "split", "drop"}:
                normalized_policy[key] = value
        if normalized_policy != policy:
            normalized["csv_export_policy"] = normalized_policy
            changed = True
            warnings.append("schema.csv_export_policy had invalid entries; normalized.")

    limits = normalized.get("csv_export_limits")
    if not isinstance(limits, dict):
        normalized["csv_export_limits"] = dict(CSV_EXPORT_LIMIT_DEFAULTS)
        changed = True
        warnings.append("schema.csv_export_limits is invalid or missing; replaced with defaults.")
    else:
        normalized_limits = dict(CSV_EXPORT_LIMIT_DEFAULTS)
        for key, default_value in CSV_EXPORT_LIMIT_DEFAULTS.items():
            try:
                value = int(limits.get(key, default_value))
            except Exception:
                value = default_value
            normalized_limits[key] = max(0, value)
        if normalized_limits != limits:
            normalized["csv_export_limits"] = normalized_limits
            changed = True
            warnings.append("schema.csv_export_limits had invalid entries; normalized.")

    columns = normalized.get("columns")
    if not isinstance(columns, dict):
        normalized["columns"] = {}
        changed = True
        warnings.append("schema.columns is invalid; reset to empty object.")
        return normalized, changed, warnings

    for name, meta in list(columns.items()):
        if not isinstance(meta, dict):
            columns[name] = {"kind": "unknown"}
            changed = True
            warnings.append(f"schema.columns.{name} is not an object; reset with descriptor metadata only.")
            continue

        had_obsolete = False
        for obsolete_key in ("flatten", "csv_export"):
            if obsolete_key in meta:
                meta.pop(obsolete_key, None)
                had_obsolete = True
        if had_obsolete:
            changed = True
            warnings.append(f"schema.columns.{name} had obsolete export metadata; simplified.")

        if "keep" not in meta or "name" not in meta:
            meta.update(_descriptor_keep_config(name, meta, normalized))
            changed = True
            warnings.append(f"schema.columns.{name} missing keep/name; generated from current policy.")

    return normalized, changed, warnings


def _normalize_schema_version(raw_version: Any) -> tuple[int, bool, str | None]:
    if raw_version is None:
        return (
            SCHEMA_VERSION_CURRENT,
            True,
            f"schema.version missing; set to {SCHEMA_VERSION_CURRENT}.",
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
                f"schema.version {parsed} is below supported "
                f"{SCHEMA_VERSION_MIN_COMPATIBLE}; set to {SCHEMA_VERSION_CURRENT}."
            ),
        )

    if parsed > SCHEMA_VERSION_CURRENT:
        return (
            parsed,
            changed,
            (
                f"schema.version {parsed} is newer than supported {SCHEMA_VERSION_CURRENT}; "
                "normalizing known container metadata only."
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

    if "num_keys" in current:
        new_num_keys = current["num_keys"]
        old_num_keys = merged.get("num_keys")
        if old_num_keys is None or old_num_keys != new_num_keys:
            merged["num_keys"] = new_num_keys
            changed = True

    if "length" in current:
        new_length = current["length"]
        old_length = merged.get("length")
        if old_length is None:
            merged["length"] = new_length
            changed = True
        elif old_length != new_length and old_length != "variable":
            merged["length"] = "variable"
            changed = True

    for obsolete_key in ("flatten", "csv_export"):
        if obsolete_key in merged:
            merged.pop(obsolete_key, None)
            changed = True

    return merged, changed


def _csv_export_policy(schema: dict[str, Any]) -> dict[str, str]:
    policy = dict(CSV_EXPORT_POLICY_DEFAULTS)
    raw = schema.get("csv_export_policy", {}) if isinstance(schema, dict) else {}
    if isinstance(raw, dict):
        for key, value in raw.items():
            if key in policy and value in {"keep", "split", "drop"}:
                policy[key] = value
    return policy


def _csv_export_limits(schema: dict[str, Any]) -> dict[str, int]:
    limits = dict(CSV_EXPORT_LIMIT_DEFAULTS)
    raw = schema.get("csv_export_limits", {}) if isinstance(schema, dict) else {}
    if isinstance(raw, dict):
        for key, default in CSV_EXPORT_LIMIT_DEFAULTS.items():
            try:
                limits[key] = max(0, int(raw.get(key, default)))
            except Exception:
                limits[key] = default
    return limits


def _column_suffix(value: Any) -> str:
    suffix = re.sub(r"[^0-9A-Za-z_]+", "_", str(value)).strip("_")
    return suffix or "value"


def _scalar_type_group(value: Any) -> str:
    if isinstance(value, (np.integer, np.floating, sympy.Float, sympy.Integer)):
        return "number"
    if isinstance(value, bool) or _is_sympy_boolean(value):
        return "bool"
    if isinstance(value, (int, float)):
        return "number"
    if isinstance(value, str):
        return "str"
    if value is None:
        return "null"
    return _infer_kind(value)


def _is_flat_sequence(values: list[Any]) -> bool:
    return all(not isinstance(item, (dict, list, tuple, set, np.ndarray)) for item in values)


def _is_mixed_flat_values(values: list[Any]) -> bool:
    groups = {_scalar_type_group(item) for item in values if item is not None}
    return len(groups) > 1


def _dict_depth(value: Any) -> int:
    if not isinstance(value, dict) or not value:
        return 0
    child_depths = []
    for item in value.values():
        if isinstance(item, dict):
            child_depths.append(1 + _dict_depth(item))
        elif isinstance(item, (list, tuple, set, np.ndarray)):
            child_depths.append(2)
        else:
            child_depths.append(1)
    return max(child_depths, default=0)


def _is_flat_dict(value: dict[Any, Any]) -> bool:
    return all(not isinstance(item, (dict, list, tuple, set, np.ndarray)) for item in value.values())


def _is_path_like_field(name: str, value: Any) -> bool:
    if name in PATH_LIKE_FIELD_NAMES:
        return True
    if not isinstance(value, str):
        return False
    lowered = name.lower()
    return lowered.endswith("_path") or lowered.endswith("_file")


def _policy_decision(policy: dict[str, str], reason_key: str) -> str:
    return policy.get(reason_key, CSV_EXPORT_POLICY_DEFAULTS.get(reason_key, "drop"))


def _value_dtype(value: Any) -> str:
    if isinstance(value, np.ndarray):
        return str(value.dtype)
    if isinstance(value, np.generic):
        return str(value.dtype)
    if isinstance(value, (sympy.Float, sympy.Integer)):
        return type(value).__name__
    if _is_sympy_boolean(value):
        return "bool"
    if value is None:
        return "null"
    return type(value).__name__


def _field_keep_config(keep: bool, name: Any, value: Any, kind: str | None = None) -> dict[str, Any]:
    config: dict[str, Any] = {
        "keep": bool(keep),
        "name": name,
    }
    if kind in {"dict", "list", "ndarray", "mixed"}:
        config["kind"] = kind
    else:
        config["dtype"] = _value_dtype(value)
    return config


def _subfield_config(column_name: str, value: Any, keep: bool = True) -> dict[str, Any]:
    return {
        "name": column_name,
        "keep": bool(keep),
        "dtype": _value_dtype(value),
    }


def _default_keep_config(field: str, desc: dict[str, Any], value: Any, schema: dict[str, Any]) -> dict[str, Any]:
    policy = _csv_export_policy(schema)
    limits = _csv_export_limits(schema)
    kind = str(desc.get("kind", "unknown"))

    if kind == "str" and _is_path_like_field(field, value):
        decision = _policy_decision(policy, "path_like_string")
        return _field_keep_config(decision == "keep", field, value)

    if kind in {"bool", "int", "float", "str", "null", "numpy_scalar", "sympy_scalar"}:
        decision = _policy_decision(policy, "scalar")
        return _field_keep_config(decision == "keep", field, value)

    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            return _field_keep_config(False, field, value, kind="ndarray")
        length = int(value.shape[0])
        desc["length"] = length
        keep = 0 < length <= limits["max_list_len"] and _policy_decision(policy, "small_list") == "split"
        values = value.tolist()
        names = {
            str(idx): _subfield_config(f"{field}_{idx}", values[idx])
            for idx in range(length)
        }
        return _field_keep_config(keep, names, value, kind="ndarray")

    if isinstance(value, (list, tuple, set)):
        seq = list(value)
        desc["length"] = len(seq)
        keep = (
            0 < len(seq) <= limits["max_list_len"]
            and _is_flat_sequence(seq)
            and not _is_mixed_flat_values(seq)
            and _policy_decision(policy, "small_list") == "split"
        )
        names = {
            str(idx): _subfield_config(f"{field}_{idx}", item)
            for idx, item in enumerate(seq)
        }
        return _field_keep_config(keep, names, value, kind="list")

    if isinstance(value, dict):
        num_keys = len(value)
        desc["num_keys"] = num_keys
        keep = (
            0 < num_keys <= limits["max_dict_keys"]
            and _dict_depth(value) <= limits["max_nested_depth"]
            and _is_flat_dict(value)
            and not _is_mixed_flat_values(list(value.values()))
            and _policy_decision(policy, "small_flat_dict") == "split"
        )
        names = {
            str(key): _subfield_config(f"{field}_{_column_suffix(key)}", item)
            for key, item in value.items()
        }
        return _field_keep_config(keep, names, value, kind="dict")

    if kind == "mixed":
        return _field_keep_config(False, field, value, kind="mixed")

    return _field_keep_config(False, field, value)


def _descriptor_keep_config(field: str, meta: dict[str, Any], schema: dict[str, Any]) -> dict[str, Any]:
    policy = _csv_export_policy(schema)
    limits = _csv_export_limits(schema)
    kind = str(meta.get("kind", "unknown"))
    dtype = str(meta.get("dtype", kind if kind != "unknown" else "object"))

    if kind == "str" and _is_path_like_field(field, ""):
        return {
            "dtype": "str",
            "keep": _policy_decision(policy, "path_like_string") == "keep",
            "name": field,
        }

    if kind in {"bool", "int", "float", "str", "null", "numpy_scalar", "sympy_scalar", "unknown"}:
        return {
            "dtype": dtype,
            "keep": _policy_decision(policy, "scalar") == "keep",
            "name": field,
        }

    if kind in {"ndarray", "list"}:
        length = meta.get("length")
        if not isinstance(length, int):
            shape = meta.get("shape")
            length = shape[0] if isinstance(shape, list) and shape and isinstance(shape[0], int) else 0
        keep = 0 < length <= limits["max_list_len"] and _policy_decision(policy, "small_list") == "split"
        element_dtype = str(meta.get("element_kind") or meta.get("dtype") or "object")
        names = {
            str(idx): {
                "name": f"{field}_{idx}",
                "keep": True,
                "dtype": element_dtype,
            }
            for idx in range(max(0, int(length)))
        }
        return {
            "kind": kind,
            "keep": keep,
            "name": names,
        }

    if kind == "dict":
        num_keys = meta.get("num_keys")
        keys = meta.get("keys_preview", [])
        if not isinstance(num_keys, int):
            num_keys = len(keys) if isinstance(keys, list) else 0
        keep = (
            0 < num_keys <= limits["max_dict_keys"]
            and isinstance(keys, list)
            and len(keys) >= num_keys
            and _policy_decision(policy, "small_flat_dict") == "split"
        )
        names = {
            str(key): {
                "name": f"{field}_{_column_suffix(key)}",
                "keep": True,
                "dtype": "object",
            }
            for key in (keys if isinstance(keys, list) else [])
        }
        return {
            "kind": "dict",
            "keep": keep,
            "name": names,
        }

    return {
        "dtype": dtype,
        "keep": False,
        "name": field,
    }


def _ensure_export_for_descriptor(
    name: str,
    descriptor: dict[str, Any],
    value: Any,
    schema: dict[str, Any],
) -> tuple[dict[str, Any], bool]:
    updated = dict(descriptor)
    changed = False
    for obsolete_key in ("flatten", "csv_export"):
        if obsolete_key in updated:
            updated.pop(obsolete_key, None)
            changed = True
    if "keep" not in updated or "name" not in updated:
        keep_config = _default_keep_config(name, updated, value, schema)
        updated.update(keep_config)
        changed = True
    if updated.get("kind") not in {"dict", "list", "ndarray", "mixed"} and "kind" in updated:
        updated.pop("kind", None)
        changed = True
    return updated, changed


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
            current, _ = _ensure_export_for_descriptor(key, current, value, schema)
            columns[key] = current
            changed = True
            continue

        merged, merged_changed = _merge_descriptor(existing, current)
        merged, export_changed = _ensure_export_for_descriptor(key, merged, value, schema)
        if merged_changed:
            columns[key] = merged
            changed = True
        elif export_changed:
            columns[key] = merged
            changed = True

    if changed:
        schema["updated_at"] = utc_now_iso()
    return changed


def _json_cell(value: Any) -> str:
    return json.dumps(make_json_compatible(value), ensure_ascii=False, separators=(",", ":"))


def validate_csv_export_schema(schema: dict[str, Any]) -> None:
    """Validate the simplified CSV export schema."""
    columns = schema.get("columns")
    if not isinstance(columns, dict):
        raise ValueError("schema.columns must be an object")

    for field, meta in columns.items():
        if not isinstance(meta, dict):
            raise ValueError(f"schema.columns.{field} must be an object")
        if not isinstance(meta.get("keep"), bool):
            raise ValueError(f"schema.columns.{field}.keep must be boolean")
        name = meta.get("name")
        if not isinstance(name, (str, dict)):
            raise ValueError(f"schema.columns.{field}.name must be a string or object")
        if isinstance(name, str) and meta["keep"] and not name:
            raise ValueError(f"schema.columns.{field}.name must be non-empty when keep is true")
        if isinstance(name, dict):
            for subkey, submeta in name.items():
                if not isinstance(submeta, dict):
                    raise ValueError(f"schema.columns.{field}.name.{subkey} must be an object")
                if not isinstance(submeta.get("keep"), bool):
                    raise ValueError(f"schema.columns.{field}.name.{subkey}.keep must be boolean")
                subname = submeta.get("name")
                if not isinstance(subname, str) or not subname:
                    raise ValueError(f"schema.columns.{field}.name.{subkey}.name must be a non-empty string")


def _field_export_config(schema: dict[str, Any], field: str) -> dict[str, Any]:
    columns = schema.get("columns", {}) if isinstance(schema, dict) else {}
    meta = columns.get(field) if isinstance(columns, dict) else None
    if not isinstance(meta, dict):
        raise ValueError(f"schema.columns.{field} must exist for CSV export")
    return meta


def _split_values_for_field(value: Any, name_map: dict[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {}
    if isinstance(value, np.ndarray):
        value = value.tolist()

    if isinstance(value, (list, tuple)):
        seq = list(value)
        for raw_idx, submeta in name_map.items():
            if not isinstance(submeta, dict) or not submeta.get("keep", False):
                continue
            try:
                idx = int(raw_idx)
            except Exception:
                continue
            if 0 <= idx < len(seq):
                row[submeta["name"]] = seq[idx]
        return row

    if isinstance(value, dict):
        for raw_key, submeta in name_map.items():
            if not isinstance(submeta, dict) or not submeta.get("keep", False):
                continue
            row[submeta["name"]] = value.get(str(raw_key), value.get(raw_key))
        return row

    return row


def _flatten_record_for_csv_internal(
    record: dict[str, Any],
    schema: dict[str, Any],
    populate_name_map: bool,
) -> tuple[dict[str, Any], bool]:
    """Apply per-column flatten policies to a record for CSV export."""
    if not isinstance(record, dict):
        return {}, False

    row: dict[str, Any] = {}
    normalized_record = make_json_compatible(record)

    validate_csv_export_schema(schema)
    columns = schema.get("columns", {})

    for key, value in normalized_record.items():
        field = str(key)
        meta = _field_export_config(schema, field)
        keep = meta["keep"]
        name = meta["name"]

        if not keep:
            continue

        if isinstance(name, str):
            if isinstance(value, (dict, list)):
                row[name] = _json_cell(value)
            else:
                row[name] = value
            continue

        if isinstance(name, dict):
            row.update(_split_values_for_field(value, name))

    return row, False


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


def csv_export_fieldnames_from_schema(schema: dict[str, Any]) -> list[str]:
    """Return default CSV columns in schema order."""
    validate_csv_export_schema(schema)
    fieldnames: list[str] = []
    seen: set[str] = set()
    columns = schema.get("columns", {})
    for meta in columns.values():
        if not isinstance(meta, dict) or not meta.get("keep", False):
            continue
        name = meta.get("name")
        columns_to_add: list[str] = []
        if isinstance(name, str):
            columns_to_add.append(name)
        elif isinstance(name, dict):
            columns_to_add.extend(
                submeta["name"]
                for submeta in name.values()
                if isinstance(submeta, dict) and submeta.get("keep", False)
            )
        for column in columns_to_add:
            if column not in seen:
                seen.add(column)
                fieldnames.append(column)
    return fieldnames


def format_csv_export_report(schema: dict[str, Any]) -> str:
    validate_csv_export_schema(schema)
    kept_scalar: list[str] = []
    split_fields: list[str] = []
    json_fields: list[str] = []
    dropped: list[str] = []
    columns = schema.get("columns", {})
    for field, meta in columns.items():
        if not meta.get("keep", False):
            dropped.append(field)
            continue
        name = meta.get("name")
        if isinstance(name, dict):
            split_fields.append(field)
        elif meta.get("kind") in {"dict", "list", "ndarray"}:
            json_fields.append(field)
        else:
            kept_scalar.append(field)

    def _summary(values: list[str]) -> str:
        if not values:
            return "none"
        preview = ", ".join(values[:12])
        extra = len(values) - 12
        return f"{preview} (+{extra} more)" if extra > 0 else preview

    lines = [
        "[CSV Export]",
        f"kept scalar fields: {_summary(kept_scalar)}",
        f"split fields: {_summary(split_fields)}",
        f"json fields: {_summary(json_fields)}",
        f"dropped fields: {len(dropped)}",
    ]
    if dropped:
        lines.append("dropped examples:")
        lines.extend(f"  - {field}" for field in dropped[:5])
    return "\n".join(lines)


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
