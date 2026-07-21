#!/usr/bin/env python3
"""Projected flat hep:feedback records (D13.8).

Wire default::

    {"uuid": "...", "logL": -12.34}

Unusable points carry ``logL = -inf`` (set by likelihood / projection).
Extra optimizer fields are top-level siblings, never a nested observables bag.
Sample lifecycle ``status`` stays on the archive path only.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

# Wire key for total log-likelihood (maps from sample observables["LogL"]).
WIRE_LOGL_KEY = "logL"

_INF_STRINGS = {
    "-inf": float("-inf"),
    "-infinity": float("-inf"),
    "+inf": float("inf"),
    "+infinity": float("inf"),
    "inf": float("inf"),
    "infinity": float("inf"),
    "nan": float("nan"),
}


def encode_feedback_float(value: float) -> float | str:
    """JSON-friendly encoding for non-finite floats (msgpack may keep native floats)."""
    x = float(value)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "-inf" if x < 0 else "+inf"
    return x


def decode_feedback_float(value: Any, *, default: float = float("-inf")) -> float:
    """Decode a feedback scalar to float (−∞ default when unusable)."""
    if value is None:
        return float(default)
    if isinstance(value, (float, int)) and not isinstance(value, bool):
        x = float(value)
        if math.isnan(x):
            return float(default)
        return x
    if isinstance(value, str):
        text = value.strip().lower().replace("∞", "inf")
        if text in _INF_STRINGS:
            x = _INF_STRINGS[text]
            if math.isnan(x):
                return float(default)
            return x
        try:
            x = float(text)
        except ValueError:
            return float(default)
        if math.isnan(x):
            return float(default)
        return x
    try:
        x = float(value)
    except (TypeError, ValueError):
        return float(default)
    if math.isnan(x):
        return float(default)
    return x


def extract_logl_total(observables: Mapping[str, Any] | None) -> float | None:
    """Read total LogL from sample observables (LogL or sum of LogL* terms)."""
    if not isinstance(observables, Mapping):
        return None
    if "LogL" in observables:
        try:
            x = float(observables["LogL"])
        except (TypeError, ValueError):
            return None
        if math.isnan(x):
            return None
        return x
    terms: list[float] = []
    for key, value in observables.items():
        name = str(key)
        if name == "LogL" or not name.startswith("LogL"):
            continue
        try:
            terms.append(float(value))
        except (TypeError, ValueError):
            continue
    if not terms:
        return None
    total = float(sum(terms))
    if math.isnan(total):
        return None
    return total


def _is_feedback_scalar(value: Any) -> bool:
    if isinstance(value, bool):
        return True
    if isinstance(value, (int, float)):
        return True
    if isinstance(value, str):
        return len(value) <= 512
    return False


def build_feedback_record(
    sample: Any,
    *,
    spec: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build a flat feedback dict for ``hep:feedback``."""
    policy = dict(spec or {})
    mode = str(policy.get("mode") or "minimal").strip().lower()
    if mode == "all":
        mode = "all_flat"
    include_logl = bool(policy.get("include_logl", True))

    uuid = str(getattr(sample, "uuid", "") or "")
    if not uuid and isinstance(sample, Mapping):
        uuid = str(sample.get("uuid") or "")
    out: dict[str, Any] = {"uuid": uuid}

    obs: Mapping[str, Any]
    params: Mapping[str, Any]
    if isinstance(sample, Mapping):
        obs = dict(sample.get("observables") or {})
        params = dict(sample.get("params") or {})
    else:
        obs = dict(getattr(sample, "observables", None) or {})
        params = dict(getattr(sample, "params", None) or {})
    bag = {**params, **obs}

    if include_logl:
        logl = extract_logl_total(obs)
        out[WIRE_LOGL_KEY] = encode_feedback_float(
            float(logl) if logl is not None else float("-inf")
        )

    if mode == "minimal":
        return out

    if mode == "fields":
        names = [str(n) for n in (policy.get("fields") or [])]
    elif mode == "all_flat":
        names = [
            str(k)
            for k in bag.keys()
            if str(k) not in {"uuid", WIRE_LOGL_KEY, "LogL", "status"}
        ]
    else:
        raise ValueError(
            f"unknown feedback_return mode: {mode!r} "
            "(expected minimal|fields|all_flat)"
        )

    for name in names:
        if name in {"uuid", WIRE_LOGL_KEY, "status", "observables"}:
            continue
        if name not in bag:
            continue
        value = bag[name]
        if not _is_feedback_scalar(value):
            continue
        key = WIRE_LOGL_KEY if name == "LogL" else name
        if key == WIRE_LOGL_KEY and WIRE_LOGL_KEY in out:
            continue
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            out[key] = encode_feedback_float(value)
        else:
            out[key] = value
    return out


def normalize_feedback_record(record: Mapping[str, Any] | None) -> dict[str, Any]:
    """Normalize a pulled feedback record to flat floats (legacy nested tolerated)."""
    if not isinstance(record, Mapping):
        return {}
    raw = dict(record)
    out: dict[str, Any] = {"uuid": str(raw.get("uuid") or "")}

    # Legacy nested bag: promote LogL / scalars if flat logL absent.
    nested = raw.get("observables")
    if isinstance(nested, Mapping):
        if WIRE_LOGL_KEY not in raw and "LogL" not in raw:
            logl = extract_logl_total(nested)
            if logl is not None:
                out[WIRE_LOGL_KEY] = float(logl)
        for key, value in nested.items():
            name = str(key)
            if name in {"uuid", "status", "observables"}:
                continue
            wire = WIRE_LOGL_KEY if name == "LogL" else name
            if wire not in raw and wire not in out and _is_feedback_scalar(value):
                out[wire] = value

    for key, value in raw.items():
        if key in {"status", "observables"}:
            continue
        if key == "uuid":
            continue
        if key == "LogL":
            out.setdefault(WIRE_LOGL_KEY, value)
            continue
        out[key] = value

    if WIRE_LOGL_KEY in out:
        out[WIRE_LOGL_KEY] = decode_feedback_float(out[WIRE_LOGL_KEY])
    return out


def feedback_logl(record: Mapping[str, Any] | None) -> float:
    """Return logL from a (possibly legacy) feedback record; default −∞."""
    norm = normalize_feedback_record(record)
    if WIRE_LOGL_KEY not in norm:
        return float("-inf")
    return decode_feedback_float(norm[WIRE_LOGL_KEY])


def is_unusable_logl(logl: float) -> bool:
    """True when logL is non-finite (nan or ±inf) — treat as failed proposal."""
    x = float(logl)
    return math.isnan(x) or math.isinf(x)


def resolve_feedback_return(
    config: Mapping[str, Any] | None,
    *,
    sampler: Any = None,
    method: str | None = None,
) -> dict[str, Any]:
    """Resolve worker_config['feedback_return'] from YAML / sampler / defaults."""
    cfg = dict(config or {})
    sampling = cfg.get("Sampling") if isinstance(cfg.get("Sampling"), Mapping) else {}
    sampling = dict(sampling or {})

    explicit = (
        sampling.get("FeedbackReturn")
        or sampling.get("feedback_return")
        or cfg.get("FeedbackReturn")
        or cfg.get("feedback_return")
    )
    if isinstance(explicit, Mapping) and explicit:
        return _normalize_spec(dict(explicit))

    if sampler is not None and hasattr(sampler, "feedback_return_spec"):
        try:
            spec = sampler.feedback_return_spec()
        except Exception:
            spec = None
        if isinstance(spec, Mapping) and spec:
            return _normalize_spec(dict(spec))

    method_name = str(
        method
        or getattr(sampler, "method", "")
        or sampling.get("Method")
        or sampling.get("method")
        or ""
    ).strip()
    return _normalize_spec(_default_spec_for_method(method_name, sampling=sampling, cfg=cfg))


def _normalize_spec(raw: Mapping[str, Any]) -> dict[str, Any]:
    mode = str(raw.get("mode") or "minimal").strip().lower()
    if mode == "all":
        mode = "all_flat"
    if mode not in {"minimal", "fields", "all_flat"}:
        raise ValueError(
            f"FeedbackReturn.mode must be minimal|fields|all_flat, got {mode!r}"
        )
    fields = [str(x) for x in (raw.get("fields") or [])]
    include_logl = bool(raw.get("include_logl", True))
    if mode == "fields" and not fields and not include_logl:
        raise ValueError(
            "FeedbackReturn mode=fields requires non-empty fields "
            "or include_logl: true"
        )
    return {
        "mode": mode,
        "fields": fields,
        "include_logl": include_logl,
    }


def _default_spec_for_method(
    method: str,
    *,
    sampling: Mapping[str, Any],
    cfg: Mapping[str, Any],
) -> dict[str, Any]:
    name = str(method or "").strip()
    lower = name.lower()
    # AdaptiveBridson: request target expression symbols as flat fields.
    if lower in {"adaptivebridson", "adaptive_bridson"}:
        fields = _target_expression_fields(sampling, cfg)
        return {
            "mode": "fields" if fields else "minimal",
            "fields": fields,
            "include_logl": True,
        }
    # Everything feedback-driven defaults to minimal uuid+logL.
    return {"mode": "minimal", "fields": [], "include_logl": True}


_SYMBOL_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\b")


def _target_expression_fields(
    sampling: Mapping[str, Any],
    cfg: Mapping[str, Any],
) -> list[str]:
    block = (
        sampling.get("AdaptiveBridson")
        or sampling.get("adaptive_bridson")
        or cfg.get("AdaptiveBridson")
        or {}
    )
    if not isinstance(block, Mapping):
        return []
    expr = str(
        block.get("target_expression")
        or block.get("TargetExpression")
        or block.get("target")
        or ""
    ).strip()
    if not expr:
        return []
    # Lightweight free-symbol scrape (no sympy dependency here).
    reserved = {
        "and",
        "or",
        "not",
        "if",
        "else",
        "True",
        "False",
        "None",
        "abs",
        "min",
        "max",
        "sin",
        "cos",
        "tan",
        "exp",
        "log",
        "sqrt",
        "pi",
        "e",
    }
    names: list[str] = []
    seen: set[str] = set()
    for match in _SYMBOL_RE.finditer(expr):
        token = match.group(1)
        if token in reserved or token in seen:
            continue
        seen.add(token)
        # Wire uses logL; request LogL from bag maps to logL in builder.
        names.append(token)
    return names


__all__ = [
    "WIRE_LOGL_KEY",
    "build_feedback_record",
    "decode_feedback_float",
    "encode_feedback_float",
    "extract_logl_total",
    "feedback_logl",
    "is_unusable_logl",
    "normalize_feedback_record",
    "resolve_feedback_return",
]
