#!/usr/bin/env python3
"""Diagnostic code -> Jarvis man target index (DESIGN_YAML_MAN_2.0 D24.10).

Every JV2-* code used in the runtime should resolve to a man page so agents can
close the validate -> man -> edit -> validate loop without reading docs/.
"""

from __future__ import annotations

from typing import Any

# kind: sampler | yaml | calculator | operas | tokens | center
# method: only for sampler pages
# path: YAML path for yaml kind (optional)
# tokens: optional list for resolve_man_request

_CODE_INDEX: dict[str, dict[str, Any]] = {
    # Schema
    "JV2-SCH-001": {"kind": "yaml", "path": None},
    "JV2-SCH-002": {"kind": "calculator", "topic": "execution"},
    "JV2-SCH-003": {"kind": "yaml", "path": None},
    "JV2-SCH-004": {"kind": "yaml", "path": None},
    # Method / Bounds
    "JV2-MTH-001": {"kind": "yaml", "path": "Sampling.Bounds"},
    "JV2-MTH-002": {"kind": "sampler"},
    "JV2-MTH-003": {"kind": "sampler"},
    "JV2-MTH-010": {"kind": "sampler", "method": "Bridson"},
    "JV2-MTH-011": {"kind": "sampler", "method": "Bridson"},
    "JV2-MTH-012": {"kind": "sampler", "method": "Bridson"},
    "JV2-MTH-013": {"kind": "sampler", "method": "Bridson"},
    "JV2-MTH-020": {"kind": "sampler", "method": "Random"},
    "JV2-MTH-021": {"kind": "sampler", "method": "Random"},
    "JV2-MTH-030": {"kind": "sampler", "method": "CSV"},
    "JV2-MTH-031": {"kind": "sampler", "method": "CSV"},
    "JV2-MTH-040": {"kind": "sampler"},
    "JV2-MTH-041": {"kind": "sampler", "method": "AdaptiveBridson"},
    "JV2-MTH-042": {"kind": "sampler", "method": "AdaptiveBridson"},
    "JV2-MTH-050": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-051": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-052": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-053": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-054": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-055": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-MTH-056": {"kind": "sampler", "method": "ToyMCMC"},
    "JV2-BND-000": {"kind": "yaml", "path": "Sampling.Bounds"},
    "JV2-BND-001": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-010": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-012": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-013": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-020": {"kind": "sampler", "method": "MultiNest"},
    "JV2-BND-021": {"kind": "sampler", "method": "MultiNest"},
    "JV2-BND-030": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-031": {"kind": "sampler", "method": "Dynesty"},
    "JV2-BND-032": {"kind": "sampler", "method": "Dynesty"},
    # Variables
    "JV2-VAR-001": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-002": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-010": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-011": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-012": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-020": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-021": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-030": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-031": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-032": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-033": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-034": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-035": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-036": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-037": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-038": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-039": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-040": {"kind": "yaml", "path": "Sampling.Variables"},
    "JV2-VAR-041": {"kind": "yaml", "path": "Sampling.Variables"},
    # Mapper
    "JV2-MAP-001": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-002": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-003": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-004": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-005": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-006": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-007": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-010": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-050": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-051": {"kind": "yaml", "path": "Sampling"},
    "JV2-MAP-052": {"kind": "yaml", "path": "Sampling"},
    # Operas
    "JV2-OPR-001": {"kind": "operas"},
    "JV2-OPR-002": {"kind": "operas"},
    # Calculator modes / pools
    "JV2-MOD-001": {"kind": "calculator", "topic": "module"},
    "JV2-MOD-002": {"kind": "calculator", "topic": "execution"},
    "JV2-MOD-003": {"kind": "calculator", "topic": "modes"},
    "JV2-MOD-004": {"kind": "calculator", "topic": "execution"},
    "JV2-MOD-005": {"kind": "calculator", "topic": "module"},
    "JV2-MOD-006": {"kind": "calculator", "topic": "modes"},
    "JV2-MOD-007": {"kind": "calculator", "topic": "pools"},
    "JV2-MOD-008": {"kind": "calculator", "topic": "pools"},
    "JV2-MOD-009": {"kind": "calculator", "topic": "pools"},
    # Archiver / Env (operational)
    "JV2-ARC-001": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-002": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-010": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-011": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-012": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-013": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-014": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ARC-015": {"kind": "yaml", "path": "EnvReqs.V2.archiver"},
    "JV2-ENV-001": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-002": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-003": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-010": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-011": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-012": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-019": {"kind": "yaml", "path": "EnvReqs.V2.worker"},
    "JV2-ENV-020": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-021": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-022": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-030": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-031": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-032": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-033": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-034": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-035": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-036": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-037": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-040": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-041": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-042": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-043": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-050": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-051": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-052": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-053": {"kind": "yaml", "path": "EnvReqs"},
    "JV2-ENV-060": {"kind": "yaml", "path": "EnvReqs.V2.check_modules"},
    "JV2-ENV-061": {"kind": "yaml", "path": "EnvReqs.V2.check_modules"},
    "JV2-ENV-062": {"kind": "yaml", "path": "EnvReqs.V2.check_modules"},
    "JV2-ENV-063": {"kind": "yaml", "path": "EnvReqs.V2.check_modules"},
    "JV2-ENV-064": {"kind": "yaml", "path": "EnvReqs.V2.check_modules"},
    # Load / encoding / yaml parse
    "JV2-LOAD-001": {"kind": "center"},
    "JV2-LOAD-002": {"kind": "center"},
    "JV2-ENC-001": {"kind": "center"},
    "JV2-YAML-001": {"kind": "center"},
    "JV2-YAML-002": {"kind": "center"},
}


def known_diagnostic_codes() -> frozenset[str]:
    return frozenset(_CODE_INDEX)


def man_target_for_code(code: str) -> dict[str, Any] | None:
    """Return the index entry for *code*, or None if unregistered."""
    return _CODE_INDEX.get(str(code or "").strip().upper())


def man_command_for(code: str, path: str | None = None) -> str:
    """Return an executable ``Jarvis man ...`` hint for validate JSON."""
    entry = man_target_for_code(code)
    if entry is None:
        # Fall back to path-based yaml man when possible.
        p = str(path or "").strip()
        if p:
            if p.startswith("$."):
                p = p[2:]
            elif p.startswith("$"):
                p = p[1:].lstrip(".")
            return f"Run: Jarvis man yaml.{p}"
        return "Run: Jarvis man"
    kind = entry.get("kind")
    if kind == "sampler":
        method = entry.get("method")
        if method:
            return f"Run: Jarvis man sampler.{method}"
        return "Run: Jarvis man sampler"
    if kind == "calculator":
        topic = entry.get("topic")
        if topic:
            return f"Run: Jarvis man calculator.{topic}"
        return "Run: Jarvis man calculator"
    if kind == "operas":
        return "Run: Jarvis man operas"
    if kind == "tokens":
        return "Run: Jarvis man tokens"
    if kind == "yaml":
        p = entry.get("path") or path
        if p:
            p = str(p)
            if p.startswith("$."):
                p = p[2:]
            elif p.startswith("$"):
                p = p[1:].lstrip(".")
            return f"Run: Jarvis man yaml.{p}"
        if path:
            p = str(path)
            if p.startswith("$."):
                p = p[2:]
            elif p.startswith("$"):
                p = p[1:].lstrip(".")
            return f"Run: Jarvis man yaml.{p}"
        return "Run: Jarvis man yaml"
    return "Run: Jarvis man"


def resolve_tokens_for_code(code: str, path: str | None = None) -> list[str]:
    """Argv tokens for ``resolve_man_request`` / ``dispatch_man`` (dotted topic)."""
    entry = man_target_for_code(code)
    if entry is None:
        if path:
            p = str(path)
            if p.startswith("$."):
                p = p[2:]
            return [f"yaml.{p}"]
        return []
    kind = entry.get("kind")
    if kind == "sampler":
        method = entry.get("method")
        return [f"sampler.{method}"] if method else ["sampler"]
    if kind == "calculator":
        topic = entry.get("topic")
        return [f"calculator.{topic}"] if topic else ["calculator"]
    if kind == "operas":
        return ["operas"]
    if kind == "tokens":
        return ["tokens"]
    if kind == "yaml":
        p = entry.get("path") or path
        if not p:
            return ["yaml"]
        p = str(p)
        if p.startswith("$."):
            p = p[2:]
        return [f"yaml.{p}"]
    return []


__all__ = [
    "known_diagnostic_codes",
    "man_command_for",
    "man_target_for_code",
    "resolve_tokens_for_code",
]
