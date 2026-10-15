"""Export sampler sources as project assets; fetching never executes or loads them."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import yaml
from jarvis_hep_sampler import SamplerRegistry
from jarvis_hep_sampler.metadata import parse_metadata

from jarvishep2.sampler_catalog import (
    SOURCE_USER,
    registered_plugin_record,
    reserved_method_names,
)
from jarvishep2.sampling.plugin_identity import (
    environment_versions,
    read_binding,
    source_path,
)

INDEX = ".jarvis-pack/samplers.json"


def sampler_assets(root: str, included: list[str], profile: str) -> dict[str, bytes]:
    """Deduplicate by full hash, including multiple revisions of one method."""
    project = Path(root).resolve()
    assets: dict[str, bytes] = {}
    entries: dict[tuple[str, str], dict[str, Any]] = {}
    registry = SamplerRegistry()

    def add(meta, path, digest, environment, used_by):
        if meta["method"] in reserved_method_names():
            return
        payload = Path(path).read_bytes()
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError(f"sampler source hash mismatch: {path}")
        method = parse_metadata(meta).method
        rel = f".jarvis-pack/samplers/{method}-{digest}.py"
        assets[rel] = payload
        key = (method, digest)
        entry = entries.setdefault(
            key,
            {
                "method": method,
                "source_sha256": digest,
                "file": rel,
                "metadata": dict(meta),
                "environment": environment,
                "used_by": [],
            },
        )
        if used_by not in entry["used_by"]:
            entry["used_by"].append(used_by)

    # Existing scans always keep their binding, including when the registry changed.
    if profile in {"repro", "full"}:
        for rel in included:
            if not rel.endswith("/sampler_identity.json"):
                continue
            raw = json.loads((project / rel).read_text())
            cfg = {
                "task_root": root,
                "task_result_dir": str((project / rel).parent.parent),
                "Sampling": {"Method": raw["metadata"]["method"]},
            }
            data = read_binding(cfg)
            add(
                data["metadata"],
                source_path(cfg, data),
                data["source_sha256"],
                data["environment"],
                rel,
            )

    for rel in included:
        if not rel.endswith((".yaml", ".yml")):
            continue
        raw = yaml.safe_load((project / rel).read_text(encoding="utf-8"))
        if not isinstance(raw, dict) or not isinstance(raw.get("Sampling"), dict):
            continue
        method = str(raw["Sampling"].get("Method") or "")
        if not method or method in reserved_method_names():
            continue
        scan = str((raw.get("Scan") or {}).get("name") or "default")
        cfg = {
            "task_root": root,
            "task_result_dir": str(project / "outputs" / scan),
            "scan_name": scan,
            "Sampling": raw["Sampling"],
        }
        bound = read_binding(cfg) if profile in {"repro", "full"} else None
        if bound:
            add(
                bound["metadata"],
                source_path(cfg, bound),
                bound["source_sha256"],
                bound["environment"],
                rel,
            )
            continue
        record = registered_plugin_record(method)
        if record is None or record.source != SOURCE_USER:
            raise ValueError(
                f"cannot package sampler {method!r} used by {rel}; load its source first"
            )
        entry = registry.get(method)
        add(
            entry.plugin,
            registry.snapshot_path(entry),
            entry.source_sha256,
            environment_versions(),
            rel,
        )
    if entries:
        index = {
            "format": "jarvis-hep-sampler-assets",
            "version": 1,
            "entries": [entries[key] for key in sorted(entries)],
        }
        assets[INDEX] = (json.dumps(index, sort_keys=True, indent=2) + "\n").encode()
    return assets


def inspect_packaged_samplers(root: str) -> list[str]:
    """Verify fetched source bytes and return actionable registration instructions."""
    project = Path(root).resolve()
    index = project / INDEX
    if not index.exists():
        return []
    data = json.loads(index.read_text())
    if data.get("format") != "jarvis-hep-sampler-assets" or data.get("version") != 1:
        raise ValueError("unsupported packaged sampler index")
    notices = []
    registry = SamplerRegistry()
    for entry in data["entries"]:
        meta = parse_metadata(entry["metadata"])
        path = (project / entry["file"]).resolve()
        if not path.is_relative_to(project / ".jarvis-pack" / "samplers"):
            raise ValueError("packaged sampler path escapes its asset directory")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != entry["source_sha256"]:
            raise ValueError(f"packaged sampler source hash mismatch: {meta.method}")
        current = registry.get(meta.method)
        if (
            current is not None
            and current.source_sha256 == digest
            and not registry.problems(current)
        ):
            notices.append(
                f"{meta.method} {meta.version} ({digest[:12]}): identical source already loaded"
            )
        else:
            reason = "not loaded" if current is None else "different source is loaded"
            import shlex

            notices.append(
                f"{meta.method} {meta.version} ({digest[:12]}): {reason}; "
                f"Jarvis sampler load {shlex.quote(str(path))}"
            )
    if notices:
        notices.append(
            "Install the sampler's imported packages/resources in the same environment; "
            "the archive includes entry source, not its dependency environment."
        )
    return notices
