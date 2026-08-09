#!/usr/bin/env python3
"""Jarvis man: YAML writing manual for humans and coding agents (DESIGN_YAML_MAN_2.0).

Focus: help write a **valid task-card YAML** (paths, keys, examples that validate).
Runtime behaviour of Portal adapters and Operas operators is documented by their
own CLIs (``Jarvis portal …`` / ``Jarvis operas …``); man points there instead of
duplicating it.

Output is always English. Man only renders surfaces the validator knows about (K1).
"""

from __future__ import annotations

import argparse
import json
import sys
from difflib import get_close_matches
from pathlib import Path
from typing import Any, Mapping

from rich import box
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.text import Text

from jarvishep2.task_schema import (
    resolve_schema_ref,
    schema_by_id,
    schema_manifest,
)

# Absolute schema ids used as bases for fragment-only ``#/$defs/...`` refs.
_CALC_SCHEMA_ID = "https://jarvis-hep.org/schema/v2/core/calculators.json"
_OPERAS_SCHEMA_ID = "https://jarvis-hep.org/schema/v2/core/operas.json"
_COMMON_SCHEMA_ID = "https://jarvis-hep.org/schema/v2/core/common.json"

_CALC_MODULE_FIELD_HELP: dict[str, str] = {
    "name": "Module id (no dots). Used in required_modules and mode references.",
    "path": "Install / pack path for this calculator (tokens allowed, quote &J/...).",
    "source": "Optional source tree copied or linked into path.",
    "deps_source": "Optional dependency source path (V1 layout).",
    "clone_shadow": "When true, use a shadow clone of the pack for this sample.",
    "required_modules": "Other Calculators.Modules names that must run first.",
    "installation": "One-time install commands (string or {cmd,cwd}).",
    "initialization": "Per-run init commands before execution.",
    "timeout": "Wall-clock timeout seconds for this module.",
    "symlink_name": "Optional symlink name inside the sample workdir.",
    "env_setup": "Environment setup commands or maps before commands run.",
    "selection": "Optional expression gate; false skips this module.",
    "make_paraller": "Legacy concurrency hint (spelling preserved); prefer Pools.",
    "execution": "Run block: path, commands, input (list), output (list) (Portal IO).",
    "modes": "Named alternate execution profiles (multimode); exclusive with execution.",
}

_OPERAS_FIELD_HELP: dict[str, str] = {
    "name": "Module id for the execution plan / logs.",
    "operator": (
        "Registered Jarvis-Operas full name. Must return a mapping (dict), "
        "not a bare scalar. Inspect signatures with: Jarvis operas info NAME"
    ),
    "call_mode": "call (sync) or acall (async wrapper).",
    "selection": "Optional expression; false skips this Opera for the sample.",
    "required_modules": "Calculator/Opera names that must finish first.",
    "kwargs": "Static keyword args merged into the operator call.",
    "timeout": "Optional timeout seconds (alias: timeout_sec).",
    "timeout_sec": "Optional timeout seconds (alias: timeout).",
    "input": (
        "List of {name, expression} or {name, entry}. expression is sympy over "
        "params/observables; entry copies an existing observable by name."
    ),
    "output": (
        "List of {name, entry} bindings from the operator mapping into observables. "
        "Use entry (not expression) for output columns."
    ),
}

_EXECUTION_EXAMPLE = """\
execution:
  path: "&J/calculators/echo"
  commands:
    - python3 -c "import json, pathlib; pathlib.Path('@Sdir/out.json').write_text(json.dumps({'y': 1.0}))"
  input:
    - name: in
      path: input.json
      type: JSON
      actions:
        - type: Dump
          variables:
            - {name: x, expression: x}
  output:
    - name: out
      path: '@Sdir/out.json'
      type: JSON
      variables:
        - {name: y, entry: y}"""

_OPERAS_EXAMPLE = """\
Operas:
  Modules:
    - name: eggbox
      operator: helper.eggbox2d
      input:
        - {name: x, expression: x}
        - {name: y, expression: y}
      output:
        - {name: z, entry: z}"""

_ENVREQS_TOP_HELP: dict[str, str] = {
    "Check_default_dependencies": (
        "Project scaffold hook. When required is true, Jarvis loads and deep-merges "
        "EnvReqs.V2 from default_yaml_path before task overrides."
    ),
    "Python": (
        "Python requirement block. Usually inherited from the project default "
        "environment YAML; task values override it."
    ),
    "CERN_ROOT": (
        "CERN ROOT requirement block. Usually inherited from the project default "
        "environment YAML; set required true only when the task needs ROOT."
    ),
    "V2": (
        "User-facing V2 runtime overlay. Jarvis project supplies defaults; task "
        "values override only the settings that need changing."
    ),
    "Runtime": (
        "Legacy V1 default-runtime reference. Do not add a top-level Runtime block "
        "to a V2 task card."
    ),
}

_ENVREQS_V2_HELP: dict[str, str] = {
    "workers": "Worker count. Project default: 4; Jarvis check forces 1 at runtime.",
    "batch_size": "Runtime batch size. Project default: 256.",
    "sample_directory": (
        "Sample bucket layout: enabled, limit, width, pack, start_bucket. "
        "Configure it only under EnvReqs.V2; project defaults are inherited there."
    ),
    "cleanup": "Cleanup policy: strategy direct or mv_to_staging, plus staging_dir.",
    "archiver": (
        "Archiver policy: mode, batch_size, flush_interval_sec, max_hdf5_bytes, "
        "delete_after_archive, handoff, and pack_buckets."
    ),
    "redis": (
        "Per-project Redis host/port/db. Jarvis may select a free port and write it "
        "back to the project default YAML."
    ),
    "factory": "Factory monitor and watchdog policy.",
    "worker": "Worker policy: force_serial_layers and sample_artifacts.",
    "check_modules": (
        "Jarvis check settings: data, n_samples, timeout_sec. Check overrides "
        "workers and archive layout while it runs."
    ),
    "checkpoint_heartbeat_sec": (
        "Checkpoint heartbeat interval in seconds; minimum 1. Project runtime "
        "default: 30."
    ),
}

_ENVREQS_NESTED_PAGES: dict[str, tuple[str, dict[str, str], str]] = {
    "check_default_dependencies": (
        "Project default-environment loader. Jarvis project create writes "
        "deps/environment_default.yaml and task cards normally point to it here.",
        {
            "required": "Enable loading the project default environment YAML.",
            "default_yaml_path": "Path to the default YAML, normally &J/deps/environment_default.yaml.",
        },
        "EnvReqs:\n  Check_default_dependencies:\n    required: true\n    default_yaml_path: \"&J/deps/environment_default.yaml\"",
    ),
    "python": (
        "V1-compatible Python requirement block. Values are commonly inherited from the project default YAML.",
        {
            "version": "Minimum Python version, for example >=3.10.",
            "Dependencies": "Python package requirement entries (list).",
        },
        "Python:\n  version: \">=3.10\"",
    ),
    "cern_root": (
        "V1-compatible CERN ROOT requirement block. The default project policy leaves ROOT optional.",
        {
            "required": "Whether root-config and the requested ROOT features must be available.",
            "version": "Minimum ROOT version, for example >=6.20.",
            "path": "Optional ROOT prefix used by @{ROOT path}.",
            "get_path_command": "Command used to discover the ROOT prefix, normally root-config --prefix.",
            "Dependencies": "ROOT feature checks such as minuit2 (list).",
        },
        "CERN_ROOT:\n  required: false\n  version: \">=6.20\"",
    ),
    "runtime": (
        "Legacy V1 compatibility reference. V2 task cards should use EnvReqs.V2 instead of a top-level Runtime block.",
        {"default_runtime_settings": "Path to a legacy runtime-default YAML file."},
        "EnvReqs:\n  Runtime:\n    default_runtime_settings: \"&J/deps/runtime_default.yaml\"",
    ),
    "sample_directory": (
        "Sample bucket layout inherited from the project default V2 YAML.",
        {
            "enabled": "Enable SAMPLE bucket directories.",
            "limit": "Maximum samples per bucket.",
            "width": "Zero-padding width for bucket names.",
            "pack": "Pack sealed buckets into tar.gz archives.",
            "start_bucket": "First bucket number.",
        },
        "sample_directory:\n  enabled: true\n  limit: 200\n  width: 6\n  pack: true\n  start_bucket: 1",
    ),
    "cleanup": (
        "Cleanup policy inherited from EnvReqs.V2 in the project default YAML.",
        {
            "strategy": "direct or mv_to_staging.",
            "staging_dir": "Optional staging directory when using mv_to_staging.",
        },
        "cleanup:\n  strategy: direct",
    ),
    "archiver": (
        "Archive policy inherited from EnvReqs.V2 in the project default YAML.",
        {
            "mode": "thread or process.",
            "batch_size": "Number of samples handled per archive batch.",
            "flush_interval_sec": "Maximum archive flush interval in seconds.",
            "max_hdf5_bytes": "Seal an HDF5 shard after this many bytes.",
            "strategy": "move or copy for the archive handoff.",
            "delete_after_archive": "Delete source sample files after archiving.",
            "handoff": "direct or staging.",
            "pack_buckets": "Pack sealed SAMPLE buckets into tar.gz archives.",
        },
        "archiver:\n  mode: process\n  batch_size: 200\n  pack_buckets: true",
    ),
    "redis": (
        "Per-project Redis settings. Jarvis can rewrite port in the default environment YAML when 6379 is occupied.",
        {
            "host": "Redis host, normally 127.0.0.1.",
            "port": "Requested Redis TCP port, normally 6379.",
            "db": "Redis database number, normally 0.",
        },
        "redis:\n  host: 127.0.0.1\n  port: 6379\n  db: 0",
    ),
    "factory": (
        "Factory monitor and watchdog policy.",
        {
            "monitor_hz": "Factory monitor frequency; project default: 120.",
            "monitor": "Legacy monitor mapping; hz or update_hz supplies the frequency.",
            "watchdog": "Watchdog policy mapping.",
        },
        "factory:\n  monitor_hz: 120\n  watchdog:\n    enabled: true",
    ),
    "worker": (
        "Per-worker execution policy.",
        {
            "force_serial_layers": "Serialize dependency layers when true.",
            "sample_artifacts": "auto, always, or never.",
        },
        "worker:\n  force_serial_layers: false\n  sample_artifacts: auto",
    ),
    "check_modules": (
        "Fixed-point Calculator smoke settings. Jarvis check forces its own worker and archive policy while running.",
        {
            "data": "CSV path, normally &J/data/check_modules_points.csv.",
            "n_samples": "Number of smoke samples when data is not used.",
            "timeout_sec": "Maximum archive wait in seconds.",
        },
        "check_modules:\n  data: \"&J/data/check_modules_points.csv\"\n  n_samples: 10\n  timeout_sec: 120",
    ),
    "watchdog": (
        "Watchdog policy nested under EnvReqs.V2.factory.",
        {
            "enabled": "Enable stale-worker detection.",
            "stale_sec": "Seconds before a worker is considered stale.",
            "poll_interval_sec": "Watchdog polling interval in seconds.",
            "max_sample_retries": "Maximum retries for one sample.",
        },
        "watchdog:\n  enabled: true\n  stale_sec: 30",
    ),
}

_ENVREQS_FIELD_TYPES: dict[str, str] = {
    "required": "boolean",
    "default_yaml_path": "string",
    "default_runtime_settings": "string",
    "version": "string",
    "Dependencies": "list",
    "path": "string",
    "get_path_command": "string",
    "enabled": "boolean",
    "limit": "integer",
    "width": "integer",
    "start_bucket": "integer",
    "strategy": "enum",
    "staging_dir": "string",
    "mode": "enum",
    "batch_size": "integer",
    "flush_interval_sec": "number",
    "max_hdf5_bytes": "integer",
    "delete_after_archive": "boolean",
    "handoff": "enum",
    "pack_buckets": "boolean",
    "host": "string",
    "port": "integer",
    "db": "integer",
    "monitor_hz": "number",
    "monitor": "mapping",
    "watchdog": "mapping",
    "force_serial_layers": "boolean",
    "sample_artifacts": "enum",
    "data": "string",
    "n_samples": "integer",
    "timeout_sec": "number",
    "stale_sec": "number",
    "poll_interval_sec": "number",
    "max_sample_retries": "integer",
}

# Keys-table navigation markers (human output). ▸ = deeper man page; · = leaf field.
_NAV_EXPAND = "\u25b8"  # ▸
_NAV_LEAF = "\u00b7"  # ·
_NAV_EXPAND_STYLE = "bold cyan"
_NAV_LEAF_STYLE = "dim"
_NAV_LEGEND = (
    f"{_NAV_EXPAND} cyan = deeper Jarvis man page (see Open next / man field in --json); "
    f"{_NAV_LEAF} dim = leaf field (Description only — not a subcommand)"
)
_LIST_NAV_LEGEND = (
    f"{_NAV_EXPAND} cyan = deeper Jarvis man page (see Open next / man field in --json); "
    f"{_NAV_LEAF} dim = listed item without a deeper page"
)


def _clean_man_prose(value: Any) -> str:
    """Use list wording in man prose instead of JSONPath-style brackets."""
    text = str(value or "")
    for old, new in (
        ("Calculators.Modules[]", "Calculators.Modules list"),
        ("Operas.Modules[]", "Operas.Modules list"),
        ("execution.input[]", "execution.input list"),
        ("execution.output[]", "execution.output list"),
        ("Variables[].distribution", "Variables list distribution"),
        ("Variables[]", "Variables list"),
        ("commands[]", "commands list"),
        ("modes[]", "modes list"),
    ):
        text = text.replace(old, new)
    return text

_DOMAINS = frozenset({"yaml", "sampler", "calculator", "operas", "tokens", "example"})
_ROOT_BLOCKS = ("Scan", "Sampling", "Calculators", "Operas", "EnvReqs", "LibDeps")
_CALC_TOPICS = ("module", "execution", "modes", "pools")
_EXAMPLE_CATALOG: dict[str, str] = {
    "bridson": "bin/quickstart_bridson_operas.yaml",
    "csv": "bin/quickstart_csv_operas.yaml",
    "dynesty": "bin/sampling/Sampling_Dynesty_Simple.yaml",
    "dynesty-full": "bin/sampling/Sampling_Dynesty_Full.yaml",
    "multinest": "bin/sampling/Sampling_MultiNest_Simple.yaml",
    "multinest-full": "bin/sampling/Sampling_MultiNest_Full.yaml",
    "calculator": "bin/quickstart_calculator.yaml",
}

# Portal owns the adapter implementation, while Jarvis man owns the task-card
# writing surface. Keep the direction-specific action vocabulary here so a user
# does not have to infer input actions from an output variables table.
_PORTAL_IO_GUIDANCE: dict[str, dict[str, dict[str, Any]]] = {
    "JSON": {
        "input": {
            "summary": "Input writes JSON values; output only reads selected values.",
            "fields": {
                "actions": (
                    "Input-only action list. Supported action: Dump. Each Dump "
                    "variable requires name; expression computes a value from the "
                    "sample namespace; entry is an optional dotted JSON path."
                ),
            },
            "actions": [
                "Dump: write each variables item into the JSON object; expression is optional, and entry targets a dotted nested path.",
                "Dump without expression reads data[name]; Dump with expression evaluates the expression and returns that named value as an observable.",
            ],
            "action_rows": [
                {
                    "name": "Dump",
                    "fields": "variables list",
                    "description": "Each item: name required; expression optional; entry optional dotted JSON path.",
                }
            ],
        },
        "output": {
            "summary": "Output reads selected values from a JSON object; it does not run input actions.",
            "fields": {
                "variables": "Output bindings: name is required; entry is an optional dotted JSON path and defaults to name.",
            },
            "example": """name: output
path: output.json
type: JSON
variables:
  - {name: result,   entry: result}
  - {name: mass,     entry: observables.mass}
  - {name: status,   entry: status}
  - {name: best_fit, entry: results.best_fit}
  - {name: errors,   entry: results.errors}
  - {name: labels,   entry: labels}""",
        },
    },
    "CSV": {},
    "DAT": {},
    "TSV": {},
    "Wolfram": {
        "input": {
            "summary": "Input writes Wolfram Language Association values; output only reads selected values.",
            "fields": {
                "actions": (
                    "Input-only action list. Supported action: Dump. Each variable "
                    "requires name; expression is optional; entry is an optional "
                    "dotted Association path, with numeric segments indexing lists."
                ),
            },
            "actions": [
                "Dump: write values into the Association; expression computes a value, and entry selects a dotted Association/list path.",
                "Dump without expression reads data[name]; Dump with expression returns the named computed value as an observable.",
            ],
            "action_rows": [
                {
                    "name": "Dump",
                    "fields": "variables list",
                    "description": "Each item: name required; expression optional; entry optional dotted Association path with numeric list indexes.",
                }
            ],
        },
        "output": {
            "summary": "Output reads selected values from a Wolfram Language Association; it does not run input actions.",
            "fields": {
                "variables": "Output bindings: name is required; entry is an optional dotted Association path with numeric list indexes.",
            },
            "extra_fields": {"variables": "list"},
        },
    },
    "Text": {
        "input": {
            "summary": "Text is input-only; the adapter performs ordered placeholder replacement or whole-file copy.",
            "fields": {
                "actions": "Input-only action list. Supported actions are Replace and File; actions run in order, and File is terminal.",
                "encoding": "Text encoding passed to the adapter; default: utf-8.",
            },
            "extra_fields": {"actions": "list", "encoding": "string"},
            "actions": [
                "Replace: variables require name and placeholder; expression is optional and values are written in scientific notation.",
                "File: source names a runtime file-path value in the sample data and copies that file byte-for-byte; it ends the action sequence.",
            ],
            "action_rows": [
                {
                    "name": "Replace",
                    "fields": "variables list",
                    "description": "Each item: name and placeholder required; expression optional; values use scientific notation.",
                },
                {
                    "name": "File",
                    "fields": "source string",
                    "description": "Copy the runtime source file byte-for-byte; terminal action.",
                },
            ],
        },
    },
    "SLHA": {
        "input": {
            "summary": "Input supports text-template replacement, structured SLHA edits, and whole-file copy; output reads blocks and decays.",
            "fields": {
                "actions": "Input action list. Supported actions: Replace, SLHA, and File; actions run in order, and File is terminal.",
                "encoding": "SLHA text encoding passed to the adapter; default: utf-8.",
            },
            "extra_fields": {"actions": "list", "encoding": "string"},
            "actions": [
                "Replace: variables use name/value or expression plus placeholder; values are formatted in scientific notation without parsing the whole file.",
                "SLHA: variables use block and entry, plus value/name/expression; entry may be a scalar or a list for a multi-index block; this performs a structured pyslha block edit.",
                "File: source names a runtime file-path value and replaces the destination contents; it ends the action sequence.",
            ],
            "action_rows": [
                {
                    "name": "Replace",
                    "fields": "variables list of dicts",
                    "description": "Template replacement: each dict needs placeholder; use value/name or expression; no full SLHA parse.",
                },
                {
                    "name": "SLHA",
                    "fields": "variables list of dicts",
                    "description": "Structured edit: each dict needs block and entry; use value/name or expression.",
                },
                {
                    "name": "File",
                    "fields": "source string",
                    "description": "Copy the runtime source file into the destination; terminal action.",
                },
            ],
        },
        "output": {
            "summary": "Output reads selected SLHA block entries and DECAY widths/branching ratios; it does not run input actions.",
            "fields": {
                "variables": "Bindings require name, block, and entry. DECAY entries use an integer parent id for width or a list [parent, daughters...] for branching ratio.",
                "encoding": "SLHA text encoding passed to the adapter; default: utf-8.",
            },
            "extra_fields": {"variables": "list", "encoding": "string"},
        },
    },
    "xSLHA": {
        "output": {
            "summary": "Output-only extended SLHA read-back; input writing and input actions are not supported.",
            "fields": {
                "variables": "Bindings require name, block, and entry; DECAY/DECAY1L entries support widths or branching ratios.",
            },
            "extra_fields": {"variables": "list"},
        },
    },
}

_PORTAL_IO_GUIDANCE.update({
    fmt: {
        "input": {
            "summary": f"Input writes {fmt} tabular data with Dump actions; output reads selected cells or columns.",
            "fields": {
                "header": "Whether the first row is a header; default: true.",
                "columns": "Optional list used to seed column names when the file is missing or empty.",
                "actions": (
                    "Input-only action list. Supported action: Dump. Each variable "
                    "requires name and column; expression is optional; row is optional "
                    "and defaults to 0."
                ),
            },
            "extra_fields": {"header": "boolean", "columns": "list", "actions": "list"},
            "actions": [
                "Dump: write each variables item to one table cell; column is required, row defaults to 0, and expression optionally computes the value.",
            ],
            "action_rows": [
                {
                    "name": "Dump",
                    "fields": "variables list",
                    "description": "Each item: name and column required; row defaults to 0; expression optional.",
                }
            ],
        },
        "output": {
            "summary": f"Output reads {fmt} table cells or columns; it does not run input actions.",
            "fields": {
                "header": "Whether the first row is a header; default: true.",
                "columns": "Optional list used to name columns when the file is missing or empty.",
                "variables": "Each binding requires name; column defaults to name, row is optional, and omitting row returns the whole column.",
            },
            "extra_fields": {"header": "boolean", "columns": "list", "variables": "list"},
        },
    }
    for fmt in ("CSV", "DAT", "TSV")
})

_PORTAL_IO_GUIDANCE["DAT"]["input"]["fields"] = {
    **_PORTAL_IO_GUIDANCE["DAT"]["input"]["fields"],
    "comment": "Comment prefix ignored while reading; default: #.",
}
_PORTAL_IO_GUIDANCE["DAT"]["input"]["extra_fields"] = {"header": "boolean", "columns": "list", "actions": "list", "comment": "string"}
_PORTAL_IO_GUIDANCE["DAT"]["output"]["fields"] = {
    **_PORTAL_IO_GUIDANCE["DAT"]["output"]["fields"],
    "comment": "Comment prefix ignored while reading; default: #.",
}
_PORTAL_IO_GUIDANCE["DAT"]["output"]["extra_fields"] = {"header": "boolean", "columns": "list", "variables": "list", "comment": "string"}

_PORTAL_IO_GUIDANCE["CSV"]["input"]["example"] = """\
name: input
path: input.csv
type: CSV
header: true
actions:
  - type: Dump
    variables:
      - name: mass
        expression: \"x * 125\"
        column: mass"""
_PORTAL_IO_GUIDANCE["DAT"]["input"]["example"] = """\
name: input
path: input.dat
type: DAT
header: true
actions:
  - type: Dump
    variables:
      - name: mass
        expression: \"x * 125\"
        column: mass"""
_PORTAL_IO_GUIDANCE["TSV"]["input"]["example"] = """\
name: input
path: input.tsv
type: TSV
header: true
actions:
  - type: Dump
    variables:
      - name: mass
        expression: \"x * 125\"
        column: mass"""
_PORTAL_IO_GUIDANCE["CSV"]["output"]["example"] = """\
name: output
path: output.csv
type: CSV
header: true
variables:
  - {name: chi2,       column: chi2,          row: 0}
  - {name: best_mass,  column: mass,          row: 0}
  - {name: samples,    column: samples}
  - {name: likelihood, column: log_likelihood, row: 0}
  - {name: point,      column: point}
  - {name: weights,    column: weight}"""
_PORTAL_IO_GUIDANCE["DAT"]["output"]["example"] = _PORTAL_IO_GUIDANCE["CSV"]["output"]["example"].replace("output.csv", "output.dat").replace("type: CSV", "type: DAT")
_PORTAL_IO_GUIDANCE["TSV"]["output"]["example"] = _PORTAL_IO_GUIDANCE["CSV"]["output"]["example"].replace("output.csv", "output.tsv").replace("type: CSV", "type: TSV")
_PORTAL_IO_GUIDANCE["Wolfram"]["input"]["example"] = """\
name: input
path: input.wl
type: Wolfram
actions:
  - type: Dump
    variables:
      - name: mass
        expression: \"x * 125\"
        entry: SMParameters.HiggsMass"""
_PORTAL_IO_GUIDANCE["Wolfram"]["output"]["example"] = """\
name: output
path: output.wl
type: Wolfram
variables:
  - {name: BR_bb,   entry: Channels.0.BR}
  - {name: mh1,     entry: Masses.0}
  - {name: status,  entry: Status}
  - {name: central, entry: Results.Central}
  - {name: errors,  entry: Results.Errors}
  - {name: points,  entry: Scan.0.Point}"""
_PORTAL_IO_GUIDANCE["Text"]["input"]["example"] = """\
name: input
path: hdecay.in
type: Text
actions:
  - type: Replace
    variables:
      - name: HiggsMass
        placeholder: \">>>HIGGS_MASS<<<\"
  - type: File
    source: upstream_output"""
_PORTAL_IO_GUIDANCE["SLHA"]["input"]["example"] = """\
name: input
path: input.slha
type: SLHA
actions:
  - type: Replace
    variables:
      - {name: M1, placeholder: \">>>>>M1<<<<<\"}
      - {name: M2, placeholder: \">>>>>M2<<<<<\", expression: \"2000\"}
      - {name: Mu, placeholder: \">>>>>Mu<<<<<\", expression: \"100\"}
      - {name: Tb, placeholder: \">>>>>Tb<<<<<\", expression: \"7\"}
  - type: SLHA
    variables:
      - {name: tanbeta, block: MINPAR, entry: 3,     expression: \"7\"}
      - {name: m0,      block: MINPAR, entry: 1,     value: 1000.0}
      - {name: mix,     block: NMIX,   entry: [1, 1], value: 0.99}"""
_PORTAL_IO_GUIDANCE["SLHA"]["input"]["examples"] = [
    _PORTAL_IO_GUIDANCE["SLHA"]["input"]["example"],
    """name: input
path: input.slha
type: SLHA
actions:
  - type: File
    source: flexop
save: false""",
]
_PORTAL_IO_GUIDANCE["SLHA"]["output"]["example"] = """\
name: output
path: output.slha
type: SLHA
variables:
  - {name: mh1,          block: MASS,  entry: 25}
  - {name: width,        block: DECAY, entry: 1000024}
  - {name: br,           block: DECAY, entry: [1000024, 23, 1000022]}
  - {name: m_neutralino, block: MASS,  entry: 1000022}
  - {name: m_chargino,   block: MASS,  entry: 1000024}
  - {name: width2,       block: DECAY, entry: 1000024}
  - {name: n11,          block: NMIX,   entry: [1, 1]}
  - {name: n12,          block: NMIX,   entry: [1, 2]}
  - {name: n22,          block: NMIX,   entry: [2, 2]}"""
_PORTAL_IO_GUIDANCE["xSLHA"]["output"]["example"] = """\
name: output
path: spectrum.xslha
type: xSLHA
variables:
  - {name: mH,         block: MASS,  entry: 25}
  - {name: width,      block: DECAY, entry: 25}
  - {name: neutralino, block: MASS,  entry: 1000022}
  - {name: chargino,   block: MASS,  entry: 1000024}
  - {name: higgs,      block: MASS,  entry: 25}"""

# Diagnostic codes most relevant to sampler Bounds (subset; D24.10 expands fully).
_METHOD_DIAGNOSTICS: dict[str, list[str]] = {
    "Bridson": ["JV2-MTH-010", "JV2-MTH-011", "JV2-MTH-012", "JV2-MTH-013"],
    "Random": ["JV2-MTH-001", "JV2-MTH-020", "JV2-MTH-021", "JV2-SCH-001"],
    "Grid": ["JV2-VAR-001", "JV2-SCH-001"],
    "CSV": ["JV2-MTH-030", "JV2-MTH-031", "JV2-SCH-001"],
    "AdaptiveBridson": ["JV2-MTH-041", "JV2-MTH-042", "JV2-SCH-001"],
    "Dynesty": ["JV2-BND-001", "JV2-BND-012", "JV2-BND-030", "JV2-SCH-001"],
    "MultiNest": ["JV2-BND-001", "JV2-BND-012", "JV2-BND-020", "JV2-SCH-001"],
}


def _console(*, force_terminal: bool | None = None) -> Console:
    return Console(
        stderr=False,
        force_terminal=force_terminal,
        highlight=False,
    )


def _bullet_panel_text(value: Any) -> str:
    """Render informational panel content as one bullet per non-empty line."""
    if isinstance(value, str):
        lines = value.splitlines()
    elif isinstance(value, (list, tuple)):
        lines = [str(item) for item in value]
    else:
        lines = str(value or "").splitlines()
    rendered: list[str] = []
    for line in lines:
        text = line.strip()
        if not text:
            continue
        rendered.append(text if text.startswith("•") else f"• {text}")
    return "\n".join(rendered)


def normalize_yaml_path(raw: str) -> str:
    """Normalize user path input to ``$.A.B.C`` form.

    Lists are represented by their field name only.  For example, the CLI
    topic for ``Calculators.Modules`` is ``yaml.Calculators.Modules``.
    """
    text = str(raw or "").strip()
    if "[" in text or "]" in text:
        raise ValueError(
            "array brackets are not supported in Jarvis man paths; "
            "use the field name and read its list type"
        )
    if text.startswith("$."):
        text = text[2:]
    elif text.startswith("$"):
        text = text[1:].lstrip(".")
    parts = [p for p in text.split(".") if p]
    return "$." + ".".join(parts) if parts else "$"


def _schema_status(schema: Mapping[str, Any]) -> str:
    status = str(schema.get("x-jarvis-status") or "stable").strip().lower()
    return "unstable" if status == "unstable" else "stable"


def _follow_pointer(doc: Mapping[str, Any], pointer: str) -> Any:
    """Resolve a JSON Pointer fragment (without leading ``#``)."""
    node: Any = doc
    text = str(pointer or "").lstrip("#")
    if text.startswith("/"):
        text = text[1:]
    if not text:
        return doc
    for part in text.split("/"):
        part = part.replace("~1", "/").replace("~0", "~")
        if not isinstance(node, Mapping) or part not in node:
            return None
        node = node[part]
    return node


def _expand_schema(
    schema: Mapping[str, Any] | None,
    *,
    base: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Expand a top-level ``$ref``, including fragment-only ``#/$defs/...``.

    Absolute refs use the global schema catalog. Fragment-only refs need *base*
    (the owning schema document). Sibling keys beside ``$ref`` are merged in.
    """
    if not isinstance(schema, Mapping):
        return {}
    ref = schema.get("$ref")
    if not isinstance(ref, str) or not ref:
        return dict(schema)

    target: Any = None
    if ref.startswith("#/"):
        if not isinstance(base, Mapping):
            # Fall back to stock resolve (leaves fragment unresolved).
            return resolve_schema_ref(schema)
        target = _follow_pointer(base, ref)
    elif "#/" in ref:
        doc_id, frag = ref.split("#/", 1)
        try:
            doc = schema_by_id(doc_id)
        except Exception:
            return resolve_schema_ref(schema)
        target = _follow_pointer(doc, "/" + frag)
        # Nested fragment targets may themselves use #/ relative to *doc*.
        if isinstance(target, Mapping):
            target = _expand_schema(target, base=doc)
    else:
        return resolve_schema_ref(schema)

    if not isinstance(target, Mapping):
        return dict(schema)
    merged = dict(target)
    for key, value in schema.items():
        if key == "$ref":
            continue
        merged[key] = value
    # If the target is still a bare fragment ref, expand once more against base.
    if isinstance(merged.get("$ref"), str) and str(merged["$ref"]).startswith("#/"):
        return _expand_schema(merged, base=base)
    return merged


def _collect_properties(
    schema: Mapping[str, Any],
    *,
    base: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Merge ``properties`` across ``allOf`` / ``$ref`` expansions (shallow)."""
    props: dict[str, Any] = {}
    node = _expand_schema(schema, base=base)
    if isinstance(node.get("properties"), Mapping):
        props.update(dict(node["properties"]))
    for item in node.get("allOf") or []:
        if isinstance(item, Mapping):
            props.update(_collect_properties(item, base=base))
    return props


def _required_keys(
    schema: Mapping[str, Any],
    *,
    base: Mapping[str, Any] | None = None,
) -> set[str]:
    node = _expand_schema(schema, base=base)
    required: set[str] = set()
    for name in node.get("required") or []:
        required.add(str(name))
    for item in node.get("allOf") or []:
        if isinstance(item, Mapping):
            required |= _required_keys(item, base=base)
    # anyOf required alternatives (e.g. Point number | point_number)
    for item in node.get("anyOf") or []:
        if isinstance(item, Mapping):
            for name in item.get("required") or []:
                required.add(str(name))
    return required


def _type_label(
    schema: Mapping[str, Any],
    *,
    base: Mapping[str, Any] | None = None,
) -> str:
    ref = str(schema.get("$ref") or "")
    if ref.endswith("positiveNumeric"):
        return "number >0"
    if ref.endswith("integerish"):
        return "integer"
    if ref.endswith("numeric"):
        return "number"
    if ref.endswith("nonemptyString"):
        return "string"
    node = _expand_schema(schema, base=base)
    if "const" in node:
        return f"const {node['const']!r}"
    if "enum" in node:
        return "enum"
    t = node.get("type")
    if isinstance(t, list):
        return "list" if "array" in t else "|".join(str(x) for x in t)
    if t:
        return "list" if t == "array" else str(t)
    if node.get("oneOf") or node.get("anyOf"):
        # Numeric unions (common.json) after full resolve.
        if node.get("x-jarvis-numeric") or any(
            isinstance(item, Mapping) and item.get("x-jarvis-numeric")
            for item in (node.get("oneOf") or [])
        ):
            return "number"
        return "union"
    return "any"


def _key_entries(
    schema: Mapping[str, Any],
    *,
    base: Mapping[str, Any] | None = None,
    fallback_help: Mapping[str, str] | None = None,
    nav_map: Mapping[str, str] | None = None,
) -> list[dict[str, Any]]:
    props = _collect_properties(schema, base=base)
    required = _required_keys(schema, base=base)
    help_map = dict(fallback_help or {})
    keys: list[dict[str, Any]] = []
    for name, raw in props.items():
        raw_map = raw if isinstance(raw, Mapping) else {}
        node = _expand_schema(raw_map, base=base)
        desc = _clean_man_prose(raw_map.get("description") or node.get("description") or "")
        if not desc.strip():
            desc = str(help_map.get(name) or "")
        keys.append(
            {
                "name": name,
                "type": _type_label(raw_map, base=base),
                "required": name in required,
                "default": node.get("default"),
                "aliases": [],
                "description": desc,
                "minimum": node.get("minimum"),
                "exclusive": bool(node.get("exclusiveMinimum") is not None)
                or "positiveNumeric" in str(raw_map.get("$ref") or ""),
                "nav": False,
            }
        )
    return _annotate_nav(keys, nav_map)


def _annotate_nav(
    keys: list[dict[str, Any]],
    nav_map: Mapping[str, str] | None = None,
) -> list[dict[str, Any]]:
    """Mark keys that open a deeper ``Jarvis man`` page.

    *nav_map* maps key name to a full man
    command such as ``Jarvis man calculator.execution``. Leaf keys get
    ``nav: false`` and no ``man`` field.
    """
    mapping = {str(k): str(v) for k, v in dict(nav_map or {}).items()}
    out: list[dict[str, Any]] = []
    for key in keys:
        item = dict(key)
        name = str(item.get("name") or "")
        # Also allow lookup without input./output. prefix for dual-direction tables.
        short = name.split(".", 1)[-1] if "." in name else name
        cmd = mapping.get(name) or mapping.get(short)
        if cmd:
            item["nav"] = True
            item["man"] = cmd
        else:
            item["nav"] = False
            item.pop("man", None)
        out.append(item)
    return out


def _bounds_schema_for_method(method: str) -> tuple[dict[str, Any], str]:
    manifest = schema_manifest()
    methods = dict(manifest.get("sampling_methods") or {})
    # Case-insensitive method match
    uri = methods.get(method)
    if uri is None:
        for name, target in methods.items():
            if name.casefold() == method.casefold():
                method = name
                uri = target
                break
    if uri is None:
        raise KeyError(method)
    root = schema_by_id(str(uri))
    props = _collect_properties(root)
    bounds = props.get("Bounds")
    if isinstance(bounds, Mapping):
        return resolve_schema_ref(bounds), method
    return {}, method


def _method_page(method: str) -> dict[str, Any]:
    from jarvishep2.distributor import Distributor

    manifest = schema_manifest()
    methods = dict(manifest.get("sampling_methods") or {})
    uri = None
    canonical = method
    for name, target in methods.items():
        if name.casefold() == method.casefold():
            canonical = name
            uri = target
            break
    if uri is None:
        raise KeyError(method)
    root = schema_by_id(str(uri))
    status = _schema_status(root)
    bounds_schema, _ = _bounds_schema_for_method(canonical)
    keys = _key_entries(bounds_schema) if status == "stable" else []
    try:
        from jarvishep2.distributor import STATELESS_METHODS, Distributor as Dist

        resume = Dist.get_resume_status(canonical)
        is_stateless = canonical in STATELESS_METHODS
    except Exception:
        resume = "unknown"
        is_stateless = False

    summary = _clean_man_prose(root.get("description") or f"Sampling.Method = {canonical}")
    if status == "unstable":
        summary = (
            f"STATUS: not finalised — {canonical} Bounds knobs are not locked yet. "
            "Jarvis validate does not reject unknown Bounds keys for this method."
        )
    page = {
        "path": "$.Sampling.Bounds",
        "context": {"method": canonical},
        "zone": str(bounds_schema.get("x-jarvis-zone") or root.get("x-jarvis-zone") or "delegated"),
        "status": status,
        "summary": summary,
        "keys": keys,
        "examples": [str(root.get("x-jarvis-example"))] if root.get("x-jarvis-example") else [],
        "diagnostics": list(_METHOD_DIAGNOSTICS.get(canonical, [])),
        "see_also": [
            "Jarvis man sampler",
            "Jarvis man yaml.Sampling.Variables",
            f"Jarvis man example.{canonical.casefold()}",
        ],
        "further_reading": None,
        "capabilities": {
            "stateless": is_stateless,
            "resume": resume,
        },
    }
    return page


def _sampler_index() -> dict[str, Any]:
    from jarvishep2.distributor import STATELESS_METHODS, Distributor

    manifest = schema_manifest()
    rows = []
    for name in sorted(manifest.get("sampling_methods") or {}):
        root = schema_by_id(str(manifest["sampling_methods"][name]))
        status = _schema_status(root)
        rows.append(
            {
                "method": name,
                "status": status,
                "stateless": name in STATELESS_METHODS,
                "resume": Distributor.get_resume_status(name),
                "summary": str(root.get("description") or ""),
            }
        )
    return {
        "path": "$.Sampling.Method",
        "context": {},
        "zone": "closed",
        "status": "stable",
        "summary": "Sampling methods registered in the schema catalog.",
        "keys": [],
        "methods": rows,
        "examples": [],
        "diagnostics": [],
        "see_also": ["Jarvis man sampler.Bridson", "Jarvis man yaml.Sampling"],
        "further_reading": None,
    }


def _root_overview() -> dict[str, Any]:
    root = schema_by_id(str(schema_manifest()["root"]))
    props = _collect_properties(root)
    required_blocks = {str(name) for name in root.get("required") or []}
    one_of_blocks = {"Calculators", "Operas"}
    blocks = []
    for name in _ROOT_BLOCKS:
        node = props.get(name)
        if not isinstance(node, Mapping):
            blocks.append({"name": name, "zone": "?", "description": ""})
            continue
        resolved = resolve_schema_ref(node)
        blocks.append(
            {
                "name": name,
                "zone": str(resolved.get("x-jarvis-zone") or node.get("x-jarvis-zone") or "?"),
                "description": _clean_man_prose(resolved.get("description") or ""),
                "requirement": "one of Calculators / Operas" if name in one_of_blocks else None,
            }
        )
    nav_for_block = {
        "Scan": "Jarvis man yaml.Scan",
        "Sampling": "Jarvis man sampler",
        "Calculators": "Jarvis man calculator",
        "Operas": "Jarvis man operas",
        "EnvReqs": "Jarvis man yaml.EnvReqs",
        "LibDeps": "Jarvis man yaml.LibDeps",
    }
    return {
        "path": "$",
        "context": {},
        "zone": str(root.get("x-jarvis-zone") or "closed"),
        "status": "stable",
        "summary": (
            "• Blocks: Jarvis-HEP V2 task-card top-level blocks.\n"
            "• Required: Scan, Sampling, and EnvReqs.\n"
            "• Conditional: at least one of Calculators or Operas must appear.\n"
            "• Optional: LibDeps.\n"
            f"• Keys: {_NAV_EXPAND} means a deeper man page; {_NAV_LEAF} means a leaf."
        ),
        "keys": _annotate_nav(
            [
                {
                    "name": b["name"],
                    "type": "object",
                    "required": b["name"] in required_blocks,
                    **({"requirement": b["requirement"]} if b.get("requirement") else {}),
                    "default": None,
                    "aliases": [],
                    "description": (
                        f"zone={b['zone']}. {b['description']} "
                        "At least one of Calculators or Operas must appear."
                        if b.get("requirement")
                        else f"zone={b['zone']}. {b['description']}"
                    ).strip(),
                }
                for b in blocks
            ],
            nav_for_block,
        ),
        "examples": [],
        "diagnostics": ["JV2-SCH-001"],
        "see_also": [
            "Jarvis man sampler",
            "Jarvis man calculator",
            "Jarvis man operas",
            "Jarvis man tokens",
        ],
        "further_reading": None,
        "blocks": blocks,
        "nav_legend": _NAV_LEGEND,
    }


def _default_environment_yaml() -> str:
    path = Path(__file__).with_name("project_template") / "deps" / "environment_default.yaml"
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def _default_v2_yaml() -> str:
    text = _default_environment_yaml()
    lines = text.splitlines()
    start = next((i for i, line in enumerate(lines) if line.strip() == "V2:"), None)
    if start is None:
        return ""
    return "EnvReqs:\n" + "\n".join(lines[start:])


def _env_key(
    name: str,
    field_type: str,
    description: str,
    *,
    required: bool = False,
    default: Any = None,
    man: str | None = None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "name": name,
        "type": field_type,
        "required": required,
        "default": default,
        "aliases": [],
        "description": description,
        "nav": bool(man),
    }
    if man:
        item["man"] = man
    return item


def _envreqs_nested_page(
    path: str,
    title: str,
    fields: Mapping[str, str],
    example: str,
    *,
    diagnostic: str = "JV2-ENV-001",
) -> dict[str, Any]:
    keys = []
    for name, description in fields.items():
        field_type = _ENVREQS_FIELD_TYPES.get(name, "any")
        keys.append(_env_key(name, field_type, description))
    return {
        "path": path,
        "context": {},
        "zone": "delegated",
        "status": "stable",
        "summary": title,
        "keys": keys,
        "examples": [example],
        "diagnostics": [diagnostic],
        "see_also": ["Jarvis man yaml.EnvReqs", "Jarvis validate TASK.yaml"],
        "further_reading": None,
    }


def _envreqs_v2_page() -> dict[str, Any]:
    nav = {
        name: f"Jarvis man yaml.EnvReqs.V2.{name}"
        for name in (
            "sample_directory",
            "cleanup",
            "archiver",
            "redis",
            "factory",
            "worker",
            "check_modules",
        )
    }
    keys = [
        _env_key("workers", "integer", _ENVREQS_V2_HELP["workers"], default=4),
        _env_key("batch_size", "integer", _ENVREQS_V2_HELP["batch_size"], default=256),
        _env_key("sample_directory", "mapping", _ENVREQS_V2_HELP["sample_directory"], man=nav["sample_directory"]),
        _env_key("cleanup", "mapping", _ENVREQS_V2_HELP["cleanup"], man=nav["cleanup"]),
        _env_key("archiver", "mapping", _ENVREQS_V2_HELP["archiver"], man=nav["archiver"]),
        _env_key("redis", "mapping", _ENVREQS_V2_HELP["redis"], man=nav["redis"]),
        _env_key("factory", "mapping", _ENVREQS_V2_HELP["factory"], man=nav["factory"]),
        _env_key("worker", "mapping", _ENVREQS_V2_HELP["worker"], man=nav["worker"]),
        _env_key("check_modules", "mapping", _ENVREQS_V2_HELP["check_modules"], man=nav["check_modules"]),
        _env_key(
            "checkpoint_heartbeat_sec",
            "number",
            _ENVREQS_V2_HELP["checkpoint_heartbeat_sec"],
            default=30,
        ),
    ]
    overlay = (
        "EnvReqs:\n"
        "  # Jarvis project supplies the rest from deps/environment_default.yaml.\n"
        "  V2:\n"
        "    check_modules:\n"
        "      n_samples: 10\n"
        "      timeout_sec: 120"
    )
    default_section = _default_v2_yaml()
    examples = [overlay]
    if default_section:
        examples.append(default_section)
    return {
        "path": "$.EnvReqs.V2",
        "context": {},
        "zone": "delegated",
        "status": "stable",
        "summary": (
            "• Role: V2 runtime settings are an overlay, not normally a hand-written full block.\n"
            "• Defaults: Jarvis project create writes deps/environment_default.yaml.\n"
            "• Merge: the loader deep-merges that EnvReqs.V2 section, then applies task values.\n"
            "• Check mode: Jarvis check temporarily forces one worker and an unpacked smoke layout."
        ),
        "keys": keys,
        "examples": examples,
        "default_yaml": default_section,
        "diagnostics": ["JV2-ENV-001"],
        "see_also": [
            "Jarvis man yaml.EnvReqs",
            "Jarvis man example.calculator",
            "Jarvis project create NAME",
            "Jarvis validate TASK.yaml",
        ],
        "further_reading": {
            "command": "Jarvis project create NAME",
            "topic": "Creates deps/environment_default.yaml with the project V2 policy",
        },
    }


def _envreqs_page(parts: list[str]) -> dict[str, Any]:
    lowered = [str(part).casefold() for part in parts]
    if not lowered:
        keys = [
            _env_key(
                "Check_default_dependencies",
                "mapping",
                _ENVREQS_TOP_HELP["Check_default_dependencies"],
                man="Jarvis man yaml.EnvReqs.Check_default_dependencies",
            ),
            _env_key("Python", "mapping", _ENVREQS_TOP_HELP["Python"], man="Jarvis man yaml.EnvReqs.Python"),
            _env_key("CERN_ROOT", "mapping", _ENVREQS_TOP_HELP["CERN_ROOT"], man="Jarvis man yaml.EnvReqs.CERN_ROOT"),
            _env_key("V2", "mapping", _ENVREQS_TOP_HELP["V2"], man="Jarvis man yaml.EnvReqs.V2"),
            _env_key("Runtime", "mapping", _ENVREQS_TOP_HELP["Runtime"], man="Jarvis man yaml.EnvReqs.Runtime"),
        ]
        task_overlay = (
            "EnvReqs:\n"
            "  Check_default_dependencies:\n"
            "    required: true\n"
            "    default_yaml_path: \"&J/deps/environment_default.yaml\"\n"
            "  V2:\n"
            "    check_modules:\n"
            "      n_samples: 10\n"
            "      timeout_sec: 120"
        )
        examples = [task_overlay]
        default_text = _default_environment_yaml()
        if default_text:
            examples.append(default_text)
        return {
            "path": "$.EnvReqs",
            "context": {},
            "zone": "delegated",
            "status": "stable",
            "summary": (
                "• Ownership: EnvReqs is delegated to project/runtime policy.\n"
                "• Defaults: Jarvis project create writes deps/environment_default.yaml.\n"
                "• Loading: Check_default_dependencies loads the file and deep-merges V2.\n"
                "• Override: task values override the inherited defaults.\n"
                "• Inheritance: Python and CERN_ROOT are also commonly inherited.\n"
                "• V2 rule: do not add a top-level Runtime block to a V2 task card."
            ),
            "keys": keys,
            "examples": examples,
            "default_yaml": default_text,
            "diagnostics": ["JV2-ENV-001"],
            "see_also": [
                "Jarvis man yaml.EnvReqs.V2",
                "Jarvis man yaml.EnvReqs.Check_default_dependencies",
                "Jarvis project create NAME",
                "Jarvis validate TASK.yaml",
            ],
            "further_reading": {
                "command": "Jarvis project create NAME",
                "topic": "Creates the default environment YAML used by EnvReqs",
            },
        }

    if lowered[0] == "v2":
        if len(lowered) == 1:
            return _envreqs_v2_page()
        nested_key = lowered[-1]
        if nested_key in _ENVREQS_NESTED_PAGES:
            title, fields, example = _ENVREQS_NESTED_PAGES[nested_key]
            return _envreqs_nested_page(
                "$.EnvReqs.V2." + ".".join(parts[1:]),
                title,
                fields,
                example,
            )

    nested_key = lowered[0]
    if len(lowered) == 1 and nested_key in _ENVREQS_NESTED_PAGES:
        title, fields, example = _ENVREQS_NESTED_PAGES[nested_key]
        return _envreqs_nested_page(
            "$.EnvReqs." + parts[0], title, fields, example
        )
    raise KeyError("EnvReqs." + ".".join(parts))


def _path_page(path: str) -> dict[str, Any]:
    """Resolve a YAML path against the root card schema (best-effort)."""
    norm = normalize_yaml_path(path)
    parts = [p for p in norm[2:].split(".") if p] if norm.startswith("$.") else []
    if parts and parts[0].casefold() == "envreqs":
        return _envreqs_page(parts[1:])
    # Special case Sampling.Bounds without method context
    if len(parts) >= 2 and parts[0].casefold() == "sampling" and parts[1].casefold() == "bounds":
        return {
            "path": "$.Sampling.Bounds",
            "context": {},
            "zone": "delegated",
            "status": "stable",
            "summary": (
                "• Scope: method-specific sampler knobs.\n"
                "• Shape: depends on Sampling.Method.\n"
                "• Details: use Jarvis man sampler.<Method> for the closed key list."
            ),
            "keys": [],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man sampler", "Jarvis man sampler.Bridson"],
            "further_reading": None,
        }
    if len(parts) >= 2 and parts[0].casefold() == "sampling" and parts[1].casefold() == "variables":
        from jarvishep2.contracts.variables import DISTRIBUTION_TYPES, PARAMS_REQUIRED

        keys = []
        for dtype in sorted(DISTRIBUTION_TYPES):
            req = sorted(PARAMS_REQUIRED.get(dtype, ()))
            keys.append(
                {
                    "name": f"distribution.type={dtype}",
                    "type": "object",
                    "required": True,
                    "default": None,
                    "aliases": [],
                    "description": "Required parameters: " + ", ".join(req),
                }
            )
        return {
            "path": "$.Sampling.Variables",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Scope: list of sampling variables.\n"
                "• Item shape: each item has name + distribution."
            ),
            "list_title": "Distributions",
            "list_rows": keys,
            "keys": [],
            "examples": [
                "Variables:\n  - name: x\n    distribution:\n      type: Flat\n      parameters: {min: 0.0, max: 1.0}"
            ],
            "diagnostics": ["JV2-VAR-001"],
            "see_also": ["Jarvis man sampler.Bridson"],
            "further_reading": None,
        }

    # Route well-known Calculator / Operas paths to curated pages (richer examples).
    if parts and parts[0].casefold() == "calculators":
        return _calculator_yaml_path_page(parts)
    if parts and parts[0].casefold() == "operas":
        if len(parts) == 1:
            return _operas_page()
        if parts[1].casefold() == "modules":
            return _operas_page()

    # Walk root properties for first-level block
    root = schema_by_id(str(schema_manifest()["root"]))
    props = _collect_properties(root)
    if not parts:
        return _root_overview()
    head = parts[0]
    match = None
    for name in props:
        if name.casefold() == head.casefold():
            match = name
            break
    if match is None:
        raise KeyError(path)

    # Owning document for fragment-only $ref (e.g. Calculators.json #/$defs/module).
    raw_block = props[match] if isinstance(props.get(match), Mapping) else {}
    base_doc: Mapping[str, Any] | None = None
    block_ref = str(raw_block.get("$ref") or "")
    if block_ref and not block_ref.startswith("#") and "#/" not in block_ref:
        try:
            base_doc = schema_by_id(block_ref)
        except Exception:
            base_doc = None
    elif match.casefold() == "calculators":
        base_doc = schema_by_id(_CALC_SCHEMA_ID)
    elif match.casefold() == "operas":
        base_doc = schema_by_id(_OPERAS_SCHEMA_ID)

    cursor = _expand_schema(raw_block, base=base_doc)
    if base_doc is None and isinstance(cursor, Mapping) and cursor.get("$defs"):
        base_doc = cursor
    trail = [match]
    for part in parts[1:]:
        bare = part
        cursor = _expand_schema(cursor, base=base_doc)
        sub = _collect_properties(cursor, base=base_doc)
        items = cursor.get("items")
        if bare.casefold() == "modules" and "Modules" in sub:
            modules_node = _expand_schema(sub["Modules"], base=base_doc)
            trail.append("Modules")
            items_schema = modules_node.get("items")
            if isinstance(items_schema, Mapping):
                cursor = _expand_schema(items_schema, base=base_doc)
            else:
                cursor = modules_node
            continue
        # items schema for list fields (the path uses the field name only)
        if isinstance(items, Mapping) and not sub:
            cursor = _expand_schema(items, base=base_doc)
            sub = _collect_properties(cursor, base=base_doc)
        hit = None
        for name in sub:
            if name.casefold() == bare.casefold() or name.casefold() == part.casefold():
                hit = name
                break
        if hit is None and isinstance(items, Mapping):
            cursor = _expand_schema(items, base=base_doc)
            sub = _collect_properties(cursor, base=base_doc)
            for name in sub:
                if name.casefold() == bare.casefold():
                    hit = name
                    break
        if hit is None:
            # $defs on the owning document
            defs: Mapping[str, Any] = {}
            if isinstance(base_doc, Mapping) and isinstance(base_doc.get("$defs"), Mapping):
                defs = base_doc["$defs"]  # type: ignore[assignment]
            elif isinstance(cursor.get("$defs"), Mapping):
                defs = cursor["$defs"]  # type: ignore[assignment]
            for name, child in defs.items():
                if name.casefold() == bare.casefold():
                    hit = name
                    cursor = _expand_schema(child, base=base_doc)
                    trail.append(name)
                    break
            else:
                raise KeyError(path)
            continue
        cursor = _expand_schema(sub[hit], base=base_doc)
        trail.append(hit)

    fallback = None
    nav_map = None
    if match.casefold() == "calculators":
        fallback = _CALC_MODULE_FIELD_HELP
        # Same expand targets as calculator module / execution pages.
        nav_map = {
            "execution": "Jarvis man calculator.execution",
            "modes": "Jarvis man calculator.modes",
            "Modules": "Jarvis man calculator.module",
            "input": "Jarvis man calculator.execution.input --type JSON",
            "output": "Jarvis man calculator.execution.output --type JSON",
        }
    elif match.casefold() == "operas":
        fallback = _OPERAS_FIELD_HELP
    keys = _key_entries(
        cursor, base=base_doc, fallback_help=fallback, nav_map=nav_map
    )
    see = ["Jarvis man yaml"]
    further = None
    if match.casefold() == "calculators":
        see = [
            "Jarvis man calculator",
            "Jarvis man calculator.execution",
            "Jarvis portal man",
        ]
        further = {
            "command": "Jarvis portal man",
            "topic": "Portal adapter runtime (IO behaviour); man covers YAML keys only",
        }
    elif match.casefold() == "operas":
        see = ["Jarvis man operas", "Jarvis operas list", "Jarvis operas info NAME"]
        further = {
            "command": "Jarvis operas info helper.eggbox2d",
            "topic": "Operator signature and return shape (must be a mapping)",
        }
    return {
        "path": "$." + ".".join(trail),
        "context": {},
        "zone": str(cursor.get("x-jarvis-zone") or "delegated"),
        "status": _schema_status(cursor),
        "summary": _clean_man_prose(
            cursor.get("description") or f"YAML path $.{'.'.join(trail)}"
        ),
        "keys": keys,
        "examples": [str(cursor["x-jarvis-example"])] if cursor.get("x-jarvis-example") else [],
        "diagnostics": [],
        "see_also": see,
        "further_reading": further,
    }


def _calculator_yaml_path_page(parts: list[str]) -> dict[str, Any]:
    """Map Calculators.* paths onto calculator man pages (stable for agents)."""
    if len(parts) == 1:
        return _calculator_page(None)
    second = parts[1].casefold()
    if second == "pools":
        return _calculator_page("pools")
    if second == "modules":
        if len(parts) == 2:
            return _calculator_page("module")
        third = parts[2].casefold()
        if third == "execution":
            if len(parts) >= 4:
                direction = parts[3].casefold()
                if direction in {"input", "output"}:
                    # Optional format in a later segment is not required.
                    return _calculator_page(
                        "execution",
                        io_type=None,
                        direction=direction if direction in {"input", "output"} else None,
                    )
            return _calculator_page("execution")
        if third == "modes":
            return _calculator_page("modes")
        if third in _CALC_MODULE_FIELD_HELP or third in {
            "name",
            "path",
            "source",
            "installation",
            "required_modules",
        }:
            # Field of module — still show full module vocabulary.
            return _calculator_page("module")
        return _calculator_page("module")
    # Unknown under Calculators: fall through via generic error
    raise KeyError(".".join(parts))


def _io_format_page(direction: str, fmt: str) -> dict[str, Any]:
    direction = str(direction).strip().lower()
    if direction not in {"input", "output"}:
        raise KeyError(f"direction must be input|output, got {direction!r}")
    manifest = schema_manifest()
    io_map = dict((manifest.get("io") or {}).get(direction) or {})
    # Runtime format list (A10): Portal may expose formats beyond the catalog.
    runtime_names: set[str] = set()
    try:
        from jarvishep2.io_portal import available_io_formats

        runtime_names = {str(x) for x in (available_io_formats(direction) or [])}
    except Exception:
        runtime_names = set(io_map)

    uri = None
    canonical = fmt
    for name, target in io_map.items():
        if name.casefold() == fmt.casefold():
            canonical = name
            uri = target
            break
    if uri is None:
        for name in runtime_names:
            if name.casefold() == fmt.casefold():
                canonical = name
                break

    known_directions = _PORTAL_IO_GUIDANCE.get(canonical, {})
    if known_directions and direction not in known_directions and uri is None:
        raise KeyError(f"Portal format {canonical} does not support {direction}.")

    keys: list[dict[str, Any]] = []
    example = None
    summary = f"{direction} format {canonical}"
    guidance = _PORTAL_IO_GUIDANCE.get(canonical, {}).get(direction, {})
    further = {
        "command": f"Jarvis portal man {canonical}",
        "topic": (
            "Portal runtime adapter behaviour for this format. "
            "Jarvis man only lists YAML keys needed to write the card."
        ),
    }
    status = "stable"
    if uri is None:
        summary = (
            f"Format {canonical!r} is available from Portal at runtime, but this Jarvis "
            "build has no bundled field schema for it (field manual missing in this build)."
        )
    else:
        schema = schema_by_id(str(uri))
        _IO_FALLBACK = {
            "name": "Logical name of this IO slot.",
            "path": "File path relative to execution.path (quote tokens like \"&J/...\").",
            "type": f"Portal format name ({canonical}).",
            "save": "When true, keep the file after the sample finishes.",
            "actions": "Ordered write actions (input). Dump + variables with expression.",
            "variables": (
                "Output bindings: {name, entry} map file payload keys into observables. "
                "Input Dump actions use variables with expression instead."
            ),
        }
        keys = _key_entries(schema, fallback_help=_IO_FALLBACK)
        for key in keys:
            if not key.get("description"):
                key["description"] = _IO_FALLBACK.get(str(key.get("name")), "")
        example = schema.get("x-jarvis-example")
        summary = _clean_man_prose(schema.get("description") or summary)
        # Ensure example tokens that start with & are shown quoted in man prose when possible
        if example and "path: &" in str(example):
            example = str(example).replace("path: &", 'path: "&')
            # naive close — only if single-token path
            # leave schema examples as-is if complex
    field_help = guidance.get("fields") or {}
    for key in keys:
        field_name = str(key.get("name") or "")
        if field_name in field_help:
            key["description"] = str(field_help[field_name])
    existing_names = {str(key.get("name") or "") for key in keys}
    for field_name, field_type in (guidance.get("extra_fields") or {}).items():
        if field_name in existing_names:
            continue
        keys.append(
            {
                "name": field_name,
                "type": field_type,
                "required": False,
                "default": None,
                "aliases": [],
                "description": str(field_help.get(field_name) or ""),
                "nav": False,
            }
        )
    if guidance.get("summary"):
        summary = f"{summary} {guidance['summary']}"
    if guidance.get("example"):
        example = str(guidance["example"])
    other = "output" if direction == "input" else "input"
    if uri is None:
        summary = f"• Format: {summary}"
    else:
        summary = f"• Schema: {summary}"
    if guidance.get("summary"):
        # The schema and Portal adapter explain different layers of the same
        # YAML entry; keep both as separate bullets in the human card.
        schema_summary = _clean_man_prose(
            schema.get("description")
            if uri is not None
            else str(summary).removeprefix("• Format: ")
        )
        adapter_summary = str(guidance["summary"])
        summary = (
            f"• Schema: {schema_summary or summary}\n"
            f"• Adapter: {adapter_summary}"
        )
    summary += f"\n• Other direction: use --direction {other}."
    page_examples = [str(example)] if example else []
    if guidance.get("examples"):
        page_examples = [str(item) for item in guidance["examples"]]
    return {
        "path": f"$.Calculators.Modules.execution.{direction}",
        "context": {"type": canonical, "direction": direction},
        "zone": "delegated",
        "status": status,
        "summary": summary,
        "keys": keys,
        "examples": page_examples,
        "diagnostics": ["JV2-SCH-002"],
        "see_also": [
            "Jarvis man calculator.execution",
            f"Jarvis man calculator.execution.{other} --type {canonical}",
            f"Jarvis portal man {canonical}",
        ],
        "further_reading": further,
        "runtime_formats": sorted(runtime_names),
        "action_help": list(guidance.get("actions") or []),
        "action_rows": list(guidance.get("action_rows") or []),
    }


def _io_format_both_directions(fmt: str) -> dict[str, Any]:
    """When --type is set without --direction, expose both input and output keys."""
    pages: list[dict[str, Any]] = []
    for direction in ("input", "output"):
        try:
            pages.append(_io_format_page(direction, fmt))
        except Exception:
            continue
    if not pages:
        raise KeyError(fmt)
    if len(pages) == 1:
        return pages[0]
    keys: list[dict[str, Any]] = []
    examples: list[str] = []
    action_help: list[str] = []
    action_rows: list[dict[str, Any]] = []
    for page in pages:
        direction = str((page.get("context") or {}).get("direction") or "")
        for key in page.get("keys") or []:
            item = dict(key)
            item["name"] = f"{direction}.{key.get('name')}"
            keys.append(item)
        examples.extend(str(x) for x in (page.get("examples") or []))
        action_help.extend(
            f"{direction}: {line}"
            for line in (page.get("action_help") or [])
        )
        for row in page.get("action_rows") or []:
            item = dict(row)
            item["name"] = f"{direction}.{row.get('name')}"
            action_rows.append(item)
    canonical = str((pages[0].get("context") or {}).get("type") or fmt)
    return {
        "path": "$.Calculators.Modules.execution",
        "context": {"type": canonical, "direction": "both"},
        "zone": "delegated",
        "status": "stable",
        "summary": (
            f"• Format: Portal {canonical}.\n"
            "• Vocabularies: input and output YAML keys differ.\n"
            "• Keys: combined rows are prefixed with input./output.\n"
            f"• Details: prefer Jarvis man calculator.execution.input|output --type {canonical}."
        ),
        "keys": keys,
        "examples": examples,
        "action_help": action_help,
        "action_rows": action_rows,
        "diagnostics": ["JV2-SCH-002"],
        "see_also": [
            f"Jarvis man calculator.execution.input --type {canonical}",
            f"Jarvis man calculator.execution.output --type {canonical}",
            f"Jarvis portal man {canonical}",
        ],
        "further_reading": {
            "command": f"Jarvis portal man {canonical}",
            "topic": "Portal runtime adapter (not a substitute for YAML field lists)",
        },
    }


def _calculator_page(
    topic: str | None,
    *,
    io_type: str | None = None,
    direction: str | None = None,
) -> dict[str, Any]:
    if topic is None:
        return {
            "path": "$.Calculators",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Flow: write inputs through Portal IO -> run commands -> read outputs.\n"
                "• Optional stage: Operas modules can consume the resulting observables.\n"
                "• Manual: Jarvis man lists YAML keys; use Jarvis portal man for adapter runtime."
            ),
            "list_title": "Topics",
            "list_rows": _annotate_nav(
                [
                    {
                        "name": t,
                        "type": "topic",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": {
                            "module": "One Calculators.Modules list entry (name, path, execution).",
                            "execution": "Run block: path, commands, input list, output list.",
                            "modes": "Multimode alternative execution profiles.",
                            "pools": "Named concurrency pools.",
                        }.get(t, ""),
                    }
                    for t in _CALC_TOPICS
                ],
                {
                    "module": "Jarvis man calculator.module",
                    "execution": "Jarvis man calculator.execution",
                    "modes": "Jarvis man calculator.modes",
                    "pools": "Jarvis man calculator.pools",
                },
            ),
            "keys": [],
            "examples": [],
            "diagnostics": [],
            "see_also": [
                "Jarvis man calculator.module",
                "Jarvis man calculator.execution",
                "Jarvis man calculator.execution.input --type JSON",
                "Jarvis man calculator.execution.output --type JSON",
                "Jarvis man calculator.modes",
                "Jarvis man calculator.pools",
                "Jarvis man operas",
                "Jarvis man tokens",
                "Jarvis portal man",
            ],
            "further_reading": {
                "command": "Jarvis portal man",
                "topic": "Portal IO adapter runtime (formats, read/write behaviour)",
            },
            "chain": [
                "1. Calculators.Modules list entries declare name/path/execution",
                "2. execution.input list writes files via Portal formats (JSON, SLHA, ...)",
                "3. execution.commands list runs the external program",
                "4. execution.output list reads files back into observables",
                "5. Operas.Modules list entries consume observables (optional)",
            ],
        }
    topic_l = topic.casefold()
    if topic_l == "execution":
        if io_type:
            if direction:
                return _io_format_page(str(direction), io_type)
            return _io_format_both_directions(io_type)
        # List formats from runtime + catalog
        try:
            from jarvishep2.io_portal import available_io_formats

            input_fmts = sorted(available_io_formats("input") or [])
            output_fmts = sorted(available_io_formats("output") or [])
        except Exception:
            manifest = schema_manifest()
            input_fmts = sorted((manifest.get("io") or {}).get("input") or {})
            output_fmts = sorted((manifest.get("io") or {}).get("output") or {})
        if direction in {"input", "output"}:
            formats = input_fmts if direction == "input" else output_fmts
            return {
                "path": f"$.Calculators.Modules.execution.{direction}",
                "context": {"direction": direction},
                "zone": "delegated",
                "status": "stable",
                "summary": (
                    f"• Purpose: execution.{direction} is a list of Portal file specifications.\n"
                    "• Formats: supported formats are listed in the Formats table below.\n"
                    f"• Details: run Jarvis man calculator.execution.{direction} --type FORMAT."
                ),
                "keys": [
                    {
                        "name": "format",
                        "type": "choice",
                        "required": True,
                        "default": None,
                        "aliases": [],
                        "description": "Portal format for each list item; choose from the Formats table below.",
                        "nav": False,
                    }
                ],
                "examples": [],
                "format_rows": [
                    {
                        "format": fmt,
                        "command": f"Jarvis man calculator.execution.{direction} --type {fmt}",
                        "description": (
                            f"Open the {fmt} Keys and Actions tables for this {direction} list."
                        ),
                    }
                    for fmt in formats
                ],
                "diagnostics": ["JV2-SCH-002"],
                "see_also": [
                    f"Jarvis man calculator.execution.{direction} --type JSON",
                    "Jarvis man calculator.execution",
                    "Jarvis portal man",
                ],
                "further_reading": {
                    "command": "Jarvis portal man",
                    "topic": "Portal runtime format behaviour",
                },
                "input_formats": input_fmts,
                "output_formats": output_fmts,
            }
        return {
            "path": "$.Calculators.Modules.execution",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Fields: path, commands, input list, and output list.\n"
                "• Paths: quote token paths, for example path: \"&J/calculators/echo\".\n"
                "• Portal: Use --type FMT --direction input|output for Portal field tables."
            ),
            "keys": _annotate_nav(
                [
                    {
                        "name": "path",
                        "type": "string",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": (
                            "Working directory for commands. Quote &J/... tokens so YAML "
                            "does not treat them as anchors."
                        ),
                    },
                    {
                        "name": "commands",
                        "type": "list",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": "Shell commands or {cmd,cwd} objects to run.",
                    },
                    {
                        "name": "input",
                        "type": "list",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": (
                            "Portal input file specs (Dump + expression). Formats: "
                            + ", ".join(input_fmts)
                        ),
                    },
                    {
                        "name": "output",
                        "type": "list",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": (
                            "Portal output file specs (variables + entry, not Load actions). "
                            "Formats: " + ", ".join(output_fmts)
                        ),
                    },
                ],
                {
                    "input": "Jarvis man calculator.execution.input --type JSON",
                    "output": "Jarvis man calculator.execution.output --type JSON",
                },
            ),
            "examples": [_EXECUTION_EXAMPLE],
            "diagnostics": [],
            "see_also": [
                "Jarvis man calculator.execution.input --type JSON",
                "Jarvis man calculator.execution.output --type JSON",
                "Jarvis man example.calculator",
                "Jarvis portal man",
                "Jarvis man operas",
            ],
            "further_reading": {
                "command": "Jarvis portal man JSON",
                "topic": "How Portal reads/writes JSON at runtime",
            },
            "input_formats": input_fmts,
            "output_formats": output_fmts,
        }
    if topic_l == "module":
        calc = schema_by_id(_CALC_SCHEMA_ID)
        fields = (calc.get("$defs") or {}).get("moduleFields") or {}
        return {
            "path": "$.Calculators.Modules",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Scope: one calculator module entry in the Calculators.Modules list.\n"
                "• Card: YAML keys for the module.\n"
                "• Runtime: install/run behaviour follows path + execution."
            ),
            "keys": _key_entries(
                fields if isinstance(fields, Mapping) else {},
                base=calc,
                fallback_help=_CALC_MODULE_FIELD_HELP,
                nav_map={
                    "execution": "Jarvis man calculator.execution",
                    "modes": "Jarvis man calculator.modes",
                },
            ),
            "examples": [
                "Modules:\n"
                "  - name: Echo\n"
                '    path: "&J/calculators/echo"\n'
                "    execution:\n"
                '      path: "&J/calculators/echo"\n'
                "      commands:\n"
                "        - python3 -c \"print(1)\"\n"
                "      # input and output are optional lists"
            ],
            "diagnostics": [],
            "see_also": [
                "Jarvis man calculator.execution",
                "Jarvis man calculator.modes",
                "Jarvis man example.calculator",
            ],
            "further_reading": {
                "command": "Jarvis portal man",
                "topic": "Portal IO formats used under execution.input/output",
            },
        }
    if topic_l == "modes":
        return {
            "path": "$.Calculators.Modules.modes",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Pack: multimode modules share one physical pack.\n"
                "• Execution: each mode has its own execution block.\n"
                "• Reference: use dot-qualified mode names."
            ),
            "keys": [
                {
                    "name": "modes",
                    "type": "list",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Alternative execution profiles; mutually exclusive with top-level execution.",
                }
            ],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man calculator.module"],
            "further_reading": None,
        }
    if topic_l == "pools":
        return {
            "path": "$.Calculators.Pools",
            "context": {},
            "zone": "closed",
            "status": "stable",
            "summary": (
                "• Role: named integer pools provide concurrent calculator capacity.\n"
                "• Value: each pool name maps to its maximum concurrent slot count."
            ),
            "keys": [
                {
                    "name": "<pool-name>",
                    "type": "integer",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Max concurrent slots for modules assigned to this pool.",
                }
            ],
            "examples": ["Pools:\n  susy: 4"],
            "diagnostics": [],
            "see_also": ["Jarvis man calculator"],
            "further_reading": None,
        }
    raise KeyError(topic)


def _operas_page() -> dict[str, Any]:
    operas = schema_by_id(_OPERAS_SCHEMA_ID)
    modules = (operas.get("properties") or {}).get("Modules")
    item: Mapping[str, Any] = {}
    if isinstance(modules, Mapping) and isinstance(modules.get("items"), Mapping):
        item = modules["items"]
    # Expand operasIO item refs for input/output when present
    keys = _key_entries(item, base=operas, fallback_help=_OPERAS_FIELD_HELP)
    # If schema leaves input/output as opaque, still ensure descriptions.
    for key in keys:
        if not key.get("description"):
            key["description"] = _OPERAS_FIELD_HELP.get(str(key.get("name")), "")
    # operator row points at Operas CLI (not another man page) — leave nav false;
    # leaf fields only on this page.
    return {
        "path": "$.Operas.Modules",
        "context": {},
        "zone": str(item.get("x-jarvis-zone") or "closed"),
        "status": "stable",
        "summary": (
            "• Role: Operas.Modules list entries write derived observables from in-process operators.\n"
            "• Contract: operator must return a mapping (dict).\n"
            "• Catalog: Jarvis man documents YAML shape; use Jarvis operas list/info for operators."
        ),
        "keys": keys,
        "examples": [_OPERAS_EXAMPLE],
        "diagnostics": ["JV2-OPR-002"],
        "see_also": [
            "Jarvis man calculator",
            "Jarvis man tokens",
            "Jarvis operas list",
            "Jarvis operas info helper.eggbox2d",
        ],
        "further_reading": {
            "command": "Jarvis operas info helper.eggbox2d",
            "topic": "Operator signature, flags, and return shape (mapping required)",
        },
    }


def _tokens_page() -> dict[str, Any]:
    keys = [
        {"name": "&J/...", "type": "token", "required": False, "default": None, "aliases": [], "description": "Project root (standalone project or cwd)."},
        {"name": "@Sdir", "type": "token", "required": False, "default": None, "aliases": [], "description": "Per-sample working directory under SAMPLE/."},
        {"name": "@PackID", "type": "token", "required": False, "default": None, "aliases": [], "description": "Physical calculator pack id for multimode affinity."},
        {"name": "@SampleID", "type": "token", "required": False, "default": None, "aliases": [], "description": "Current sample uuid."},
        {"name": "${LibDeps:NAME}", "type": "token", "required": False, "default": None, "aliases": [], "description": "Path from LibDeps-built shared libraries."},
        {"name": "${Scan:KEY}", "type": "token", "required": False, "default": None, "aliases": [], "description": "Value from the Scan block (unknown keys may fall back silently)."},
        {"name": "@{ROOT path}", "type": "token", "required": False, "default": None, "aliases": [], "description": "ROOT file path helper where supported."},
    ]
    return {
        "path": "$",
        "context": {"domain": "tokens"},
        "zone": "open",
        "status": "stable",
        "summary": "Path and expression tokens used in task YAML.",
        "list_title": "Tokens",
        "list_rows": keys,
        "keys": [],
        "examples": ["path: &J/calculators/runtime\ncommand: ${LibDeps:flexiblesusy}/bin/Spectrum"],
        "diagnostics": [],
        "see_also": ["Jarvis man calculator.execution", "Jarvis man yaml.LibDeps"],
        "further_reading": None,
    }


def _example_page(name: str | None) -> dict[str, Any]:
    root = Path(__file__).with_name("project_template")
    if name is None:
        return {
            "path": "$",
            "context": {"domain": "example"},
            "zone": "open",
            "status": "stable",
            "summary": "Minimal task cards shipped with the package (same files as project create).",
            "list_title": "Examples",
            "list_rows": _annotate_nav(
                [
                    {
                        "name": key,
                        "type": "example",
                        "required": False,
                        "default": None,
                        "aliases": [],
                        "description": rel,
                    }
                    for key, rel in _EXAMPLE_CATALOG.items()
                ],
                {
                    key: f"Jarvis man example.{key}"
                    for key in _EXAMPLE_CATALOG
                },
            ),
            "keys": [],
            "examples": [],
            "diagnostics": [],
            "see_also": ["Jarvis man example.bridson"],
            "further_reading": None,
        }
    key = name.casefold()
    rel = _EXAMPLE_CATALOG.get(key)
    if rel is None:
        raise KeyError(name)
    path = root / rel
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    return {
        "path": "$",
        "context": {"example": key},
        "zone": "open",
        "status": "stable",
        "summary": f"Shipped template card: {rel}",
        "keys": [],
        "examples": [text],
        "diagnostics": [],
        "see_also": ["Jarvis man yaml", "Jarvis validate"],
        "further_reading": None,
        "body": text,
    }


def _center_page() -> dict[str, Any]:
    return {
        "path": "$",
        "context": {"domain": "center"},
        "zone": "open",
        "status": "stable",
        "summary": (
            "• Purpose: write a correct task-card YAML.\n"
            "• Portal runtime: use Jarvis portal.\n"
            "• Operator catalogs: use Jarvis operas.\n"
            f"• Keys: {_NAV_EXPAND} (cyan) means a deeper man page exists; "
            f"{_NAV_LEAF} (dim) means a leaf field — read Description only.\n"
            "• Topic syntax: use one dotted path only; do not split segments with spaces.\n"
            "• Example: Jarvis man calculator.execution.input --type JSON.\n"
            "• Navigation: Jarvis man calculator.module opens calculator.execution for the execution row.\n"
            "• Leaf paths have no subcommand, for example calculator.module.path is not a topic."
        ),
        "list_title": "Domains",
        "list_rows": _annotate_nav(
            [
                {
                    "name": "yaml",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Task-card structure and path lookup",
                },
                {
                    "name": "sampler",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Sampling.Method Bounds knobs and capabilities",
                },
                {
                    "name": "calculator",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Calculators.Modules / execution / modes / pools",
                },
                {
                    "name": "operas",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Operas.Modules list YAML shape",
                },
                {
                    "name": "tokens",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "&J / @Sdir / ${LibDeps:} / ${Scan:}",
                },
                {
                    "name": "example",
                    "type": "domain",
                    "required": False,
                    "default": None,
                    "aliases": [],
                    "description": "Emit a minimal runnable card",
                },
            ],
            {
                "yaml": "Jarvis man yaml",
                "sampler": "Jarvis man sampler",
                "calculator": "Jarvis man calculator",
                "operas": "Jarvis man operas",
                "tokens": "Jarvis man tokens",
                "example": "Jarvis man example",
            },
        ),
        "keys": [],
        "examples": [],
        "diagnostics": [],
        "see_also": [
            "Jarvis man yaml",
            "Jarvis man sampler.Bridson",
            "Jarvis man calculator.execution.input --type JSON",
            "Jarvis man calculator.execution.output --type JSON",
            "Jarvis man operas",
            "Jarvis validate TASK.yaml",
            "Jarvis portal man",
            "Jarvis operas list",
        ],
        "further_reading": {
            "command": "Jarvis portal man  |  Jarvis operas info NAME",
            "topic": "Runtime behaviour and operator catalogs (outside YAML field lists)",
        },
        "nav_legend": _NAV_LEGEND,
    }


def _coerce_man_topic(args: list[str]) -> str:
    """Return the single dotted topic path from positional man args.

    Only one positional topic is allowed (plus flags like ``--type``)::

        Jarvis man calculator.execution.input --type JSON

    Space-separated topics are rejected — use dots instead of spaces.
    """
    parts = [str(a).strip() for a in args if str(a).strip()]
    if not parts:
        return ""
    if any("[" in part or "]" in part for part in parts):
        raise KeyError(
            "array brackets are no longer supported; use the list field name "
            "without brackets, for example: Jarvis man yaml.Calculators.Modules.execution"
        )
    if len(parts) > 1:
        suggested = ".".join(p.lstrip(".") for p in parts)
        # Use a plain string so CLI prints without KeyError's extra quotes.
        raise KeyError(
            f"use a single dotted path (not spaces): Jarvis man {suggested}"
        )
    return parts[0].lstrip(".")


def resolve_man_request(
    tokens: list[str],
    *,
    io_type: str | None = None,
    direction: str | None = None,
    code: str | None = None,
) -> dict[str, Any]:
    """Resolve argv tokens after ``man`` into a page dict.

    Topic grammar (exactly one dotted path; no space-separated segments)::

        Jarvis man
        Jarvis man yaml
        Jarvis man yaml.Calculators.Modules.execution
        Jarvis man sampler
        Jarvis man sampler.Bridson
        Jarvis man calculator
        Jarvis man calculator.module
        Jarvis man calculator.execution
        Jarvis man calculator.execution.input --type JSON
        Jarvis man calculator.execution.output --type JSON
        Jarvis man operas
        Jarvis man tokens
        Jarvis man example.calculator
    """
    if code:
        from jarvishep2.man_codes import man_target_for_code, resolve_tokens_for_code

        entry = man_target_for_code(str(code))
        if entry is None:
            raise KeyError(code)
        tokens = resolve_tokens_for_code(str(code))
        # Fall through to normal token resolution.

    args = [str(t) for t in tokens if str(t).strip()]
    if not args:
        return _center_page()

    topic = _coerce_man_topic(args)
    # Bare $.… path → yaml walker
    if topic.startswith("$"):
        return _path_page(topic)

    segments = [p for p in topic.split(".") if p]
    if not segments:
        return _center_page()

    head = segments[0].casefold()
    rest = segments[1:]

    # Optional redundant yaml. prefix before another domain (A4).
    if (
        head == "yaml"
        and rest
        and rest[0].casefold() in _DOMAINS - {"yaml"}
    ):
        head = rest[0].casefold()
        rest = rest[1:]

    if head == "yaml":
        if not rest:
            return _root_overview()
        try:
            return _path_page(".".join(rest))
        except KeyError:
            choices = list(_ROOT_BLOCKS) + [
                "Bounds",
                "Variables",
                "Modules",
                "Method",
                "execution",
                "input",
                "output",
            ]
            suggestions = get_close_matches(rest[0], choices, n=5, cutoff=0.4)
            raise KeyError(
                f"{rest[0]} (did you mean: {', '.join(suggestions) or 'nothing close'})"
            ) from None
    if head == "sampler":
        if not rest:
            return _sampler_index()
        return _method_page(rest[0])
    if head == "calculator":
        topic_name = rest[0] if rest else None
        path_direction = direction
        path_type = io_type
        # calculator.execution.input.<FMT>
        if topic_name and topic_name.casefold() == "execution" and len(rest) >= 2:
            if rest[1].casefold() in {"input", "output"}:
                path_direction = path_direction or rest[1].casefold()
                if len(rest) >= 3 and not path_type:
                    path_type = rest[2]
            elif not path_type:
                # calculator.execution.JSON → type only
                path_type = rest[1]
        return _calculator_page(
            topic_name, io_type=path_type, direction=path_direction
        )
    if head == "operas":
        return _operas_page()
    if head == "tokens":
        return _tokens_page()
    if head == "example":
        return _example_page(rest[0] if rest else None)

    # Bare card path (Sampling.Bounds, Calculators.Modules, …)
    try:
        return _path_page(topic)
    except KeyError:
        choices = (
            list(_DOMAINS)
            + list(_ROOT_BLOCKS)
            + list(schema_manifest().get("sampling_methods") or {})
            + [
                "calculator.module",
                "calculator.execution",
                "calculator.execution.input",
                "calculator.execution.output",
                "sampler.Bridson",
            ]
        )
        suggestions = get_close_matches(topic, choices, n=5, cutoff=0.4)
        raise KeyError(
            f"{topic} (did you mean: {', '.join(suggestions) or 'nothing close'})"
        ) from None


def _print_page(page: Mapping[str, Any], *, as_json: bool) -> None:
    if as_json:
        # Stable public keys only (+ optional extras for indexes)
        payload = {
            "path": page.get("path"),
            "context": page.get("context") or {},
            "zone": page.get("zone"),
            "status": page.get("status"),
            "summary": page.get("summary"),
            "keys": page.get("keys") or [],
            "examples": page.get("examples") or [],
            "diagnostics": page.get("diagnostics") or [],
            "see_also": page.get("see_also") or [],
            "further_reading": page.get("further_reading"),
        }
        for extra in (
            "methods",
            "blocks",
            "chain",
            "capabilities",
            "input_formats",
            "output_formats",
            "runtime_formats",
            "format_rows",
            "list_title",
            "list_rows",
            "action_help",
            "action_rows",
            "body",
            "default_yaml",
        ):
            if extra in page:
                payload[extra] = page[extra]
        json.dump(payload, sys.stdout, indent=2, sort_keys=True)
        sys.stdout.write("\n")
        return

    console = _console()
    title = str(page.get("path") or "Jarvis man")
    ctx = page.get("context") or {}
    if ctx.get("method"):
        title = f"{title} · {ctx['method']}"
    if ctx.get("type"):
        title = f"{title} · type={ctx['type']}"
    if ctx.get("direction") and ctx.get("direction") != "both":
        title = f"{title} · {ctx['direction']}"
    status = str(page.get("status") or "stable")
    summary = str(page.get("summary") or "")
    if status == "unstable":
        summary = f"[bold red]STATUS: not finalised[/]\n{summary}"
    console.print(
        Panel(
            _bullet_panel_text(summary),
            title=title,
            box=box.ROUNDED,
            border_style="cyan",
        )
    )

    format_rows = list(page.get("format_rows") or [])
    if format_rows:
        table = Table(title=f"Formats · {title}", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("Format")
        table.add_column("Description", overflow="fold")
        for row in format_rows:
            table.add_row(
                str(row.get("format") or "—"),
                str(row.get("description") or "—"),
            )
        console.print(table)

    list_rows = list(page.get("list_rows") or [])
    if list_rows:
        list_title = str(page.get("list_title") or "List")
        table = Table(title=f"{list_title} · {title}", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("", justify="center", width=1, no_wrap=True)
        table.add_column("Name")
        table.add_column("Type")
        table.add_column("Description", overflow="fold")
        for row in list_rows:
            navigable = bool(row.get("nav"))
            mark = Text(_NAV_EXPAND, style=_NAV_EXPAND_STYLE) if navigable else Text(_NAV_LEAF, style=_NAV_LEAF_STYLE)
            name_cell = Text(
                str(row.get("name") or ""),
                style=_NAV_EXPAND_STYLE if navigable else "bold",
            )
            table.add_row(
                mark,
                name_cell,
                str(row.get("type") or "—"),
                str(row.get("description") or "—"),
            )
        console.print(table)
        if any(row.get("nav") for row in list_rows):
            console.print(
                Text(str(page.get("nav_legend") or _LIST_NAV_LEGEND), style="dim")
            )
        deeper = [str(row.get("man")) for row in list_rows if row.get("nav") and row.get("man")]
        if deeper:
            console.print(
                Panel(
                    _bullet_panel_text(deeper),
                    title=f"{_NAV_EXPAND} Open next",
                    box=box.ROUNDED,
                    border_style="cyan",
                )
            )

    if page.get("methods"):
        table = Table(title="Methods", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("", justify="center", width=1, no_wrap=True)
        table.add_column("Method")
        table.add_column("Status")
        table.add_column("Stateless")
        table.add_column("Resume")
        table.add_column("Summary", overflow="fold")
        for row in page["methods"]:
            method_name = str(row.get("method") or "")
            table.add_row(
                Text(_NAV_EXPAND, style=_NAV_EXPAND_STYLE),
                Text(method_name, style=_NAV_EXPAND_STYLE),
                str(row.get("status")),
                "yes" if row.get("stateless") else "no",
                str(row.get("resume")),
                str(row.get("summary") or ""),
            )
        console.print(table)
        console.print(
            Text(
                f"{_NAV_EXPAND} cyan = Jarvis man sampler.<Method>",
                style="dim",
            )
        )

    if page.get("chain"):
        console.print(
            Panel(
                _bullet_panel_text(page["chain"]),
                title="Calculator chain",
                box=box.ROUNDED,
            )
        )

    if page.get("capabilities"):
        cap = page["capabilities"]
        console.print(
            Panel(
                _bullet_panel_text(
                    f"stateless: {cap.get('stateless')}\nresume: {cap.get('resume')}"
                ),
                title="Capabilities",
                box=box.ROUNDED,
            )
        )

    keys = list(page.get("keys") or [])
    if keys:
        table = Table(title=f"Keys · {title}", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("", justify="center", width=1, no_wrap=True)
        table.add_column("Name")
        table.add_column("Type")
        table.add_column("Required")
        table.add_column("Description", overflow="fold")
        for key in keys:
            navigable = bool(key.get("nav"))
            if navigable:
                mark = Text(_NAV_EXPAND, style=_NAV_EXPAND_STYLE)
                name_cell = Text(str(key.get("name") or ""), style=_NAV_EXPAND_STYLE)
            else:
                mark = Text(_NAV_LEAF, style=_NAV_LEAF_STYLE)
                name_cell = Text(str(key.get("name") or ""), style="bold")
            table.add_row(
                mark,
                name_cell,
                str(key.get("type")),
                str(key.get("requirement") or ("yes" if key.get("required") else "optional")),
                str(key.get("description") or "—"),
            )
        console.print(table)
        # Legend whenever a Keys table is shown (center page also explains in summary).
        legend = str(page.get("nav_legend") or _NAV_LEGEND)
        console.print(Text(legend, style="dim"))
        # List follow-up commands for ▸ rows (agents + humans).
        deeper = [
            str(k.get("man"))
            for k in keys
            if k.get("nav") and k.get("man")
        ]
        if deeper:
            console.print(
                Panel(
                    _bullet_panel_text(deeper),
                    title=f"{_NAV_EXPAND} Open next",
                    box=box.ROUNDED,
                    border_style="cyan",
                )
            )
    elif status == "unstable":
        console.print(
            Panel(
                _bullet_panel_text(
                    "No Keys table: Bounds vocabulary is not finalised for this method."
                ),
                title=f"Keys · {title}",
                border_style="red",
                box=box.ROUNDED,
            )
        )

    action_rows = list(page.get("action_rows") or [])
    if action_rows:
        table = Table(title=f"Actions · {title}", box=box.SIMPLE_HEAVY, show_header=True)
        table.add_column("Action")
        table.add_column("Fields")
        table.add_column("Description", overflow="fold")
        for row in action_rows:
            table.add_row(
                str(row.get("name") or "—"),
                str(row.get("fields") or "—"),
                str(row.get("description") or "—"),
            )
        console.print(table)

    examples = list(page.get("examples") or [])
    if page.get("body"):
        examples = [str(page["body"])]
    for index, body in enumerate(examples, start=1):
        example_title = "YAML" if len(examples) == 1 else f"YAML example {index}"
        console.print(
            Panel(body, title=example_title, box=box.ROUNDED, border_style="green")
        )

    if page.get("default_yaml"):
        console.print(
            Panel(
                str(page["default_yaml"]),
                title="Project default YAML",
                box=box.ROUNDED,
                border_style="yellow",
            )
        )

    diags = list(page.get("diagnostics") or [])
    if diags:
        console.print(
            Panel(
                _bullet_panel_text(diags),
                title="Diagnostics that fire here",
                box=box.ROUNDED,
            )
        )

    see = list(page.get("see_also") or [])
    further = page.get("further_reading")
    if further:
        see.append(f"Further reading: {further.get('command')} ({further.get('topic')})")
    if see:
        console.print(
            Panel(
                _bullet_panel_text(see),
                title="See also",
                box=box.ROUNDED,
            )
        )


def _print_man_help() -> None:
    """Render ``Jarvis man -h`` in the same card language as man pages."""
    console = _console()
    console.print(
        Panel(
            _bullet_panel_text(
                "• Purpose: write a correct task-card YAML.\n"
                "• Topic: use one dotted path; list fields use their field name alone.\n"
                "• Domains: yaml, sampler, calculator, operas, tokens, and example.\n"
                "• Coding agents: use --json for structured keys, requirements, examples, and actions.\n"
                "• Runtime: use Jarvis portal for Portal adapter behaviour and Jarvis operas for operator catalogs."
            ),
            title="Jarvis man",
            box=box.ROUNDED,
            border_style="cyan",
        )
    )

    usage = Table(box=box.SIMPLE_HEAVY, show_header=False, expand=True)
    usage.add_column("Command", style="bold")
    usage.add_column("Meaning", overflow="fold")
    usage.add_row(
        "Jarvis man TOPIC ...",
        "Open the center page or one dotted YAML/manual topic.",
    )
    console.print(Panel(usage, title="Usage", box=box.ROUNDED, border_style="cyan"))

    options = Table(box=box.SIMPLE_HEAVY, show_header=True, expand=True)
    options.add_column("Option", style="bold", no_wrap=True)
    options.add_column("Description", overflow="fold")
    for option, description in (
        ("-h, --help", "Show this card."),
        ("--json", "Emit a machine-readable page; preferred for coding agents."),
        ("--code JV2-XXX-NNN", "Jump from a diagnostic code to its relevant man page."),
        ("--type FMT", "Open a Portal format field table, for example --type JSON."),
        (
            "--direction input|output",
            "Select Portal input or output vocabulary; prefer the direction in the dotted topic.",
        ),
    ):
        options.add_row(option, description)
    console.print(Panel(options, title="Options", box=box.ROUNDED, border_style="cyan"))

    examples = [
        "Jarvis man",
        "Jarvis man calculator.module",
        "Jarvis man calculator.execution.input --type JSON",
        "Jarvis man calculator.execution.input --type JSON --json  # coding agent",
        "Jarvis man sampler.Bridson",
        "Jarvis man yaml.Calculators.Modules.execution",
        "Jarvis man example.calculator",
    ]
    console.print(
        Panel(
            _bullet_panel_text(examples),
            title="Examples · one dotted path",
            box=box.ROUNDED,
            border_style="green",
        )
    )


def build_man_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="Jarvis man",
        description=(
            "YAML writing manuals (English). Focus: correct task-card keys and examples. "
            "Topic path uses dots: Jarvis man calculator.execution.input --type JSON. "
            "List fields use their field name without array suffixes. "
            "Domains: yaml, sampler, calculator, operas, tokens, example. "
            "Runtime: Jarvis portal … / Jarvis operas …"
        ),
        epilog=(
            "Examples (one dotted path; do not split with spaces):\n"
            "  Jarvis man\n"
            "  Jarvis man calculator.module\n"
            "  Jarvis man calculator.execution.output --type JSON\n"
            "  Jarvis man sampler.Bridson\n"
            "  Jarvis man yaml.Calculators.Modules.execution\n"
            "  Jarvis man example.calculator\n"
            "List fields use the field name alone; their Keys type is list.\n"
            "Not: Jarvis man calculator execution input"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "tokens",
        nargs="*",
        help="Single dotted topic path, e.g. calculator.execution.input (no spaces)",
    )
    parser.add_argument("--json", action="store_true", help="Machine-readable page")
    parser.add_argument("--code", metavar="JV2-XXX-NNN", help="Jump from a diagnostic code")
    parser.add_argument(
        "--type",
        dest="io_type",
        metavar="FMT",
        help="Portal format for calculator.execution field tables (e.g. JSON)",
    )
    parser.add_argument(
        "--direction",
        choices=("input", "output"),
        default=None,
        help="Portal direction (prefer calculator.execution.input|output in the topic path)",
    )
    return parser


def dispatch_man(argv: list[str] | None = None) -> int:
    raw_argv = list(argv or [])
    if any(token in {"-h", "--help"} for token in raw_argv):
        _print_man_help()
        return 0
    parser = build_man_parser()
    args = parser.parse_args(raw_argv)
    try:
        page = resolve_man_request(
            list(args.tokens or []),
            io_type=getattr(args, "io_type", None),
            direction=getattr(args, "direction", None),
            code=getattr(args, "code", None),
        )
    except KeyError as exc:
        print(f"Unknown man topic: {exc}", file=sys.stderr)
        print("Try: Jarvis man   or   Jarvis man sampler.Bridson", file=sys.stderr)
        return 2
    except Exception as exc:  # pragma: no cover - defensive
        print(f"Jarvis man failed: {exc}", file=sys.stderr)
        return 1
    _print_page(page, as_json=bool(args.json))
    return 0


__all__ = [
    "build_man_parser",
    "dispatch_man",
    "normalize_yaml_path",
    "resolve_man_request",
]
