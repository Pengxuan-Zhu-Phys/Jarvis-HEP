# Task-card JSON Schema design

## Goal

Jarvis2 task cards remain YAML documents, but their stable user-facing
structure is validated by a versioned JSON Schema before the existing Python
contracts run. This gives editors useful completion and type diagnostics while
preserving the V2 runtime's richer semantic checks.

The schema set uses JSON Schema Draft 2020-12. `jarvishep2/schema/manifest.json`
is the central, data-only index: it defines the root card schema, every bundled
schema resource, and the sampler/I-O-format dispatch table. The validator only
loads files named there; it never fetches a remote schema.

## Validation layers

```text
YAML parse + V1/default normalisation
  -> ASCII task-card gate: keys and string values only
  -> JSON Schema: card structure, types, ownership zones
    -> Python contracts: cross-field, numerical, path, runtime semantics
```

JSON Schema errors are converted into the normal `ValidationReport` as
`JV2-SCH-*` diagnostics. Therefore `Jarvis2 validate`, `Jarvis2 run`, and
`Jarvis2 check` present one coherent report.

Every diagnostic also supplies a modification suggestion and, where its schema
owns a safe minimal fragment, a YAML example. See
`docs/validation-diagnostics.md` for the rendering and authoring contract.

## Schema scope

The schema is deliberately strict for the common, user-authored interfaces:

| Block | Schema responsibility | Python responsibility |
| --- | --- | --- |
| `Scan` | name/path scalar shapes and sample-directory structure | resolved output paths and project layout |
| `Sampling` | method enum, variables and distribution object shapes, method branches | value relations such as `min < max` and method-specific policy |
| `Calculators.Modules` | module/execution shape, command and I/O list shapes | command resolution, module installation, execution graph |
| Calculator `input` / `output` | registered format, common fields, and JSON-specific layouts | Portal adapter availability and file contents |
| `Operas.Modules` | name/operator/call-mode/input/output structure | operator discovery, signatures, expression evaluation |

Each object schema has an `x-jarvis-zone` owner marker, self-tested when the
catalog loads:

| Zone | Meaning |
| --- | --- |
| `closed` | Jarvis2 owns the interface; an unknown key is an error. |
| `delegated` | A downstream component owns its keywords; Jarvis2 validates only its own declared envelope. |
| `open` | An identified, un-migrated surface; nested keys pass through with one warning and an explicit reason. |

`Scan`, `Sampling`, `Calculators`, Calculator modules, and Opera modules are
closed interfaces. `Calculators.path`, `Modules[].deps_source`, and
`Operas.Modules[].selection` are documented V1 compatibility fields rather
than accidental exceptions. `EnvReqs` sibling V1 blocks, Portal I/O payloads,
Opera `kwargs`, and dynesty pass-through blocks are delegated. The legacy
RLTPMCMC `Control`, `Reward`, `PPO`, and `Diagnostics` blocks are temporarily
open; their presence emits `JV2-SCH-004`.

## Task-card text encoding

Task-card keys and every string value must be ASCII. This is enforced before
schema validation as `JV2-ENC-001`, so a non-ASCII key receives an honest
encoding diagnostic rather than a misleading unknown-key error. Each issue
identifies the parsed YAML path and one-based character positions. This also
protects `Scan.name` before it can become an output directory, tar member, or
HDF5 attribute.

Chinese (and every other non-ASCII script) remains fully supported in YAML
comments because comments are removed by the YAML parser:

```yaml
# Chinese explanation is allowed here.
Scan:
  name: ascii-scan
```

Human-readable diagnostics escape non-ASCII input as `\uXXXX` and use an
ASCII summary table, so table alignment stays deterministic even while
explaining a rejected key.

## File composition

The layout follows the successful V1 model while removing its runtime `$ref`
patching:

```text
schema/
  manifest.json                    # central composition/dispatch table
  task-card-v2.schema.json          # root card surface
  core/                             # Scan, common Sampling, Calculators, Operas
  sampling/methods/<method>.json    # one file per supported sampling method
  io/input/<format>.json            # one file per input format
  io/output/<format>.json           # one file per output format
```

Validation first applies the root schema, then looks up `Sampling.Method` and
each calculator I/O `type` in `manifest.json` and validates that object against
the declared local file. Adding a sampler or Portal format is therefore an
additive change: add its schema file, list it under `schema_files`, add its
dispatch entry, and add fixtures. The generic Python loader does not change.

## I/O format design

Every calculator I/O item supplies the common fields:

```yaml
name: oupjson
path: "@Sdir/output.json"
type: JSON
save: false
```

`type` is owned by the live Portal registry, not by the manifest. The names
are case-sensitive after surrounding whitespace is normalised. A manifest
entry is optional enrichment: if Portal supports a format without a bundled
Jarvis2 schema it is accepted and delegated to Portal; if Portal does not
support the format, `JV2-SCH-002` lists the formats valid for that direction.

JSON has the first strict format-specific layout:

* JSON input: `actions` with `Dump` actions and variables containing
  `name`/`expression`.
* JSON output: `variables` containing `name`/`entry`.
* Text is an input-only Portal format and requires `type: Text`.

Bundled formats have their own input and/or output schema file. Unknown format
names are `JV2-SCH-002` structural errors. Adding a Portal adapter does not
require a synchronous schema release; add a local schema later when Jarvis2
can usefully validate format-specific fields.

## Compatibility policy

The loader still accepts V1-compatible task cards and applies its existing
normalisation before schema validation. The schema only closes blocks whose V2
surface is already stable. V1 compatibility fields accepted by the current
runtime are explicitly represented where necessary, rather than silently
discarded.

The Python contracts remain the authority for semantics. A rule must not be
duplicated in JSON Schema when it needs filesystem access, expression parsing,
runtime registry discovery, or a numerical relationship between fields.

## Adding a format

1. Add `io/input/<format>.json` and/or `io/output/<format>.json`.
2. List the files in `manifest.json.schema_files` and optionally map the
   Portal format name in `manifest.json.io` for detailed local validation.
3. Add positive and negative card fixtures.
4. Confirm the matching Jarvis-Portal adapter is registered.
5. Keep content-dependent validation in Python/Portal, not JSON Schema.

## Rollout

Schema validation is enabled by default for loaded task cards. It emits
deterministically sorted diagnostics containing the JSON path, schema message,
allowed keys, and (for close spellings) a “Did you mean?” rename. Multiple
unknown keys each receive their own spelling suggestion; very large allowed-key
lists are compacted. Numeric YAML strings such as `1e-5` remain compatible with
V1 and emit `JV2-SCH-003`, suggesting canonical numeric notation; invalid
numeric strings explain the accepted forms (`0.05` or `1.0e-5`) directly.
Existing Python
`JV2-*` diagnostics continue to run, so a card can report both structural and
semantic corrections in one command. The corpus regression test covers all
example `*/bin/*.yaml` cards and allows only the explicitly unimplemented
sampler methods to remain runtime errors.
