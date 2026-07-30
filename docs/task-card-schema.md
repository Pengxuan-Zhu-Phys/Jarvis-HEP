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
    -> JSON Schema: card structure, types, closed user interfaces
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
| Calculator `input` / `output` | common fields and format-specific JSON layouts | Portal adapter availability and file contents |
| `Operas.Modules` | name/operator/call-mode/input/output structure | operator discovery, signatures, expression evaluation |

`EnvReqs`, project metadata, and integration extension blocks stay open at this
layer because they are normalised from V1/default files or have runtime-owned
semantics.

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

Every calculator I/O item first matches `ioCommon`:

```yaml
name: oupjson
path: "@Sdir/output.json"
type: JSON
save: false
```

`type` is the discriminator. Inputs and outputs use separate `oneOf` branches
because their allowed fields differ. JSON is implemented as the first strict
format:

* JSON input: `actions` with `Dump` actions and variables containing
  `name`/`expression`.
* JSON output: `variables` containing `name`/`entry`.
* Text input/output: a separate lightweight branch preserving the current
  `name`/`type`/optional `path` and `save` interface.

CSV, SLHA, XML, text, and future Portal adapters each get their own input and
output schema files. Unknown format names are `JV2-SCH-002` structural errors;
a format must be added to the manifest and its Portal adapter together.

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
2. List the files in `manifest.json.schema_files` and map the format names in
   `manifest.json.io`.
3. Add positive and negative card fixtures.
4. Confirm the matching Jarvis-Portal adapter is registered.
5. Keep content-dependent validation in Python/Portal, not JSON Schema.

## Rollout

Schema validation is enabled by default for loaded task cards. It emits
deterministic `JV2-SCH-001` diagnostics containing the JSON path and the
schema's error message. Existing Python `JV2-*` diagnostics continue to run,
so a card can report both structural and semantic corrections in one command.
