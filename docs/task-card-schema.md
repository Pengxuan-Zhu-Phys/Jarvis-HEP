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
YAML parse
  -> ASCII task-card gate on the original user document: keys and string values only
  -> V1/default normalisation and runtime metadata injection
  -> JSON Schema: card structure, types, ownership zones
    -> Python contracts: cross-field, numerical, path, runtime semantics
```

JSON Schema errors are converted into the normal `ValidationReport` as
`JV2-SCH-*` diagnostics. Therefore `Jarvis2 validate`, `Jarvis2 run`, and
`Jarvis2 check` present one coherent report.

`EnvReqs.V2.checkpoint_heartbeat_sec` controls the stateless-sampler resume
checkpoint interval (default: `30`).  A heartbeat only requests a checkpoint;
the sampler captures its small generator state at the next complete batch.
The Archiver publishes an HDF5-durable `sample_index` prefix, so the stored
state is never newer than the contiguous on-disk prefix.  Reducing this value
reduces the maximum replay window after an interruption without adding a
per-sample checkpoint cost.

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

## Calculator modes

A calculator with independent program modes is written once as a parent and
expanded during task-card loading.  The runtime sees ordinary dotted module
names for dependency and I/O purposes, but their mutable calculator pack is
shared: every mode of one parent uses the same physical `@PackID` pool.

```yaml
Calculators:
  Modules:
    - name: Prospino
      path: "&J/calculators/runtime/Prospino/@PackID"
      clone_shadow: true
      modes:
        - name: ng
          installation:                 # config ng, then make if needed
            - "./configure ng && make"
          execution:
            commands: ["./prospino ng"]
            output:
              - {name: prospino_ng, path: output.dat, type: DAT}
        - name: ns
          installation:                 # config ns, then make if needed
            - "./configure ns && make"
          execution:
            commands: ["./prospino ns"]
            output:
              - {name: prospino_ns, path: output.dat, type: DAT}

    - name: Analysis
      required_modules: [Prospino]       # depends on Prospino.ng and Prospino.ns
      execution: {commands: ["./analyse"], output: []}

    - name: NGOnly
      required_modules: [Prospino.ng]    # depends on this mode only
      execution: {commands: ["./ng-only"], output: []}
```

There is no `@Mode` token and no `${mode_dir}` token.  All modes use the
parent path, for example `.../Prospino/001`; the acquired pack is rebuilt for
the requested mode only when its successful install stamp is for another mode
or has an outdated fingerprint/epoch.  The old stamp is removed before a
rebuild and written only after success, so an interrupted rebuild never looks
healthy.  `clone_shadow: true` and a parent path containing `@PackID` are
required for a multi-mode calculator.

The parent owns the shared physical settings: `source`, `deps_source`,
`clone_shadow`, `make_paraller`, and `symlink_name` cannot be overridden by a
mode. A mode may override `env_setup`, `timeout`, `selection`, and
`required_modules`. It may also specify an execution directory: the precedence
is `mode.execution.path` > `mode.path` > parent `path`. Parent `path` remains
the shared PackID installation anchor, so a mode path does not create a second
PackID pool. The parent `installation` is the shared base build: it
runs only when a PackID is first created, explicitly reinstalled, or its
parent fingerprint changes. A later mode switch runs only that mode's optional
`installation`. Parent and mode `initialization` both run for every sample
call, in that order. A parent with `modes` must not contain `execution`: every
mode must define its own execution block.

The dot is reserved as the module/mode separator.  Thus a bare name in
`required_modules` means every mode of that module, while `Module.mode` means
only that one mode; list several qualified names to depend on several modes.
Mode names, resulting module names, runtime paths, and output `name` values
must be unambiguous. `Calculators.Pools` names the physical parent pool:
write `Prospino: 4`, not `Prospino.ng` or `Prospino.ns`.

Redis keeps free shared PackIDs in mode-affinity lists. A Worker first tries a
pack already built for its requested mode, then an unassigned pack. When that
mode's warm pack is busy, the Worker waits for it for a bounded interval before
borrowing another mode's pack and rebuilding it. A mode whose `selection` is
false is removed before Redis acquisition, so it neither occupies nor relabels
a PackID. Within one sample,
modes of the same parent run serially because they alias the same physical
pool; independent calculator parents can still run concurrently. Before each
such serial group, the Worker greedily chooses a pending mode with the most
free warm packs, reducing avoidable rebuilds without changing physics order
between declared dependency layers.

Pool capacity determines whether affinity can settle. As an operational rule,
prefer `Calculators.Pools.<parent> >= number of modes + number of Workers` when
resources allow. A pool exactly equal to the mode count can still thrash under
Worker contention; a pool smaller than the mode count must rebuild because all
modes cannot be resident simultaneously. `JV2-MOD-009` reports undersized pools
as a warning rather than rejecting the task.

Each object schema has an `x-jarvis-zone` owner marker, self-tested when the
catalog loads:

| Zone | Meaning |
| --- | --- |
| `closed` | Jarvis2 owns the interface; an unknown key is an error. |
| `delegated` | A downstream component owns its keywords; Jarvis2 validates only its own declared envelope. |
| `open` | An identified, un-migrated surface; nested keys pass through with one warning and an explicit reason. |

The task-card root is closed: `Scan`, `Sampling`, `Calculators`, `Operas`,
`EnvReqs`, and `LibDeps` are the complete top-level vocabulary. A misspelling
such as `Calculater` is an error before the scan can start. `LibDeps` is an
explicit closed schema rather than a root-level exception. Likelihood
expressions belong in `Sampling.LogLikelihood`; top-level `Likelihood`,
`Mapper`, and `project_name` are not V2 task-card interfaces. Remove top-level
`Utils` and migrate
`Utils.interpolations_1D` to Jarvis-Operas with `interp1.*` for custom curves
or `dmddxe.*` for built-in direct-detection limits. `Calculators.path`,
`Modules[].deps_source`, and `Operas.Modules[].selection` are documented V1
compatibility fields rather than accidental exceptions. `EnvReqs` sibling V1
blocks, Portal I/O payloads,
Opera `kwargs`, and dynesty pass-through blocks are delegated. The legacy
RLTPMCMC `Control`, `Reward`, `PPO`, and `Diagnostics` blocks are temporarily
open; their presence emits `JV2-SCH-004`.

## LibDeps

`LibDeps` is a closed, V1-compatible declaration for shared native libraries.
`Modules` are installed once by the control process after environment preflight
and before Redis or Workers start. `required_modules` defines the build order;
independent modules may build concurrently up to `make_paraller`.

```yaml
LibDeps:
  path: "&J/deps/library"
  make_paraller: 4
  Modules:
    - name: ExampleLibrary
      required_modules: []
      installed: false             # accepted for V1 cards; V2 ignores it
      installation:
        path: "${LibDeps:path}/ExampleLibrary"
        source: "&J/deps/source/example.tar.gz"
        commands:
          - "cd ${LibDeps:path}"
          - "mkdir -p ${path}"
          - "make -j${LibDeps:make_paraller}"
```

The supported static tokens are `${LibDeps:path}`, `${LibDeps:make_paraller}`,
`${LibDeps:<module name>}`, and `@{ROOT path}`. The ROOT token requires
`EnvReqs.CERN_ROOT.path` or `get_path_command`. Each successful module writes
`<installation path>/.jarvis_install_stamp.json`; the graph control file is
`<LibDeps.path>/jarvis_install.json`. Matching stamps are reused. Set the
control file's `reinstall` value to `true` to rebuild all libraries on the next
run; V2 deliberately has no interactive reinstall prompt.

Use `Jarvis2 run TASK.yaml --skip-library-installation` (or `check`) only when
every declared installation path already exists. The command emits a warning,
does not build anything, and fails before Redis if a path is missing. Module
command output is recorded under `logs/<scan>/library-<module>.log`.

## Task-card text encoding

Task-card keys and every string value must be ASCII. `JV2-ENC-001` scans the
original parsed user document, before normalisation and runtime metadata
injection. Thus generated values such as `scan_name` and `task_result_dir` are
never presented as something the user must edit. Each issue identifies the
parsed YAML path and one-based character positions. This also protects
`Scan.name` before it can become an output directory, tar member, or HDF5
attribute.

Encoding and schema/contract validation still run in the same report. The
redundant unknown-key schema diagnostic for the exact non-ASCII key is removed,
but unrelated unknown keys in that object retain their normal diagnostic and
spelling suggestion. A card with both a non-ASCII key and `naem`, for example,
reports `JV2-ENC-001` and the suggestion to rename `naem` to `name` together.

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
