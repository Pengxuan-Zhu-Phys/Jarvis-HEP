# Changelog

User guides and reference documentation are available in
[Jarvis-Books/Docs/Jarvis-HEP](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/README.md).

## 2.1.1 — 2026-10-27
### Changed

- Remove the bundled `LatinHypercube` implementation and built-in catalog
  entry. The same algorithm is available as an explicitly loaded external
  sampler in Jarvis-HEP-Sampler's
  [`examples/LatinHypercube.py`](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP-Sampler/blob/main/examples/LatinHypercube.py).
- Require `Jarvis-HEP-Sampler>=1.2.0` without an upper bound. Default
  environment cards and runtime dependency checks use the same minimum.
- Require `Jarvis-Operas>=1.5.1` for the public unit-contract API.

### Added

- Support optional Opera input/output `unit` bindings through registered
  Jarvis-Operas contracts. Inputs declare the source unit; outputs select the
  stored unit. `unit_policy: strict` requires confirmed input units, while
  `canonical` records assumptions for inputs without a source declaration.
- Preserve canonical results and record identities when projecting output units.
  Store units and per-call audit in HDF5 and CSV, including when
  `store_samples: false`. Downstream unit-enabled Operas consume recorded sources.
- Accept JSON Pointer entries in unit-enabled modules alongside dot paths;
  output `entry: ""` preserves the complete canonical result. Unit diagnostics
  identify the affected binding and retain the underlying JO error code.

- Run external feedback samplers through the SDK's generational ask/tell
  interface. `Jarvis sampler load SliceMCMC.py` followed by `Jarvis run` now
  evaluates proposals in Workers and returns LogLikelihood to the plugin.
  The SDK's `SliceMCMC.eggbox.task.yaml` demonstrates the Opera path with a
  20,000 global trial limit.
- Save posterior transitions in `DATABASE/chain.csv`, separately from all
  trial evaluations in `samples.csv`. Sampler summaries distinguish proposal
  budgets, actual Worker evaluations, selection rejects, failures and retained draws.
- Restore feedback plugins from the scan's own source and complete durable
  barriers. Recovery replays unfinished proposals deterministically, deduplicates
  archived evaluations and truncates chain output to its checkpointed prefix.
- Versioned plugin feedback distinguishes valid `-inf` log density from missing
  or invalid scores and calculator failures. Generational evaluations reject
  boolean, NaN and positive-infinite scores; existing built-in feedback stays
  compatible.

- Download the PDF manual with `Jarvis download-manual`, with support for a
  specific GitHub Release, destination path and JSON output. Existing files
  require `--force`; failed downloads preserve the previous copy.
- Show installed versions of all five Jarvis packages below the `Jarvis -v`
  banner. `--version --json` provides the same inventory as structured data.
  Text output respects `NO_COLOR` and omits color codes when redirected.

### Fixed

- Format tiny elapsed times correctly in progress records instead of interpreting
  scientific-notation digits as milliseconds.

- Emit budget progress and periodic algorithm summaries for external feedback
  samplers. Logs distinguish proposals, Worker evaluations, retained transitions,
  failures and durable recovery progress, including resumed runs and early stops.
- Save generational checkpoints at the configured heartbeat instead of at every
  256 evaluations. Archive/chain safety barriers still bound replay memory;
  checkpoint records show the actual saved generation, trial and chain cursors.
- Flush and fsync checkpoint temporary files before atomic replacement, then sync
  the parent directory on POSIX. Serialization or write failures preserve the
  previous checkpoint and clean up the temporary file.

- Keep Redis sampling statistics updating when host metrics are unavailable.
  If Redis fails, retain the previous statistics with a stale indicator while
  host metrics continue to refresh. Each source recovers independently, and
  repeated failures do not flood the log.

### Compatibility and upgrade notes

- `unit` remains optional; cards without unit requests need no migration.
  Function unit contracts are required when requesting conversion. Write
  dimensionless units as the string `"1"`, rather than the YAML number `1`.
- Unit checkpoint resume supports Grid, Random, CSV and Bridson. Changes to
  bindings, unit policy, function contracts or JO version reject incompatible
  resume state. Other samplers support fresh unit runs and reject unit resume.

- Native generational plugins require SDK 1.2.0 with host contract 3. Upgrade
  both HEP and SDK before loading SliceMCMC; plugin API and stream protocol stay 1.
- Use `Sampling.Bounds.max_evaluations` for the global proposal ceiling, including
  initialization, burn-in and selection rejects. It is not a retained-chain target.
  `Sampling.LogLikelihood` is required, and `FeedbackReturn.include_logl` cannot
  be disabled. Failure handling belongs to the plugin; SliceMCMC halts by default.

- Upgrade the SDK to 1.2.0 or newer, obtain `examples/LatinHypercube.py` and
  its matching `.task.yaml` card from Jarvis-HEP-Sampler, and load the file
  with `jsampler load` or `Jarvis sampler load` before starting a new scan
  using `Sampling.Method: LatinHypercube`. An unloaded method is no longer
  listed by `Jarvis sampler list` or accepted in a new task card.
- Keep existing scans' method names and checkpoints unchanged. Scans with a
  valid source binding resume from their saved source with `--resume`, even
  when the former package file or external registration is absent. Moving
  the example does not change the algorithm or the `centered` option.

## 2.1.0 — 2026-10-18

Requires `Jarvis-HEP-Sampler>=1.1.0` (no upper bound) and `packaging>=23.0`;
pip installs both with Jarvis-HEP.

### Added

- External sampler plug-ins through Jarvis-HEP-Sampler. `Jarvis sampler
  init|check|load|list|info|remove` is the same command as `jsampler`; a loaded
  plug-in is used as `Sampling.Method` like a built-in method.
- `LatinHypercube`, shipped as an official plug-in.
- A scan keeps the sampler source it started with. Reloading or removing a
  registration affects later scans only; `--resume` uses the bound source.
  Project packages include the sampler source each scan is bound to.
- `Jarvis man sampler.<Method>` and the `Jarvis man` sampler index include
  registered plug-ins, read from the registry without running plug-in code.
  Plug-in problems use the `JV2-PLG-*` codes.
- `Jarvis agent` describes the `Jarvis sampler` commands and their effects.

### Changed

- `--resume` follows the checkpoint only. Without a checkpoint, or with
  `checkpoint.enabled: false`, it warns and starts from scratch. Changing
  `Sampling.Bounds` or the selection invalidates the checkpoint.
- Bridson writes its point set once to `points.npy` and resumes from it
  without regenerating the points.
- Running a scan again over an existing `DATABASE/` recomputes every point.
  At the end, rows beyond the new scan size and duplicate indices are removed;
  shard files keep their names and unchanged shards are not rewritten.
- `Jarvis run` and `Jarvis check` stop before starting Redis when a registered
  external sampler uses a built-in method name. The error gives the exact
  `Jarvis sampler remove <name> --yes` command.
- Grid and CSV checkpoints record a fingerprint of the point set or the CSV
  file. If it changed, or the checkpoint comes from an earlier release,
  `--resume` discards the checkpoint and starts from scratch.
- An unknown `Sampling.Method` suggests the correctly cased name or
  `Jarvis sampler load`.
- `EnvReqs` Python dependency versions accept `>` as well as `>=`. The default
  environment cards require `Jarvis-HEP-Sampler>=1.1.0` and `packaging>=23.0`.

### Fixed

- `--resume` after an interrupted fresh run no longer counts rows left from
  the previous scan as done.

## 2.0.17 — 2026-10-15

First PyPI release since 2.0.14. It includes the changes prepared as 2.0.15
and 2.0.16 that do not need the Sampler SDK, and has no Jarvis-HEP-Sampler
dependency.

### Added

- `Jarvis redis info|check|install|start`, each with its own help and JSON.
  `install` and `start` print the steps for the detected system; normal scans
  still manage their own Redis.
- `Jarvis agent` for coding agents: CLI contract, search, bundled schema and
  templates, `inspect` for effective card values and their sources, and
  `request` for a feedback email draft (never sent).
- Validation JSON includes manual and navigation hints.

### Changed

- Explicit duplicate YAML keys fail with `JV2-YAML-003`, the field path and
  both locations. Legal merge overrides still work.
- Unknown manual paths and conflicting selectors fail with guidance.

### Fixed

- An unknown Operas function such as `nosuch.f(z)` fails validation with
  `JV2-EXP-001` and the YAML path.
- Worker start-up errors reach the control process; the scan stops and
  cleans up instead of waiting.
- `Jarvis monitor --once --json` prints one JSON object, reports unknown counts
  as `null`, and also returns JSON on errors.
- Passing a scan card to `plot` explains that a plot card is needed and points
  to `Jarvis gen-plot-yaml`.
- Conflicting CLI intent flags are reported before or after the subcommand.
- Scan CPU on macOS and Windows is an estimate from observed process deltas
  and no longer freezes when processes exit.
- Nested samplers report `engine returned` instead of claiming convergence
  when the engine gives no reason.

## Earlier releases

See [GitHub Releases](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases)
and [PyPI](https://pypi.org/project/Jarvis-HEP/#history).
