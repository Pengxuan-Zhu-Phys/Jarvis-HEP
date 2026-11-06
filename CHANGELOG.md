# Changelog

User guides and reference documentation are available in
[Jarvis-Books/Docs/Jarvis-HEP](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/README.md).

## Unreleased

### Fixed

- Align generated parameter-plot limits with finite sampling bounds, and use
  logarithmic axes for Log and Log-Normal variables, including direct Mapper aliases.
- Preserve canonical module metadata when using legacy import paths, avoiding
  relative-import deprecation warnings without changing module identity.
- Prevent watchdog recovery from combining different tasks' brief idle transitions
  into one timeout. Recheck task ownership before stopping a Worker and publish
  heartbeats in order so an older idle update cannot overwrite busy.
- Report a refused unit resume (`JV2-UNIT-004`) as CRITICAL, and say that the
  scan needs a new name or a moved-away DATABASE: a fresh run under the same name
  keeps the old rows, so it cannot resume either.
- Resume Dynesty and MultiNest from the engine snapshot in `runtime.pickle`.
  The snapshot is taken between dynesty iterations and does not wait for DATABASE.
- Confirm checkpoint progress with per-sample archive receipts that the
  Archiver posts after each HDF5 fsync. MCMC checkpoints still advance only
  after those receipts. Grid, Random, Bridson, CSV, AdaptiveBridson, Dynesty
  and MultiNest resume from the landed checkpoint file alone.
- Report a stalled or stopped Archiver at the end of a scan with the number of
  rows still missing, instead of waiting for the full timeout.
- Keep the sample index on the failure row that the watchdog writes after the
  last retry.

### Changed

- Allow 30 seconds for the same idle in-flight task, expire Worker heartbeats
  after 60 seconds, and check Workers every 2 seconds by default. Existing
  watchdog settings remain configurable.
- Dynesty and MultiNest `--resume` reads only `runtime.pickle`. A missing file
  is an error. `nested_engine.pkl` is a tools copy written by the checkpoint
  process; a fresh scan removes the old copy, and resume never reads it.
- Reject Dynesty and MultiNest checkpoints from earlier releases; the scan
  starts fresh and the log gives the reason.
- DATABASE is append-only. A fresh scan appends to an existing DATABASE, and
  `--resume` evaluates the samples after the last checkpoint again, so their
  rows can appear twice. The Archiver no longer reads or rewrites existing rows,
  and the end-of-scan cleanup of old rows is removed.
- Resume of Grid, Random, Bridson, CSV, AdaptiveBridson, Dynesty and MultiNest
  reads the landed `runtime.pickle` only and never scans DATABASE.

## 2.1.2

### Fixed

- Preserve nuisance pass-condition failures as failed samples; reject unknown
  nuisance methods, target modes and failure policies before execution.
- Reject non-positive Log prior minima, and report numeric range violations
  directly instead of claiming an integer or float is not a number.
- Keep run and performance counters consistent across nested evaluations and
  Worker retries; use the sampling method name in run summaries.
- Use pipes for FileOperation IPC so killed Workers leave no named semaphore
  queues, and report a stopped FileOperation child promptly.
- Resolve project-local Python modules and `${LibDeps:path}` installation paths;
  reject unknown module dependencies and cycles.
- Keep a blank `LogL` column in CSV exports when every evaluation failed.
- Keep DATABASE archiving fast in long fresh Dynesty and MultiNest scans. Each
  batch no longer re-reads all stored rows to look for duplicates.
- Show generational plugin budgets and a valid calculator-pool example in `man`;
  remove internal milestone text and unrelated nested-method smoke messages.
- Keep only the encrypted archive after a successful `project pack --encrypt`.
- Include the bundled Dynesty licence, author credits and provenance in wheels
  and source distributions.
- Give manifest packages the same project root directory, checksums and
  `.jarvis-pack` records as profile packages. Resolve relative archive outputs
  beside the manifest and preserve file/directory distinctions in exclusions.
- Exclude `.redis/` and `.jarvis2/` runtime directories from share and repro
  packages, including manifest packages with directory includes.
- Show relevant archiver, sample-directory and module-check examples in
  validation errors, and clarify which default environment sections are loaded.

### Changed

- Require `Jarvis-Operas>=1.5.2` and `JarvisPLOT>=2.1.14` so upgrades include
  the matching operator and plotting fixes. Sampler remains `>=1.2.0` with no upper bound.
- Write one `images/<scan>/plot.yaml` per scan with the sample plot and any
  available level-set overlay or Dynesty/MultiNest diagnostics. The manual
  `gen-plot-yaml` command produces the same card as scan finalization.

### Added

- Support `Jarvis --refs --json` for structured framework and sampler citations.
- Accept `Jarvis -V` as an alias for `-v` and `--version`, including JSON output.
  Help continues to advertise only the lowercase `-v` short option.

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
