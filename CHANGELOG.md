# Changelog

Developer documentation lives in
[Jarvis-Books/Docs/Jarvis-HEP](https://github.com/Pengxuan-Zhu-Phys/Jarvis-Books/blob/main/Docs/Jarvis-HEP/README.md).

## Unreleased

Requires `Jarvis-HEP-Sampler>=1.0.3` (no upper bound). This source also needs
an SDK build with host contract 2, which is not yet in a public SDK release.

### Added

- External sampler plug-ins through Jarvis-HEP-Sampler. `Jarvis sampler
  init|check|load|list|info|remove` is the same command as `jsampler`; a loaded
  plug-in is used as `Sampling.Method` like a built-in method.
- `LatinHypercube`, shipped as an official plug-in.
- Project packages include the sampler source each scan is bound to.

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
