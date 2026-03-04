# Jarvis-HEP Maintenance Context

Last updated: 2026-03-03
Audience: maintainers and future feature implementers.
Goal: preserve behavior while enabling safe iteration.

## 1. Project Snapshot

Jarvis-HEP is a YAML-driven orchestration framework for HEP scans:

- Builds a runtime context from one task YAML.
- Executes parameter generation and module workflows.
- Computes likelihood and nuisance pass conditions.
- Writes observables to HDF5, then converts to CSV with schema rules.

Current release state around this document:

- `master` and `v2` are aligned.
- `v1.5` tag exists as pre-PyPI publication baseline.
- `v1.5` and earlier are GitHub-only release versions.
- PyPI publication starts from `v1.6.0`.
- PyPI entry point is present (`Jarvis`) via `jarvishep.client:main`.
- v1.6.5 sampler hardening P0 has been implemented and smoke-tested (see `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`).

## 2. Runtime Architecture (Top-Down)

Entry and modes:

- package entry point: `jarvishep.client:main` (command: `Jarvis`).
- module entry point: `python -m jarvishep`.
- runtime dispatch implementation is in `jarvishep/client.py`.

Core orchestration:

- `jarvishep/core.py` performs init order:
  - parse args
  - infer runtime paths and project directories
  - configure logger
  - load and validate config
  - build workflow, factories, samplers
  - execute scan/convert/plot/monitor mode

Runtime path resolution:

- `jarvishep/base.py` is the path contract center.
- `&SRC` resolves to the installed `jarvishep` package root.
- `&J` resolves to runtime task root.
- no legacy `src` fallback is kept in `decode_path`; cards must use current package paths.

Config and module definitions:

- `jarvishep/config.py` loads YAML and normalizes optional sections.
- It resolves library/calculator/operas command paths and placeholders.
- It validates against sampler schema.

Execution graph:

- `jarvishep/workflow.py` resolves dependencies by module IO.
- `jarvishep/moduleManager.py` executes modules layer-by-layer.
- `jarvishep/modulePool.py` manages calculator instances and installation state.

Data and output:

- `jarvishep/hdf5writer.py` writes JSON-serialized records into HDF5 datasets.
- `jarvishep/observable_io.py` handles:
  - value normalization (`make_json_compatible`)
  - schema inference and merge
  - flatten policy for CSV export (`scalar/json/split/drop`)
  - split-column `name_map` aliases

Sampling and likelihood:

- `jarvishep/Sampling/*` hosts different samplers and nuisance support.
- `jarvishep/Module/likelihood.py` computes final `LogL`.
- `jarvishep/Module/nuisance_passCondition.py` computes boolean gates.

## 3. Hard Invariants (Do Not Break)

- Invariant A: runtime/source decoupling.
  - No feature may require running from repo root.
  - No new hardcoded legacy launcher assumptions or implicit local source-tree assumptions.

- Invariant B: observable compatibility.
  - Observables may include scalar, string, list, dict, ndarray-like values.
  - Serialization path must keep backward compatibility for existing scalar workflows.

- Invariant C: Likelihood/PassCondition output contracts.
  - `Likelihood` expression output must be scalar-castable (`float(...)`).
  - `PassCondition` expression output must be boolean-castable (`bool(...)`).
  - Array inputs are allowed only if expression reduces to scalar/bool output.

- Invariant D: schema compatibility.
  - Existing schema fields and modes must remain supported.
  - Invalid user-edited schema should normalize with warnings, not crash by default.

- Invariant E: failure behavior.
  - Failures should be explicit in logs.
  - No silent data loss in write/convert path.

## 4. Engineering Principles

- Principle 1: preserve user workflows first.
  - Keep backward compatibility unless a deliberate breaking change is announced.

- Principle 2: explicit contracts over hidden behavior.
  - Data shape/type expectations should be documented and tested.

- Principle 3: deterministic IO.
  - Schema and output conversion should be reproducible from saved artifacts.

- Principle 4: isolate side effects.
  - Keep path resolution, conversion, and module execution boundaries clear.

- Principle 5: incremental migration.
  - Add adapters/fallbacks before removing legacy behavior.

## 5. Explicit “Do Not” Rules

- Do not reintroduce coupling to repository layout (legacy launcher path or relative source-tree assumptions).
- Do not remove `&J`/`&SRC` semantics without a full migration path.
- Do not force structured observables back to scalar-only handling.
- Do not change flatten mode names or behavior silently.
- Do not bypass schema writeback when normalization modifies user schema.
- Do not hide broad exceptions without logging context.
- Do not add CLI options without updating:
  - `jarvishep/card/argparser.json`
  - help text behavior
  - README usage section
- Do not break command entry point `Jarvis`.
- Do not reintroduce runtime dependency installation flags or code paths; dependency resolution belongs to `pip install Jarvis-HEP` and `pyproject.toml`.
- Do not mix release/tag version intent with unrelated feature commits.

## 6. Risk Register From Current Review

- Risk R1: version mismatch.
  - `pyproject.toml` must align with PyPI policy (`>=1.6.0`) while preserving GitHub `v1.5` baseline.

- Risk R2: sampler-path test coverage is still shallow.
  - High-churn sampler code needs stronger contract and regression tests.

- Risk R3: singleton lifecycle in `ModuleManager`.
  - Hidden global state may leak across runs in long-lived processes.

- Risk R4: same-layer output collisions.
  - Multiple modules writing same observable key can cause last-write-wins ambiguity.

- Risk R5: full-memory CSV conversion.
  - conversion currently collects all rows before write, risky for very large snapshots.

- Risk R6: sampler/factory interface drift.
  - `WorkerFactory.submit_task` call signatures diverged across samplers and can break runtime.
  - Status: mitigated in v1.6.5 (`submit_task(sample_info)` unified, legacy interface removed).

- Risk R7: selection API contract ambiguity.
  - `super().evaluate_selection(...)` calls exist without a clearly guaranteed base implementation.
  - Status: mitigated in v1.6.5 (base contract restored for `evaluate_selection` / `check_evaluation`).

- Risk R8: silent module failure handling.
  - module execution path catches exceptions and continues, which can hide data-quality failures.
  - Status: mitigated in v1.6.5 (fail-fast default policy, explicit continue policy via configuration).

- Risk R9: nested parallelism and busy-wait overhead.
  - multiple thread pools + polling loops increase CPU overhead and scheduling contention.
  - Status: partially mitigated in v1.6.5 (hot-path busy-wait loops removed; nested execution simplification remains a v1.7.0 focus).

- Risk R10: numeric edge/precision hazards.
  - Grid endpoint transforms can hit infinities; Dynesty precision/mapping choices can corrupt outputs.
  - Status: partially mitigated in v1.6.5 (Dynesty precision/mapping fixed; Grid endpoint boundary guard still pending).

## 7. Safe Change Workflow

- Step 1: identify affected invariants in section 3.
- Step 2: patch smallest surface area first.
- Step 3: run minimal checks:
  - unit tests under `tests/`
  - `python -m build`
  - `Jarvis --help`
  - one smoke YAML run (or convert mode smoke)
- Step 4: update docs (`README`, schema docs, release notes) when behavior changes.
- Step 5: if compatibility changed, add migration note and explicit warning.

## 8. Release Workflow Notes

- Use `jhrel` as packaging helper, but keep release notes maintained in `docs/releases/`.
- Keep GitHub tag, release notes, and package version in sync with policy:
  - GitHub-only historical line: `<=1.5.x`
  - PyPI line: `>=1.6.0`
- Before PyPI publish:
  - dry build
  - TestPyPI verification
  - clean env install check
  - CLI smoke check

## 9. File Map For Frequent Maintenance

- Runtime entry: `Jarvis`, `jarvishep/client.py`, `jarvishep/__main__.py`
- Path resolution: `jarvishep/base.py`
- Config and schema validation: `jarvishep/config.py`, `jarvishep/card/schema/*`
- Workflow execution: `jarvishep/workflow.py`, `jarvishep/moduleManager.py`, `jarvishep/modulePool.py`
- Likelihood/nuisance: `jarvishep/Module/likelihood.py`, `jarvishep/Module/nuisance_*`
- Data output: `jarvishep/hdf5writer.py`, `jarvishep/observable_io.py`, `jarvishep/utils.py`
- Release and packaging: `pyproject.toml`, `jhrel`, `docs/releases/*`

## 10. Sampler Hardening Status (2026-03-03)

- Implemented in v1.6.5:
  - `submit_task(sample_info)` contract unified across samplers.
  - Base sampler contract restored for `evaluate_selection` / `check_evaluation`.
  - DNN correctness fixes (gating/import path) and JSON serialization hardening.
  - `Sample.close()` lifecycle hardening in sampler paths.
  - Module failure policy made explicit (fail-fast default).
  - Busy-wait futures polling removed in sampler hot paths.
  - Dynesty float16 downgrade removed; `samples_u` mapping fixed.
- Remaining follow-up:
  - deeper DNN torch-chain regression coverage beyond smoke tests.
  - broader concurrency simplification in v1.7.0 blueprint track.
  - SAMPLE directory archive/compression pipeline (non-blocking background process + YAML switch) tracked in `docs/releases/v1.6.5_SAMPLE_ARCHIVE_BLUEPRINT.md`.
- Residual P1 breakdown:
  - closure record for previously unresolved v1.6.5 P1 items is tracked in `docs/releases/v1.6.5_REMAINING_P1_BREAKDOWN.md`.
- Tracking source:
  - implementation backlog and residuals are maintained in `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`.

## 11. New-Chat Bootstrap Checklist

- Read:
  - `docs/TODO.md`
  - `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`
  - `docs/.internal/AGENT_PRIVATE_NOTES.md`
- Start from one P0 item at a time, each with:
  - code patch
  - tests
  - docs update
- Run gates after each meaningful patch:
  - `python -m unittest discover -s tests -p 'test_*.py'`
  - `python -m build`
  - `Jarvis --help`

## 12. SAMPLE Archive Compression Plan

Objective:

- Reduce file/inode pressure caused by huge `SAMPLE` bucket directories while keeping sampling throughput stable.

Execution model:

- Archive must run in a background independent process (not in Factory/sampler hot path).
- Main process behavior is enqueue-only and non-blocking.
- Archive target remains bucket-level (`SAMPLE/000001 -> SAMPLE_ARCHIVE/000001.tar.gz`), active bucket excluded.

YAML control key (top-level `Directory_Setting`):

- `archive_samples` (bool, default `true`)
  - `true`: enable archive dispatch.
  - `false`: disable archive dispatch.

Implementation notes:

- `ConfigLoader` now normalizes `Directory_Setting.archive_samples` default to `true`.
- `Core` runtime sample context now carries `sample.archive_samples` for downstream archive scheduler use.
- Archive format default: `tar.gz` for cross-machine compatibility.
- Safety path: temp archive write -> verify -> atomic rename -> optional prune.
