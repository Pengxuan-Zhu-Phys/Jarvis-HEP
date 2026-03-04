# Jarvis-HEP Todo And Roadmap

Last updated: 2026-03-03
Scope: PyPI line maintenance, sampler hardening, release stability.
Active release track: `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`
Planned next release track: `docs/releases/v1.7.0_CONCURRENCY_SIMPLIFICATION_BLUEPRINT.md`
Previous release tracks:
- `docs/releases/v1.6.4_TASKLIST.md` (completed)
- `docs/releases/v1.6.3_TASKLIST.md` (completed)
- `docs/releases/v1.6.2_TASKLIST.md` (completed)
- `docs/releases/v1.6.2_DEPENDENCY_ALIGNMENT_TASKLIST.md` (completed)

## Version Policy

- `v1.5` and earlier are GitHub-only release lines.
- PyPI release line starts from `v1.6.0`.
- Do not publish any `<=1.5.x` package to PyPI.

## A. Completed Milestones

- [DONE] Decoupled runtime task root from source root.
- [DONE] Added installable entry point (`Jarvis`) via `pyproject.toml`.
- [DONE] Switched CLI entry target from legacy package entry to `jarvishep.client:main` and removed legacy repo-root launcher script.
- [DONE] Added `jhrel` release helper script for build/upload/reinstall verification workflow.
- [DONE] Added schema-driven observable export with array/list support in HDF5 flow.
- [DONE] Added CSV flatten rules with `split/json/drop/scalar` modes.
- [DONE] Added `name_map` support for user-controlled flattened column names.
- [DONE] Added schema warning normalization and persistence path.
- [DONE] Added unit tests for observable schema IO core behavior.

## A1. Active Focus (Current Cycle)

- [ACTIVE] Sampler robustness/correctness track is managed in `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`.
- [ACTIVE] Prioritize all sampler P0 items before performance-only optimizations.
- [ACTIVE] Keep release discipline: behavior fix + tests + docs in the same cycle.
- [STATUS] v1.6.5 sampler hardening P0/P1 items are complete (as of 2026-03-03); details and closure notes are tracked in `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`.
- [ACTIVE] SAMPLE archive-pressure mitigation tracked in `docs/MAINTENANCE_CONTEXT.md` section 12 (background process, non-blocking, YAML `Directory_Setting.archive_samples` default on); detailed draft also in `docs/releases/v1.6.5_SAMPLE_ARCHIVE_BLUEPRINT.md`.

## A2. Planned Focus (Next Cycle: v1.7.0 Concurrency Simplification)

- [PLANNED] Track docs: `docs/releases/v1.7.0_CONCURRENCY_SIMPLIFICATION_BLUEPRINT.md`, `docs/adr/ADR_001_concurrency_model.md`.
- [DEPENDENCY] v1.7.0 implementation starts only after v1.6.5 sampler hardening is accepted (v1.6.5 P0 complete; remaining v1.6.5 P1 items tracked/triaged).

### v1.7.0 P0

- [P0] [Owner: TBD] Introduce bounded sample scheduler abstraction (single long-lived executor + `max_pending` backpressure).
  Acceptance: scheduler enforces in-flight window and blocks/awaits when full without short-timeout polling loops.
  Dependency: requires v1.6.5 submit-contract/failure-policy stabilization.

- [P0] [Owner: TBD] Remove nested execution submission in `modulePool` runtime path.
  Acceptance: no runtime `executor.submit(...).result()` nesting inside per-sample module execution path.
  Dependency: requires v1.6.5 module execution correctness and failure handling checks.

- [P0] [Owner: TBD] Standardize sampler collection loops on bounded wait/as_completed behavior.
  Acceptance: no busy-wait loops in first-party sampler hot paths; all exceptions consumed and surfaced with sample UUID context.

### v1.7.0 P1

- [P1] [Owner: TBD] Add deterministic scheduling diagnostics (submission/completion counters, pending window, exception summary, saved paths).
  Acceptance: logs include start/end + submitted/completed/failed/cancelled + explicit output path messages.

- [P1] [Owner: TBD] Add resource-stability hooks (FD checkpoints + idle-spin warning hints).
  Acceptance: periodic FD growth checks at run start/milestones/shutdown, with clear warnings on abnormal drift.

- [P1] [Owner: TBD] Align dynesty/TPMCMC special paths with bounded scheduling policy.
  Acceptance: documented and tested policy for dynesty nested pooling interactions and TPMCMC exchange pause behavior.

### v1.7.0 P2

- [P2] [Owner: TBD] Add short micro-bench + regression checks for concurrency simplification.
  Acceptance: throughput regression target and idle CPU target are measured in short controlled runs (no long benchmark jobs).

- [P2] [Owner: TBD] Evaluate optional guarded module-level parallel mode (off by default).
  Acceptance: explicit opt-in, bounded, and compatibility-tested; default remains sample-level parallelism only.

## B. P0 (Before PyPI Production Release)

- [P0] Align package versioning with release policy.
  Acceptance: `pyproject.toml` version matches the intended tag/release version.
  Notes: first PyPI publication target is `1.6.0`; keep GitHub `v1.5` as pre-PyPI baseline.

- [P0] Finalize PyPI metadata and dependency boundaries.
  Acceptance: `pip install Jarvis-HEP` from clean env works; `Jarvis --help` succeeds; dependency list is justified and minimal.

- [P0] Create reproducible release checklist (build/test/publish/verify).
  Acceptance: one documented checklist covers TestPyPI and PyPI paths with explicit rollback plan.

- [P0] Create GitHub release note template and changelog routine.
  Acceptance: every tag has corresponding release note content source file under `docs/releases/`.

- [P0] Clean package payload.
  Acceptance: wheel/sdist do not include transient artifacts (`dist/`, `*.egg-info`, `__pycache__`, `.DS_Store`, local results).

## C. P1 (Stability And Engineering Quality)

- [P1-1] Expand automated tests around CLI modes.
  Acceptance: test coverage for `--mkproject`, `--convert`, and argument conflict behavior.

- [P1-0] Complete package namespace migration from `src` to `jarvishep`.
  Acceptance: runtime no longer depends on implicit `src` top-level imports or `sys.path` insertion hacks.

- [P1-2] Expand tests for execution path.
  Acceptance: smoke tests for `Core -> Workflow -> ModuleManager` with minimal fake modules.

- [P1-3] Add tests for `Likelihood/PassCondition` scalar contract.
  Acceptance: explicit cases for scalar success and array-output failure behavior.

- [P1-4] Add CI pipeline.
  Acceptance: PR/branch checks run unit tests + package build + basic import smoke checks.

- [P1-5] Add schema migration guardrails.
  Acceptance: schema `version` bump policy and migration helper exist, with backward-compatible defaults.

- [P1-6] Improve CSV conversion scalability.
  Acceptance: conversion can stream/iterate for larger HDF5 snapshots without loading all rows into memory at once.

## D. P2 (Architecture Refactor Candidates)

- [P2-1] Remove singleton state from `ModuleManager`.
  Acceptance: manager lifecycle is explicit and isolated per run, with no cross-run hidden state.

- [P2-2] Formalize module output conflict policy.
  Acceptance: if two modules write same observable key in one layer, behavior is deterministic and documented.

- [P2-3] Separate old experimental code paths.
  Acceptance: legacy code under `Sampling/Old` and embedded third-party sources are isolated with clear maintenance status.

- [P2-4] Add typed contracts for observables and module IO.
  Acceptance: key runtime interfaces have type hints + runtime validation where needed.

## E. Documentation Backlog

- [DOC] Add a user-facing flattening cookbook with practical examples (array, dict, mixed).
- [DOC] Add a “from YAML to output files” execution lifecycle diagram.
- [DOC] Add maintainer release playbook linked from `README.md`.

## F. Working Rules For Updating This File

- Keep completed items in section A; do not delete historical milestones.
- New tasks must include acceptance criteria.
- P0 items are blockers for production release.
- Update date at top on each maintenance pass.
