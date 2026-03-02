# Jarvis-HEP Todo And Roadmap

Last updated: 2026-03-02
Scope: PyPI publication prep, release hardening, maintainability upgrades.
Active release track: `docs/releases/v1.6.1_TASKLIST.md`

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
