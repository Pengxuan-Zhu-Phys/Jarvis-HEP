# Agent Private Notes (Internal)

Last updated: 2026-03-02
Purpose: keep implementation memory compact between maintenance turns.

## Current Priorities

- Keep decoupling intact (`task_root` vs `src_root`).
- Keep observable structured data support intact end-to-end.
- Keep release process reproducible and low-risk.
- Enforce version policy: GitHub-only `<=1.5.x`, PyPI starts from `1.6.0`.

## High-Value Quick Checks Before Any Non-Trivial Change

- `python -m unittest discover -s tests -p 'test_*.py'`
- `python -m build`
- `Jarvis --help`
- one quick convert smoke (`--convert`) on a small snapshot

## Behavioral Contracts To Remember

- `Likelihood`: output must be scalar-castable.
- `PassCondition`: output must be bool-castable.
- array/list/dict observables are legal for storage/export path.
- CSV flatten behavior is schema-controlled; respect user `name_map`.

## Local Risk Heuristics

- If edit touches `Base.decode_path`, re-check all marker semantics.
- If edit touches `observable_io`, run conversion tests immediately.
- If edit touches `moduleManager` or `modulePool`, reason about concurrency and shared state.
- Avoid sweeping refactor in same patch as release/versioning changes.

## Pending Technical Debt To Revisit

- version sync between tag/release/package metadata.
- broader test matrix beyond schema IO.
- reduce singleton/global runtime state.
- streaming conversion for very large HDF5 files.
