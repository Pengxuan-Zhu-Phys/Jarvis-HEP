# Agent Private Notes (Internal)

Last updated: 2026-03-03
Purpose: keep implementation memory compact between maintenance turns.

## Current Priorities

- Close out v1.6.5 sampler hardening release docs and patch artifact delivery.
- Keep residual P1 items explicit (Grid boundary guard, deeper DNN regression tests).
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
- Grid endpoint `ppf/logit` boundary protection.
- full torch-chain DNN regression coverage.

## Next Chat Handoff (Sampler Cycle)

- Primary tracking doc:
  - `docs/releases/v1.6.5_SAMPLER_HARDENING_TASKLIST.md`

- v1.6.5 completed items:
  1. unified `WorkerFactory.submit_task(sample_info)` call contract and removed legacy interface
  2. restored base `evaluate_selection` / `check_evaluation` sampler contract
  3. fixed DNN gating/import issues and hardened DNN JSON run-info serialization
  4. enforced `Sample.close()` lifecycle in sampler paths
  5. set module failure policy to fail-fast by default with explicit continue option
  6. removed hot-path busy-wait futures polling
  7. fixed Dynesty precision (`float64`) and `samples_u` source mapping
  8. added sampler smoke/regression tests and generated patch artifact

- v1.6.5 residual follow-up:
  1. add deeper DNN torch-chain regression test path (beyond lightweight smoke)
  2. keep tracking CLI banner assertion drift separately from sampler hardening
  3. v1.6.5 P1 closure record (nested executor de-nesting + Grid endpoint guard) is tracked in `docs/releases/v1.6.5_REMAINING_P1_BREAKDOWN.md`

- Hot files likely touched first:
  - `jarvishep/factory.py`
  - `jarvishep/core.py`
  - `jarvishep/Sampling/randoms.py`
  - `jarvishep/Sampling/grid.py`
  - `jarvishep/Sampling/mcmc.py`
  - `jarvishep/Sampling/tpmcmc.py`
  - `jarvishep/Sampling/dnn.py`
  - `jarvishep/Sampling/dynesty.py`
  - `jarvishep/moduleManager.py`
  - `jarvishep/sample.py`
  - `tests/test_sampler_hardening_smoke.py`

- Non-negotiable checks per patch:
  - `python -m unittest discover -s tests -p 'test_*.py'`
  - `python -m build`
  - `Jarvis --help`
  - `python tests/test_sampler_hardening_smoke.py`
  - `python tests/test_execution_smoke.py`

- Avoid:
  - changing release/version metadata in same commit as sampler bugfixes
  - broad refactors before P0 contract bugs are closed
