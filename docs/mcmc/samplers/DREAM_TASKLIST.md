# DREAM Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Upgrade `DREAMLite` implementation to full `DREAM` sampler while keeping `DREAMLite` backward-compatible entry.

## Task Checklist

- [x] Add `DREAMChain` engine (`DE + snooker + archive + adaptive crossover weights`).
- [x] Add `DREAM` sampler class (`Sampling.Method: "DREAM"`).
- [x] Keep `DREAMLite` as compatibility profile on top of `DREAM` runtime.
- [x] Add Distributor route for `DREAM` and keep `DREAMLite` route.
- [x] Add schema mapping (`DREAM_schema.json`).
- [x] Add EggBox Operas quick YAML for `DREAM`.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_dream.py`
  - `jarvishep/Sampling/dream.py`
  - `jarvishep/card/schema/DREAM_schema.json`
  - `bin/EggBox/Example_DREAM_Operas.yaml`
- Modified:
  - `jarvishep/Sampling/dream_lite.py`
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`
  - `tests/test_mcmc_state_machine.py`
  - `tests/test_sampler_hardening_smoke.py`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "dream or demcmc or distributor"` -> `17 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "dream_sampler or dream_lite or demcmc"` -> `3 passed`
