# DREAM-lite Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `DREAMLite` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `DREAMLiteChain` engine (`DE` + `snooker` mixed moves).
- [x] Add `DREAMLite` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "DREAMLite"`).
- [x] Add sampler schema entry (`DREAMLite_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_dream_lite.py`
  - `jarvishep/Sampling/dream_lite.py`
  - `jarvishep/card/schema/DREAMLite_schema.json`
  - `bin/EggBox/Example_DREAMLite_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "dream or demcmc or distributor"` -> `8 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "dream_lite or demcmc or dram"` -> `3 passed`
