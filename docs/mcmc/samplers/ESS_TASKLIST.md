# ESS Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `ESS` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `ESSChain` engine.
- [x] Add `ESS` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "ESS"`).
- [x] Add sampler schema entry (`ESS_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_ess.py`
  - `jarvishep/Sampling/ess.py`
  - `jarvishep/card/schema/ESS_schema.json`
  - `bin/EggBox/Example_ESS_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "ess or slice or distributor"` -> `12 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "ess_sampler or slicemcmc or ptensemble"` -> `3 passed`
