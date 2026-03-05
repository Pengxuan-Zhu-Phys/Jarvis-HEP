# EnsembleMCMC Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `EnsembleMCMC` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `EnsembleChain` stretch-move engine.
- [x] Add `EnsembleMCMC` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "EnsembleMCMC"`).
- [x] Add sampler schema entry (`Ensemble_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_ensemble.py`
  - `jarvishep/Sampling/ensemblemcmc.py`
  - `jarvishep/card/schema/Ensemble_schema.json`
  - `bin/EggBox/Example_EnsembleMCMC_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "ensemble or dream or demcmc or distributor"` -> `10 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "ensemblemcmc or dream_lite or demcmc"` -> `3 passed`
