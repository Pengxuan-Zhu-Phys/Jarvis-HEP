# SliceMCMC Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `SliceMCMC` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `SliceChain` engine.
- [x] Add `SliceMCMC` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "SliceMCMC"`).
- [x] Add sampler schema entry (`SliceMCMC_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_slice.py`
  - `jarvishep/Sampling/slicemcmc.py`
  - `jarvishep/card/schema/SliceMCMC_schema.json`
  - `bin/EggBox/Example_SliceMCMC_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "slice or ptensemble or distributor"` -> `11 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "slicemcmc or ptensemble or ensemblemcmc"` -> `3 passed`
