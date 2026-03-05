# DEMCMC Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `DEMCMC` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `DEMCMCChain` population proposal engine.
- [x] Add `DEMCMC` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "DEMCMC"`).
- [x] Add sampler schema entry (`DEMCMC_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/engine_demcmc.py`
  - `jarvishep/Sampling/demcmc.py`
  - `jarvishep/card/schema/DEMCMC_schema.json`
  - `bin/EggBox/Example_DEMCMC_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "demcmc or distributor or state_machine_smoke"` -> `9 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "demcmc or ammcmc or robustam or dram"` -> `4 passed`
