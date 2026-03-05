# HMC Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `HMC` sampler class with placeholder leapfrog dynamics on shared runtime.

## Task Checklist

- [x] Add `HMCChain` engine.
- [x] Add `HMC` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "HMC"`).
- [x] Add sampler schema entry (`HMC_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "mala or hmc or nuts or distributor"` -> `16 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "mala_sampler or hmc_sampler or nuts_sampler or ess_sampler"` -> `4 passed`
