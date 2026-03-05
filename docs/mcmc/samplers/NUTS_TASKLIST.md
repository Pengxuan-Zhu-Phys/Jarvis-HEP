# NUTS Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `NUTS` sampler class with placeholder stochastic tree-depth proposals on shared runtime.

## Task Checklist

- [x] Add `NUTSChain` engine.
- [x] Add `NUTS` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "NUTS"`).
- [x] Add sampler schema entry (`NUTS_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "mala or hmc or nuts or distributor"` -> `16 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "mala_sampler or hmc_sampler or nuts_sampler or ess_sampler"` -> `4 passed`
