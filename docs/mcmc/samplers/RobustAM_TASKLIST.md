# RobustAM Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `RobustAM` sampler class using shared MCMC state-machine runtime with AM adaptation + heavy-tail/global-jump proposal mixture.

## Task Checklist

- [x] Add `RobustAMChain` engine.
- [x] Add `RobustAM` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "RobustAM"`).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/robustam_chain.py`
  - `jarvishep/Sampling/robustam.py`
  - `bin/EggBox/Example_RobustAM_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `tests/test_mcmc_state_machine.py`
  - `tests/test_sampler_hardening_smoke.py`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py` -> `11 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "mcmc or tpmcmc or ammcmc or robustam or distributor"` -> `5 passed`
  - `pytest -q` -> `89 passed`
