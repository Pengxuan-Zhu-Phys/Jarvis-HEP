# AMMCMC Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add an independent `AMMCMC` sampler class on top of shared `MCMCStateMachineBase`, with EggBox Operas runnable YAML and minimal regression tests.

## Task Checklist

- [x] Add AMMCMC engine (`AMMCMCChain`) for adaptive covariance proposals.
- [x] Add AMMCMC sampler class and wire method name.
- [x] Add Distributor route for `Sampling.Method: "AMMCMC"`.
- [x] Add EggBox Operas quick example YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/ammcmc_chain.py`
  - `jarvishep/Sampling/ammcmc.py`
  - `bin/EggBox/Example_AMMCMC_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `tests/test_mcmc_state_machine.py`
  - `tests/test_sampler_hardening_smoke.py`

## Validation Snapshot

- Executed checks:
  - `pytest -q tests/test_mcmc_state_machine.py` -> `9 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "mcmc or tpmcmc or ammcmc or distributor"` -> `4 passed`
  - `pytest -q` -> `86 passed`

## Notes

- AMMCMC keeps Factory and lifecycle invariants:
  - `submit_task(sample_info)` only
  - event-driven wait
  - `sample.close()` finally path
  - bounded in-flight submissions
