# DRAM Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `DRAM` sampler class on shared state-machine runtime.
Current delivery is DRAM-lite for v1.7.0 phase: adaptive MH with rejection-triggered proposal scaling (`dr_scale_factors`) as delayed-rejection behavior placeholder.

## Task Checklist

- [x] Add `DRAMChain` engine.
- [x] Add `DRAM` sampler class.
- [x] Add Distributor method route (`Sampling.Method: "DRAM"`).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/dram_chain.py`
  - `jarvishep/Sampling/dram.py`
  - `bin/EggBox/Example_DRAM_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `tests/test_mcmc_state_machine.py`
  - `tests/test_sampler_hardening_smoke.py`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py` -> `13 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "mcmc or tpmcmc or ammcmc or robustam or dram or distributor"` -> `6 passed`
  - `pytest -q` -> `92 passed`
