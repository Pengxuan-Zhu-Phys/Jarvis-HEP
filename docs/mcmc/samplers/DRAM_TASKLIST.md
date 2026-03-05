# DRAM Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `DRAM` sampler class on shared state-machine runtime.
This delivery upgrades DRAM-lite to full delayed-rejection runtime:
- parallel multi-stage state-machine path for staged proposals
- stage-1 Metropolis accept/reject
- stage-2 delayed-rejection acceptance correction
- adaptive covariance inherited from AM core

## Task Checklist

- [x] Add parallel multi-stage runtime base for staged MCMC proposals.
- [x] Upgrade `DRAMChain` to staged delayed-rejection flow.
- [x] Upgrade `DRAM` sampler class to multi-stage runtime.
- [x] Add Distributor method route (`Sampling.Method: "DRAM"`).
- [x] Add EggBox Operas quick YAML.
- [x] Add stage-transition unit test + smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/Source/MCMC/state_machine_multistage_base.py`
  - `jarvishep/Sampling/Source/MCMC/dram_chain.py`
  - `jarvishep/Sampling/dram.py`
  - `bin/EggBox/Example_DRAM_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `tests/test_mcmc_state_machine.py`
  - `tests/test_sampler_hardening_smoke.py`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py` -> `14 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "dram or ammcmc or robustam or mcmc or tpmcmc"` -> `5 passed`
