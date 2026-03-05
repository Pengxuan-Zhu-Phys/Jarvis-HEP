# PT-Ensemble Sampler Task List (Jarvis-HEP)

Last updated: 2026-03-05  
Status: Completed

## Scope

Add independent `PTEnsemble` sampler class on shared state-machine runtime.

## Task Checklist

- [x] Add `PTEnsemble` sampler class (tempered + stretch-move walkers).
- [x] Reuse existing tempering exchange hooks from `TPMCMC`.
- [x] Add Distributor method route (`Sampling.Method: "PTEnsemble"`).
- [x] Add sampler schema entry (`PTEnsemble_schema.json` + preference mapping).
- [x] Add EggBox Operas quick YAML.
- [x] Add minimal unit/smoke tests.
- [x] Git add + commit (no push).

## Files Added/Modified

- Added:
  - `jarvishep/Sampling/pt_ensemble.py`
  - `jarvishep/card/schema/PTEnsemble_schema.json`
  - `bin/EggBox/Example_PTEnsemble_Operas.yaml`
- Modified:
  - `jarvishep/distributor.py`
  - `jarvishep/card/preference.json`

## Validation Snapshot

- Executed:
  - `pytest -q tests/test_mcmc_state_machine.py -k "ptensemble or ensemble or distributor"` -> `10 passed`
  - `pytest -q tests/test_sampler_hardening_smoke.py -k "ptensemble or ensemblemcmc or dream_lite"` -> `3 passed`
