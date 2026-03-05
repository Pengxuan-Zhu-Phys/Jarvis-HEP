# Tasks: MCMC All-in-One (Documentation-First Roadmap)

Last updated: 2026-03-05  
Canonical design baseline: `docs/MCMC_STATE_MACHINE_DESIGN.md`

## Phase 1: Docs (this PR)

- [x] Task P1-DOC-01
  - Objective: update canonical design doc with all-in-one roadmap section and DoD template.
  - Files to add/modify: `docs/MCMC_STATE_MACHINE_DESIGN.md`
  - YAML schema keys to add: none (docs-only phase)
  - Minimal runnable example YAML (later): `bin/EggBox/Example_ToyMCMC_Toy.yaml`
  - Minimal test case name: `test_mcmc_doc_contracts_reference`
  - DoD checklist: canonical section 11 exists; phased roadmap + catalog + DoD template documented; no code changes.

- [x] Task P1-DOC-02
  - Objective: publish PRD, architecture context, algorithm catalog, test plan, and this task file.
  - Files to add/modify:
    - `docs/mcmc/PRD_MCMC_AllInOne.md`
    - `docs/mcmc/ARCH_MCMC_StateMachine.md`
    - `docs/mcmc/ALGO_CATALOG.md`
    - `docs/mcmc/TEST_PLAN.md`
    - `docs/mcmc/TASKS_MCMC_AllInOne.md`
  - YAML schema keys to add: none (docs-only phase)
  - Minimal runnable example YAML (later): `bin/EggBox/Example_TPMCMC_Toy.yaml`
  - Minimal test case name: `test_mcmc_docs_presence`
  - DoD checklist: all docs exist; cross-links to canonical doc included; terminology consistent.

## Phase 2.0: Infra Polish

- [ ] Task P2.0-INFRA-01
  - Objective: formalize engine protocol typing/validation and runtime adapter checks.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/state_machine_base.py`
    - `jarvishep/Sampling/Source/MCMC/chain_runtime.py`
    - `jarvishep/Sampling/Source/MCMC/mcmc_chain.py`
    - `tests/test_mcmc_state_machine.py`
  - YAML schema keys to add: none
  - Minimal runnable example YAML (later): `bin/EggBox/Example_ToyMCMC_Toy.yaml`
  - Minimal test case name: `test_engine_protocol_contract`
  - DoD checklist: protocol checks in place; invalid engine contract fails explicitly; no deadlock regression.

- [ ] Task P2.0-INFRA-02
  - Objective: standardize `ChainHistory` event metadata schema and metrics field naming.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/chain_history.py`
    - `jarvishep/Sampling/Source/MCMC/metrics_bus.py`
    - `docs/MCMC_STATE_MACHINE_DESIGN.md`
    - `tests/test_mcmc_state_machine.py`
  - YAML schema keys to add: none
  - Minimal runnable example YAML (later): `bin/EggBox/Example_ToyMCMC_Toy.yaml`
  - Minimal test case name: `test_chain_history_meta_standard`
  - DoD checklist: event schema stable; metrics names consistent; backward compatibility note documented.

- [ ] Task P2.0-INFRA-03
  - Objective: normalize shared MCMC schema keys and defaults across samplers.
  - Files to add/modify:
    - `jarvishep/card/schema/mcmc_schema.json`
    - `jarvishep/card/schema/tpmcmc_schema.json`
    - `jarvishep/config.py`
    - `docs/mcmc/PRD_MCMC_AllInOne.md`
  - YAML schema keys to add:
    - shared normalization for `num_chains`, `num_iters`, `selection`, proposal-scale naming
  - Minimal runnable example YAML (later): `bin/EggBox/Example_TPMCMC_Toy.yaml`
  - Minimal test case name: `test_mcmc_schema_normalization`
  - DoD checklist: normalized keys accepted; invalid combinations rejected with clear diagnostics.

## Phase 2.1: AM / RobustAM / DRAM

- [x] Task P2.1-AM-01
  - Objective: add AMMCMC engine and thin sampler entry.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/ammcmc_chain.py`
    - `jarvishep/Sampling/ammcmc.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/ammcmc_schema.json`
  - YAML schema keys to add:
    - `adapt.enabled`, `adapt.start_iter`, `adapt.window`, `adapt.scale`, `adapt.eps`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_AMMCMC_Operas.yaml`
  - Minimal test case name: `test_ammcmc_sampler_smoke_closes_samples`
  - DoD checklist: runs on Gaussian toy; acceptance and covariance sanity metrics available; no deadlock.
  - Tracking doc: `docs/mcmc/samplers/AMMCMC_TASKLIST.md`

- [x] Task P2.1-RAM-02
  - Objective: add RobustAM engine (AM + global/heavy-tail jumps).
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/robustam_chain.py`
    - `jarvishep/Sampling/robustam.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/robustam_schema.json`
  - YAML schema keys to add:
    - AM keys + `global_jump_prob`, `heavy_tail_df`, `global_scale`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_RobustAM_Operas.yaml`
  - Minimal test case name: `test_robustam_sampler_smoke_closes_samples`
  - DoD checklist: improved mode-jump behavior on mixture toy; bounded in-flight preserved.
  - Tracking doc: `docs/mcmc/samplers/RobustAM_TASKLIST.md`

- [x] Task P2.1-DRAM-03
  - Objective: add DRAM engine with delayed rejection stages.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/state_machine_multistage_base.py`
    - `jarvishep/Sampling/Source/MCMC/dram_chain.py`
    - `jarvishep/Sampling/dram.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/dram_schema.json`
  - YAML schema keys to add:
    - AM keys + `dr_steps`, `dr_scale_factors`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_DRAM_Operas.yaml`
  - Minimal test case name: `test_dram_sampler_smoke_closes_samples`
  - DoD checklist: full stage-2 delayed-rejection acceptance correction implemented; no lifecycle regression; async fail-path safe.
  - Tracking doc: `docs/mcmc/samplers/DRAM_TASKLIST.md`

## Phase 2.2: DEMCMC / DREAM-lite / Ensemble / PT-Ensemble

- [x] Task P2.2-DE-01
  - Objective: add DEMCMC population engine.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/engine_demcmc.py`
    - `jarvishep/Sampling/demcmc.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/DEMCMC_schema.json`
  - YAML schema keys to add:
    - `de_gamma`, `de_noise`, `de_crossover`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_DEMCMC_Operas.yaml`
  - Minimal test case name: `test_demcmc_mixture_toy`
  - DoD checklist: population diversity metrics available; no queue exhaustion.
  - Tracking doc: `docs/mcmc/samplers/DEMCMC_TASKLIST.md`

- [x] Task P2.2-DREAM-02
  - Objective: add DREAM-lite enhanced DE variant.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/engine_dream_lite.py`
    - `jarvishep/Sampling/dream_lite.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/DREAMLite_schema.json`
  - YAML schema keys to add:
    - DEMCMC keys + `dream_snooker_prob`, `dream_archive_size`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_DREAMLite_Operas.yaml`
  - Minimal test case name: `test_dream_lite_periodic_ridge_toy`
  - DoD checklist: move-type metrics logged; async stability verified.
  - Tracking doc: `docs/mcmc/samplers/DREAMLite_TASKLIST.md`

- [x] Task P2.2-ENS-03
  - Objective: add EnsembleMCMC and optional PT-Ensemble.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/engine_ensemble.py`
    - `jarvishep/Sampling/ensemblemcmc.py`
    - `jarvishep/Sampling/pt_ensemble.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/Ensemble_schema.json`
    - `jarvishep/card/schema/PTEnsemble_schema.json`
  - YAML schema keys to add:
    - `stretch_a`, `move_mix`, and optional tempering keys
  - Minimal runnable example YAML (created):
    - `bin/EggBox/Example_EnsembleMCMC_Operas.yaml`
    - `bin/EggBox/Example_PTEnsemble_Operas.yaml`
  - Minimal test case name: `test_ensemble_affine_toy`
  - DoD checklist: walker metrics + optional swap metrics validated; no deadlock.
  - Tracking docs:
    - `docs/mcmc/samplers/EnsembleMCMC_TASKLIST.md`
    - `docs/mcmc/samplers/PTEnsemble_TASKLIST.md`

## Phase 2.3: Slice / ESS

- [x] Task P2.3-SLICE-01
  - Objective: add SliceMCMC (univariate and random-direction modes).
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/engine_slice.py`
    - `jarvishep/Sampling/slicemcmc.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/SliceMCMC_schema.json`
  - YAML schema keys to add:
    - `slice.mode`, `slice.width`, `slice.max_steps_out`, `slice.max_shrink`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_SliceMCMC_Operas.yaml`
  - Minimal test case name: `test_slice_banana_toy`
  - DoD checklist: bracket-step metrics available; bounded in-flight preserved.
  - Tracking doc: `docs/mcmc/samplers/SliceMCMC_TASKLIST.md`

- [x] Task P2.3-ESS-02
  - Objective: add ESS with prior-structure validation.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/engine_ess.py`
    - `jarvishep/Sampling/ess.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/ESS_schema.json`
  - YAML schema keys to add:
    - `ess.enabled`, `ess.prior_cov`, `ess.angle_bracket_init`
  - Minimal runnable example YAML (created): `bin/EggBox/Example_ESS_Operas.yaml`
  - Minimal test case name: `test_ess_gaussian_prior_toy`
  - DoD checklist: prior checks fail clearly when unsupported; no lifecycle regressions.
  - Tracking doc: `docs/mcmc/samplers/ESS_TASKLIST.md`

## Phase 2.4: Gradient-provider Hooks + MALA/HMC/NUTS Placeholder

- [ ] Task P2.4-GRAD-01
  - Objective: define gradient-provider interface and fallback policy.
  - Files to add/modify:
    - `jarvishep/Sampling/Source/MCMC/gradient_provider.py`
    - `jarvishep/Sampling/Source/MCMC/state_machine_base.py`
    - `docs/MCMC_STATE_MACHINE_DESIGN.md`
    - `jarvishep/card/schema/gradient_hooks_schema.json`
  - YAML schema keys to add:
    - `grad.provider`, `grad.clip`, `grad.timeout`, `grad.fallback`
  - Minimal runnable example YAML (later): `bin/EggBox/Example_GradHooks_Toy.yaml`
  - Minimal test case name: `test_gradient_provider_contract`
  - DoD checklist: interface only, no forced gradient dependency, clear fallback behavior.

- [x] Task P2.4-MALAHMC-02
  - Objective: add placeholder sampler entries for MALA/HMC/NUTS on top of gradient interface.
  - Files to add/modify:
    - `jarvishep/Sampling/mala.py`
    - `jarvishep/Sampling/hmc.py`
    - `jarvishep/Sampling/nuts.py`
    - `jarvishep/distributor.py`
    - `jarvishep/card/schema/MALA_schema.json`
    - `jarvishep/card/schema/HMC_schema.json`
    - `jarvishep/card/schema/NUTS_schema.json`
  - YAML schema keys to add:
    - `mala.step_size`, `hmc.step_size`, `hmc.leapfrog_steps`, `nuts.max_depth`
  - Minimal runnable example YAML (created):
    - `bin/EggBox/Example_MALA_Operas.yaml`
    - `bin/EggBox/Example_HMC_Operas.yaml`
    - `bin/EggBox/Example_NUTS_Operas.yaml`
  - Minimal test case name: `test_mala_hmc_nuts_placeholder_configs`
  - DoD checklist: config/schema and interfaces compile/run in placeholder mode; no deadlock; lifecycle invariants preserved.
  - Tracking docs:
    - `docs/mcmc/samplers/MALA_TASKLIST.md`
    - `docs/mcmc/samplers/HMC_TASKLIST.md`
    - `docs/mcmc/samplers/NUTS_TASKLIST.md`
