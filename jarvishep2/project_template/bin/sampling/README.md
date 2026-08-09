# Nested sampling YAML templates

Copy a `Sampling_*.yaml` **Sampling:** block into your task card. Replace
`Variables` / `LogLikelihood` with your model.

## Method = engine (no `Bounds.dynamic`)

| Method | Engine | Use when |
|--------|--------|----------|
| **Dynesty** | always `DynamicNestedSampler` | evidence + dynamic batches / ESS stop |
| **MultiNest** | always static `NestedSampler` | classic static nested (not Fortran MultiNest) |

Validation rejects `Bounds.dynamic` / `Dynamic`. Switch Method instead of a flag.

## Files

| File | Level | Notes |
|------|--------|--------|
| `Sampling_Dynesty_Simple.yaml` | everyday | few knobs |
| `Sampling_Dynesty_Full.yaml` | full | Dynamic `run_nested` + constructor surface |
| `Sampling_MultiNest_Simple.yaml` | everyday | few knobs |
| `Sampling_MultiNest_Full.yaml` | full | static `run_nested` + constructor only |

## Shared everyday surface (Simple)

```yaml
Sampling:
  Method: Dynesty   # or MultiNest
  Variables: [...]
  Bounds:
    seed: 21
    nlive: 100
    dlogz: 0.5        # Dynesty → dlogz_init; MultiNest → static dlogz
    run_nested:
      print_progress: true
  LogLikelihood: [...]
```

Optional proposal knobs under `Bounds` (both methods): `bound`, `sample`, `walks`, …

## Evidence keys

| Method | Prefer | Meaning |
|--------|--------|---------|
| MultiNest | `Bounds.dlogz` | static stop |
| Dynesty | `Bounds.dlogz` | remapped to `run_nested.dlogz_init` for the initial phase |

Full Dynesty cards may also set `run_nested.dlogz_init` / `nlive_batch` / `n_effective` / …

## Do not put in YAML

HEP injects: `loglikelihood`, `prior_transform`, `ndim`, `pool`, `rstate`.  
Also omit: `print_func`, `use_pool`, `logl_args` / `ptform_args` (not user surface).

## Products

| Method | Nested CSV |
|--------|------------|
| Dynesty | `DATABASE/dynesty_result.csv` |
| MultiNest | `DATABASE/multinest_result.csv` |

Plus `samples.hdf5` / `sampler_summary.json` for both.
