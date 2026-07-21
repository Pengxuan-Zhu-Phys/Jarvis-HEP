# Nested sampling YAML templates (Dynesty / MultiNest)

Copy a `Sampling_*.yaml` **Sampling:** block into your task card (or start from the
full fragment). Physics (`Variables` / `LogLikelihood`) are placeholders — replace
with your model.

| File | Method | Mode | Use when |
|------|--------|------|----------|
| `Sampling_Dynesty_Simple.yaml` | Dynesty | Dynamic (default) | laptop / smoke, few knobs |
| `Sampling_Dynesty_Full_Static.yaml` | Dynesty | Static NestedSampler | full constructor + static `run_nested` surface |
| `Sampling_Dynesty_Full_Dynamic.yaml` | Dynesty | DynamicNestedSampler | full dynamic batch / stop API |
| `Sampling_MultiNest_Simple.yaml` | MultiNest | Static only | smoke / production with few knobs |
| `Sampling_MultiNest_Full.yaml` | MultiNest | Static only | full static surface (same engine as Dynesty static) |

## Notes

- **MultiNest** in V2 is **not** Fortran MultiNest / pymultinest — it is static
  vendored `NestedSampler` on Redis Workers (same UUID + feedback path as Dynesty).
- **HEP injects** (do not set in YAML): `loglikelihood`, `prior_transform`, `ndim`,
  `pool`, `rstate`.
- Constructor kwargs: flat under `Bounds` and/or nested `Bounds.sampler:`.
- Run kwargs: always under `Bounds.run_nested:`.
- Evidence threshold: static uses `dlogz`; dynamic uses `dlogz_init` (a top-level
  `dlogz` is remapped to `dlogz_init` when `dynamic: true`).
- Products: `DATABASE/dynesty_result.csv` or `multinest_result.csv`, plus
  `samples.hdf5` / `sampler_summary.json`.
