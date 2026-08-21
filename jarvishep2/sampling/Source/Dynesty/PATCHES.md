# Patches vs stock dynesty 3.1.0

Vendored from PyPI `dynesty==3.1.0`, then modified for the **Jarvis Sample UUID channel**
(V1 route A — same idea as `Jarvis-HEP` `Sampling/Source/Dynesty`).

Upgrade history:

| Vendor | Notes |
|--------|--------|
| 3.0.0 | Initial D13.5 vendor + UUID / logging patches |
| **3.1.0** | Upstream refresh (2026-07); re-applied Jarvis patches; print/utils API follows 3.1 |

## UUID channel contract

1. User `prior_transform(u)` may return `np.append(coords, uuid_str)`.
2. Live points store coords in `live_v` and uuid in `live_uid` (`U36`).
3. `loglikelihood` is called with the **full** vector including trailing uuid
   (so `NestedLikelihoodBridge` can stamp `Sample.uuid`).
4. Dead-point records store `uid`; results expose `samples_uid`.

## Files touched

| File | Change |
|------|--------|
| `jarvis_uuid.py` | **new** — split/join helpers |
| `jarvis_logging.py` | **new** — Jarvis logger bridge (progress/warnings) |
| `sampler.py` | live init split; `live_uid`; dead-point uid; `samples_uid`; progress via Jarvis logger |
| `dynamicsampler.py` | progress via `emit_progress` + logger bind; UUID unpack of `_initialize_live_points` (u,v,uid,logl,blobs) in `sample_initial` + batch live init |
| `utils.py` | `IteratorResult`/`Short` `uidstar`; `RunRecord['uid']`; `samples_uid`; `print_fn` → Jarvis logger |
| `pool.py` | companion pools live in `jarvishep2.Sampling.redis_evaluation_pool` |
| `__init__.py` | vendored `__version__ = "3.1.0+jarvishep2"`; install warnings bridge |

## Upstream 3.1 highlights absorbed

- `_parse_pool_queue` moved into `utils` (used by static/dynamic factories)
- Progress / checkpoint helpers refined in `utils` (we still prefer Jarvis logger over tqdm)
- Minor pool/queue_size and sampler docstring/f-string cleanups

## MultiNest

V2 `MultiNest` is **not** Fortran MultiNest; it is static `NestedSampler` on this same
vendored tree (`multinest_sampler.py` forces `dynamic=False`). Updating this package
updates both Dynesty and MultiNest engines.
