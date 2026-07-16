# Patches vs stock dynesty 3.0.0

Vendored from PyPI `dynesty==3.0.0`, then modified for the **Jarvis Sample UUID channel**
(V1 route A — same idea as `Jarvis-HEP` `Sampling/Source/Dynesty`).

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
| `dynamicsampler.py` | progress via `emit_progress` + logger bind |
| `utils.py` | `RunRecord['uid']`; `samples_uid`; IteratorResult `uidstar`; V1-style `print_fn` → logger |
| `pool.py` | companion pools live in `jarvishep2.Sampling.redis_evaluation_pool` |

## Not fully ported from V1 2.1.3 fork

DynamicNestedSampler batch paths in `dynamicsampler.py` still mostly stock;
`uidstar` defaults keep them importable. Full dynamic-batch uuid parity can be
tightened under D13.6 if needed. Static `NestedSampler` + Dynamic main path
init go through the patched live-point code.
