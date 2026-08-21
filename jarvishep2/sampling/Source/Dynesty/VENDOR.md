# Vendored dynesty

- **Upstream version**: 3.1.0 (PyPI wheel `dynesty-3.1.0-py3-none-any.whl`)
- **Purpose**: Jarvis-HEP V2 nested sampling (route A — in-tree fork with UUID channel)
- **Used by**: `Sampling.Method: Dynesty` and `Sampling.Method: MultiNest` (static NestedSampler)
- **License**: MIT (see LICENSE / AUTHORS)
- **Patches**: See `PATCHES.md` (Jarvis Sample UUID channel + logging bridge)
- **Pool**: `jarvishep2.Sampling.redis_evaluation_pool.RedisEvaluationPool`

Import:

```python
from jarvishep2.Sampling.Source.Dynesty.py import dynesty
# dynesty.__version__ == "3.1.0+jarvishep2" when not installed as a site package
```

Never import from V1 `jarvishep` or from a system-wide `dynesty` package for science runs
(the UUID channel would be missing).
