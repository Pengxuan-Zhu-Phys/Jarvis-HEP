# Vendored dynesty

- **Upstream version**: 3.0.0 (PyPI wheel `dynesty-3.0.0-py3-none-any.whl`)
- **Purpose**: Jarvis-HEP V2 D13.5 (route A — in-tree fork with UUID channel)
- **License**: MIT (see LICENSE / AUTHORS)
- **Patches**: See `PATCHES.md` (Jarvis Sample UUID channel)
- **Pool**: `jarvishep2.Sampling.redis_evaluation_pool.RedisEvaluationPool`

Import:

```python
from jarvishep2.Sampling.Source.Dynesty.py import dynesty
```

Never import from V1 `jarvishep`.
