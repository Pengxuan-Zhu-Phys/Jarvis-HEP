"""Latin hypercube design: one point per stratum in every variable.

Shipped with Jarvis-HEP as ``Sampling.Method: LatinHypercube``. It uses only
the public Jarvis-HEP-Sampler API, so it runs through the same path as a user
plug-in. The same file is the SDK's ``examples/LatinHypercube.py``.
"""

import numpy as np

PLUGIN = {
    "api_version": 1,
    "method": "LatinHypercube",
    "version": "0.1.0",
    "profile": "independent",
    "description": "Latin hypercube design with point_number strata per variable.",
}

OPTIONS_SCHEMA = {
    "type": "object",
    "properties": {
        "centered": {
            "type": "boolean",
            "description": "Put each point at the centre of its stratum.",
        },
    },
    "additionalProperties": False,
}


def prepare(setup):
    rng = setup.rng()
    perms = np.stack([rng.permutation(setup.point_number) for _ in range(setup.dimension)])
    return {"perms": perms, "centered": bool(setup.options.get("centered", False))}


def draw(context, rng, plan):
    strata = plan["perms"][:, context.candidate_index]
    offset = 0.5 if plan["centered"] else rng.random(context.dimension)
    return (strata + offset) / context.point_number


def settings_rows(setup, plan):
    return [("centered", plan["centered"])]
