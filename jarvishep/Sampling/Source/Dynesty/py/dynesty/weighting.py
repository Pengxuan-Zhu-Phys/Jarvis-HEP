#!/usr/bin/env python
from __future__ import annotations

import warnings

import numpy as np
from scipy.special import logsumexp


def compute_weights(results):
    """Derive evidence and posterior weights."""
    logl = results.logl
    logz = results.logz
    logvol = results.logvol
    logwt = results.logwt
    samples_n = results.samples_n

    if np.ptp(logz) == 0:
        warnings.warn(
            """The calculation of weights is seeing same
logz values associated with all the samples. It may mean somethings is
wrong with your likelihood."""
        )
        zweight = np.ones(len(logl)) / len(logl)
    else:
        logz_remain = logl[-1] + logvol[-1]
        logz_tot = np.logaddexp(logz[-1], logz_remain)
        lzones = np.ones_like(logz)
        logzin = logsumexp(
            [lzones * logz_tot, logz],
            axis=0,
            b=[lzones, -lzones],
        )
        logzweight = logzin - np.log(samples_n)
        logzweight -= logsumexp(logzweight)
        zweight = np.exp(logzweight)

    pweight = np.exp(logwt - logz[-1])
    pweight /= np.sum(pweight)
    return zweight, pweight
