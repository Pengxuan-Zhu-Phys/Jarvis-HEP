#!/usr/bin/env python3 

import time, sys, os

# basic numeric setup
import numpy as np
from numpy import linalg

# plotting
import matplotlib
from matplotlib import pyplot as plt

# seed the random number generator
rstate = np.random.default_rng(4232)

import dynesty
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt 
# Define the dimensionality of our problem.

# Define our 3-D correlated multivariate normal log-likelihood.
ndim = 3  # number of dimensions
C = np.identity(ndim)  # set covariance to identity matrix
C[C==0] = 0.95  # set off-diagonal terms (strongly correlated)
Cinv = linalg.inv(C)  # precision matrix
lnorm = -0.5 * (np.log(2 * np.pi) * ndim + np.log(linalg.det(C)))  # ln(normalization)

# 3-D correlated multivariate normal log-likelihood
def loglikelihood(x):
    """Multivariate normal log-likelihood."""
    return -0.5 * np.dot(x, np.dot(Cinv, x)) + lnorm

# prior transform
def prior_transform(u):
    """Transforms our unit cube samples `u` to a flat prior between -10. and 10. in each variable."""
    return 10. * (2. * u - 1.)

# sampler = dynesty.NestedSampler(loglikelihood, prior_transform, ndim,
#                                 bound='single', nlive=1000, rstate=rstate)

# # sample from the distribution
# sampler.run_nested(dlogz=0.01)

# # grab our results
# res = sampler.results
# print(res.niter , res.nlive)



# Sample from our distribution.
dsampler = dynesty.DynamicNestedSampler(loglikelihood, prior_transform, ndim, 
                                        bound='single', sample='unif',rstate=rstate)
dsampler.reset()
# dsampler.run_nested(nlive_init=50, nlive_batch=50,
#                     maxiter=14275, use_stop=False, 
#                     wt_kwargs={'pfrac': 0.0})
# dres_z = dsampler.results

dsampler.run_nested(nlive_init=50, nlive_batch=50,
                    maxiter=14275, use_stop=False, 
                    wt_kwargs={'pfrac': 1.0})
dres_p = dsampler.results

# fig, axes = plt.subplots(2, 5, figsize=(12, 3))
# axes = axes.reshape((2, 5))
# fg, ax = dyplot.cornerpoints(dres_p, cmap='viridis', truths=np.zeros(ndim),
#                              fig=(fig, axes[:, 0:2]))

# fig, axes = dyplot.traceplot(dres_z, truths=np.zeros(ndim), show_titles=True, trace_cmap='plasma',
#                              title_kwargs={'fontsize': 28, 'y': 1.05}, quantiles=None,
#                              fig=plt.subplots(3, 2, figsize=(12, 10)))
# fig.tight_layout()

fig, axes = dyplot.traceplot(dres_p, truths=np.zeros(ndim), show_titles=True, trace_cmap='viridis',
                             title_kwargs={'fontsize': 28, 'y': 1.05}, quantiles=None,
                             fig=plt.subplots(3, 2, figsize=(12, 10)))
fig.tight_layout()

# plot evidence-oriented dynamic run (right)
# fg, ax = dyplot.cornerpoints(dres_z, cmap='inferno', truths=np.zeros(ndim),
#                              fig=(fig, axes[:, 3:5]))

plt.show()