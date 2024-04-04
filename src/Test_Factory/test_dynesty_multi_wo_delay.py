#!/usr/bin/env python3 

import time, sys, os
sys.path.append("/home/buding/Jarvis-HEP/src/Sampling/dynesty/py")
import numpy as np
from numpy import linalg

import matplotlib
from matplotlib import pyplot as plt

rstate = np.random.default_rng(736109)

# import dynesty
from dynesty import *
from dynesty import plotting as dyplot
from dynesty.results import print_fn
import matplotlib.pyplot as plt 
import math

import random
import pandas as pd 
from subprocess import Popen
from uuid import uuid4

tmax = 5.0 * np.pi
r = 1.3  # radius
w = 0.1  # width
c1 = np.array([-2.5, 0.])  # center of shell 1
c2 = np.array([2.5, 0.])  # center of shell 2
const = math.log(1. / math.sqrt(2. * math.pi * w**2))  # normalization constant

def loglike(x, **kwarg):
    t = 2.0 * tmax * x - tmax
    if t[0] > t[1]:
        return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0
    else:
        return 0.
    # print("\tMain fun Call for loglike, .... ", x, kwarg, ret)

# define the prior transform
def prior_transform(x):
    return x


# def logcirc(theta, c):
#     d = np.sqrt(np.sum((theta - c)**2, axis=-1))  # |theta - c|
#     return const - (d - r)**2 / (2. * w**2)

# # log-likelihood of two shells
# def loglike(theta):
#     return np.logaddexp(logcirc(theta, c1), logcirc(theta, c2))

# # our prior transform
# def prior_transform(x):
#     return 12. * x - 6.



# tstart = time.time()
# import dynesty.pool as dypool
# with dypool.Pool(2, loglike, prior_transform) as pool:
#     # The important thing that we provide the loglikelihood/prior transform from 
#     # the pool    
    # psampler = dynesty.DynamicNestedSampler(
        # pool.loglike, pool.prior_transform, 
        # ndim=2, bound='multi', sample='unif', rstate=rstate,
        # pool=pool)
    # psampler.run_nested(dlogz_init=0.01, nlive_init=320, nlive_batch=600,
                    #  wt_kwargs={'pfrac': 1.0}, stop_kwargs={'pfrac': 1.0}, checkpoint_file="./test.sav", checkpoint_every=1)
#     # psampler = dynesty.DynamicNestedSampler.restore("./test.sav", pool=pool)
#     # psampler.run_nested(resume=True, checkpoint_file="./test.sav", checkpoint_every=1)
# pres = psampler.results
# # print(dres[''])
dsampler = dynesty.DynamicNestedSampler(
    loglike, prior_transform, 
    ndim=2, bound='multi', sample='unif', rstate=rstate
    )
dsampler.run_nested(
    dlogz_init=0.01, 
    nlive_init=500, 
    nlive_batch=500,
    wt_kwargs={'pfrac': 0.0}, 
    stop_kwargs={'pfrac': 0.0}, 
    print_progress=True,
    checkpoint_file="test.sav", checkpoint_every=10
    )
# dsampler = dsampler.restore("test.sav")
# dsampler.run_nested(resume=True)

pres = dsampler.results


from copy import deepcopy
dres = deepcopy(pres)
dres = dict(dres)
print(dres.keys())
for kk, vv in dres.items():
    # print(kk, vv.shape)
    if type(vv) is np.ndarray:
        print(kk, vv.shape)
        print(kk, vv[0])

# # dlogz_final = 0.1
# ncall = dsampler.ncall 
# nit = dsampler.it
# for it, results in enumerate(dsampler.sample_initial(dlogz=dlogz_final)):
#     # split up our results
#     (worst, ustar, vstar, loglstar, logvol, logwt, logz, logzvar,
#      h, nc, worst_it, boundidx, bounditer, eff, delta_logz, blob) = results
#     # add number of function calls
#     ncall += nc
#     nit += 1
#     # print results
#     # print_fn(results, nit, ncall, dlogz=dlogz_final)
#     # print(worst)
#     # print
#     # df = 
#     # print(dsampler.results.samples)
#     print("\n", dsampler.results.samples.shape, dsampler.results.samples_id[-1], "\n")
#     time.sleep(5)
#     # print(dsampler.results.samples_id[-1])
# for kk, vv in dict(dsampler.results).items():
#     print(kk, vv)
# # add the remaining live points back into our final results 
# # (they are removed from our set of dead points each time we start sampling)
# for it2, results in enumerate(dsampler.sampler.add_live_points()):
#     # split up results
#     (worst, ustar, vstar, loglstar, logvol, logwt, logz, logzvar,
#      h, nc, worst_it, boundidx, bounditer, eff, delta_logz,blob) = results
#     # print results
#     print_fn(results, nit, ncall, add_live_it=it2+1, dlogz=dlogz_final)

# dres_z = dsampler.results
# print(dres_z)

# fig, axes = dyplot.cornerplot(pres, quantiles=None, color='darkviolet',
#                               span=([-6, 6], [-6, 6]),
#                               fig=plt.subplots(2, 2, figsize=(10, 10)))
# fig, axes = dyplot.cornerplot(pres, quantiles=None, color='darkviolet',
#                               span=([0, 1], [0, 1]),
#                               fig=plt.subplots(2, 2, figsize=(10, 10)))
# plt.show()
# # plt.cla()
# lnz_truth = -2 * np.log(10. * 0.999999426697)

# fig, axes = dyplot.runplot(pres, color='navy', logplot=True,
#                            lnz_truth=lnz_truth, truth_color='black')
# plt.show()

# fig, axes = dyplot.traceplot(pres, truths=np.zeros(2), show_titles=True, trace_cmap="plasma", quantiles=None)
# plt.show()

# print(dict(pres).items())
fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.1, 0.1, 0.88, 0.88])
ax.plot([0., 1.], [0., 1.], '-')

# fg, ax = dyplot.cornerpoints(pres, cmap='plasma', truths=np.zeros(ndim),
#                              fig=(fig, axes[:, 0:2]))
fg, ax = dyplot.cornerpoints(pres, cmap='plasma', fig=(fig, ax), plot_kwargs={"marker": 'o', "s": 6})
ax[0, 0].set_title('No Bound', fontsize=26)
ax[0, 0].set_xlim([0, 1])
ax[0, 0].set_ylim([0, 1])
plt.show()