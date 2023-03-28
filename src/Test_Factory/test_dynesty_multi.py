#!/usr/bin/env python3 

import time, sys, os
sys.path.append("/home/buding/Jarvis-HEP/src/Sampling/dynesty/py")
import numpy as np
from numpy import linalg

import matplotlib
from matplotlib import pyplot as plt

rstate = np.random.default_rng()

# import dynesty
from dynesty import *
from dynesty import plotting as dyplot
from dynesty.results import print_fn
import matplotlib.pyplot as plt 
from asyncio import sleep
import random
import pandas as pd 
from subprocess import Popen

tmax = 5.0 * np.pi
def loglike(x,  **kwarg):
    # print(x[:-1])
    # tsleep = 0.05 * random.random()
    # tstart = time.time()
    # print("LL called for {} at {}".format(x, kwarg))
    # sleep(1)
    # p = Popen("sleep {:.2f}".format(tsleep), shell=True)
    # # time.sleep(random.random())
    # while p.poll() is None:
    #     time.sleep(0.0001)
    t = 2.0 * tmax * x - tmax
    # print("\tLL ended for {} at {:.2f}".format(x, time.time()-tstart))
    if t[0] > t[1]:
        return (2.0 + np.cos(t[0] / 2.0) * np.cos(t[1] / 2.0)) ** 5.0
    else:
        return -1.0
# define the prior transform
def prior_transform(x):
    return x

# import dynesty.pool as dypool
# with dypool.Pool(16, loglike, prior_transform) as pool:
# #     # The important thing that we provide the loglikelihood/prior transform from 
# #     # the pool    
#     psampler = dynesty.DynamicNestedSampler(
#         pool.loglike, pool.prior_transform, 
#         ndim=2, bound='multi', sample='unif', rstate=rstate,
#         pool=pool)
#     psampler.run_nested(dlogz_init=0.05, nlive_init=500, nlive_batch=500,
#                      wt_kwargs={'pfrac': 0.0}, stop_kwargs={'pfrac': 0.0}, 
#                      checkpoint_file="./test.sav", checkpoint_every=1)
# #     psampler = dynesty.DynamicNestedSampler(
# #         pool.loglike, pool.prior_transform, 
# #         ndim=2,
# #         pool=pool
# #     )
# #     # psampler = dynesty.DynamicNestedSampler.restore("./test.sav", pool=pool)
# #     psampler.restore("./test.sav", pool=pool)
# #     # psampler = psampler.restore("./test.sav", pool=pool)
# #     psampler.run_nested(resume=True, checkpoint_file="./test.sav", checkpoint_every=1)
# pres = psampler.results
# print(pres)
# dres = dict(pres)
# dres = dres.pop("bound")
# # print(dict(pres).keys(), dres.keys())
# for kk, vv in dict(pres).items():
#     print("{}, \n\t{}".format(kk, vv))
# import json 
# with open("./dynesty.json", 'w') as f1:
#     json.dump(dres, f1, indent=4)
    

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
    checkpoint_file="./test.sav", checkpoint_every=1
    )

# dsampler = dynesty.DynamicNestedSampler.restore("./test.sav")
# dsampler.run_nested(
#     dlogz_init=0.01, 
#     nlive_init=500, 
#     nlive_batch=500,
#     wt_kwargs={'pfrac': 0.0}, 
#     stop_kwargs={'pfrac': 0.0},
#     checkpoint_file="./test.sav", checkpoint_every=1
#     )
pres = dsampler.results
# print(dict(pres))
# knl = []
# from copy import deepcopy
# kil = deepcopy(list(dict(pres).keys()))
# dd = {}


# for kk, vv in dict(pres).items():
#     if type(vv) not in  [int, float]:
#         print(kk, type(vv[0]), vv[0])
#         dd[kk] = vv
#     else:
#         knl.append(kk)
#         kil.remove(kk)
#         # print("Not iterable key {}".format(kk))
# print("not iterable keys {}".format(knl))
# print("iterable keys {}".format(kil))
# df = pd.DataFrame(dd)
# print(df)

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
#                               span=[[0, 1], [0, 1]],
#                               fig=plt.subplots(2, 2, figsize=(10, 10)))
# plt.show()

fig, axes = dyplot.cornerplot(pres, quantiles=None, color='darkviolet',
                              span=([0, 1], [0, 1]),
                              fig=plt.subplots(2, 2, figsize=(10, 10)))
plt.show()

lnz_truth = -2 * np.log(10. * 0.999999426697)
print(lnz_truth)
fig, axes = dyplot.runplot(pres, color='navy', logplot=True,
                           lnz_truth=lnz_truth, truth_color='black')
plt.show()

fig, axes = dyplot.traceplot(pres, truths=np.zeros(2), show_titles=True, trace_cmap="plasma", quantiles=None)
plt.show()


fig = plt.figure(figsize=(6, 6))
ax = fig.add_axes([0.1, 0.1, 0.88, 0.88])
ax.plot([0., 1.], [0., 1.], '-')
fg, ax = dyplot.cornerpoints(pres, cmap='plasma', fig=(fig, ax), plot_kwargs={"marker": 'o', "s": 6})
ax[0, 0].set_title('No Bound', fontsize=26)
ax[0, 0].set_xlim([0, 1])
ax[0, 0].set_ylim([0, 1])
plt.show()