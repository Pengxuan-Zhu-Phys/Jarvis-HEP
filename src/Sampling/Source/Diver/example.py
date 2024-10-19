#!/usr/bin/env python3 

import numpy as np
import math
from scipy.optimize import differential_evolution


# Define the parameter bounds and settings
param_dim = 2
NP = 10
numgen = 15
numciv = 1
nDerived = 0
path = 'example_f/output/example'
Cr = 0.9
tol = 1e-3
lamb = 0.8
F = 0.6
lowerbounds = np.array([-50.0] * param_dim)
upperbounds = np.array([50.0] * param_dim)
boundranges = upperbounds - lowerbounds


# Define the target functions
def constant(params):
    return 0.0 * params[0]


def step(params):
    if params[0] > 0.0:
        return 0.0
    else:
        return 1.0


def linear(params):
    if params[0] > 0.0:
        return params[0]
    else:
        return 0.0


def gauss(params):
    return np.sum(params**2)


def manygauss(params):
    discrete_indices = [0, 2, 4]  # Python uses 0-based indexing
    val = 0.0
    for i in range(len(params)):
        if i not in discrete_indices:
            val += params[i]**2
    return val + 1.0


def spikygauss(params):
    val = np.sum(params**2 + 1e32 * np.sin(params)**2)
    return val


def rosenbrock(params):
    return (1 - params[0])**2 + 100 * (params[1] - params[0]**2)**2 + 1


def rastrigin(params):
    A = 10
    return A * len(params) + np.sum(params**2 - A * np.cos(2 * math.pi * params))


def eggcarton(params):
    return -(2 + np.cos(0.5 * params[0]) * np.cos(0.5 * params[1]))**5 + 1


# Define the prior function (flat prior)
def flatprior(X):
    return 1.0 / np.prod(boundranges)


# Helper function to convert bounds to the format for scipy's differential evolution
def get_bounds(lowerbounds, upperbounds):
    return [(low, high) for low, high in zip(lowerbounds, upperbounds)]


# Define the optimization procedure using DEOptimizer
def example_f():
    optimizer = DEOptimizer()
    result = optimizer.diver(
        func=rosenbrock,  # Use the rosenbrock function or any other function like gauss, manygauss, etc.
        lowerbounds=lowerbounds,
        upperbounds=upperbounds,
        NP=NP,
        maxgen=maxgen,
        maxciv=maxciv,
        F=F,
        Cr=Cr,
        lamb=lambda_,
        doBayesian=False,
        path=path,
        verbose=3,
        prior=flatprior,
        lambdajDE=True,
        discard_unfit_points=True
    )
    
    print("\nMinimum found:", result)

# Run the example function
if __name__ == "__main__":
    example_f()
