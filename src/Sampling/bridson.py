#!/usr/bin/env python3 
from calendar import day_name
from lib2to3.pgen2.token import RPAR
import logging
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import pandas as pd 
from numpy.lib.function_base import meshgrid
from pandas.core.series import Series
from sympy import sympify
from sympy.geometry import parabola
import time 
from Func_lib import decode_path_from_file
from random import randint
from Sampling.sampler import SamplingVirtial
import json
from scipy.special import gammainc

class Bridson(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Bridson"

    def load_schema_file(self):
        self.schema = self.path['BridsonSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        self.active_list = [np.random.uniform(np.zeros(self.dims.size), self.dims)]
        self.samples = []
        self.grid = self.initialize_grid()
        return self

    def __next__(self):
        if not self.active_list:
            raise StopIteration
        # Sampling logic here, modifying to yield a sample at a time
        point = self.active_list.pop()
        new_points = self.generate_points_around(point)
        for new_point in new_points:
            if self.is_valid_point(new_point):
                self.active_list.append(new_point)
                self.samples.append(new_point)
                return new_point
        raise StopIteration

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._radius = self.config['Sampling']['Radius']
        self._k = self.config['Sampling']['MaxAttempt']
        self._dimensions = len(self.vars)

    def initialize(self):
        self.logger.warning("Initializing the Bridson Sampling")
        if self._dimensions < 2 or self._dimensions >= 5:
            self.logger.error("Bridson Sampling Algorithm only support 2d to 4d parameter space.")
            sys.exit(2)
        try:
            t0 = time.time()
            dims = np.array([var.parameters["length"] for var in self.vars])
            self._P = Bridson_sampling(
                dims=dims, radius=self._radius, k=self._k, hypersphere_sample=hypersphere_surface_sample
            )
            self.info["NSamples"] = self._P.shape[0]
            self.info["t0"]       = time.time() - t0 
            self.logger.info("Bridson Sampler obtains {} samples in {:.2f} sec".format(self.info['NSamples'], self.info['t0']))
        except:
            self.logger.error("Bridson Sampler meets error when trying to scan the parameter space.")
            sys.exit(2)






# Can be optimized for Bridson algorithm by excluding all points within the r/2 sphere
def hypersphere_volume_sample(center,radius,k=1):
    ndim = center.size
    x = np.random.normal(size=(k, ndim))
    ssq = np.sum(x**2,axis=1)
    fr = radius*gammainc(ndim/2,ssq/2)**(1/ndim)/np.sqrt(ssq)
    frtiled = np.tile(fr.reshape(k,1),(1,ndim))
    p = center + np.multiply(x,frtiled)
    return p


# Uniform sampling on the sphere's surface
def hypersphere_surface_sample(center,radius,k=1):
    ndim = center.size
    vec = np.random.standard_normal(size=(k, ndim))
    vec /= np.linalg.norm(vec, axis=1)[:,None]
    p = center + np.multiply(vec, radius)
    return p


def squared_distance(p0, p1):
    return np.sum(np.square(p0-p1))

def Bridson_sampling(dims=np.array([1.0,1.0]), radius=0.05, k=30, hypersphere_sample=hypersphere_volume_sample):
    # References: Fast Poisson Disk Sampling in Arbitrary Dimensions
    #             Robert Bridson, SIGGRAPH, 2007

    ndim=dims.size

    # size of the sphere from which the samples are drawn relative to the size of a disc (radius)
    sample_factor = 2
    if hypersphere_sample == hypersphere_volume_sample:
        sample_factor = 2
        
    # for the surface sampler, all new points are almost exactly 1 radius away from at least one existing sample
    # eps to avoid rejection
    if hypersphere_sample == hypersphere_surface_sample:
        eps = 0.001
        sample_factor = 1 + eps
    
    def in_limits(p):
        return np.all(np.zeros(ndim) <= p) and np.all(p < dims)

    # Check if there are samples closer than "squared_radius" to the candidate "p"
    def in_neighborhood(p, n=2):
        indices = (p / cellsize).astype(int)
        indmin = np.maximum(indices - n, np.zeros(ndim, dtype=int))
        indmax = np.minimum(indices + n + 1, gridsize)
        
        # Check if the center cell is empty
        if not np.isnan(P[tuple(indices)][0]):
            return True
        a = []
        for i in range(ndim):
            a.append(slice(indmin[i], indmax[i]))
        with np.errstate(invalid='ignore'):
            if np.any(np.sum(np.square(p - P[tuple(a)]), axis=ndim) < squared_radius):
                return True

    def add_point(p):
        points.append(p)
        indices = (p/cellsize).astype(int)
        P[tuple(indices)] = p

    cellsize = radius/np.sqrt(ndim)
    gridsize = (np.ceil(dims/cellsize)).astype(int)

    # Squared radius because we'll compare squared distance
    squared_radius = radius*radius

    # Positions of cells
    P = np.empty(np.append(gridsize, ndim), dtype=np.float32) #n-dim value for each grid cell
    # Initialise empty cells with NaNs
    P.fill(np.nan)

    points = []
    add_point(np.random.uniform(np.zeros(ndim), dims))
    while len(points):
        i = np.random.randint(len(points))
        p = points[i]
        del points[i]
        Q = hypersphere_sample(np.array(p), radius * sample_factor, k)
        for q in Q:
            if in_limits(q) and not in_neighborhood(q):
                add_point(q)
    return P[~np.isnan(P).any(axis=ndim)]



