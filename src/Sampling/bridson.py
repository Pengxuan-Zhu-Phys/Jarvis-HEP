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
from random import randint
from Sampling.sampler import SamplingVirtial
import json
from scipy.special import gammainc
from sample import Sample
import concurrent.futures
from utils import format_duration
class Bridson(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Bridson"
        self._P     = None
        self._index = 0 
        self.tasks  = {}
        self.info   = {}
        self.barinfo = {}

    def load_schema_file(self):
        self.schema = self.path['BridsonSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.set_bucket_alloc()
        self.init_generator()

    def __iter__(self):
        if self._P is None:
            self.initialize()  # ensure the _P is generated before iteration 
        self._index = 0  # Ensure the index starting from 0
        return self

    def __next__(self):
        if self._P is None or not isinstance(self._P, np.ndarray):
            raise StopIteration  # Stop iteration, if _P is not defined or _P is not np.array
        if self._index < len(self._P):
            result = self._P[self._index]
            result = self.map_point_into_distribution(result)
            self._index += 1
            self.progress_bar()
            return result
        else:
            raise StopIteration

    def next_sample(self):
        return self.__next__()

    def progress_bar(self): 
        if self.barinfo == {}:
            self.barinfo = {
                "total":    len(self._P), 
                "t0":       time.time(), 
                "permille":  0
            }
            self.logger.warning(
                "{}‰ of {}/{} samples submited in {}".format(
                    self.barinfo['permille'], 
                    self._index, 
                    self.barinfo['total'], 
                    format_duration(time.time() - self.barinfo['t0'])
                    )
                )
        else: 
            permille = int(self._index / self.barinfo["total"] * 1000)
            if permille != self.barinfo["permille"]: 
                self.barinfo.update({"permille": permille})
                self.logger.warning(
                    "{}‰ of {}/{} samples submited in {}".format(
                        self.barinfo['permille'], 
                        self._index, 
                        self.barinfo['total'], 
                        format_duration(time.time() - self.barinfo['t0']))
                    )
        

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii]/self.vars[ii]._parameters['length'])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        super().init_generator()
        self.load_variable()
        self._radius = self.config['Sampling']['Radius']
        self._k = self.config['Sampling']['MaxAttempt']
        self._dimensions = len(self.vars)
        self.total_core = self.config['Sampling'].get("MaxWorker", self.total_core)
        self.load_nuisance_sampler()            
            
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
            self._index           = 0
            self.logger.info("Bridson Sampler obtains {} samples in {:.2f} sec".format(self.info['NSamples'], self.info['t0']))
        except:
            self.logger.error("Bridson Sampler meets error when trying to scan the parameter space.")
            sys.exit(2)

    def run_nested(self):
        if not self._with_nuisance: 
            self.run_wo_nuisance()
        else:
            self.run_w_nuisance()
     
    def run_w_nuisance(self):
        total_cores = self.total_core
        from copy import deepcopy
        exhausted = False 
        while True: 
            while not exhausted and len(self.tasks) < total_cores: 
                try: 
                    param = self.next_sample()
                    nuis_param = self.nuisance_sampler
                except StopIteration: 
                    exhausted = True
                    break 
                
                if param is None: 
                    exhausted = True 
                    break 
                
                sample = Sample(param)
                sconfig = deepcopy(self.info['sample'])
                sconfig['save_dir'] = self.bucket_alloc.next_bucket_dir()
                sample.set_config(sconfig)

                sample.start()
                future = self.factory.submit_task(sample.info)
                self.tasks[future] = sample
        
            done, _ = concurrent.futures.wait(
                list(self.tasks.keys()), 
                timeout=0.01, 
                return_when=concurrent.futures.FIRST_COMPLETED
            )
            
            for future in done: 
                sample = self.tasks.pop(future) 

                try:
                    future.result()
                except Exception as e: 
                    self.logger.error(f"Sample error info -> {e}")
                
    
                self.nuisance_sampler.renew_sample_info(sample.info)
                if sample.info['status'] == "Accept": 
                    sample.close()
                else: 
                    sample.record()
                    sample.combine_nuisance_card()
                    future = self.factory.submit_task(sample.info)
                    self.tasks[future] = sample 

            if exhausted and not self.tasks: 
                break
                
     
        
    def run_wo_nuisance(self):
        total_cores = self.total_core
        from copy import deepcopy
        exhausted = False 
        while True:
            while not exhausted and len(self.tasks) < total_cores: 
                try: 
                    param = self.next_sample()
                except StopIteration: 
                    exhausted = True 
                    break 
                
                if param is None: 
                    exhausted = True 
                    break 
                sample = Sample(param)
                sconfig = deepcopy(self.info['sample'])
                sconfig['save_dir'] = self.bucket_alloc.next_bucket_dir()
                sample.set_config(sconfig)
                
                future = self.factory.submit_task(sample.info)
                self.tasks[future] = sample

            # futures = [f for f, _ in self.tasks]
            done, _ = concurrent.futures.wait(
                list(self.tasks.keys()), 
                timeout=0.01, 
                return_when=concurrent.futures.FIRST_COMPLETED
            )
            
            for future in done: 
                sample = self.tasks.pop(future) 
                try: 
                    future.result() 
                except Exception as e: 
                    self.logger.error(f"Sample failed: {e}")
                finally: 
                    sample.close()
          
            if exhausted and not self.tasks: 
                break
            

    def finalize(self):
        pass

     


    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Bridson sampler")

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



