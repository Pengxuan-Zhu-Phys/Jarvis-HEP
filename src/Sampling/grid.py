#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import pandas as pd 
from numpy import meshgrid
from pandas.core.series import Series
from sympy import sympify
from sympy.geometry import parabola
import time 
from random import randint
# from Sampling.sampler import SamplingVirtial
from sampler import SamplingVirtial
import json
from scipy.special import gammainc
from sample import Sample
import concurrent.futures
import itertools


class Grid(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Grid"
        self._P     = None
        self._index = 0 
        self.tasks  = []
        self.info   = {}

    def load_schema_file(self):
        self.schema = self.path['GridSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
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
            return result
        else:
            # raise StopIteration
            return None

    def next_sample(self):
        return self.__next__()

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._dimensions = len(self.vars)

    def initialize(self):
        self.logger.warning("Initializing the Grid Sampling")
        try:
            t0 = time.time()
            dims = np.array([var.parameters["num"] for var in self.vars])
            self._P = grid_sampling(dimensions=dims)
            self.info["NSamples"] = self._P.shape[0]
            self.info["t0"]       = time.time() - t0 
            self._index           = 0
            self.logger.info("Grid Sampler obtains {} samples in {:.2f} sec".format(self.info['NSamples'], self.info['t0']))
        except:
            self.logger.error("Grid Sampler meets error when trying scan the parameter space.")
            sys.exit(2)

    def run_nested(self):
        total_cores = os.cpu_count()
        from copy import deepcopy

        while True:
            while len(self.tasks) < total_cores: 
                try: 
                    param = self.next_sample()
                    if param is not None: 
                        sample = Sample(param)
                        sample.set_config(deepcopy(self.info['sample']))
                        future = self.factory.submit_task(sample.params, sample.info)
                        self.tasks.append(future)
                    else: 
                        break
                except StopIteration:
                    break
            
            done, _ = concurrent.futures.wait(self.tasks, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED)
            # Remove completed futures
            self.tasks = [f for f in self.tasks if f not in done]
            # Exit loop if no tasks are pending and no more samples
            if not self.tasks and not param:
                break  

    def finalize(self):
        pass


    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Grid sampler")

def grid_sampling(dimensions):
    """
    Generate grid samples based on dimensions and dimension-specific number of steps.

    Args:
        dimensions (list of tuple): A list of tuples where each tuple specifies the number of step for each dimension.

    Returns:
        numpy.ndarray: A 2D array where each row represents a grid sample.
    """

    # Generate grid points for each dimension based on the number of steps
    grid_ranges = [
        np.linspace(0., 1., steps) 
        for steps in dimensions
    ]
    # Generate the Cartesian product of all grid points
    grid_samples = np.array(list(itertools.product(*grid_ranges)))
    return grid_samples
