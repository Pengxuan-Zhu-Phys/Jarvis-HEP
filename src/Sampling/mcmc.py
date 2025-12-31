#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
# from mpmath.functions.functions import re
import numpy as np
import time 
from Sampling.sampler import SamplingVirtial, BoolConversionError
import json
from sample import Sample
import concurrent.futures
import sympy as sp 


class MCMC(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "MCMC"
        self._P     = None
        self._index = None
        self.tasks  = []
        self.info   = {}
        self._selectionexp = None
        self.future_to_sample = {}
        self.chain_next = []

    def load_schema_file(self):
        self.schema = self.path['MCMCSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        if self._index is None:
            self.initialize()  # ensure the _P is generated before iteration 
        self._index = 0  # Ensure the index starting from 0
        return self

    def __next__(self):# Stop iteration, if _P is not defined or _P is not np.array
        while self.chain_next: 
            cid = self.chain_next.pop(0)
            is_selection = False 
            while not is_selection: 
                proposal = next(self._chains[cid])
                param = self.map_point_into_distribution(proposal)
                if self._selectionexp: 
                    is_selection = self.evaluate_selection(self._selectionexp, param)
                else: 
                    is_selection = True
            return param, cid 
        else:
            # raise StopIteration
            return None, None

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
        smp = self.config["Sampling"]['Bounds']
        self._nchains = smp['num_chains']
        self._niters  = smp["num_iters"]
        self._proposal_scales = smp["proposal_scale"]
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the MCMC Sampling")
        self._index = 0 
        try: 
            t0 = time.time() 
            from Source.MCMC.mcmc_chain import MCMCChain
            self._chains = [
                MCMCChain(np.random.rand(self._dimensions), self._proposal_scales, self._niters)
                for ii in range(self._nchains)
            ]
            self.info["t0"] = time.time() - t0 
            self.logger.info("MCMC Sampler initialized in {:.2f} sec".format(self.info['t0']))
            self._chain_iters = [0 for ii in range(self._nchains)]
            self.chain_next = [ii for ii in range(self._nchains)]
        except: 
            self.logger.error("MCMC Sampler meets error when trying scan the parameter space.")
            sys.exit(2)

    def run_nested(self):
        from copy import deepcopy
        while True:
            while len(self.tasks) < self._nchains: 
                try: 
                    param, cid = self.next_sample()
                    if param is not None: 
                        sample = Sample(param)
                        sample.set_config(deepcopy(self.info['sample']))
                        sample.info['chain_id'] = cid 
                        future = self.factory.submit_task(sample.params, sample.info)
                        self.tasks.append(future)
                        self.future_to_sample[future] = sample
                    else: 
                        break
                except StopIteration:
                    break
            
            done, _ = concurrent.futures.wait(self.tasks, timeout=0.01, return_when=concurrent.futures.FIRST_COMPLETED)
            
            self.tasks = [f for f in self.tasks if f not in done]
            for future in done:
                try: 
                    result = future.result() 
                    sample = self.future_to_sample.pop(future, None)
                    if sample: 
                        cid = sample.info['chain_id']
                        self.chain_next.append(cid)
                        self._chains[cid].update(sample.info['observables']['LogL'])
                        self._chain_iters[cid] += 1
                except Exception as e: 
                    self.logger.error(f"Error processing sample: {e}")
            if all(lgh >= self._niters for lgh in self._chain_iters):
                break 

    def finalize(self):
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for MCMC sampler")

    def evaluate_selection(self, expression, variables):
        return super().evaluate_selection(expression, variables)