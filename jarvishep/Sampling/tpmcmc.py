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
# from jarvishep.Sampling.sampler import SamplingVirtial
from jarvishep.Sampling.sampler import SamplingVirtial
import json
from scipy.special import gammainc
from jarvishep.sample import Sample
import concurrent.futures
from copy import deepcopy

class TPMCMC(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "TPMCMC"
        self._P     = None
        self._index = 0 
        self.tasks  = set()
        self.info   = {}
        self.status = "init"
        self.future_to_sample = {}
        self.chain_next = []
        

    def load_schema_file(self):
        self.schema = self.path['TPMCMCSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        self.status = "run"
        self._index = 0  # Ensure the index starting from 0
        return self

    def _refresh_chain_queue_if_needed(self):
        if not hasattr(self, "_exchange_iterval"):
            return
        if not hasattr(self, "_chain_length"):
            return
        if not hasattr(self, "_exchange_cycle"):
            return
        if not hasattr(self, "_chain_iters"):
            return
        if not hasattr(self, "_nchains"):
            return
        if not all(x == self._exchange_iterval for x in self._chain_length):
            return
        idxs = [next(self._exchange_cycle) for ii in range(self._nchains)]
        for pair in list(zip(idxs[::2], idxs[1::2])):
            self.exchange_chain(pair)
        self.chain_next = [ii for ii in range(self._nchains)]
        self._chain_length = [0 for ii in range(self._nchains)]
        self.logger.warning(f"Each chain iters {self._chain_iters[0]} samples")
        time.sleep(5)

    def __next__(self):
        self._refresh_chain_queue_if_needed()
        # print("In next(): self.chain_next -> ", self.chain_next)
        while self.chain_next:
            cid = self.chain_next.pop(0)
            if self._chain_length[cid] >= self._exchange_iterval:
                continue
            # print("In next: cid -> ", cid)
            is_seletion = False 
            while not is_seletion: 
                proposal = next(self._chains[cid])
                # print("In next(): next proposal -> ", proposal)
                param = self.map_point_into_distribution(proposal)
                if self._selectionexp: 
                    is_seletion = self.evaluate_selection(self._selectionexp, param)
                else: 
                    is_seletion = True
                    # print("In next(): param -> ", is_seletion, param)
            # print("In next: return the samples")
            return param, cid 
        raise StopIteration
            

    def next_sample(self):
        return self.__next__()

    def exchange_chain(self, pair):
        from math import exp 
        likelihood1 = exp(self._chains[pair[0]].last_loglikelihood)
        likelihood2 = exp(self._chains[pair[1]].last_loglikelihood)
        acceptance_prob = likelihood2 / likelihood1
        # Swap the last states and their likelihoods
        if np.random.rand() < acceptance_prob:
            self._chains[pair[0]].param, self._chains[pair[1]].param = self._chains[pair[1]].param, self._chains[pair[0]].param
            self._chains[pair[0]].last_loglikelihood, self._chains[pair[1]].last_loglikelihood = self._chains[pair[1]].last_loglikelihood, self._chains[pair[0]].last_loglikelihood

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
        smp = self.config["Sampling"]["Bounds"]
        self._nchains   = smp["num_chains"]
        self._niters    = smp["num_iters"]
        self._exchange_iterval = smp["exchange_interval"]
        self._proposal_scales  = smp["proposal_scales"]
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]
        self._total_cores = self._nchains

    def initialize(self):
        self.logger.warning("Initializing the temporal parallel MCMC (TPMCMC) Sampling")
        try:
            t0 = time.time()
            from jarvishep.Sampling.Source.MCMC.mcmc_chain import MCMCChain
            from jarvishep.Sampling.Source.MCMC.tpmcmc_data import TPMCMCData
            self._chains = [
                MCMCChain(np.random.rand(self._dimensions), self._proposal_scales[ii], self._niters)
                for ii in range(self._nchains)
            ]
            # Initialize a counter to track iterations between exchanges.
            self._exchange_counter = 0 
            self.data = TPMCMCData()
            self._global_iter       = 0 
            self.info["t0"]         = time.time() - t0 
            self.logger.info("TPMCMC Sampler initialized in {:.2f} sec".format(self.info['t0']))
            self._chain_length = [0 for ii in range(self._nchains)]
            self._chain_iters = [0 for ii in range(self._nchains)]
            from itertools import cycle
            self._exchange_cycle = cycle(range(self._nchains))
            self.chain_next = [ii for ii in range(self._nchains)]
        except:
            self.logger.error("TPMCMC Sampler meets error when trying scan the parameter space.")
            sys.exit(2)

    def run_nested(self):
        self.tasks = set()
        self.future_to_sample = {}
        exhausted = False
        while (not exhausted) or self.tasks:
            self._refresh_chain_queue_if_needed()
            while not exhausted and len(self.tasks) < self._total_cores:
                try: 
                    param, cid = self.next_sample()
                except StopIteration:
                    if not self.tasks:
                        exhausted = True
                    break

                sample = Sample(param)
                sample.set_config(deepcopy(self.info['sample']))
                sample.info["chain_id"] = cid
                future = self.factory.submit_task(sample.info)
                self.tasks.add(future)
                self.future_to_sample[future] = sample

            if all(lgh >= self._niters for lgh in self._chain_iters):
                exhausted = True

            if not self.tasks:
                continue

            done, _ = concurrent.futures.wait(
                self.tasks,
                return_when=concurrent.futures.FIRST_COMPLETED,
            )
            self.tasks.difference_update(done)
            for future in done: 
                sample = self.future_to_sample.pop(future, None)
                try: 
                    future.result()
                    if sample:
                        cid = sample.info["chain_id"]
                        self._chains[cid].update(sample.info['observables']['LogL'])
                        self._chain_length[cid] += 1
                        self._chain_iters[cid] += 1
                        if self._chain_iters[cid] < self._niters:
                            self.chain_next.append(cid)
                except Exception as exc:
                    suuid = sample.uuid if sample else "UNKNOWN"
                    self.logger.error(f"[WorkerFactory] future exception consumed: uuid={suuid} error={exc}")
                    raise
                finally:
                    if sample is not None:
                        sample.close()
            if all(lgh >= self._niters for lgh in self._chain_iters):
                exhausted = True

    def finalize(self):
        pass


    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for TPMCMC sampler")
    
