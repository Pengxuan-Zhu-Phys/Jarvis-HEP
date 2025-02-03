#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import numpy as np
import time 
from Sampling.sampler import SamplingVirtial, BoolConversionError
import json
from sample import Sample
import concurrent.futures
import itertools
import sympy as sp 
from test.test_imaplib import RemoteIMAPTest

# from _pytest.mark import param


class RandomS(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Random"
        self._P     = None
        self._index = None
        self.tasks  = []
        self.info   = {}
        self._selectionexp = None
        self.future_to_sample = {}

    def load_schema_file(self):
        self.schema = self.path['RandomSchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def __iter__(self):
        if self._index is None:
            self.initialize()  # ensure the _P is generated before iteration 
        self._index = 0  # Ensure the index starting from 0
        return self

    def __next__(self):# Stop iteration, if _P is not defined or _P is not np.array
        if self._index < self._maxp:
            is_selection = True 
            while is_selection: 
                temp    = np.random.rand(self._dimensions)
                param   = self.map_point_into_distribution(temp)
                if self._selectionexp: 
                    is_selection = self.evaluate_selection(self._selectionexp, param)
            self._index += 1
            return param
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
        self._maxp  = int(self.config['Sampling']['Point number'])
        if "selection" in self.config["Sampling"]:
            self._selectionexp = self.config["Sampling"]["selection"]

    def initialize(self):
        self.logger.warning("Initializing the Random Sampling")
        self._index = 0 
        if self._selectionexp:
            try:
                self.info["t0"]       = time.time() 
                temp    = np.random.rand(self._dimensions)
                param   = self.map_point_into_distribution(temp)
                self.evaluate_selection(self._selectionexp, param)
            except BoolConversionError:
                self.logger.error("Wrong selection condition in input YAML -> \n\t{}".format(self._selectionexp))
                sys.exit(2)
            except:
                self.logger.error("Random Sampler meets error when trying scan the parameter space.")
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
                        self.future_to_sample[future] = sample
                    else: 
                        break
                except StopIteration:
                    break
            
            done, _ = concurrent.futures.wait(self.tasks, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED)
            # Remove completed futures
            # done, _ = concurrent.futures.wait(self.tasks, timeout=0.1, return_when=concurrent.futures.FIRST_COMPLETED)
    
            # Process completed futures
            for future in done:
                try:
                    result = future.result()  # Retrieve result of completed sample

                    # Retrieve the corresponding sample instance
                    sample = self.future_to_sample.pop(future, None)

                    if sample:
                        self.logger.info(f"Sample completed with result: {result}, Sample Info: {sample.info['observables']['z']}")

                        # If needed, store results in a list or database
                        # self.completed_samples.append((sample.info, result))
                        # time.sleep(10)

                except Exception as e:
                    self.logger.error(f"Error processing sample: {e}")

            # Remove completed futures from task list
            # self.tasks = [f for f in self.tasks if f not in done]
            
            
            
            self.tasks = [f for f in self.tasks if f not in done]
            # Exit loop if no tasks are pending and no more samples
            if not self.tasks and not param:
                break  

    def finalize(self):
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Random sampler")

    def evaluate_selection(self, expression, variables):
        return super().evaluate_selection(expression, variables)