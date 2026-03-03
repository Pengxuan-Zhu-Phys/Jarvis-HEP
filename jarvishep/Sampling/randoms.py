#!/usr/bin/env python3 
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
# from mpmath.functions.functions import re
import numpy as np
import time 
from jarvishep.Sampling.sampler import SamplingVirtial, BoolConversionError
import json
from jarvishep.sample import Sample
import concurrent.futures
import sympy as sp 

# from _pytest.mark import param


class RandomS(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Random"
        self._P     = None
        self._index = None
        self.tasks  = set()
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
            if self._selectionexp:
                is_selection = False
                while not is_selection:
                    temp    = np.random.rand(self._dimensions)
                    param   = self.map_point_into_distribution(temp)
                    if self._selectionexp: 
                        is_selection = self.evaluate_selection(self._selectionexp, param)
                self._index += 1
                return param
            else: 
                temp = np.random.random(self._dimensions)
                param = self.map_point_into_distribution(temp)
                self._index += 1 
                return param 
        else:
            raise StopIteration

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
        # else: 
        #     try: 
        #         self.info['to']     = time.time() 
                

    def run_nested(self):
        total_cores = os.cpu_count() or 1
        self.tasks = set()
        self.future_to_sample = {}
        exhausted = False
        base_sample_cfg = self.info['sample']

        while (not exhausted) or self.tasks:
            while not exhausted and len(self.tasks) < total_cores:
                try: 
                    param = self.next_sample()
                except StopIteration:
                    exhausted = True
                    break

                sample = Sample(param)
                sample.set_config(self.build_sample_config(base_sample_cfg))
                future = self.factory.submit_task(sample.info)
                self.tasks.add(future)
                self.future_to_sample[future] = sample

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
                except Exception as exc:
                    suuid = sample.uuid if sample else "UNKNOWN"
                    self.logger.error(f"[WorkerFactory] future exception consumed: uuid={suuid} error={exc}")
                    raise
                finally:
                    if sample is not None:
                        sample.close()

    def finalize(self):
        pass

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Random sampler")
