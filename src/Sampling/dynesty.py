#!/usr/bin/env python3 

from lib2to3.pgen2.token import RPAR
import os, sys
import numpy as np
from Sampling.sampler import SamplingVirtial
from sample import Sample
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor
from uuid import uuid4

def prior_transform(u):
    uuid = str(uuid4())
    ret = np.append(u, [uuid])
    return ret
class Dynesty(SamplingVirtial):
    def __init__(self) -> None:
        super().__init__()
        self.load_schema_file()
        self.method = "Dynesty"
        self.sampler = None
        # self.max_workers = os.cpu_count()
        self.max_workers = 2

    def load_schema_file(self):
        self.schema = self.path['DynestySchema']

    def set_config(self, config_info) -> None:
        self.config = config_info
        self.init_generator()

    def set_factory(self, factory) -> None:
        self.factory = factory
        self.logger.warning("WorkerFactory is ready for Bridson sampler")

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

    def map_point_into_distribution(self, row) -> np.ndarray:
        result = {}
        for ii in range(len(row)):
            result[self.vars[ii].name] = self.vars[ii].map_standard_random_to_distribution(row[ii]/self.vars[ii]._parameters['length'])
        return result

    def set_logger(self, logger) -> None:
        super().set_logger(logger)
        self.logger.warning("Sampling method initializaing ...")

    def init_generator(self):
        self.load_variable()
        self._nlive = self.config['Sampling']['Bounds']['nlive']
        self._rstate = np.random.default_rng(self.config['Sampling']['Bounds']['rseed'])
        self._dlogz = self.config['Sampling']['Bounds']['dlogz']
        self._dimensions = len(self.vars)

    def initialize(self):
        self.logger.warning("Initializing the Dynesty Sampling")

    def run_nested(self):
        def log_likelihood(params):
            param = params[0:-1].astype(np.float16)
            uuid = params[-1]

            pars = self.map_point_into_distribution(param)
            sample = Sample(pars)
            sample.update_uuid(uuid)
            sample.set_config(deepcopy(self.info['sample']))
            future = self.factory.submit_task(sample.params, sample.info)
            result = future.result()
            return result 
        
        # from dynesty import DynamicNestedSampler
        from Source.Dynesty.py.dynesty import DynamicNestedSampler
        with ThreadPoolExecutor(max_workers=4) as executor:
            self.sampler = DynamicNestedSampler(
                loglikelihood=log_likelihood, 
                prior_transform=prior_transform,
                ndim=self._dimensions,
                nlive=self._nlive,
                pool=executor, 
                rstate=self._rstate,
                queue_size=2*self.max_workers
            )
            self.sampler.run_nested(
                dlogz_init = self._dlogz
            )

    def finalize(self):
        print(self.sampler.results)
        pass 