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
        # for dynesty sampling method, next() method is only for testing the assembing line, not the real sampling method
        u = np.random.random(self._dimensions)
        param = self.map_point_into_distribution(u)
        return param
        

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
        print(self.vars[0]._parameters)
        self._nlive = self.config['Sampling']['Bounds']['nlive']
        self._rstate = np.random.default_rng(self.config['Sampling']['Bounds']['rseed'])
        self._dlogz = self.config['Sampling']['Bounds']['dlogz']
        self._dimensions = len(self.vars)

    def initialize(self):
        self.logger.warning("Initializing the Dynesty Sampling")

    def init_sampler_db(self):
        self.info['db'] = {
            "nested result":    os.path.join(self.info['sample']['task_result_dir'], "dynesty_result.csv")
        }

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
        
        self.init_sampler_db()

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

    def save_dynesty_results_to_csv(self):
        import pandas as pd 
        data = {
            "uuid":             self.sampler.results['samples_uid'],
            "log_weight":       self.sampler.results['logwt'],
            "log_Like":         self.sampler.results['logl'],
            "log_PriorVolume":  self.sampler.results['logvol'],
            "log_Evidence":     self.sampler.results['logz'],
            "log_Evidence_err": self.sampler.results['logzerr'],
            "samples_nlive":    self.sampler.results['samples_n'],
            "ncall":            self.sampler.results['ncall'],
            "samples_it":       self.sampler.results['samples_it'],
            "samples_id":       self.sampler.results['samples_id'],
            "information":      self.sampler.results['information']
        }

        for ii in range(self.sampler.results['samples'].shape[1]):
            data[f'samples_v[{ii}]'] = self.sampler.results['samples'][:, ii]
        for ii in range(self.sampler.results['samples_u'].shape[1]):
            data[f'samples_u[{ii}]'] = self.sampler.results['samples'][:, ii]

        df = pd.DataFrame(data)
        df.to_csv(self.info['db']['nested result'], index=False)
        self.logger.info(f"Results saved to {self.info['db']['nested result']}")
        self.logger.info(self.sampler.results.summary())
