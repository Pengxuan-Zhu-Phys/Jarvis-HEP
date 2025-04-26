#!/usr/bin/env python3self.pack 

import os, sys 
import json
import numpy as np 
import sympy as sp 
from copy import deepcopy 
from time import sleep
from base import Base
from uuid import uuid4
from inner_func import update_const, update_funcs
class Sample(Base):
    def __init__(self, params):
        self._params = params
        self._uuid = str(uuid4())
        self._likelihood = None  # Initialize likelihood with None
        self.processed = False
        self.observables = params 
        self.observables['uuid'] = self.uuid
        self._u = None

    @property
    def uuid(self):
        return self._uuid  

    @property
    def u(self):
        return self._u 

    @property
    def params(self):
        return self._params  

    @property
    def likelihood(self):
        return self._likelihood

    def update_uuid(self, uuid):
        self._uuid = uuid 
        self.observables['uuid'] = self.uuid

    @likelihood.setter
    def likelihood(self, value):
        self._likelihood = value  # 允许更新likelihood值

    def set_config(self, config):
        self.info = config
        self.create_info()

    def create_info(self):
        save_dir = self.manage_directories(self.info['sample_dirs'])
        # print(f"{self.uuid} -> save_dir is \n{save_dir}")
        self.info.update({
            "uuid": self.uuid,
            "save_dir": os.path.join(save_dir, self.uuid), 
            "run_log":  os.path.join(save_dir, self.uuid, "Sample_running.log")
        })
    
    def manage_directories(self, base_path):
        return super().manage_directories(base_path)
    
    def evaluate_output(self, outputs):
        custom_functions = update_funcs({})
        cunsom_constants = update_const({})
        result = {}
        for op in outputs: 
            if op in self.info['observables']:
                result[op] = self.info['observables'][op]
            else: 
                try:
                    symbols = {var: sp.symbols(var) for var in self.info['observables']}
                    locals_context = {**custom_functions, **cunsom_constants, **symbols}
                    expr = sp.sympify(op, locals=locals_context)
                    res = expr.subs(self.info['observables'])
                    result[op] = res
                except: 
                    raise ValueError
        return result
