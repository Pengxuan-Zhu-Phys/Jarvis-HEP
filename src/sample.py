#!/usr/bin/env python3self.pack 

import logging
import os, sys 
import json
import numpy as np 
from pandas import Series
from sympy.core import numbers as SCNum
from sympy import *
from copy import deepcopy 
from time import sleep
from base import Base
from uuid import uuid4

class Sample(Base):
    def __init__(self, params):
        self._params = params
        self._uuid = str(uuid4())
        self._likelihood = None  # 初始化likelihood为None
        # 可以添加更多属性
        self.processed = False
        self.observables = params 
        self.observables['uuid'] = self.uuid

    @property
    def uuid(self):
        return self._uuid  # 提供只读访问

    @property
    def params(self):
        return self._params  # 提供对data的只读访问，如果需要

    @property
    def likelihood(self):
        return self._likelihood

    @likelihood.setter
    def likelihood(self, value):
        self._likelihood = value  # 允许更新likelihood值

    def set_config(self, config):
        self.info = config
        self.create_info()

    def create_info(self):
        save_dir = self.manage_directories(self.info['sample_dirs'])
        print(f"{self.uuid} -> save_dir is \n{save_dir}")
        self.info.update({
            "uuid": self.uuid,
            "save_dir": os.path.join(save_dir, self.uuid), 
            "run_log":  os.path.join(save_dir, self.uuid, "Sample_running.log")
        })
    
    def manage_directories(self, base_path):
        return super().manage_directories(base_path)