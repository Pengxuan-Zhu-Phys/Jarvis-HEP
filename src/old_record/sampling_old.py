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
from old_record.Func_lib import decode_path_from_file
from random import randint
import json
from base import Base

pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 


class SamplingVirtial(Base):
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        super().__init__()
        self.info           = {}
        self.status         = 'FREE'
        self.runing_card    = None
        self.logger         = None
        self.path           = {}
        self.pack           = None
        self.samples        = {}
        self.timelock       = None 
        self.run_info       = {}
        self.manager        = None
        self.lock           = "OPEN"
        
    
    @abstractmethod
    def set_config(self, cf):
        self.scf = cf

    @abstractmethod
    def set_scan_path(self, pth):
        self.path['save dir'] = pth

    @abstractmethod
    def generate_events(self):
        pass 

    @abstractmethod
    def initialize_generator(self):
        pass

    @abstractmethod
    def subprocess_multi_core_tag(self):
        pass 

    @abstractmethod
    def afterrun_generator(self):
        pass 

    @abstractmethod
    def prerun_generator(self):
        pass 

    @abstractmethod
    def resume_generator(self):
        pass 

    @abstractmethod
    def subpstatus(self):
        pass 
    
    @abstractmethod
    def init_logger(self, cf):
        pass 
    
    @abstractmethod
    def decode_sampling_variable_setting(self, prs):
        pass 
    
    @abstractmethod
    def decode_likelihood(self, lik):
        if lik[0:4] == "CHI2":
            self.pars['likelihood'] = {
                "method":   "chisquare",
                "fcsinc":   [],
                "expression":   lik[5:].strip()
            } 
        elif lik[0:4] == "LIKI":
            self.pars['likelihood'] = {
                "method":   "likelihood",
                "fcsinc":   [],
                "expression":   lik[5:].strip()
            }
        else:
            self.pars['likelihood'] = {
                "method":   "user_defined",
                "fcsinc":   [],
                "expression":   lik
            }
        if "&FC" in self.pars['likelihood']['expression']:
            for item in self.exprs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['expression'] = self.pars['likelihood']['expression'].replace(item, self.exprs[item]['expr'])
            for item in self.fcs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['fcsinc'].append(item)


    @abstractmethod
    def decode_function(self):
        self.fcs = {}
        self.exprs = {}
        self.func = {}
        for sc in self.scf.sections():
            if "Function" == sc[0:8] and self.scf.get(sc, "method") == "expression":
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                self.exprs[name] = {
                    "expr":   self.scf.get(sc, "expression")
                }
            elif "Function" == sc[0:8] and self.scf.get(sc, "method") == "interpolate 1d":
                from pandas import read_csv
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                from old_record.Func_lib import decode_path_from_file
                data = read_csv(decode_path_from_file(self.scf.get(sc, "file")))
                self.fcs[name] = {
                    "data":      data 
                }
                if self.scf.has_option(sc, "kind"):
                    method = self.scf.get(sc, "kind").strip().lower()
                    if method not in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']:
                        method = "linear"
                else:
                    method = "linear"
                fill_value = ""
                if self.scf.has_option(sc, "fill_value"):
                    fill_value = self.scf.get(sc, "fill_value").strip().lower()
                from scipy.interpolate import interp1d
                if fill_value == "extrapolate":
                    self.fcs[name]['expr'] = interp1d(data['x'], data['y'], kind=method, fill_value=fill_value)
                else:
                    self.fcs[name]['expr'] = interp1d(data['x'], data['y'], kind=method)    
        from inner_func import update_funcs
        for kk, vv in self.fcs.items():
            self.func[kk] = vv['expr']
        self.func = update_funcs(self.func)

    @abstractmethod
    def check_generator_status(self):
        pass 
    
    @abstractmethod
    def get_empty_data_row(self):
        raw = {
            "ID":   None,
            "Status":   "Free"
        }
        for var, item in self.pars['vars'].items():
            raw[var] = None
        self.pars['emptyData'] = Series(raw)
        
    @abstractmethod
    def update_sampling_status(self, sts):
        with open(self.path['run_info'], 'r') as f1:
            self.run_info = json.loads(f1.read())
        if "sampling" not in self.run_info.keys():
            self.run_info['sampling'] = {}
        self.run_info['sampling'].update(sts)
        with open(self.path['run_info'], 'w') as f1:
            json.dump(self.run_info, f1, indent=4)

    @abstractmethod
    def update_run_status(self, block, udt):
        with open(self.path['run_info'], 'r') as f1:
            self.run_info = json.loads(f1.read())
        if block not in self.run_info.keys():
            self.run_info[block] = {}
        self.run_info[block].update(udt)
        with open(self.path['run_info'], 'w') as f1:
            json.dump(self.run_info, f1, indent=4)

    @abstractmethod
    def exit_and_errors(self, errs="Unexpected Error! Jarvis exit the current jobs!"):
        try:
            self.logger.error(errs)
        except:
            print(errs)
        sys.exit(1)

    @abstractmethod
    def hire_manager(self):
        pass 
