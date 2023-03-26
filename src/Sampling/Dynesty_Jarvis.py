#!/usr/bin/env python3 
from calendar import day_name
from lib2to3.pgen2.token import RPAR
import logging
import os, sys
sys.path.append("/home/buding/Jarvis-HEP/src/Sampling/dynesty/py")
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
from Func_lib import decode_path_from_file
from random import randint
from sampling import Sampling_method
from dynesty import * 
from dynesty.results import print_fn 


class Dynesty(Sampling_method):
    def __init__(self) -> None:
        super().__init__()
        
    def set_config(self, cf):
        return super().set_config(cf)
    
    def set_scan_path(self, pth):
        return super().set_scan_path(pth)
    
    def initialize_generator(self, cf):
        self.cf = cf 
        self.init_logger(self.cf['logging']['scanner'])
        self.set_pars()
        print(self.pars)
        self.runing_card = os.path.join(self.path['save dir'],'dynesty_run.json')
        self.status = "INIT"
        from Func_lib import get_time_lock
        self.timelock = get_time_lock(self.cf['default setting']['sampling']['TSavePoint'])
        sys.exit()

    def set_pars(self):
        if not self.scf.has_option("Sampling_Setting", "likelihood"):
            self.exit_and_errors("\nDynesty is an importance sampling algorithm, likelihood is needed!")
        self.pars = {
            "minR":     float(self.scf.get("Sampling_Setting", "min R")),
            "likelihood":   self.scf.get("Sampling_Setting", "likelihood"),
            "cubeids":  [],
            "dims":     [],
            "ndim":     0
            }
        self.decode_sampling_variable_setting(self.scf.get("Sampling_Setting", "variables"))
        self.decode_function()        
        self.decode_selection()   

    def decode_selection(self):
        if self.scf.has_option("Sampling_Setting", "selection"):
            self.pars['selection'] = self.scf.get("Sampling_Setting", "selection")
            self.pars['filter'] = "ON"
        else:
            self.pars['selection'] = "True"
            self.pars['filter'] = "OFF"            
        
    def decode_sampling_variable_setting(self, prs):
        self.pars['vars'] = {}
        nn = 0
        from math import log10 
        for line in prs.split('\n'):
            ss = line.split(",")
            if len(ss) == 4 and ss[1].strip().lower() == "flat":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("{} + ({} - {}) * cube{}".format(float(ss[2].strip()), float(ss[3].strip()), float(ss[2].strip()), nn))
                }
                self.pars['cubeids'].append("cube{}".format(nn))
                self.pars['dims'].append(1.0)
            elif len(ss) == 5 and ss[1].strip().lower() == "flat":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("{} + ({} - {}) * cube{} / {}".format(float(ss[2].strip()), float(ss[3].strip()), float(ss[2].strip()), nn, float(ss[4].strip())))
                }                
                self.pars['dims'].append(float(ss[4].strip()))
                self.pars['cubeids'].append("cube{}".format(nn))
            elif len(ss) == 4 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "log",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("10**({} + ({} - {}) * cube{})".format(log10(float(ss[2].strip())), log10(float(ss[3].strip())), log10(float(ss[2].strip())), nn))
                }
                if not (float(ss[2].strip()) > 0 and float(ss[3].strip())):
                    self.logger.error(
                        "Illegal Variable setting of {}\n\t\t-> The lower limit and high limit should be > 0 ".format(ss))
                    sys.exit(0)
                self.pars['cubeids'].append("cube{}".format(nn))
                self.pars['dims'].append(1.0)
            elif len(ss) == 5 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "log",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("10**({} + ({} - {}) * cube{} / {})".format(log10(float(ss[2].strip())), log10(float(ss[3].strip())), log10(float(ss[2].strip())), nn, float(ss[4].strip())))
                }
                if not (float(ss[2].strip()) > 0 and float(ss[3].strip())):
                    self.logger.error(
                        "Illegal Variable setting of {}\n\t\t-> The lower limit and high limit should be > 0 ".format(ss))
                    sys.exit(0)
                self.pars['cubeids'].append("cube{}".format(nn))
                self.pars['dims'].append(float(ss[4].strip()))
            else:
                self.exit_and_errors("Illegal Variable setting in: {} ".format(line))

            nn += 1
        self.pars['ndim'] = nn


    def prerun_generator(self):
        self.sampler = dynesty.DynamicNestedSampler(
            self.evalate_likelihood, 
            self.prior_transform, 
            ndim=self.pars['ndim'], 
            bound='multi', sample='unif', rstate=rstate
            )

    # def 

    # def evalate_likelihood(self, cube):
