#!/usr/bin/env python3
from    Sampling.PyneSTY.py.dynesty.results import print_fn
from    Sampling.PyneSTY.py.dynesty import DynamicNestedSampler
from    Sampling.PyneSTY.py.dynesty import pool as dypool
import  Sampling.PyneSTY.py.dynesty.utils

# from Sampling.dynesty.py.dynesty.results import print_fn
# from Sampling.dynesty.py.dynesty import DynamicNestedSampler
# from Sampling.dynesty.py.dynesty import pool as dypool
# import Sampling.dynesty.py.dynesty.utils

import dill 
from sympy import sympify
import dynesty
import time
import json
import pandas as pd
from re import S
from dynesty import *
from Sampling.sampler import Sampling_method
from random import randint
from Func_lib import decode_path_from_file
from sympy.geometry import parabola
from sympy import sympify
from pandas.core.series import Series
from numpy.lib.function_base import meshgrid
import numpy as np
from mpmath.functions.functions import re
from abc import ABCMeta, abstractmethod
from calendar import day_name
from lib2to3.pgen2.token import RPAR
import logging
import os
import sys
pwd = os.path.abspath(os.path.dirname(__file__))
# print(dynesty._DYNASTY)
# from dynesty import plotting as dyplot
from manager import SManager
import asyncio 
from copy import deepcopy

class Dynesty(Sampling_method):
    def __init__(self) -> None:
        super().__init__()
        self.manager = None
        from multiprocessing import Lock
        # self.lock = "OPEN"
        self.lock = False

    def set_config(self, cf):
        return super().set_config(cf)

    def set_scan_path(self, pth):
        return super().set_scan_path(pth)

    def initialize_generator(self, cf):
        self.cf = cf
        self.init_logger(self.cf['logging']['scanner'])
        self.logger.info("Dynesty sampling algorithm is used in this scan")
        self.set_pars()
        self.set_run_nested_pars()
        # print(self.pars)
        self.runing_card = self.path['run_info']
        self.status = "INIT"
        from Func_lib import get_time_lock
        self.timelock = get_time_lock(
            self.cf['default setting']['sampling']['TSavePoint'])
        # sys.exit()

    def set_run_nested_pars(self):
        from copy import deepcopy
        self.info['run_nested'] = deepcopy(
            self.cf['default setting']['run nested'])
        # print(self.scf.options("Sampling_Setting"))
        if "nlive_init" in self.scf.options("Sampling_Setting"):
            val = int(self.scf.get("Sampling_Setting", "nlive_init"))
            if val < 2 * self.pars['ndim']:
                self.logger.warn("Too small nlive_init => {}, \n\tusing the default value => {}".format(
                    val, self.info['run_nested']['nlive_init']))
            else:
                self.info['run_nested']["nlive_init"] = val
        if "nlive_batch" in self.scf.options("Sampling_Setting"):
            val = int(self.scf.get("Sampling_Setting", "nlive_batch"))
            if val < 10 * self.pars['ndim']:
                self.logger.warn("Too small nlive_batch => {}, \n\tusing the default value => {}".format(
                    val, self.info['run_nested']['nlive_batch']))
            else:
                self.info['run_nested']["nlive_batch"] = val
        if "dlogz_init" in self.scf.options("Sampling_Setting"):
            val = float(self.scf.get("Sampling_Setting", "dlogz_init"))
            if val >= 0. and val <= 1.:
                self.info['run_nested']["dlogz_init"] = val
            else:
                self.logger.warn("Illigal dlogz_init => {}, \n\tthe value should be a float number in range [0.0, 1.0]\n\tJarvis using the default weigth_pfrac => {}".format(
                    val, self.info['run_nested']["dlogz_init"]))
        if "weight_pfrac" in self.scf.options("Sampling_Setting"):
            val = float(self.scf.get("Sampling_Setting", "weight_pfrac"))
            if val >= 0.0 and val <= 1.0:
                self.info['run_nested']["wt_kwargs"]['pfrac'] = val
            else:
                self.logger.warn("Illigal weight_pfrac => {}, \n\tthe value should be a float number in range [0.0, 1.0]\n\tJarvis using the default weigth_pfrac => {}".format(
                    val, self.info['run_nested']["wt_kwargs"]['pfrac']))
        if "stop_pfrac" in self.scf.options("Sampling_Setting"):
            val = float(self.scf.get("Sampling_Setting", "stop_pfrac"))
            if val >= 0.0 and val <= 1.0:
                self.info['run_nested']["stop_kwargs"]['pfrac'] = val
            else:
                self.logger.warn("Illigal stop_pfrac => {}, \n\tthe value should be a float number in range [0.0, 1.0]\n\tJarvis using the default weigth_pfrac => {}".format(
                    val, self.info['run_nested']["stop_kwargs"]['pfrac']))

    def set_pars(self):
        if not self.scf.has_option("Sampling_Setting", "likelihood"):
            self.exit_and_errors(
                "\nDynesty is an importance sampling algorithm, likelihood is needed!")
        self.pars = {
            "minR":     float(self.scf.get("Sampling_Setting", "min R")),
            "likelihood":   self.scf.get("Sampling_Setting", "likelihood"),
            "cubeids":  [],
            "dims":     [],
            "ndim":     0
        }
        if self.cf['default setting']['sampling']['use_consts']:
            from inner_func import _Constant
            self.pars['constant'] = deepcopy(_Constant)
        self.decode_sampling_variable_setting(
            self.scf.get("Sampling_Setting", "variables"))
        self.decode_function()
        self.decode_likelihood(self.pars['likelihood'])
        self.decode_selection()
        # print(self.pars)
        # print(self.fcs)
        # print(self.func)

    def decode_selection(self):
        if self.scf.has_option("Sampling_Setting", "selection"):
            self.pars['selection'] = self.scf.get(
                "Sampling_Setting", "selection")
            self.pars['filter'] = "ON"
        else:
            self.pars['selection'] = "True"
            self.pars['filter'] = "OFF"

    def decode_sampling_variable_setting(self, prs):
        self.pars['vars'] = {}
        self.pars['vname'] = []
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
                self.pars['vname'].append(ss[0].strip())
                self.pars['dims'].append(1.0)
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
                self.pars['vname'].append(ss[0].strip())
                self.pars['dims'].append(1.0)

            else:
                self.exit_and_errors(
                    "Illegal Variable setting in: {} ".format(line))

            nn += 1
        self.pars['ndim'] = nn
        self.pars['vname'] = tuple(self.pars['vname'])

    def prerun_generator(self):
        # print(self.cf['default setting']['sampling']['seed'])
        self.logger.warning("Jarvis initilized the dynamic nested sampler")
        self.path['Samples_info'] = os.path.join(
            self.path['save dir'], "Samples_info")
        self.path['ruid'] = os.path.abspath(
            os.path.join(self.path['Samples_info'], "run_uids.json"))
        self.path['ckpf'] = os.path.abspath(os.path.join(
            self.path['Samples_info'], "checkpoint_dynesty.sav"))
        self.path['slive_info'] = os.path.join(
            self.path['Samples_info'], "sinfo.json")
        self.path['alldata'] = os.path.join(
            self.path['Samples_info'], "FDIR", "AllData.csv")
        self.path['prwb'] = os.path.abspath(os.path.join(
            self.path['Samples_info'], "prwb.json"
        ))
        # self.live_sample_
        # print(self.path)
        self.status = "READY"
        self.update_sampling_status({
            "method":   "Dynesty",
            "status":   self.status
        })
        self.update_run_status("output",
                               {
                                   "ruid": self.path['ruid'],
                                   "ckpf": self.path['ckpf'],
                                   "slif": self.path['slive_info'],
                                   "prwb": self.path['prwb']
                               }
                               )
        if not os.path.exists(self.path['Samples_info']):
            os.makedirs(self.path['Samples_info'])
        self.make_sinfo()

    def evalate_likelihood(self, cube, manager):
        logl = 0.0

        from sample import Sample
        # uid = kwarg['uid']
        uid = cube[-1]
        self.logger.warning("Sampling evaluting {}".format(uid))
        
        self.samples[uid] = Sample()
        self.samples[uid].id = uid
        self.samples[uid].status = "Ready"
        self.samples[uid].pack = self.pack
        self.samples[uid].func = deepcopy(self.func)
        self.samples[uid].par = self.transfer_vars_from_array_to_dict(cube[0:-1])
        self.samples[uid].path['Samples_info'] = self.path['Samples_info']
        self.samples[uid].path['slive_info'] = self.path['slive_info']
        self.samples[uid].path['scanner_run_info'] = self.path['run_info']
        self.samples[uid].path['ruid'] = self.path['ruid']
        self.samples[uid].mana = manager
        # self.samples[uid].lock = self.lock
        self.samples[uid].likelihood = self.pars['likelihood']
        if self.cf['default setting']['sampling']['use_consts']:
            self.sampler[uid].const = self.pars['constant']
        self.samples[uid].update_dirs()
        self.samples[uid].init_logger(self.cf['logging']['scanner'])

        if self.samples[uid].status is "Ready":
            self.samples[uid].start_run()
        self.samples[uid].eval_loglike()
        # time.sleep(20)
        logl = self.samples[uid].logl
        vrs = self.samples[uid].vars

        if self.samples[uid].status == "Done":
            pass

        # print(cube[0:-1], cube[1], logl)
        # uid = str(cube[-1])
        # # cube = cube[0:-1]
        # from numpy import sin 
        # if cube[0] > cube[1]:
        #     logl = (2.0 + sin(float(cube[0]) / 2.0) * sin(float(cube[1]) / 2.0)) ** 4
        # # else:
        # #     return 0.
        self.samples.pop(uid)
        return logl

    def prior_transform(self, cube):
        v = []
        cdt = self.transfer_cube_from_array_to_dict(cube)
        for kk in range(len(self.pars['vname'])):
            expr = sympify(self.pars['vars'][self.pars['vname'][kk]]['expr'])
            v.append(expr.subs(cdt))
        # time.sleep(5)
        from Func_lib import get_sample_id
        v.append(get_sample_id())
        return np.asarray(v)

    def generate_events(self):
        self.hire_manager()
        self.logger.warning("Start sampling")
        use_pool = True
        if use_pool:
            with dypool.Pool(2, self.evalate_likelihood, self.prior_transform, manager=self.manager) as pool:
                self.sampler = DynamicNestedSampler(
                    self.evalate_likelihood,
                    self.prior_transform,
                    ndim=self.pars['ndim'],
                    **self.cf['default setting']['dynesty paras'],
                    rstate=np.random.default_rng(
                        self.cf['default setting']['sampling']['seed']),
                    queue_size=2,
                    logger=self.logger,
                    pool=pool,
                    manager=self.manager
                )
                # print(pool.loglike)
                # time.sleep(10)
                self.sampler.run_nested(
                    **self.cf['default setting']['run nested']
                )
        else:
            Sampling.dynesty.py.dynesty.utils.pickle_module = dill 
            self.sampler = DynamicNestedSampler(
                self.evalate_likelihood,
                self.prior_transform,
                ndim=self.pars['ndim'],
                **self.cf['default setting']['dynesty paras'],
                rstate=np.random.default_rng(
                    self.cf['default setting']['sampling']['seed']),
                logger=self.logger
            )
            self.sampler.run_nested(
                **self.cf['default setting']['run nested'],
                checkpoint_file=self.path['ckpf']
            )
        print("Sampler after run nested", self.info)
        self.rest = self.sampler.results

    def init_logger(self, cf):
        self.logger = logging.getLogger("Dynesty")
        handler = {
            "stream":   logging.StreamHandler(),
            "ff":       logging.FileHandler(cf['logging path'])
        }
        self.logger.setLevel(logging.DEBUG)
        handler['stream'].setLevel(logging.INFO)
        handler['ff'].setLevel(logging.DEBUG)
        self.logger.addHandler(handler['stream'])
        self.logger.addHandler(handler['ff'])
        from logging import Formatter
        handler['stream'].setFormatter(
            Formatter(cf['stream_format'],  "%m/%d %H:%M:%S"))
        handler['ff'].setFormatter(Formatter(cf['file_format']))

    def transfer_cube_from_array_to_dict(self, cube):
        cdt = {}
        for ii in range(len(cube)):
            cdt['cube{}'.format(ii)] = cube[ii]
        return cdt

    def transfer_vars_from_array_to_dict(self, v):
        cdt = {}
        print(v, self.pars['vname'])
        # sys.exit()
        for ii in range(len(v)):
            cdt[self.pars['vname'][ii]] = v[ii]
        return cdt

    def make_sinfo(self):
        from Func_lib import format_PID
        sinfo = {
            "NSpack":   self.cf['default setting']['sampling']['NSpack'],
            "NLpack":   1,
            "LPid":     format_PID(1),
            "NLPp":     0,
            "RunSP":    self.path['ruid'],
            "Status":   "FREE"
        }
        with open(self.path['slive_info'], 'w') as f1:
            json.dump(sinfo, f1, indent=4)
        with open(sinfo['RunSP'], 'w') as f1:
            json.dump({}, f1, indent=4)

    def hire_manager(self):
        self.logger.warning("Hiring task manager for sampler")
        # print(self.pack)
        self.manager = SManager()
        from Func_lib import format_PID
        sinfo = {
            "NSpack":   self.cf['default setting']['sampling']['NSpack'],
            "NLpack":   1,
            "LPid":     format_PID(1),
            "NLPp":     0,
            "RunSP":    self.path['ruid'],
            "Status":   "FREE"
        }
        self.manager.sinf.update(sinfo)
        self.manager.info['pack'] = self.pack
        self.manager.path = {
            "ruid": self.path['ruid'],
            "prwb": self.path['prwb'],
            "sinf": self.path['slive_info']
        }
        self.manager.make_factory()


        time.sleep(1)
        # sys.exit()