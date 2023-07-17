#!/usr/bin/env python3self.pack 

import logging
import os, sys 
import json
import numpy as np 
from pandas import Series
from sympy.core import numbers as SCNum
from inner_func import _AllSCNum
from sympy import *
from copy import deepcopy 
from time import sleep

class Sample():
    def __init__(self) -> None:
        self.par    = None
        self.status = "init"
        self.pack   = None
        self.id     = None
        self.vars   = {}
        self.path   = {}
        self.info   = {}
        self.const  = None
        self.worker = {}
        self.children = []
        self.mother = None
        self.func   = None
        self.expr   = None
        self.likelihood = None
        self.fdir   = None
        self.mana   = None

    def update_dirs(self):
        # self.lock.acquire(timeout=2)
        # while self.lock != "OPEN":
        #     sleep(0.00001)

        self.get_fdir()
        self.path['info'] = os.path.abspath(os.path.join(
            self.path['Samples_info'], self.fdir, self.id
        ))
        self.info['RunWeb'] = {
            "status":   self.status,
            "fdir":     self.fdir,
            "path":     self.path['info']
        }

        if not os.path.exists(self.path['info']):
            os.makedirs(self.path['info'])
        self.mana.ruid[self.id] = self.info['RunWeb']
        from IOs import to_file_woblk
        to_file_woblk(dict(self.mana.ruid), self.mana.path['ruid'], method="to_json")
        # self.lock.release()
        # self.lock = "OPEN"
        # print("sample 55: ", dict(self.mana.sinf))
        # self.mana.sinf['Status'] = "FREE"



    def delete_uid_dirs(self):
        self.mana.ruid.pop(self.id)
        from IOs import to_file_woblk
        to_file_woblk(dict(self.mana.ruid), self.mana.path['ruid'], method='to_json')


    def get_fdir(self):
        # self.lock = "CLOSE"
        # while self.mana.sinf['Status'] != "FREE":
        #     print(self.id, " is waiting ...")
        #     sleep(0.00001)
            
        self.mana.sinf['Status'] = self.id
        if self.mana.sinf['NLPp'] < self.mana.sinf['NSpack']:
            self.mana.sinf['NLPp'] += 1
            self.fdir = deepcopy(self.mana.sinf['LPid'])
        else:
            from Func_lib import format_PID
            self.mana.sinf['NLPp'] = 0
            self.mana.sinf['NLpack'] += 1 
            self.mana.sinf['LPid'] = format_PID(self.mana.sinf['NLpack'])
        from IOs import to_file_woblk
        to_file_woblk(dict(self.mana.sinf), self.mana.path['sinf'], method="to_json")
        # self.mana.sinf['Lock'] = "OPEN"

    def get_par(self, pars):
        pass
    
    def init_logger(self, cf):
        self.logger = logging.getLogger("Sample@id: {} ".format(self.id))
        self.path['logpath'] = os.path.join(self.path['info'], "running.log")
        self.path['datapath'] = os.path.join(self.path['info'], "data.csv")
        self.path['tempath'] = os.path.join(self.path['info'], "temp")
        self.handler = {
            "stream":   logging.StreamHandler(),
            "ff":       logging.FileHandler(self.path['logpath'])
        }
        self.logger.setLevel(logging.DEBUG)
        self.handler['stream'].setLevel(logging.WARNING)
        self.handler['ff'].setLevel(logging.DEBUG)
        self.logger.addHandler(self.handler['stream'])
        self.logger.addHandler(self.handler['ff'])
        from logging import Formatter
        self.handler['stream'].setFormatter(Formatter(cf['stream_format'],  "%m/%d %H:%M:%S"))
        self.handler['ff'].setFormatter(Formatter(cf['file_format']))
        self.logger.warning("Sample created in the disk")
        self.logger.info("\nSample info: {}\n".format(self.par))
        self.logger.info("Current status is {}".format(self.status))
        # print("sample 108 =>", self.mana.sinf)

        # self.mana.pack['TestFunction']['workers'].update({"001": self.id})
        # print("sample 108 =>", self.mana.pack['TestFunction']['workers'])
        # import time 
        # time.sleep(4)
    
    def close_logger(self):
        self.handler['ff'].close()
        import pandas as pd 
        self.vars = pd.Series(self.vars)
        self.vars.to_csv(self.path['datapath'])
        if os.path.exists(self.path['tempath']):
            from shutil import rmtree
            rmtree(self.path['tempath'])
    
    def start_run(self):
        self.status = "Running"
        self.par['Status'] = "Running"
        self.info['current_depth'] = 0
        from copy import deepcopy
        self.vars = dict(deepcopy(self.par))
        self.vars.pop("Status")
        if not os.path.exists(self.path['tempath']):
            os.makedirs(self.path['tempath'])
        with open(self.path['scanner_run_info'], 'r') as f1:
            temp = json.loads(f1.read())
            self.path['scanner_logging_path'] = temp['log']['scanner_logging_path']
            self.info['log'] = {
                "scanner_format":   temp['sp']['cf']['logging']['scanner']['file_format'],
                "stream_format":    temp['sp']['cf']['logging']['package']['stream_format'],
                "file_format":      temp['sp']['cf']['logging']['package']['file_format']
            }
        self.run_next_layer_function(0)
        
    def run_next_layer_function(self, dps):
        pkgs = self.pack['tree'].layer[dps]
        from program import program
        for pname in pkgs:
            if "\n@" in pname:
                pkg, psect = pname.split("\n@")
                self.worker[pkg] = program()
                self.worker[pkg].set_logger_setting({
                    "name": "\n  {} @ {}".format(pkg, self.id),
                    "scanner_logging_path": self.path['scanner_logging_path'],
                    "ff_logging_path":      self.path['logpath'],
                    "stream_format":        self.info['log']['stream_format'],
                    "scanner_format":       self.info['log']['scanner_format'],
                    "file_format":          self.info['log']['file_format']
                })
                from copy import deepcopy
                print("Sample 157:", self.id, pname, self.mana)
                self.worker[pkg].suid   = deepcopy(self.id)
                self.worker[pkg].config = deepcopy(self.pack['include'][psect])
                self.worker[pkg].config['name'] = deepcopy(pkg)
                self.worker[pkg].vars   = dict(deepcopy(self.vars))
                self.worker[pkg].path   = deepcopy(self.path)
                self.worker[pkg].mana   = self.mana
                self.worker[pkg].config['paraller number'] = deepcopy(self.pack['config']['paraller number'])
                self.worker[pkg].init()            
            else:
                pkg = pname
                self.worker[pkg] = program()
                self.worker[pkg].set_logger_setting({
                    "name": "\n  {} @ {}".format(pkg, self.id),
                    "scanner_logging_path": self.path['scanner_logging_path'],
                    "ff_logging_path":      self.path['logpath'],
                    "stream_format":        self.info['log']['stream_format'],
                    "scanner_format":       self.info['log']['scanner_format'],
                    "file_format":          self.info['log']['file_format']
                })
                from copy import deepcopy
                print("Sample 157:", self.id, pname, self.mana)

                self.worker[pkg].suid = deepcopy(self.id)
                self.worker[pkg].config = deepcopy(self.pack['include'][pkg])
                self.worker[pkg].config['name'] = deepcopy(pkg)
                self.worker[pkg].vars   = dict(deepcopy(self.vars))
                self.worker[pkg].path   = deepcopy(self.path)
                self.worker[pkg].mana   = self.mana
                self.worker[pkg].config['paraller number'] = deepcopy(self.pack['config']['paraller number'])
                self.worker[pkg].init()
        # sys.exit()
    
    def update_variable(self, vars):
        for kk, vv in vars.items():
            if kk not in self.vars.keys():
                self.vars[kk] = vv
                
    def evaluate_likelihood(self):
        # print(self.likelihood)
        if self.likelihood is not None:
            from sympy import sympify
            from inner_func import update_funcs
            try:
                self.func = update_funcs(self.func)
                expr = sympify(self.likelihood)
                expr = expr.subs(self.expr)
                expr = expr.subs(self.vars)
                expr = eval(str(expr), self.func)
                self.likelihood = expr 
                self.vars['Likelihood'] = expr 
                self.logger.warning("Likelihood is evaluating as -> {}".format(expr))
                # self.delete_uid_dirs()
                self.status = "Done"
            except:
                self.logger.error("Likelihood expression can not be evaluated")
                # self.delete_uid_dirs()
                self.status = "Stoped"
        else:
            self.status = "Done"
        
    def eval_loglike(self):
        '''
        Evaluation loglikelihood for external sampling algorithm (like the dynesty in this version)

        Parameters 
        -----------
        parmeters in self:
        self.pars: dict
            all parameters of the sample, updating in the evluation progress.

        self.func: dict 
            A function list.
        
        self.pack: dict 
            A list of HEP package from external

        Returns
        ___________
        This function will not return a value. Please not using this function in your code. 
        -----------
        self.status: str
            updating the status of the sample.
            "Ready": Sample is created in the disk, before call the package.
            "Running": Sample is calling the HEP package.
            "Erroring": Find error when calling the HEP package.
            "Stoped": Force stopping all package in calculation, and resetted all package called. 
            "Finish": End calling the package, all calculation of HEP package is done. 

        self.logl: float
            The value of log(likelihood) after self.status is "Finish", else this value will return -np.inf
        '''
        self.logl = -np.inf
        from time import sleep 
        while self.status not in ["Done", "Finish", "Stoped"]:
            self.update_status()
        if self.status == "Finish":
            try:
                self.logger.warning("Calculation Finished")
                from sympy import sympify
                expr = sympify(self.likelihood['expression'])
                if self.expr is not None:
                    expr = expr.subs(self.expr)
                if self.vars is not None:
                    expr = expr.subs(self.vars)
                if self.const is not None:
                    expr = expr.subs(self.const)
                if not isinstance(expr.evalf(), _AllSCNum):
                    expr = expr.subs(self.func)
                    expr = eval(str(expr), self.func)
                if isinstance(expr.evalf(), _AllSCNum):
                    self.logl = float(expr)
                self.vars.update({"LogL": self.logl, "ID": self.id})
                self.logger.warning("The output variables are summarized as followed: \n\n{}\n".format(Series(self.vars)))
                self.status = "Done"
            except:
                self.logger.error("Errors in evaluation loglike, updating the sample status into Stoped.")
                self.logger.info("log(likelihood) => {},\n\texpression => {},\n\tvariables => {},\n\tconstants => {},\n\tlogl => {}".format(self.likelihood['expression'], self.expr, self.vars, self.logl))
                self.status = "Stoped"
            self.delete_uid_dirs()
        self.close_logger()

    def update_status(self):
        if self.status == "Erroring":
            errt = False
            for pkg, worker in self.worker.items():
                worker.update_status()
                if not worker.status == "done":
                    errt = True
            if not errt:
                self.status = "Stoped"
                self.close_logger()
        elif self.status == "Running":
            nxtt = False
            errt = False
            for pkg, worker in self.worker.items():
                worker.update_status()
                if not worker.status == "done":
                    nxtt = True
                if worker.status == "done":
                    self.update_variable(worker.vars)
                if worker.status == "error":
                    errt = True
            if errt:
                self.logger.error("Calculation stopped by error!")
                self.status = "Erroring"
                for pkg, worker in self.worker.items():
                    worker.stop_calculation()
            if not nxtt and self.status != "Erroring":
                if self.info['current_depth'] < len(self.pack['tree'].layer) - 1:
                    self.info['current_depth'] += 1
                    self.run_next_layer_function(self.info['current_depth'])
                elif self.info['current_depth'] == len(self.pack['tree'].layer) - 1:
                    self.status = "Finish"
        elif self.status == "Finish":
            self.logger.warning("Calculation Finished")
            self.logger.warning("The output variables are summarized as followed: \n\n{}\n".format(Series(self.vars)))
            self.evaluate_likelihood()                
            self.close_logger()
        
        if self.status != self.par['Status']:
            self.par['Status'] = self.status
            self.logger.info("Sample status update into  {}".format(self.status))
        
