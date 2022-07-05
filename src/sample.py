#!/usr/bin/env python3
import imp
import logging
import os, sys 
import json


class Sample():
    def __init__(self) -> None:
        self.par = None
        self.status = "init"
        self.pack = None
        self.id = None
        self.vars = {}
        self.path = {}
        self.info = {}
        self.worker = {}
        self.children = []
        self.mother = None
        self.func = None
        self.expr = None
        self.likelihood = None

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
    
    def close_logger(self):
        self.handler['ff'].close()
        # self.vars['ID'] = str(self.vars['ID'])
        # with open(self.path['datapath'], 'w') as f1:
            # json.dump(self.vars, f1, indent=2)
        import pandas as pd 
        self.vars = pd.Series(self.vars)
        self.vars.to_csv(self.path['datapath'])
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
        for pkg in pkgs:
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
            self.worker[pkg].config = deepcopy(self.pack['include'][pkg])
            self.worker[pkg].config['name'] = deepcopy(pkg)
            self.worker[pkg].vars   = dict(deepcopy(self.vars))
            self.worker[pkg].path   = deepcopy(self.path)
            self.worker[pkg].config['paraller number'] = deepcopy(self.pack['config']['paraller number'])
            self.worker[pkg].init()
    
    def update_variable(self, vars):
        for kk, vv in vars.items():
            if kk not in self.vars.keys():
                self.vars[kk] = vv
                
    def evaluate_likelihood(self):
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
                self.status = "Done"
            except:
                self.logger.error("Likelihood expression can not be evaluated")
                self.status = "Stoped"
        else:
            self.status = "Done"
        
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
            from pandas import Series
            self.logger.warning("The output variables are summarized as followed: \n\n{}\n".format(Series(self.vars)))
            self.evaluate_likelihood()                
            self.close_logger()
        
        if self.status != self.par['Status']:
            self.par['Status'] = self.status
            self.logger.info("Sample status update into  {}".format(self.status))
        