#!/usr/bin/env python3 
from calendar import day_name
from lib2to3.pgen2.token import RPAR
import logging
import os, sys
from re import S 
from abc import ABCMeta, abstractmethod
from mpmath.functions.functions import re
import pandas as pd 
from numpy.lib.function_base import meshgrid
from pandas.core.series import Series
from sympy import sympify
from sympy.geometry import parabola
import time 
from Func_lib import decode_path_from_file



class Sampling_method():
    __metaclass__ = ABCMeta
    def __init__(self) -> None:
        self.info           = {}
        self.status         = 'FREE'
        self.runing_card    = None
        self.logger         = None
        self.path           = {}
        self.pack           = None
        self.samples        = {}
        self.timelock       = None 
        
    
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
        for sc in self.scf.sections():
            if "Function" == sc[0:8] and self.scf.get(sc, "method") == "expression":
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                self.exprs[name] = {
                    "expr":   self.scf.get(sc, "expression")
                }
            elif "Function" == sc[0:8] and self.scf.get(sc, "method") == "interpolate 1d":
                from pandas import read_csv
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                from Func_lib import decode_path_from_file
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
        
    
class Grid(Sampling_method):
    def __init__(self) -> None:
        super().__init__()
        
    def set_config(self, cf):
        return super().set_config(cf)
    
    def set_scan_path(self, pth):
        return super().set_scan_path(pth)
    
    def initialize_generator(self, cf):
        self.cf = cf 
        self.init_logger(self.cf['logging']['scanner'])
        self.logger.info("Grid Sampling Algorithm is used in this scan")
        self.set_pars()
        self.runing_card = os.path.join(self.path['save dir'], "grid_run.json")
        self.status = "INIT"
        from Func_lib import get_time_lock
        self.timelock = get_time_lock(self.cf['default setting']['sampling']['TSavePoint'])
        
    def decode_sampling_variable_setting(self, prs):
        self.pars['vars'] = {}
        from math import log10
        for line in prs.split('\n'):
            ss = line.split(',')
            if len(ss) == 4 and ss[1].strip().lower() == "flat":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "Nbin":     self.cf['default setting']['sampling']['BinNum']
                }
            elif len(ss) == 5 and ss[1].strip().lower() == "flat":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "Nbin":     int(ss[4].strip())
                }
            elif len(ss) == 4 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "log",
                    "min":      log10(float(ss[2].strip())),
                    "max":      log10(float(ss[3].strip())),
                    "Nbin":     self.cf['default setting']['sampling']['BinNum']
                }
            elif len(ss) == 5 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "log",
                    "min":      log10(float(ss[2].strip())),
                    "max":      log10(float(ss[3].strip())),
                    "Nbin":     int(ss[4].strip())
                }
            elif len(ss) == 4 and ss[1].strip().lower() == "count":
                    self.pars['vars'][ss[0].strip()] = {
                    "prior":    "count",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "step":     self.cf['default setting']['sampling']['CountStep']
                }
            elif len(ss) == 5 and ss[1].strip().lower() == "count":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "count",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "step":     float(ss[4].strip())
                }
             
    def resume_generator(self, rerun=None):
        self.get_empty_data_row()
        self.status = "RESUME"
        self.path['Samples_info']   = os.path.join(self.path['save dir'], "Samples_info")
        self.path['datadir']        = os.path.join(self.path["save dir"], "GridData.csv")
        self.path['datapuredir']    = os.path.join(self.path['save dir'], "GridData_pure.csv")
        self.path['dataalldir']     = os.path.join(self.path['save dir'], "AllData.csv")
        self.path['pointdir']       = os.path.join(self.path['save dir'], "Point.csv")     
        self.pars['points']         = pd.read_csv(self.path['pointdir'])
        self.combine_data()
        self.check_samples_status_number(True)              
        self.logger.info("Grid reset running samples to free status ...")
        self.reset_running_sample_to_free()
        self.check_samples_status_number(True)   
        if rerun is not None:
            self.clear_rerun_sample_data(rerun)
        
        if os.path.exists(self.path['datadir']):
            from shutil import move, copyfile
            tpath = os.path.join(os.path.dirname(self.path['datadir']), ".Running", "{}_{}".format(os.path.basename(self.path['datadir']), self.timelock ))
            move(self.path['datadir'], tpath)
        if os.path.exists(self.path['datapuredir']):
            from shutil import move, copyfile
            tpath = os.path.join(os.path.dirname(self.path['datapuredir']), ".Running", "{}_{}".format(os.path.basename(self.path['datapuredir']), self.timelock ))
            move(self.path['datapuredir'], tpath)
        if os.path.exists(self.path['dataalldir']):
            from shutil import move, copyfile
            tpath = os.path.join(os.path.dirname(self.path['dataalldir']), ".Running", "{}_{}".format(os.path.basename(self.path['dataalldir']), self.timelock ))
            move(self.path['dataalldir'], tpath)            
            
        self.info['nrunningcore']   = 0
        self.pars['Data'] = []
        self.pars['AllData'] = []
        self.status = "READY"         

    def clear_rerun_sample_data(self, rerun):
        self.logger.info("Cleanning the sample data with rerun tag: {}".format(rerun))
        allpoints = pd.read_csv(self.path["dataalldir"])
        gridpoint = pd.read_csv(self.path['datadir'])
        purepoint = pd.read_csv(self.path['datapuredir'])
        for tag in rerun:
            ids = self.pars['points'].loc[self.pars['points']["Status"] == tag]
            discard = "|".join(list(ids.ID))
            if discard:
                allpoints = allpoints[ allpoints['ID'].str.contains(discard) == False ]
                gridpoint = gridpoint[ gridpoint['ID'].str.contains(discard) == False ]
                purepoint = purepoint[ purepoint['ID'].str.contains(discard) == False ]
            
        allpoints.to_csv(self.path['dataalldir'], index=False)
        gridpoint.to_csv(self.path['datadir'], index=False)
        purepoint.to_csv(self.path['datapuredir'], index=False)

    def reset_running_sample_to_free(self):
        while self.info['nsample']['runn'] > 0:
            sp = self.pars['points'].loc[self.pars['points']['Status'] == "Running"].iloc[0]
            spinfopath = os.path.join(self.path['Samples_info'], str(sp.ID))
            from shutil import rmtree
            rmtree(spinfopath)
            sp.Status = "Free" 
            self.pars['points'].iloc[self.pars['points']['ID'] == sp.ID] = sp 
            self.check_samples_status_number(False)
        while self.info['nsample']['redy'] > 0:
            sp = self.pars['points'].loc[self.pars['points']['Status'] == "Ready"].iloc[0]
            spinfopath = os.path.join(self.path['Samples_info'], str(sp.ID))
            from shutil import rmtree
            rmtree(spinfopath)
            sp.Status = "Free" 
            self.pars['points'].iloc[self.pars['points']['ID'] == sp.ID] = sp 
            self.check_samples_status_number(False)
        while self.info['nsample']['fini'] > 0:
            sp = self.pars['points'].loc[self.pars['points']['Status'] == "Finish"].iloc[0]
            spinfopath = os.path.join(self.path['Samples_info'], str(sp.ID))
            from shutil import rmtree
            rmtree(spinfopath)
            sp.Status = "Free" 
            self.pars['points'].iloc[self.pars['points']['ID'] == sp.ID] = sp 
            self.check_samples_status_number(False)
        
        
        self.pars['points'].to_csv(self.path['pointdir'], index=False)
        # print(len(os.listdir(self.path['Samples_info'])))
        # for ff in os.listdir(self.path['Samples_info']):
            # print(self.pars['points'].loc[self.pars['points']['ID'] == ff])
          
    def prerun_generator(self):
        self.get_empty_data_row()
        data = []
        name_list = []
        from numpy import arange, linspace, logspace, vstack, meshgrid
        for var, item in self.pars['vars'].items():
            if item['prior'] == "count":
                data.append(arange(item['min'], item['max'], item['step']))
                name_list.append(var)
            elif item['prior'] == "flat":
                data.append(linspace(item['min'], item['max'], item['Nbin']))
                name_list.append(var)
            elif item['prior'] == "log":
                data.append(logspace(item['min'], item['max'], item['Nbin']))
                name_list.append(var)
        res = vstack(meshgrid(*data)).reshape(len(data), -1).T
        
        points = []
        from copy import deepcopy
        from Func_lib import get_sample_id
        for item in res:
            cube = deepcopy(self.pars['emptyData'])
            cube['ID'] = get_sample_id()
            for ii, name in enumerate(name_list):
                cube[name] = item[ii]
            if self.pars['selection']['tag']:
                if self.pars['selection']['expr'].subs(dict(cube)):
                    points.append(cube)
            else:
                points.append(cube)
        self.pars["points"] = pd.DataFrame(points)
        self.pars["Data"] = []
        self.pars['AllData'] = []
        # print(self.pars['points'].iloc[0].ID, self.pars['points'].shape)
        self.status = "READY"
        self.check_samples_status_number(output=True)
        self.info['nrunningcore']   = 0
        self.path['Samples_info']   = os.path.join(self.path['save dir'], "Samples_info")
        self.path['datadir']        = os.path.join(self.path["save dir"], "GridData.csv")
        self.path['datapuredir']    = os.path.join(self.path['save dir'], "GridData_pure.csv")
        self.path['dataalldir']     = os.path.join(self.path['save dir'], "AllData.csv")
        self.path['pointdir']       = os.path.join(self.path['save dir'], "Point.csv")
        self.pars['points'].to_csv(self.path['pointdir'], index=False)
        self.record_points()
        if not os.path.exists(self.path['Samples_info']):
            os.makedirs(self.path['Samples_info'])
        
    def set_pars(self):
        pars = {
                "variable": self.scf.get("Sampling_Setting", "variables")
        }
        self.pars = {
            "selection":   {
                "tag":  False,
                "expr": ""
            }
        }
        self.decode_sampling_variable_setting(pars['variable'])
        self.decode_function()
        if self.scf.has_option("Sampling_Setting", "selection"):
            from sympy import sympify
            self.pars['selection'] = {
                "tag":  True,
                "expr": sympify(self.scf.get("Sampling_Setting", "selection"))
            }
        if self.scf.has_option("Sampling_Setting", "likelihood"):
            pars['likelihood'] = self.scf.get("Sampling_Setting", "likelihood")
            self.decode_likelihood(pars['likelihood'])        

    def init_logger(self, cf):
        self.logger = logging.getLogger("Grid")
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
        handler['stream'].setFormatter(Formatter(cf['stream_format'],  "%m/%d %H:%M:%S"))
        handler['ff'].setFormatter(Formatter(cf['file_format']))

    def generate_events(self):
        while not self.status == "FINISH":
            if self.status == "READY":
                self.status = "RUNNING"
                self.get_new_sample()
            elif self.status == 'RUNNING':
                if self.info['nrunningcore'] < self.pack['config']['paraller number']:
                    self.check_running_samples_and_run_next_step()
                    if self.info['nsample']['redy'] > 0 and self.info['nsample']['runn'] < self.pack['config']['paraller number']:
                        self.ready_sample_start_run()
                    elif self.info['nsample']['free'] > 0 and self.info['nsample']['runn'] < self.pack['config']['paraller number']:
                        self.get_new_sample()
                else:
                    self.status = "FINISH"    
            self.check_generator_status()   
            self.record_points()

    def check_running_samples_and_run_next_step(self):
        doneid  = []
        stopid  = []  
        finishid= []
        for kk, smp in self.samples.items():
            if smp.status not in ['Done', "Stoped"]:
                smp.update_status()
                if smp.status == "Done":
                    doneid.append(smp.id)
                elif smp.status == "Stoped":
                    stopid.append(smp.id)                
                if smp.status == "Finish":
                    finishid.append(smp.id)
        if stopid:
            for sid in stopid:
                smp = self.pars['points'].loc[self.pars['points']['ID'] == sid].iloc[0]
                if self.samples[sid].par['Status'] != smp['Status']:
                    self.pars['points'].iloc[self.pars['points']['ID'] == sid] = self.samples[sid].par
            for sid in stopid:
                self.pars['AllData'].append(self.samples[sid].vars)
                self.samples.pop(sid)
        if doneid:
            for sid in doneid:
                smp = self.pars['points'].loc[self.pars['points']['ID'] == sid].iloc[0]
                if self.samples[sid].par['Status'] != smp['Status']:
                    self.pars['points'].iloc[self.pars['points']['ID'] == sid] = self.samples[sid].par
            for sid in doneid:
                self.pars["Data"].append(self.samples[sid].vars)
                self.pars['AllData'].append(self.samples[sid].vars)
                self.samples.pop(sid)
        if finishid:
            for sid in finishid:
                smp = self.pars['points'].loc[self.pars['points']['ID'] == sid].iloc[0]
                if self.samples[sid].par['Status'] != smp['Status']:
                    self.pars['points'].iloc[self.pars['points']['ID'] == sid] = self.samples[sid].par
        if stopid + finishid + doneid:
            self.check_samples_status_number(True)
            self.record_points()

    def record_points(self):
        from Func_lib import get_time_lock
        tl = get_time_lock(self.cf['default setting']['sampling']['TSavePoint'])
        if tl != self.timelock:
            self.logger.info("Backup running points information to disk ... ")
            self.timelock = tl 
            self.afterrun_generator()
            
    def afterrun_generator(self):
        if len(self.pars['Data']) > 0:
            self.pars['points'].to_csv(self.path['pointdir'], index=False)
            from copy import deepcopy
            self.pars['SData'] = deepcopy(self.pars['Data'])
            self.pars['SData'] = pd.DataFrame(self.pars['SData'])
            self.pars['SData'].to_csv(self.path['datadir'], index=False)
            self.pars['SAllData'] = deepcopy(self.pars['AllData'])
            self.pars['SAllData'] = pd.DataFrame(self.pars['SAllData'])
            self.pars['SAllData'].to_csv(self.path['dataalldir'], index=False)
            self.pars['PureData'] = deepcopy(self.pars['SData'])
            strout = []
            for stt in self.pack['stroutput']:
                if stt in self.pars['SData'].columns:
                    strout.append(stt)
            self.pars['PureData'] = self.pars['PureData'].drop(strout, axis=1)
            self.pars['PureData'].to_csv(self.path['datapuredir'], index=False)
 
    def combine_data_file(self, dname, pname):
        rpth  = os.path.join(self.path['save dir'], ".Running")
        cmbd  = []
        for ff in os.listdir(rpth):
            if dname in ff:
                cmbd.append(os.path.join(rpth, ff))
        
        if cmbd:
            self.logger.info("Combining data files '{}'".format(dname))
            if os.path.exists(self.path[pname]):
                cmbd.append(self.path[pname])
            gds = []
            for ff in cmbd:
                ds = pd.read_csv(ff)
                gds.append(ds)
                
    def combine_data(self):
        rpth  = os.path.join(self.path['save dir'], ".Running")
        print(rpth)
        gd    = []
        gdp   = []
        ad    = []
        for ff in os.listdir(rpth):
            fp = os.path.join(rpth, ff)
            if "GridData.csv" in ff:
                gd.append(fp)
            elif "GridData_pure.csv" in ff:
                gdp.append(fp)
            elif "AllData.csv" in ff:
                ad.append(fp)
                
        if gd + gdp + ad:        
            self.logger.info("Combining data files")
            if gd:
                if os.path.exists(self.path['datadir']):
                    gd.append(self.path['datadir'])
                gds = []
                for ff in gd:
                    ds = pd.read_csv(ff)
                    gds.append(ds)
                    os.remove(ff)
                df = pd.concat(gds)
                df.to_csv(self.path['datadir'], index=False)

            if gdp:
                if os.path.exists(self.path['datapuredir']):
                    gdp.append(self.path['datapuredir'])
                gds = []
                for ff in gdp:
                    ds = pd.read_csv(ff)
                    gds.append(ds)
                    os.remove(ff)
                df = pd.concat(gds)
                df.to_csv(self.path['datapuredir'], index=False)

            if ad:
                if os.path.exists(self.path['dataalldir']):
                    ad.append(self.path['dataalldir'])
                gds = []
                for ff in ad:
                    ds = pd.read_csv(ff)
                    gds.append(ds)
                    os.remove(ff)
                df = pd.concat(gds)
                df.to_csv(self.path['dataalldir'], index=False)
            
    def ready_sample_start_run(self):
        sampleid = []
        for kk, smp in self.samples.items():
            if smp.status == "Ready":
                sampleid.append(smp.id) 
        # print(sampleid)
        if sampleid:
            for sid in sampleid:
                self.samples[sid].start_run()
                smp = self.pars['points'].loc[self.pars['points']['ID'] == sid].iloc[0]
                if self.samples[sid].par['Status'] != smp['Status']:
                    self.pars['points'].iloc[self.pars['points']['ID'] == sid] = self.samples[sid].par
            self.check_samples_status_number(True)
            self.record_points()
        
    def helper_func_update_points_status(self, sid, status):
        self.pars['points'].iloc[self.pars['points']['ID'] == sid] = self.samples[sid].par

    def find_free_point_and_make_sample(self):
        new = self.pars['points'].loc[self.pars['points']['Status'] == "Free"].iloc[0]
        new.Status = "Ready"
        self.pars['points'].iloc[self.pars['points']['ID'] == new.ID] = new
        self.check_samples_status_number(True)        
        from sample import Sample
        self.samples[new.ID] = Sample()
        self.samples[new.ID].id = new.ID 
        self.samples[new.ID].status = "Ready"
        self.samples[new.ID].pack = self.pack         
        self.samples[new.ID].par = new
        self.samples[new.ID].path['info'] = os.path.join(self.path['Samples_info'], str(new.ID))
        
        self.samples[new.ID].path['scanner_run_info'] = self.path['run_info']
        if not os.path.exists(self.samples[new.ID].path['info']):
            os.makedirs(self.samples[new.ID].path['info'])
        else:
            from shutil import rmtree
            rmtree(self.samples[new.ID].path['info'])
            os.makedirs(self.samples[new.ID].path['info'])
        self.samples[new.ID].init_logger(self.cf['logging']['scanner'])

        # print(self.pars['emptyData'])
        
    def get_new_sample(self):
        if self.pars['points'].loc[self.pars['points']['Status'] == "Free"].shape[0]:
            self.find_free_point_and_make_sample()
            self.record_points()
        
    def check_samples_status_number(self, output=False):
        self.info['nsample'] = {
            "tot":  self.pars['points'].shape[0],
            "free": self.pars['points'].loc[self.pars['points']['Status'] == "Free"].shape[0],
            "redy": self.pars['points'].loc[self.pars['points']['Status'] == "Ready"].shape[0],
            "runn": self.pars['points'].loc[self.pars['points']['Status'] == "Running"].shape[0],
            "fini": self.pars['points'].loc[self.pars['points']['Status'] == "Finish"].shape[0],
            "stop": self.pars['points'].loc[self.pars['points']['Status'] == "Stoped"].shape[0],
            "done": self.pars['points'].loc[self.pars['points']['Status'] == "Done"].shape[0]
        }
        if output:
            self.logger.info("{} samples in total: {} Free, {} Ready, {} Running, {} Finish, {} Done, {} Stopped.".format(
                self.info['nsample']['tot'],
                self.info['nsample']['free'],
                self.info['nsample']['redy'],
                self.info['nsample']['runn'], 
                self.info['nsample']['fini'], 
                self.info['nsample']['done'],
                self.info['nsample']['stop']
            ))

    def check_generator_status(self):
        if self.info['nsample']['tot'] == self.info['nsample']['done'] + self.info['nsample']['stop']:
            self.status = "FINISH"
        elif self.info['nsample']['runn'] > 0:
            self.status = "RUNNING"
        # elif self.info['nsample']['tot'] 
        

class Possion_Disk(Sampling_method):
    def __init__(self) -> None:
        super().__init__()

    def set_config(self, cf):
        return super().set_config(cf)        

    def set_scan_path(self, pth):
        return super().set_scan_path(pth)

    def initialize_generator(self, cf):
        self.cf = cf
        self.init_logger(cf['logging']['scanner'])
        self.logger.info("Possion Disk Sampling Algorithm is used in this scan")
        self.set_pars()
        self.runing_card = os.path.join(self.path['save dir'], 'possion_run.json')
        self.status = "INIT"
        from Func_lib import get_time_lock 
        self.timelock = get_time_lock(self.cf['default setting']['sampling']['TSavePoint'])
        # print(self.path, self.runing_card)


    def set_pars(self):
        pars = {
                "variable": self.scf.get("Sampling_Setting", "variables"),
                "likelihood":   self.scf.get("Sampling_Setting", "likelihood")
        }
        self.pars = {"minR":     float(self.scf.get("Sampling_Setting", "min R"))}
        self.decode_sampling_variable_setting(pars['variable'])
        self.decode_function()
        self.decode_likelihood(pars['likelihood'])
        # print(self.pars.items())
        
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
            elif len(ss) == 4 and ss[1].strip().lower() == "log":
                self.pars['vars'][ss[0].strip()] = {
                    "prior":    "flat",
                    "min":      float(ss[2].strip()),
                    "max":      float(ss[3].strip()),
                    "expr":     sympify("10**({} + ({} - {}) * cube{})".format(float(ss[2].strip()), float(ss[3].strip()), float(ss[2].strip()), nn))
                }
                if not (float(ss[2].strip()) > 0 and float(ss[3].strip())):
                    self.logger.error("Illegal Variable setting of {}\n\t\t-> The lower limit and high limit should be > 0 ".format(ss))
                    sys.exit(0)
            else:
                self.logger.error("Illegal Variable setting in: {} ".format(line))
                sys.exit(0)
            nn += 1
        self.logger.warning(self.pars['vars'])

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
        if "&FC" in self.pars['likelihood']['expression']:
            for item in self.exprs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['expression'] = self.pars['likelihood']['expression'].replace(item, self.exprs[item]['expr'])
            for item in self.fcs:
                if item in self.pars['likelihood']['expression']:
                    self.pars['likelihood']['fcsinc'].append(item)

    def decode_function(self):
        self.fcs = {}
        self.exprs = {}
        for sc in self.scf.sections():
            if "Function" == sc[0:8] and self.scf.get(sc, "method") == "expression":
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                self.exprs[name] = {
                    "expr":   self.scf.get(sc, "expression")
                }
            elif "Function" == sc[0:8] and self.scf.get(sc, "method") == "interpolate 1d":
                from pandas import read_csv
                name = "&FC_{}".format(self.scf.get(sc, "name"))
                from Func_lib import decode_path_from_file
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

    def init_logger(self, cf):
        self.logger = logging.getLogger("Possion")
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
        handler['stream'].setFormatter(Formatter(cf['stream_format'],  "%m/%d %H:%M:%S"))
        handler['ff'].setFormatter(Formatter(cf['file_format']))

    def prerun_generator(self):
        self.get_empty_data_row()
        # print(self.pars['emptyCube'])


    def get_empty_data_row(self):
        raw = {
            "ID":   None,
            "Status":   "Free"
        }
        self.pars['ndim'] = 0
        for var, item in self.pars['vars'].items():
            raw[var] = None
            self.pars['ndim'] += 1  
        self.pars['emptyData'] = Series(raw)
        raw.pop('Status')
        self.pars['emptyCube'] = Series(raw)
