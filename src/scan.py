#!/usr/bin/env python3

import configparser
from copy import deepcopy
import json
import os
from subprocess import Popen, run
import sys
import logging 
import logging.config
import time 
from program import Pack

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.abspath(os.path.join(pwd, "Sampling"))) 

class Scan():
    def __init__(self) -> None:
        self.method = None
        from configparser import ConfigParser
        self.cf = ConfigParser()
        self.status = "init"
        self.smp = {}
        self.paths = {}
        self.op = {}
        self.rerun = None
        self.logger = {} 

    def show_logo(self):
        with open(os.path.join(jpath, "src/card/logo"), 'r') as f1:
            print("\n{}".format(f1.read()))
        

    def init_scan_by_config_file(self, cfile):
        self.paths['config'] = os.path.abspath(cfile)
        self.paths['pwd'] = os.path.dirname(self.paths['config'])
        self.cf.read(self.paths['config'])
        self.basic_check_config()
        self.init_sampling()
        self.load_scanner_logger_setting()
        self.check_sampling_mode()
        if self.mode == "new":
            self.build_library()
            self.build_package()
            self.make_generator()
            self.generator.prerun_generator()
            self.generator.pack = self.pack.pack
            self.generator.generate_events()
            self.generator.afterrun_generator()
        elif self.mode == "resume":
            self.resume_task()
            self.generator.pack = self.pack.pack
            self.generator.generate_events()
            self.generator.afterrun_generator()
            self.generator.combine_data()

            

    def resume_task(self):
        with open(self.smp['output']['libs'], 'r') as f1:
            self.logger['logger'].info("Jarvis reload library setting !")
            self.libs = json.loads(f1.read())
            
        self.pack = Pack()
        self.pack.cf = deepcopy(self.cf)
        with open(self.smp['output']['pkgs'], 'r') as f1:
            self.pack.pack = json.loads(f1.read())
        packconfig = {
                    "name":                 "Package Builder",
                    "scanner_logging_path": self.smp['log']['scanner_logging_path'],
                    "stream_format":        self.smp['sp']['cf']['logging']['scanner']['stream_format'],
                    "file_format":          self.smp['sp']['cf']['logging']['package']['file_format'],
                    "scanner_format":       self.smp['sp']['cf']['logging']['scanner']['file_format']
        }        
        self.pack.set_logger_setting(packconfig)
        self.pack.smp = self.smp['sp']['package']   
        self.pack.info['save dir'] = deepcopy(self.smp['Scan']['save dir'])
        self.pack.info['output'] = deepcopy(self.smp['output'])
        self.pack.resume_by_card()
        self.resume_generator()
        for pkg, item in self.pack.pack['include'].items():
            self.logger['logger'].info("Jarvis reloading Program -> {}\t workers".format(pkg))
            with open(item['run info'], 'r') as f1:
                runinfo = json.loads(f1.read())
            for worker, status in runinfo['workers'].items():
                if status == "running":
                    runinfo['workers'].update({
                        worker: "died"
                    })
            with open(item['run info'], 'w') as f1:
                json.dump(runinfo, f1, indent=4)

            # print(pkg, item['run info'])
            # self.load_package_setting()
            
            

        
    def load_scanner_logger_setting(self):
        self.smp['sp']['cf']['logging']['scanner']['logging config'] = self.decode_path(self.smp['sp']['cf']['logging']['scanner']['logging config'])
        # logging.config.fileConfig(self.smp['sp']['cf']['logging']['scanner']['logging config'])
        self.logger['logger'] = logging.getLogger(self.smp['sp']['cf']['logging']['scanner']['name'])
        self.smp['log'] = {
            "scanner_logging_path": os.path.join(self.smp['Scan']['save dir'], "{}.log".format(self.smp['Scan']['name']))
        }
        self.smp['sp']['cf']['logging']['scanner']['logging path'] = self.smp['log']['scanner_logging_path']
        self.logger['handler'] = {
            "stream":   logging.StreamHandler(),
            "ff":       logging.FileHandler(self.smp['log']['scanner_logging_path'])
        }
        self.logger['logger'].setLevel(logging.DEBUG)
        self.logger['handler']['stream'].setLevel(logging.INFO)
        self.logger['handler']['ff'].setLevel(logging.DEBUG)
        self.logger['logger'].addHandler(self.logger['handler']['stream'])
        self.logger['logger'].addHandler(self.logger['handler']['ff'])
        from logging import Formatter
        self.logger['handler']['stream'].setFormatter(Formatter(self.smp['sp']['cf']['logging']['scanner']['stream_format'], "%m/%d %H:%M:%S"))
        self.logger['handler']['ff'].setFormatter(Formatter(self.smp['sp']['cf']['logging']['scanner']['file_format']))

    def check_sampling_mode(self):            
        if ".Running" not in os.listdir(self.smp['Scan']['save dir']):
            with open(self.decode_path(self.scm['logo path']), 'r') as f1:
                self.logger['logger'].info("\n{}".format(f1.read()))
            self.smp['path'] = self.paths
            self.smp['path']['&J'] = jpath
            self.smp['output'] = {
                "path": self.smp['Scan']['save dir'],
                "temp": os.path.join(self.smp['Scan']['save dir'], ".Running"),
                "info": os.path.join(self.smp['Scan']['save dir'], "Run.json")
            }
            os.mkdir(self.smp['output']['temp'])
            with open(self.smp['output']['info'], 'w') as f1:
                json.dump(self.smp, f1, indent=4)
            self.mode = "new"
            self.logger['logger'].info("Jarvis RUN NEW task: {}".format(self.smp['Scan']['name']))
        elif ".Finish" not in os.listdir(self.smp['Scan']['save dir']):
            with open(os.path.join(self.smp['Scan']['save dir'], "Run.json"), 'r') as f1:
                self.smp = json.loads(f1.read())
            self.mode = "resume"
            if self.cf.has_option("Scan", "rerun tag"):
                rt = self.cf.get("Scan", "rerun tag").split(",")
                rt = [x.strip() for x in rt]
                self.rerun = rt
            self.logger['logger'].info("Jarvis RESUME the task: {}".format(self.smp['Scan']['name']))
        else:
            self.mode = "finish"
            self.logger['logger'].info("Jarvis Finished the task: {}".format(self.smp['Scan']['name']))

    def build_library(self):
        from copy import deepcopy
        from library import Library
        self.libs = Library()
        self.libs.cf = deepcopy(self.cf)
        libsconfig = {
                    "name":                 "Library Builder",
                    "scanner_logging_path": self.smp['log']['scanner_logging_path'],
                    "stream_format":        self.smp['sp']['cf']['logging']['scanner']['stream_format'],
                    "file_format":          self.smp['sp']['cf']['logging']['package']['file_format'],
                    "scanner_format":       self.smp['sp']['cf']['logging']['scanner']['file_format']
        }
        self.libs.set_logger_setting(libsconfig)
        self.libs.set_scan_info(self.smp['output']['info'])
        self.libs.load_setting()
        self.libs.build()
        self.libs.output_card()
        self.smp['sp']['library'] = {
            "info file":    self.libs.libs['info file']
        }
        self.smp['output']['libs'] = self.libs.libs['info file']
        with open(self.smp['output']['info'], 'w') as f1:
                json.dump(self.smp, f1, indent=4)
        
    def build_package(self): 
        from copy import deepcopy
        from program import Pack
        self.pack = Pack()
        self.pack.cf = deepcopy(self.cf)
        packconfig = {
                    "name":                 "Package Builder",
                    "scanner_logging_path": self.smp['log']['scanner_logging_path'],
                    "stream_format":        self.smp['sp']['cf']['logging']['scanner']['stream_format'],
                    "file_format":          self.smp['sp']['cf']['logging']['package']['file_format'],
                    "scanner_format":       self.smp['sp']['cf']['logging']['scanner']['file_format']
        }
        self.pack.set_logger_setting(packconfig)
        self.pack.set_scan_info(self.smp['output']['info'])
        self.pack.load_setting()
        self.smp['sp']['package'] = {}
        self.pack.smp = self.smp['sp']['package']
        self.pack.info['save dir'] = deepcopy(self.smp['Scan']['save dir'])
        self.pack.output_card()
        self.smp['output']['pack'] = os.path.join(self.smp['Scan']['save dir'], "packages_setting.json")        
        self.smp['output']['pkgs'] = self.pack.info['output']['pkgs']
        with open(self.smp['output']['info'], 'w') as f1:
                json.dump(self.smp, f1, indent=4)

    def init_sampling(self):
        self.smp['sp'] = {
            'cfpath': self.decode_path(self.scm['Sampling method'][self.smp['Scan']['sampling method']])
        }
        with open(self.smp['sp']['cfpath'], 'r') as f1:
            self.smp['sp']['cf'] = json.loads(f1.read())
        self.smp['Scan']['save dir'] = os.path.join(
            self.smp['Scan']['save dir'], self.smp['Scan']['name'])
        if not os.path.exists(self.smp['Scan']['save dir']) and not self.smp['Scan']['force new run']:
            os.makedirs(self.smp['Scan']['save dir'])
        elif self.smp['Scan']['force new run']:
            if os.path.exists(self.smp['Scan']['save dir']):
                from shutil import rmtree
                rmtree(self.smp['Scan']['save dir'])
            os.makedirs(self.smp['Scan']['save dir'])

    def make_generator(self):
        if self.smp['Scan']['sampling method'] == "Poisson Disk":
            from poissondisk import Possion_Disk
            self.generator = Possion_Disk()

        elif self.smp['Scan']['sampling method'] == "Grid":
            from sampling import Grid
            self.generator = Grid()
        self.generator.set_config(self.cf)
        self.generator.set_scan_path(self.smp['Scan']['save dir'])
        self.generator.path['run_info'] = self.smp['output']['info']
        self.generator.initialize_generator(self.smp['sp']['cf'])    
        if self.mode == "new":
            self.pack.generator_vars = deepcopy(self.generator.pars['vars'])
            self.pack.pack['tree'].plot_input = deepcopy(self.decode_path(self.scm['Plot_input']))
            self.pack.pack['tree'].BPlot_path = deepcopy(self.scm['BudingPLOT path'])
            self.pack.make_tree()

    def resume_generator(self):
        self.make_generator()     
        self.generator.resume_generator(self.rerun)
        
    def basic_check_config(self):
        from Func_lib import ck_sect
        # ck_sect(self.cf, "Config")
        if not self.cf.has_section("Config"):
            self.paths['scheme'] = os.path.join(
                jpath, "src/card/perference.json")
            print("\tNo scan configure scheme declared, Using default scan scheme in path:\n\t\t->\t{}".format(self.scm))
        elif not self.cf.has_option("Config", "scheme"):
            self.paths['scheme'] = os.path.join(
                jpath, "src/card/perference.json")
            print("\tNo scheme declared in Section -> Config, Using default scheme in path:\n\t\t->\t{}".format(self.scm))
        else:
            self.paths['scheme'] = self.decode_path(
                self.cf.get("Config", "scheme"))
        import json
        with open(self.paths['scheme'], 'r') as f1:
            self.scm = json.loads(f1.read())

        for ss, sv in self.scm['parser_default'].items():
            ck_sect(self.cf, ss)
            self.smp[ss] = {}
            for kk, vv in sv.items():
                if not self.cf.has_option(ss, kk):
                    print(
                        "\tNo Option ->\t{}  in Section [{ss}], Please check your parser file".format(kk, ss))
                    sys.exit(0)
                else:
                    val = self.cf.get(ss, kk)
                    if vv['type'] == "string":
                        pass
                    elif vv['type'] == "path":
                        val = self.decode_path(val)
                    elif vv['type'] == "bool":
                        val = val.lower() == 'true'
                    if vv['smpset']:
                        self.smp[ss][kk] = val

    def decode_path(self, pathdir):
        if "&pwd" in pathdir:
            pathdir = pathdir.replace("&pwd", self.paths['pwd'])
        if "&J" in pathdir:
            pathdir = pathdir.replace("&J", jpath)
        return pathdir

    def decode_command(self, cmd, pwd):
        command = {
            "path": "",
            "cmd":  ""
        }
        import re
        decode = re.compile(r'[$](.*?)[$]', re.S)
        dl = re.findall(decode, cmd)
        for dd in dl:
            dds = dd.split(':')
            cmd = cmd.replace("${}$".format(dd), self.cf.get(dds[0], dds[1]))
        if "cd " == cmd[0:3]:
            command['path'] = cmd[3:]
            command['cmd']  = 'pwd'
        else:
            command['path'] = pwd
            command['cmd']  = self.decode_path(cmd)
        return command

