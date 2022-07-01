#!/usr/bin/env python3 

from copy import deepcopy
import json
import logging
import os, sys

from sympy import true 

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class Library():
    def __init__(self) -> None:
        self.cf = None 
        self.logger = None
        self.scan_info = None
        self.rinfo = None 
        
    def set_logger_setting(self, config):
        self.log_config = config
        self.logger = logging.getLogger(config['name'])
        handler = {
            "stream":   logging.StreamHandler(),
            "scanner":  logging.FileHandler(config['scanner_logging_path']),
        }
        self.logger.setLevel(logging.DEBUG)
        handler['stream'].setLevel(logging.WARNING)
        handler['scanner'].setLevel(logging.WARNING)
        self.logger.addHandler(handler['stream'])
        self.logger.addHandler(handler['scanner'])
        from logging import Formatter
        handler['stream'].setFormatter(Formatter(config['stream_format'], "%m/%d %H:%M:%S"))
        handler['scanner'].setFormatter(Formatter(config['scanner_format']))
        self.logger.warning("Library initialize logger successful !!!\n")
        
    def load_setting(self):
        if not self.cf.has_section("Library_Settings"):
            self.logger.warning("Jarvis not found the Library setting in this run")
            self.libs = {
                "modes":    False
            }
        else:
            self.logger.warning("Jarvis loading Library packages")
            self.libs = {
                "modes":    True,
                "config":   {
                    "status":   None,
                    "root required": False
                },
                "include":  {}
            }
        if self.libs['modes']:
            from Func_lib import checking_os
            self.libs['config']['OS required'] = list(map(str.strip, self.cf.get("Library_Settings", "OS required").lower().split(',')))
            if not checking_os() in self.libs['config']['OS required']:
                self.libs['config']['status'] = False
                self.libs['config']['Current OS'] = checking_os()
                self.logger.error("Library Checking: OS requiring {}, -> No\n\tCurrent OS is {}".format(self.libs['config']['OS required'], self.libs['config']['Current OS']))
                sys.exit(0)
            else:
                self.libs['config']['status'] = True
                self.libs['config']['Current OS'] = checking_os() 
                self.logger.warning("Library Checking: OS requiring {}, -> Yes".format(self.libs['config']['OS required']))
            if not self.cf.has_option("Library_Settings", "make paraller"):
                self.libs['config']['paraller number'] = "cpu_count"
            else:
                self.libs['config']['paraller number'] = self.cf.get("Library_Settings", "make paraller")
            from Func_lib import count_cores
            self.libs['config']['paraller number'] = count_cores(self.libs['config']['paraller number'])
            self.cf.set("Library_Settings", "make paraller", str(self.libs['config']['paraller number']))
            if self.cf.has_option("Library_Settings", "root required"):
                if eval(self.cf.get("Library_Settings", "root required")):
                    self.libs['config']['root required'] = True
                    self.logger.info("Library Checking: CERN ROOT checking ...")
                else:
                    self.libs['config']['root required'] = False
                from Func_lib import check_cern_root
                if self.libs['config']['root required'] and self.cf.has_option("Library_Settings", "root path"):
                    self.libs['config']['root installed'], self.libs['config']['root path'] = check_cern_root(self.cf.get("Library_Settings", "root path"))
                    self.cf.set("Library_Settings", "root path", self.libs['config']['root path'])
                    if not self.libs['config']['root installed']: 
                        self.logger.error("\tCERN ROOT found: -> No")
                        sys.exit(0)
                    else:
                        self.logger.warning("CERN ROOT found: -> Yes\n\t=> ROOT path: -> {}".format(self.libs['config']['root path']))
            if self.cf.has_option("Library_Settings", "python version"):
                from Func_lib import check_python2_version
                pytag, self.libs['config']['python'] = check_python2_version(self.cf.get("Library_Settings", "python version"))
                if pytag:
                    self.logger.warning("Library Checking: python version requiring >= {}, -> Yes\n\t=>Current python {} found in your system".format(self.cf.get("Library_Settings", "python version"), self.libs['config']['python']))
                else:
                    self.logger.error("Libarary Checking python version requiring >= {}, -> No \n\t=>Current python {} found in your system".format(self.cf.get("Library_Settings", "python version"), self.libs['config']['python']))
                    sys.exit(0)
            if self.cf.has_option("Library_Settings", "package path"):
                self.libs['home path'] = self.decode_path(self.cf.get("Library_Settings", "package path"))
                self.cf.set("Library_Settings", "package path", self.libs['home path'])
            if self.cf.has_option("Library_Settings", "include"):
                incl = self.cf.get("Library_Settings", "include").split(',')
                self.libs['incl'] = []
                for pkg in incl:
                    self.libs['incl'].append( pkg.strip() )
                for pkg in self.libs['incl']:
                    self.load_library_package_setting(pkg)
                for pkg in self.libs['incl']:
                    self.decode_cmd_to_list(pkg)  

    def set_scan_info(self, sinfo):
        with open(sinfo, 'r') as f1:
            self.scan_info = json.loads(f1.read())          
        
    def decode_path(self, pathdir):
        if "&pwd" in pathdir:
            pathdir = pathdir.replace("&pwd", self.scan_info['path']['pwd'])
        if "&J" in pathdir:
            pathdir = pathdir.replace("&J", self.scan_info['path']['&J'])
        return pathdir

    def load_library_package_setting(self, pkg):
        if self.cf.has_section(pkg):
            rdic= dict(self.cf.items(pkg))
            rdic['required package']    = eval(rdic['required package'])
            rdic['installed']           = eval(rdic['installed'])
            rdic['install path']        = self.decode_path(rdic['install path'])
            self.cf.set(pkg, 'install path', rdic['install path'])
            rdic['source file']         = self.decode_path(rdic['source file'])
            self.cf.set(pkg, 'source file', rdic['source file'])
            self.libs['include'][pkg] = rdic

    def decode_cmd_to_list(self, pkg):
        if self.cf.has_section(pkg):
            cmds = self.libs['include'][pkg]['install cmd'].split('\n')
            res = []
            pwd = self.libs['home path']
            ii = 0
            for cmd in cmds:
                cmd = self.decode_command(cmd, pwd)
                pwd = cmd['path']
                res.append(cmd)
                ii += 1
            self.libs['include'][pkg]['install cmd'] = res

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

    def decode_layers(self):
        from tree import Tree 
        self.tree = Tree()
        for pkg in self.libs['include']:
            self.tree.addNode(pkg, self.libs['include'][pkg])
        self.tree.makeTree()
               
    def build(self):
        self.libs['status'] = "init"
        self.libs['info file'] = os.path.join(self.libs['home path'], 'Library_info.json')
        if os.path.exists(self.libs['info file']):
            self.rinfo = {
                "info file": self.libs['info file']
            } 
            with open(self.libs['info file'], 'r') as f1:
                self.libs['prev'] = json.loads(f1.read())
                self.libs['info'] = {}
            self.libs['status'] = "adding"
        else:
            self.libs['info'] = {}
            self.libs['prev'] = {}      
        # print(self.libs['prev'])
        self.decode_layers()
        if self.libs['status'] == "init":
            for layer in self.tree.layer:
                for pkg in layer:
                    self.libs['info'][pkg] = lib_package()
                    libconfig = deepcopy(self.log_config)
                    libconfig['name'] = pkg 
                    libconfig['ff_logging_path'] = os.path.join(self.libs['home path'], "Lib_{}_install.log".format(pkg))
                    self.libs['include'][pkg]["workerconfig"] = libconfig
                    self.libs['info'][pkg].set_logger_setting(libconfig)
                    self.libs['info'][pkg].set_commands(self.libs['include'][pkg]['install cmd'])
                cktg = True
                while cktg:
                    self.checking_library_package_status()
                    if self.libs['status'] == "ready":
                        cktg = False
        elif self.libs['status'] == "adding":
            self.compare_new_old_library_setting()
            for pkg in self.libs['include']:
                if pkg in self.libs['incl']:
                    self.libs['info'][pkg] = lib_package()
                    libconfig = deepcopy(self.log_config)
                    libconfig['name'] = pkg 
                    libconfig['ff_logging_path'] = os.path.join(self.libs['home path'], "Lib_{}_install.log".format(pkg))
                    self.libs['include'][pkg]["workerconfig"] = libconfig
                    self.libs['info'][pkg].set_logger_setting(libconfig)
                    self.libs['info'][pkg].set_commands(self.libs['include'][pkg]['install cmd'])
                else:
                    self.libs['include'][pkg] = self.libs['prev']['include'][pkg]
            ctag = True
            while ctag:
                self.checking_library_package_status()
                if self.libs['status'] == "ready":
                    ctag = False
        self.logger.warning("Library is ready for the scan task!\n")           

    def compare_new_old_library_setting(self):
        self.logger.warning("Jarvis find the previous Library Setting, Checking the previous Library info ...")
        if self.libs['prev']['config'] == self.libs['config']:
            self.logger.warning("Jarvis find Library Setting is same as the previous one")
            isli = []
            for pkg in self.libs['prev']['include']:
                if self.libs['prev']['include'][pkg]['installed']:
                    isli.append(pkg)
            if isli:
                self.logger.warning("Jarvis find the Packages:\n\t\t{}\n\tare installed".format(" ,".join(isli)))
                from Func_lib import wait_for_input
                ipt = wait_for_input("\tEnter 'n' to select the package should be reinstalled\n\tEnter 'y' to use the previous library:\t", 10, "y")
                if ipt == "y":
                    self.libs['incl'] = self.rm_items_from_list(isli, self.libs['incl'])
            else:
                self.logger.warning("Jarvis not find the installed packages, Jumping into next steps")
        else:
            self.logger.info("Jarvis find previous Library is not same, using init mode to install library package ?")
            from Func_lib import wait_for_input
            ipt = wait_for_input("\tEnter 'n' to select the package should be reinstalled\n\tEnter 'y' to use the previous library:\t", 10, "y")  
            isli = []
            for pkg in self.libs['prev']['include']:
                if self.libs['prev']['include'][pkg]['installed']:
                    isli.append(pkg)
            if ipt == "y":
                self.libs['incl'] = self.rm_items_from_list(isli, self.libs['incl'])
            if ipt == "n":
                self.libs['rein'] = self.ask_list(isli)
                self.libs['incl'] = self.rm_items_from_list(isli, self.libs['incl'])
                self.libs['incl'] += self.libs['rein']                
                self.logger.info("Jarvis will reinstall Library packages: {}".format(", ".join(self.libs['rein'])))
                              
    def checking_library_package_status(self):
        self.libs['status'] = "ready"
        for pkg in self.libs['info'].keys():
            if not self.libs['info'][pkg].status == "finish":
                self.libs['info'][pkg].update_status()
                self.libs['status'] = "runing"
            elif not self.libs['include'][pkg]['installed']:
                self.libs['include'][pkg]['installed'] = True
            else:
                pass 
        
    def output_card(self):
        if self.libs['status'] == "ready":
            self.libs['prev'] = {}
            self.libs['info file']
            with open(self.libs['info file'], 'w') as f1:
                self.libs.pop("info")
                json.dump(self.libs, f1, indent=4)

    def ask_list(self, li):
        atag = "-99"
        sel = []
        while not (atag == "0"):
            sreout = "\n\t========================================\n\tCode\tPackage\t\tStatus\n\t========================================\n"
            for ii in range(len(li)):
                if li[ii] not in sel:
                    sreout += "\t{}\t{}\t\tSKIP\n".format(ii+1, li[ii])
                else:
                    sreout += "\t{}\t{}\t\tSELECTED\n".format(ii+1, li[ii])
            sreout += "\t========================================\n\tEnter the code you selected\n\tEnter 0 to continue:\t "
            atag = input(sreout)
            if atag.isdigit():
                if (int(atag) > 0 and int(atag) <= len(li) and (li[int(atag)-1] not in sel)):
                    sel.append(li[int(atag) -1])
                elif (int(atag) > 0 and int(atag) <= len(li) and (li[int(atag)-1] in sel)):
                    sel.remove(li[int(atag)-1])
            else:
                print("\tillegal input, Please select the code in the table !")
        return sel

    def rm_items_from_list(self, l1, l2):
        res = []
        for it in l2:
            res.append(it.strip())
        for it in l1:
            res.remove(it)
        return res
       

class lib_package():
    def __init__(self) -> None:
        self.cmds = []
        self.status = "free"
        self.subp = None
        self.pid = None
        import logging, logging.config
        self.logger = None

    def set_commands(self, cmds):
        from copy import deepcopy
        self.cmds = deepcopy(cmds) 
    
    def set_logger_setting(self, config):
        self.logger = logging.getLogger(config['name'])
        handler = {
            "stream":   logging.StreamHandler(),
            "scanner":  logging.FileHandler(config['scanner_logging_path']),
            "ff":       logging.FileHandler(config['ff_logging_path'])
        }
        self.logger.setLevel(logging.DEBUG)
        handler['stream'].setLevel(logging.WARNING)
        handler['scanner'].setLevel(logging.WARNING)
        handler['ff'].setLevel(logging.DEBUG)
        self.logger.addHandler(handler['stream'])
        self.logger.addHandler(handler['scanner'])
        self.logger.addHandler(handler['ff'])
        from logging import Formatter
        handler['stream'].setFormatter(Formatter(config['stream_format'], "%m/%d %H:%M:%S"))
        handler['scanner'].setFormatter(Formatter(config['scanner_format']))
        handler['ff'].setFormatter(Formatter(config['file_format']))
        self.logger.warning("Logger initilized .... ")

    def update_status(self):
        if self.subp == None:
            self.status = "free"
            self.run_next_command()
            self.logger.warning("Start installation ...")
        elif self.subp.poll() == None:
            self.status = "runing"
            # print(self.subp, self.subp.poll())
        elif self.subp.poll() is not None and self.cmds:
            # print(self.subp, self.subp.poll())

            output = self.subp.stdout.read().decode()
            if output.strip() != "":
                output = output.split("\n")
                output.remove("")
                output = "\n\t".join(output)
                output += "\n"
                output = "Screen Output -> \n\t" + output
                self.logger.info(output)
            self.run_next_command()
            self.status = "runing"
        else:
            self.status = "finish"
            output = self.subp.stdout.read().decode()
            if output.strip() != "":
                output = output.split("\n")
                output.remove("")
                output = "\n\t".join(output)
                output += "\n"
                output = "Screen Output -> \n\t" + output
                self.logger.info(output)            
            self.logger.warning("Installation Finished!!!")

    def run_next_command(self):
        cmd = self.cmds.pop(0)
        self.logger.info("Package running command -> {}\n\tin path ->{}".format(cmd['cmd'], cmd['path']))
        from subprocess import Popen, PIPE, STDOUT
        self.subp = Popen(cmd['cmd'], shell=True, stdout=PIPE, stderr=STDOUT, cwd=cmd['path'])


