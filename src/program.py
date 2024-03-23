#!/usr/bin/env python3 

from curses import noecho
import logging
import os, sys 
import json

from tree import Tree

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

class Pack():
    def __init__(self) -> None:
        self.info = {}
        self.logger = None 
        self.scan_info = None 
        self.rinfo = None 
        self.smp = None
        self.generator_vars = None 
        self.mana = None
        
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
        self.logger.warning("Jarvis loading Programs ... !!!\n")
        
    def set_scan_info(self, sinfo):
        with open(sinfo, "r") as f1:
            self.scan_info = json.loads(f1.read())
        
    def load_setting(self):
        if not self.cf.has_section("Program_Settings"):
            self.logger.error("Jarvis not found the Program setting in this run")
            self.pack = {
                "modes":    False
            }
        else:
            self.logger.warning("Jarvis loading Program Setting ... ")
            self.pack = {
                "modes":    True,
                "config":   {
                    "status":   None
                },
                "include":  {},
                "stroutput":    []
            }
        if self.pack['modes']:
            if not self.cf.has_option("Program_Settings", "make paraller"):
                self.pack['config']['paraller number'] = "cpu_count"
            else:
                self.pack['config']['paraller number'] = self.cf.get("Program_Settings", "make paraller")
            from Func_lib import count_cores
            self.pack['config']['paraller number'] = count_cores(self.pack['config']['paraller number'])
            if self.cf.has_option("Program_Settings", "package path"):
                self.pack['home path'] = self.decode_path(self.cf.get("Program_Settings", "package path"))
                self.cf.set("Program_Settings", "package path", self.pack['home path'])
                
            if self.cf.has_option("Program_Settings", "include"):
                self.pack['incl'] = []
                incl = self.cf.get("Program_Settings", "include").split(',')
                for pkg in incl:
                    self.get_pkgs(pkg.strip())

    def get_pkgs(self, pkg_sect):
        if self.cf.has_section(pkg_sect):
            if "include" in self.cf.options(pkg_sect):
                self.pack['include'][pkg_sect] = self.load_multi_package(pkg_sect)
                for pp in self.pack['include'][pkg_sect]['include']:
                    self.pack['incl'].append(pp)
                    self.pack['include'][pkg_sect]['modes'][pp] = self.decode_multi_package_cmd(pp, pkg_sect)
                self.logger.warning("Jarvis loading program {}\n\t-> Mulitple mode: {}".format(pkg_sect, self.pack['include'][pkg_sect]['include']))
            else:
                self.logger.warning("Jarvis loading program {}\n\t-> Single mode".format(pkg_sect))
                self.pack['incl'].append(pkg_sect)
                self.pack['include'][pkg_sect] = self.load_package(pkg_sect)
                self.decode_package_cmd(pkg_sect)

    def decode_path(self, pathdir):
        if "&pwd" in pathdir:
            pathdir = pathdir.replace("&pwd", self.scan_info['path']['pwd'])
        if "&J" in pathdir:
            pathdir = pathdir.replace("&J", self.scan_info['path']['&J'])
        return pathdir

    def load_multi_package(self, pkg):
        if self.cf.has_section(pkg):
            rdic = dict(self.cf.items(pkg))
            rdic['required package'] = eval(rdic['required package'])
            rdic['clone shadow']     = eval(rdic['clone shadow'])
            rdic['source file']      = self.decode_path(rdic['source file'])
            rdic['workers']          = {}
            rdic['multi']            = True
            rdic['modes']            = {}
            rdic['include']          = list(map(str.strip, rdic["include"].split(",")))
            self.cf.set(pkg, "source file", rdic['source file'])
            rdic['install path']     = self.decode_path(rdic['install path'])
            if rdic['clone shadow']:
                rdic['install path'] += "/@PackID"  
            self.cf.set(pkg, 'install path', rdic['install path'])
            rdic['run info'] = self.decode_path(self.decode_command(rdic['run info'], rdic['install path'])['cmd'])
            rdic['install cmd'] = self.decode_multi_cmds(rdic['install cmd'], rdic['install path'])
            return rdic

    def load_package(self, pkg):
        if self.cf.has_section(pkg):
            rdic = dict(self.cf.items(pkg))
            rdic['required package'] = eval(rdic['required package'])
            rdic['clone shadow']     = eval(rdic['clone shadow'])
            rdic['multi']            = False
            rdic['source file']      = self.decode_path(rdic['source file'])
            rdic['workers']          = {}
            self.cf.set(pkg, "source file", rdic['source file'])
            rdic['install path']     = self.decode_path(rdic['install path'])
            if rdic['clone shadow'] and "/@PackID" != rdic['install path'][-8:]:
                rdic['install path'] += "/@PackID"  
            self.cf.set(pkg, 'install path', rdic['install path'])
            return rdic
  
    def resume_by_card(self):
        self.pack['status'] = "ready"
        self.logger.warning("Jarvis reloading program setting ... ")
        self.pack['tree'] = Tree()
        self.resume_tree()
        
    def resume_tree(self):
        self.logger.warning("Jarvis resume Program Trees ...")
        for pkg in self.pack['include']:
            self.pack['tree'].addNode(pkg, self.pack['include'][pkg])
        self.pack['tree'].makeTree()
        self.pack['tree'].pars = self.generator_vars

    def decode_multi_package_cmd(self, pkg, pkg_sect):
        if self.cf.has_section(pkg):
            rdic = dict(self.cf.items(pkg))
            rdic['command path'] = self.decode_path(self.decode_command(rdic['command path'], self.pack['include'][pkg_sect]['install path'])['cmd'])
            rdic['prerun command'] = self.decode_multi_cmds(rdic['prerun command'], rdic['command path'])
            rdic['excute command'] = self.decode_multi_cmds(rdic['excute command'], rdic['command path'])
            rdic['pkg_name'] = pkg
            rdic.update(self.decode_inputs_via_rdic(rdic))
            rdic.update(self.decode_outputs_via_rdic(rdic))
            rdic.pop('ips')
            rdic.pop('ops')
            return rdic

    def decode_package_cmd(self, pkg):
        if self.cf.has_section(pkg):
            self.pack['include'][pkg]['command path'] = self.decode_path(self.decode_command(self.pack['include'][pkg]['command path'], self.pack['include'][pkg]['install path'])['cmd'])
            self.pack['include'][pkg]['run info'] = self.decode_path(self.decode_command(self.pack['include'][pkg]['run info'], self.pack['include'][pkg]['install path'])['cmd'])

            self.decode_cmds(pkg, "install cmd", "install path")
            self.decode_cmds(pkg, "prerun command", "command path")
            self.decode_cmds(pkg, "excute command", "command path")
            
            self.decode_inputs(pkg)
            self.decode_outputs(pkg)

            self.pack['include'][pkg].pop("ips")
            self.pack['include'][pkg].pop("ops")

            self.pack['status'] = "ready"

    def output_card(self):
        if self.pack['status'] == 'ready':
            for pkg in self.pack['include']:
                self.smp[pkg] = self.pack['include'][pkg]['run info']
            self.info['output'] = {}
            self.info['output']['pack'] = os.path.join(self.info['save dir'], "packages_setting.json")
            for pkg in self.pack['include']:
                infodir = os.path.dirname(self.pack['include'][pkg]['run info'])
                if not os.path.exists(infodir):
                    os.makedirs(infodir)
                # self.pack['include'][pkg].pop("ips")
                # self.pack['include'][pkg].pop("ops")
                with open(self.pack['include'][pkg]['run info'], 'w') as f1:
                    json.dump(self.pack['include'][pkg], f1, indent=4)
            self.info['output']['pkgs'] = os.path.join(self.info['save dir'], "Program_info.json")
            with open(self.info['output']['pkgs'], 'w') as f1:
                json.dump(self.pack, f1, indent=4)
        from tree import Tree
        self.pack['tree'] =  Tree()
            
    def make_tree(self):
        self.logger.warning("Jarvis make Program Trees ...")
        for pkg in self.pack['include']:
            self.pack['tree'].addNode(pkg, self.pack['include'][pkg])
        self.pack['tree'].makeTree()
        self.pack['tree'].pars = self.generator_vars
        self.pack['tree'].makeDrawInfo()
        with open(self.info['output']['pack'], 'w') as f1:
            json.dump(self.pack['tree'].info, f1, indent=4)
        with open(self.pack['tree'].plot_input, 'r') as f1:
            inp = f1.read()
            inp = inp.replace(">>>SAVE DIR<<<", self.info['save dir'])
            inp = inp.replace(">>>INFOPATH<<<", self.info['output']['pack'])
        with open(os.path.join(self.info['save dir'], "Plot_Jarvis.ini"), 'w') as f1:
            f1.write(inp)
        # from subprocess import Popen, PIPE, STDOUT
        # subp = Popen("{} {}".format(self.pack['tree'].BPlot_path, "Plot_Jarvis.ini" ), shell=True, stdout=PIPE, stderr=STDOUT, cwd=self.info['save dir'])
        from plot import Plot 
        fig = Plot()
        fig.read_config(os.path.join(self.info['save dir'], "Plot_Jarvis.ini"))
        fig.plot()
                         
    def decode_outputs(self, pkg):                    
        if not self.pack['include'][pkg]['output variables'].strip() == '':
            self.pack['include'][pkg]['output file'] = self.decode_path(self.decode_command(self.pack['include'][pkg]['output file'], self.pack['include'][pkg]['install path'])['cmd'])
            from IOs.IOs import IOs
            self.pack['include'][pkg]['ops'] = IOs()
            self.pack['include'][pkg]['ops'].setinfos(
                finf=self.pack['include'][pkg]['output file'],
                vinf=self.pack['include'][pkg]['output variables']              
            )
            self.pack['include'][pkg]['ops'].logger = self.logger
            self.pack['include'][pkg]['ops'].pkg = pkg 
            for ipf in self.pack['include'][pkg]['ops'].files.keys():
                self.pack['include'][pkg]['ops'].files[ipf]['path'] = self.decode_path( self.decode_command( self.pack['include'][pkg]['ops'].files[ipf]['path'], self.pack['include'][pkg]['command path'] )['cmd'] )

            self.pack['include'][pkg]['ops'].type = "output"
            self.pack['include'][pkg]['ops'].build()
            self.pack['include'][pkg]['output variables'] = self.pack['include'][pkg]['ops'].vars
            self.pack['include'][pkg]['output file'] = self.pack['include'][pkg]['ops'].finfs
            self.pack['include'][pkg]['output'] = self.pack['include'][pkg]['ops'].files
        else:
            self.logger.error('\033[91m Illegal output variable in program "{}" for Jarvis, Please check your configure file\033[0m'.format(pkg))
            sys.exit(1)
 
    def decode_outputs_via_rdic(self, rdic):
        if not rdic['output variables'].strip() == "":
            rdic['output file'] = self.decode_path(self.decode_command(rdic['output file'], rdic['install cmd'])['cmd'])
            from IOs.IOs import IOs
            rdic['ops'] = IOs()
            rdic['ops'].setinfos(
                finf = rdic['output file'],
                vinf = rdic['output variables']
            )
            rdic['ops'].logger = self.logger
            rdic['ops'].pkg = rdic['pkg_name']
            for ipf in rdic['ops'].files.keys():
                rdic['ops'].files[ipf]['path'] = self.decode_path(
                    self.decode_command(rdic['ops'].files[ipf]['path'], rdic['command path'])['cmd']
                )
            
            rdic['ops'].type = "output"
            rdic['ops'].build()
            rdic['output variables'] = rdic['ops'].vars
            rdic['output file'] = rdic['ops'].finfs 
            rdic['output'] = rdic['ops'].files
            return rdic 
        else:
            self.logger.error('\033[91m Illegal output variable in program "{}" for Jarvis, Please check your configure file\033[0m'.format(rdic['pkg_name']))
            sys.exit(1)

    def decode_inputs(self, pkg):
        if not self.pack['include'][pkg]['input variables'].strip() == '':
            from IOs.IOs import IOs
            self.pack['include'][pkg]['input file'] = self.decode_path(self.decode_command(self.pack['include'][pkg]['input file'], self.pack['include'][pkg]['install path'])['cmd'])
            self.pack['include'][pkg]['ips'] = IOs()
            self.pack['include'][pkg]['ips'].setinfos(
                finf=self.pack['include'][pkg]['input file'],
                vinf=self.pack['include'][pkg]['input variables']  
            )
            self.pack['include'][pkg]['ips'].logger = self.logger
            self.pack['include'][pkg]['ips'].pkg = pkg 
            for ipf in self.pack['include'][pkg]['ips'].files.keys():
                self.pack['include'][pkg]['ips'].files[ipf]['path'] = self.decode_path( self.decode_command( self.pack['include'][pkg]['ips'].files[ipf]['path'], self.pack['include'][pkg]['command path'] )['cmd'] )
            
            self.pack['include'][pkg]['ips'].type = "input"
            self.pack['include'][pkg]['ips'].build()
            
            self.pack['include'][pkg]['input variables'] = self.pack['include'][pkg]['ips'].vars
            self.pack['include'][pkg]['input file'] = self.pack['include'][pkg]['ips'].finfs
            self.pack['include'][pkg]['input'] = self.pack['include'][pkg]['ips'].files
        else:
            self.logger.error('\033[91m Illegal input variable in program "{}" for Jarvis, Please check your configure file\033[0m'.format(pkg))
            sys.exit(1)

    def decode_inputs_via_rdic(self, rdic):
        if not rdic['input variables'].strip() == "":
            from IOs.IOs import IOs 
            rdic['input file'] = self.decode_path(self.decode_command(rdic['input file'], rdic['install cmd'])['cmd'])
            rdic['ips'] = IOs()
            rdic['ips'].setinfos(
                finf=rdic['input file'],
                vinf=rdic['input variables']
            )
            rdic['ips'].logger = self.logger 
            rdic['ips'].pkg = rdic['pkg_name'] 
            for ipf in rdic['ips'].files.keys():
                rdic['ips'].files[ipf]['path'] = self.decode_path(self.decode_command(rdic['ips'].files[ipf]['path'], rdic['command path'])['cmd'])

            rdic['ips'].type = "input"
            rdic['ips'].build()
            rdic['input variables'] = rdic['ips'].vars
            rdic['input file'] = rdic['ips'].finfs 
            rdic['input'] = rdic['ips'].files
            return rdic
        else:
            self.logger.error('\033[91m Illegal input variable in program "{}" for Jarvis, Please check your configure file\033[0m'.format(rdic['pkg_name']))
            sys.exit(1)

    def decode_multi_cmds(self, cmds, pwd):
        res = []
        ii = 0
        commands = cmds.split('\n')
        for cc in commands:
            cmd = self.decode_command(cc, pwd)
            pwd = cmd['path']
            res.append(cmd)
        return res 

    def decode_cmds(self, pkg, cmds, paths):
        commands = self.pack['include'][pkg][cmds].split('\n')
        res  = []
        ii   = 0
        pwd  = self.pack['include'][pkg][paths]
        for cmd in commands:
            cmd = self.decode_command(cmd, pwd)
            pwd = cmd['path']
            res.append(cmd)
            ii += 1
        self.pack['include'][pkg][cmds] = res 
 
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
            command['path'] = self.decode_path(cmd[3:])
            command['cmd']  = 'pwd'
        else:
            command['path'] = pwd
            command['cmd']  = self.decode_path(cmd)
        return command 
        
class program():
    def __init__(self) -> None:
        self.cmds = []
        self.suid = None
        self.status = "free"
        self.subp = None
        self.pid = None
        self.path = None
        import logging, logging.config
        self.logger = None
        self.config = None
        self.vars = None
        self.inp   = {}
        self.oup   = {}
        self.logpath = None
        self.mana = None

        
    def set_logger_setting(self, config):
        self.logger = logging.getLogger(config['name'])
        handler = {
            "stream":   logging.StreamHandler(),
            "scanner":  logging.FileHandler(config['scanner_logging_path']),
            "ff":       logging.FileHandler(config['ff_logging_path'])
        }
        self.logpath = open(config['ff_logging_path'], 'a')
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
        self.logger.warning("Initialize logger successful !!!\n")
        
    def decode_cmd(self, cmd):
        if self.config['clone shadow']:
            cmd['cmd'] = cmd['cmd'].replace("@PackID", self.packid)
            cmd['path'] = cmd['path'].replace("@PackID", self.packid)
        if not os.path.exists(cmd['path']):
            os.makedirs(cmd['path'])
        return cmd 
    
    def decode_path(self, path):
        if self.config['clone shadow']:
            path = path.replace("@PackID", self.packid)
        return path 
    
    def run_next_command(self):
        if self.status == "installing":
            cmd = self.decode_cmd(self.config['install cmd'].pop(0))
            self.logger.info("Program running install command -> {}\n\t in path -> {}".format(cmd['cmd'], cmd['path']))
        elif self.status == "prerun":
            if self.config['multi']:
                cmd = self.decode_cmd(self.config['modes'][self.config['name']]['prerun command'].pop(0))
            else:
                cmd = self.decode_cmd(self.config['prerun command'].pop(0))
            self.logger.info("Program running prerun command -> {}\n\t in path -> {}".format(cmd['cmd'], cmd['path']))
        elif self.status == "running":
            if self.config['multi']:
                cmd = self.decode_cmd(self.config['modes'][self.config['name']]['excute command'].pop(0))
            else:
                cmd = self.decode_cmd(self.config['excute command'].pop(0))                
            self.logger.info("Program excute command -> {}\n\t in path -> {}".format(cmd['cmd'], cmd['path']))
        from subprocess import Popen, PIPE, STDOUT 
        # self.subp = Popen(cmd['cmd'], shell=True, stdout=PIPE, stderr=STDOUT, cwd=cmd['path']) 
        self.subp = Popen(cmd['cmd'], shell=True, stdout=self.logpath, stderr=self.logpath, cwd=cmd['path']) 
        self.pid = self.subp.pid
                
    def clear_space(self):
        if self.config['clone shadow']:
            from shutil import rmtree
            rmtree(self.config['install path'].replace("@PackID", self.packid))
            os.makedirs(self.config['install path'].replace("@PackID", self.packid))
        else:
            from shutil import rmtree
            rmtree(self.config['install path'])
            os.makedirs(self.config['install path'])

    def get_package(self):
        with open(self.config['run info'], 'r') as f1:
            self.run_info = json.loads(f1.read())
        if self.config['clone shadow']:
            if self.run_info['workers']:
                for iid, sst in self.run_info['workers'].items():
                    if sst == "died":
                        self.packid = iid
                        self.clear_space()
                        self.status = "installing"
                        break
                    if sst == 'free':
                        self.packid = iid 
                        self.status = "prerun"
                        break
                if self.packid is None and len(self.run_info['workers'].keys()) < self.config['paraller number']:
                    self.packid = "{:03d}".format(len(self.run_info['workers'].keys()) + 1)
                    self.status = "installing"
                    self.run_next_command()
                elif self.packid is None and len(self.run_info['workers'].keys()) >= self.config['paraller number']:
                    self.status = "waiting"
                elif self.packid is not None:
                    self.run_next_command()
            else:
                self.packid = "{:03d}".format(len(self.run_info['workers'].keys()) + 1)
                self.status = 'installing'
                self.run_next_command()
            if self.status != "waiting":
                self.config['PackID'] = self.packid
                self.update_runweb_status("running")
        else:
            if self.run_info['workers'][self.packid] == "died":
                self.clear_space()
                self.status = 'installing'
                self.update_runweb_status("running")
                self.run_next_command()
            elif self.run_info['workers'][self.packid] == "free":
                self.status = "prerun"
                self.run_next_command()
            elif self.run_info['workers'][self.packid] == "running":
                self.status = "waiting"
    
    def get_package_pool(self):
        # print(self.config.keys())
        print("program 480:", self.mana.pack[self.config['name']]['workers'], "<<<<")
        if self.config['clone shadow']:
            print("program 483", self.mana.pack[self.config['name']]['workers'])
            # print(self.packid, self.config)
            if self.mana.pack[self.config['name']]['workers']:       
                if self.packid is None:         
                    for worker, stt in self.mana.pack[self.config['name']]['workers'].items():
                        if stt == "died":
                            self.packid = worker 
                            self.clear_space()
                            self.status = 'installing'
                            break
                        if stt == 'free':
                            self.packid = worker
                            self.status = "prerun"
                            break 
                    if self.packid is None and len(self.mana.pack[self.config['name']]['workers']) < self.config['paraller number']:
                        self.packid = "{:03d}".format(len(self.mana.pack[self.config['name']]['workers'].keys()) + 1)
                        self.status = "installing"
                        self.run_next_command()
                    elif self.packid is None and len(self.mana.pack[self.config['name']]['workers']) >= self.config['paraller number']:
                        self.status = "waiting"
                    elif self.packid is not None:
                        self.run_next_command()
            else:
                self.packid = "{:03d}".format(len(self.mana.pack[self.config['name']]['workers'].keys()) + 1)
                self.status = "installing"
                self.run_next_command()
            # print("program 494: packid -> {}, \t calculating sample {}".format(self.packid, self.suid))
            # print("before Getting ID", self.mana.pack[self.config['name']]['workers'])
            # self.mana.pack[self.config['name']]['workers'][self.packid] = self.suid
            if self.status != "waiting":
                self.config['PackID'] = self.packid
                self.update_mana_status()
            # print("after getting ID", self.mana.pack[self.config['name']]['workers'], self.status, self.packid)
            # import time
            # time.sleep(3)
            
        # else:


        # import time 

        # time.sleep(10)

    def check_worker_status(self, pkgname):
        for worker, status in self.mana.pack[pkgname]['workers'].items():
            if status == "free":
                return worker
        return None


                        
    def init(self):
        if self.config['clone shadow']:
            self.packid = None
        else:
            self.packid = "Naruto"
        # print(self.mana, self.config['name'])
        # import time
        # time.sleep(10)
        if self.mana is None:
            self.get_package()
        else:
            self.get_package_pool()
        
    def update_status(self):
        if self.status == "installing":
            self.update_installing()
        elif self.status == "prerun":
            self.update_prerun()
        elif self.status == "running":
            self.update_running()
        elif self.status == "finish":
            self.status = "done"
            self.read_output()
            self.logpath.close()
            self.update_runweb_status('free')
        elif self.status == "error":
            self.stop_calculation()
            self.logpath.close()
            self.status = "done"
        elif self.status == "waiting":
            self.get_package()

    def update_installing(self):
        if self.subp == None:
            self.logger.warning("Program start installing ...")
            self.run_next_command()
        elif self.subp.poll() == None:
            pass 
        elif self.subp.poll() is not None and self.config['install cmd']:
            # from Func_lib import decode_stdout 
            # output = self.subp.stdout.read().decode()
            # if output.strip() != "":
            #     output = decode_stdout(output)
            #     self.logger.info(output)
            self.subp = None
            self.run_next_command()
        else:
            # from Func_lib import decode_stdout
            # output = self.subp.stdout.read().decode()
            # if output.strip() != "":
            #     output = decode_stdout(output)
            #     self.logger.info(output)
            self.logger.warning("Installation Finished !!!")
            self.status = "prerun"  

    def update_running(self):
        if self.config['multi']:
            if self.subp == None:
                self.logger.warning("Program start running ...")
                self.run_next_command()
            elif self.subp.poll() == None:
                pass 
            elif self.subp.poll() is not None and self.config['modes'][self.config['name']]['excute command']:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.run_next_command()            
            else:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.status = "finish"        
        else:
            if self.subp == None:
                self.logger.warning("Program start running ...")
                self.run_next_command()
            elif self.subp.poll() == None:
                pass 
            elif self.subp.poll() is not None and self.config['excute command']:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.run_next_command()
            else:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.status = "finish"

    def update_prerun(self):
        if self.config['multi']:
            if self.subp == None:
                self.logger.warning("Program prepare to run ...")
                self.run_next_command()
            elif self.subp.poll() == None:
                pass 
            elif self.subp.poll() is not None and self.config['modes'][self.config['name']]['prerun command']:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.subp = None
                self.run_next_command()    
            else:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.subp = None
                self.prepare_input()
                self.status = "running"        
        else:
            if self.subp == None:
                self.logger.warning("Program prepare to run ...")
                self.run_next_command()
            elif self.subp.poll() == None:
                pass 
            elif self.subp.poll() is not None and self.config['prerun command']:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.subp = None
                self.run_next_command()
            else:
                # from Func_lib import decode_stdout
                # output = self.subp.stdout.read().decode()
                # if output.strip() != "":
                #     output = decode_stdout(output)
                #     self.logger.info(output)
                self.subp = None
                self.prepare_input()
                self.status = "running"
               
    def prepare_input(self):
        if self.config['multi']:
            ips = self.config['modes'][self.config['name']]['input'].items()
        else:
            ips = self.config['input'].items() 
        from IOs.IOs import InputsFile
        for kk, ff in ips:
            ff['path'] = self.decode_path(ff['path'])
            if os.path.exists(ff['path']):
                ff['file'] = InputsFile()
                ff['file'].file = ff['path']
                ff['file'].para = ff['vars']
                ff['file'].pkg = self.config['name']
                ff['file'].vars = self.vars 
                ff['file'].samppath = self.path['info']
                ff['file'].set_variables()
            else:
                self.logger.error("Input file {}, {} not found! Calculation will continue, but please check your output!".format(kk, ff['path']))

    def read_output(self):   
        if self.config['multi']:
            ops = self.config['modes'][self.config['name']]['output'].items()
        else:
            ops = self.config['output'].items()         
        from IOs.IOs import OutputsFile
        for kk, ff in ops:
            ff['path'] = self.decode_path(ff['path'])
            if os.path.exists(ff['path']):
                ff['file'] = OutputsFile()
                ff['file'].pkg = self.config['name']
                ff['file'].file = ff['path']
                ff['file'].para = ff['vars']
                ff['file'].samppath = self.path['info']
                ff['file'].logger = self.logger
                ff['file'].get_variables()
                self.vars.update(ff['file'].vars)
                if ff['file'].error:
                    self.status = "error"
            else:
                self.logger.error("Output file {}, {} can not be found.".format(kk, ff['path']))
                self.status = "error"
        from pandas import Series
        self.logger.warning("Calculation is done!")
        self.logger.info("The output variable is summarized as follow\n\n{}\n".format(Series(self.vars)))
                    
    def setexprvalue(self, expr):
        from sympy import sympify
        expr = sympify(expr)
        return(expr.subs(self.vars))            
                
    def update_runweb_status(self, status):
        with open(self.config['run info'], 'r') as f1:
            self.run_info = json.loads(f1.read())
        self.run_info['workers'][self.packid] = status
        with open(self.config['run info'], 'w') as f1:
            json.dump(self.run_info, f1, indent=4)
             
    def update_mana_status(self):
        self.mana.pack[self.config['name']]['workers'][self.packid] =  self.status


    def stop_calculation(self):
        self.logger.warning("Sample is force stopped by error calculation")
        if "force quit" not in self.config.keys():
            self.config['force quit'] = True
        if self.subp == None:
            self.status = 'done'
            self.update_runweb_status("free")
        elif self.subp.poll() is not None:
            self.status = "done"
            self.update_runweb_status("free")
        elif self.subp.poll() == None and self.config['force quit']:
            from Func_lib import force_kill_process
            force_kill_process(self.pid)
            self.status = 'done'
            self.update_runweb_status("free")
        else:
            from Func_lib import force_kill_process
            force_kill_process(self.pid)
            self.status = "done"
            self.update_runweb_status("died")
            
            

        