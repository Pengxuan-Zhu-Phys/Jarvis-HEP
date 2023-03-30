#!/usr/bin/env python3

from copy import deepcopy
import imp
import logging
import os
from plistlib import FMT_XML
import sys
from re import L
import json
from matplotlib.pyplot import axis
import numpy
import xslha
import pyslha
import xmltodict
import pandas as pd 

jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


class IOs():
    def __init__(self) -> None:
        self.files = {}
        self.finfs = {}
        self.vinfs = None
        self.type = None
        self.logger = None
        self.pkg = None
        self.vars = []

    def setinfos(self, finf, vinf):
        self.set_files(finf)
        self.vinfs = vinf.split("\n")

    def set_files(self, finf):
        fs = finf.split("\n")
        for fline in fs:
            fline = fline.split(",")
            self.files[fline[0].strip()] = {
                "path": fline[1].strip()
            }

    def build(self):
        if self.type == "input":
            for kk, file in self.files.items():
                file['file'] = InputsFile()
                file['file'].set(file['path'])
                file['file'].logger = self.logger
                file['file'].pkg = self.pkg
            for vinf in self.vinfs:
                vinf = list(map(str.strip, vinf.split(',')))
                if vinf[1] not in self.files.keys():
                    self.logger.error('\033[91mVariable "{}" in program "{}": input file "{}" not found\033[0m'.format(
                        vinf[0], self.pkg, vinf[1]))
                    sys.exit(1)
                else:
                    self.files[vinf[1]]['file'].addvars(vinf)
        elif self.type == "output":
            for kk, file in self.files.items():
                file['file'] = OutputsFile()
                file['file'].set(file['path'])
                file['file'].logger = self.logger
                file['file'].pkg = self.pkg
            for vinf in self.vinfs:
                vinf = list(map(str.strip, vinf.split(',')))
                if vinf[1] not in self.files.keys():
                    self.logger.error('\033[91mVariable "{}" in program "{}": output file "{}" not found\033[0m'.format(
                        vinf[0], self.pkg, vinf[1]))
                    sys.exit(1)
                else:
                    self.files[vinf[1]]['file'].addvars(vinf)
        self.collect_vars()

    def collect_vars(self):
        for ff, file in self.files.items():
            file['vars'] = file['file'].para
            self.vars += file['file'].vars
            self.finfs[ff] = file['path']
            file.pop("file")


class IOfile():
    def __init__(self) -> None:
        self.file = None
        self.vars = None
        self.type = None
        self.logger = None
        self.pkg = None
        self.para = {}
        self.samppath = None

    def set(self, file):
        self.file = file
        self.vars = []

    def addvars(self):
        pass


class InputsFile(IOfile):
    def __init__(self) -> None:
        super().__init__()

    def addvars(self, vinf):
        if vinf[2].lower() == "position":
            self.set_inp_position(vinf)
            self.set_par_position(vinf)
        elif vinf[2].lower() == "json":
            self.set_inp_json(vinf)
            self.set_par_json(vinf)
        elif vinf[2].lower() == "file":
            self.set_inp_file(vinf)
            self.set_par_file(vinf)
        elif vinf[2].lower() == "replace":
            self.set_inp_replace(vinf)
            self.set_par_replace(vinf)

    def set_inp_position(self, vinf):
        if len(vinf) != 5:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Position" method, 5 items ( Name ID, Input file ID, Method, Line number, Column number ) need to be provived.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        elif not (vinf[3].isdigit() and vinf[4].isdigit()):
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Position" method, Line number must be interger number\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "position",
                "code": [int(vinf[3]), int(vinf[4])]
            }
            self.vars.append(res)

    def set_par_position(self, vinf):
        if "position" in self.para:
            self.para['position'].append(
                {
                    "expr": vinf[0],
                    "code": [int(vinf[3]), int(vinf[4])]
                }
            )
        else:
            self.para['position'] = [{
                "expr": vinf[0],
                "code": [int(vinf[3]), int(vinf[4])]
            }]

    def set_position(self):
        from sympy import sympify
        from numpy import loadtxt, savetxt
        data = loadtxt(self.file, ndmin=2)
        for var in self.para['position']:
            expr = sympify(var['expr'], locals=self.vars)
            var['value'] = expr.subs(self.vars)
            data[var['code'][0]-1][var['code'][1]-1] = var['value']
        savetxt(self.file, data)

    def set_inp_json(self, vinf):
        if len(vinf) < 3:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Json" method, 3 items ( Name ID, Input file ID, Method ) need to be provided.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "json"
            }
            self.vars.append(res)

    def set_par_json(self, vinf):
        expr_id = 0
        if len(vinf) == 4:
            expr_id = 3
        if "json" in self.para:
            self.para['json'].append({
                "key":  vinf[expr_id],
                "expr": vinf[0],
            })
        else:
            self.para['json'] = [{
                "key":  vinf[expr_id],
                "expr": vinf[0],
            }]

    def set_json(self):
        from sympy import sympify
        vars = {}
        js = None
        for var in self.para['json']:
            expr = sympify(var['expr'], locals=self.vars)
            var['value'] = expr.subs(self.vars)
            vars[var['key']] = float(var['value'])
        with open(self.file, 'r') as f1:
            js = json.loads(f1.read())
        js.update(vars)
        with open(self.file, 'w') as f1:
            json.dump(js, f1, indent=4, cls=NpEncoder)

    def set_inp_file(self, vinf):
        if not (len(vinf) == 3 or len(vinf) == 4):
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "File" method, 3 or 4 items ( Name ID, Input file, ID, Method, "save" ) need to be provided.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "file",
                "save": False
            }
            if len(vinf) == 4 and vinf[3].lower() == "save":
                res['save'] = True
            self.vars.append(res)

    def set_par_file(self, vinf):
        par = {
            "expr": vinf[0],
            "save": len(vinf) == 4 and "save" in vinf
        }
        self.para['file'] = par

    def set_file(self):
        if self.para['file']['expr'] in self.vars:
            from sympy import sympify
            expr = sympify(self.para['file']['expr'], locals=self.vars)
            self.para['file']['value'] = expr.subs(self.vars)
            from shutil import move
            move(self.para['file']['value'], self.file)
        if self.para['file']['save']:
            from shutil import copyfile
            copyfile(self.file, os.path.join(
                self.samppath, os.path.basename(self.file)))

    def set_inp_replace(self, vinf):
        if len(vinf) != 4:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Replace" method, 4 items ( Name ID, Input file ID, Method, Match ) need to be provived.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "replace",
                "code": vinf[3]
            }
            self.vars.append(res)

    def set_par_replace(self, vinf):
        if "replace" in self.para:
            self.para['replace'].append({
                "expr": vinf[0],
                "code": vinf[3]
            })
        else:
            self.para['replace'] = [{
                "expr": vinf[0],
                "code": vinf[3]
            }]

    def set_replace(self):
        content = None
        from sympy import sympify
        with open(self.file, "r") as f1:
            content = f1.read()
        for var in self.para['replace']:
            expr = sympify(var['expr'], locals=self.vars)
            try:
                var['value'] = expr.subs(self.vars)
            except:
                var['value'] = expr
            content = content.replace(var['code'], str(var['value']))
        with open(self.file, 'w') as f1:
            f1.write(content)

    def set_variables(self):
        if "position" in self.para:
            self.set_position()
        if "json" in self.para:
            self.set_json()
        if "file" in self.para:
            self.set_file()
        if "replace" in self.para:
            self.set_replace()


class OutputsFile(IOfile):
    def __init__(self) -> None:
        super().__init__()
        self.error = False
        self.etype = None
        self.vars = {}

    def addvars(self, vinf):
        if vinf[2].lower() == "position":
            self.set_oup_position(vinf)
            self.set_par_position(vinf)
        elif vinf[2].lower() == "json":
            self.set_oup_json(vinf)
            self.set_par_json(vinf)
        elif vinf[2].lower() == "file":
            self.set_oup_file(vinf)
            self.set_par_file(vinf)
        elif vinf[2].lower() == "slha":
            self.set_oup_slha(vinf)
            self.set_par_slha(vinf)
        elif vinf[2].lower() == "xml":
            self.set_oup_xml(vinf)
            self.set_par_xml(vinf)
        elif vinf[2].lower() == "yoda":
            self.set_oup_yoda(vinf)
            self.set_par_yoda(vinf)

    def set_oup_yoda(self, vinf):
        if len(vinf) < 5:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Yoda" method, 5 items ( Name ID, Input file ID, Method, Hist Type, Hist Name ) need to be provived.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        elif not vinf[3].lower() in ['hist1d', "hist2d"]:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Yoda" method, Hist Type -> {} is not support!\033[0m'.format(vinf[0], self.pkg, vinf[3]))
            sys.exit(1)
        if vinf[3].lower() == "hist1d" and len(vinf) != 7: 
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Yoda" method, Hist1D Type needs 7 item ->\n\t\t( Name ID, Input file ID, Method, Hist type, Hist Name, Xlow, XHigh )!\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        elif vinf[3].lower() == "hist2d" and len(vinf) != 9:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Yoda" method, Hist2D Type needs 9 item ->\n\t\t( Name ID, Input file ID, Method, Hist type, Hist Name, Xlow, Xhigh, Ylow, Yhigh )!\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "yoda",
                "type": vinf[3].lower(),
                "name": vinf[4],
                "code": list(map("{:.6e}".format, list(map(float, vinf[5:]))))
            }
            self.vars.append(res)

    def set_par_yoda(self, vinf):
        if "yoda" in self.para:
            self.para['yoda'].append(
                {
                    "expr": vinf[0],
                    "type": vinf[3].lower(),
                    "name": vinf[4],
                    "code": list(map("{:.6e}".format, list(map(float, vinf[5:]))))
                }
            )
        else:
            self.para['yoda'] = [{
                "expr": vinf[0],
                "type": vinf[3].lower(),
                "name": vinf[4],
                "code": list(map("{:.6e}".format, list(map(float, vinf[5:]))))
            }]

    def get_yoda(self):
        ftem = YodaFile()
        ftem.file = self.file
        ftem.vinf = self.para['yoda']
        ftem.read()
        self.vars.update(ftem.vars)

    def set_oup_position(self, vinf):
        if len(vinf) != 5:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Position" method, 5 items ( Name ID, Input file ID, Method, Line number, Column number ) need to be provived.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        elif not (vinf[3].isdigit() and vinf[4].isdigit()):
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Position" method, Line number must be interger number\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "position",
                "code": [int(vinf[3].strip()), int(vinf[4].strip())]
            }
            self.vars.append(res)

    def set_par_position(self, vinf):
        if "position" in self.para:
            self.para['position'].append(
                {
                    "expr": vinf[0],
                    "code": [int(vinf[3]), int(vinf[4])]
                }
            )
        else:
            self.para['position'] = [{
                "expr": vinf[0],
                "code": [int(vinf[3]), int(vinf[4])]
            }]

    def get_position(self):
        from numpy import loadtxt
        data = loadtxt(self.file, ndmin=2)
        for var in self.para['position']:
            try:
                self.vars[var['expr']] = data[var['code']
                                              [0] - 1][var['code'][1] - 1]
            except:
                self.error = True
                self.etype = "ValueError"
                self.logger.error("\033[91m {}: Variable {} cannot be found in the output file {}.\033[0m".format(
                    self.etype, var['expr'], self.file))
                break
                # self.logger.error("\033[91m \033[0m")

    def set_oup_json(self, vinf):
        if len(vinf) < 4:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "Json" method, 4 items ( Name ID, Input file ID, Method, JSON Key ) need to be provided.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "json"
            }
            self.vars.append(res)

    def set_par_json(self, vinf):
        if "json" in self.para:
            self.para['json'].append({
                "key":  vinf[0],
                "expr": vinf[3]
            })
        else:
            self.para['json'] = [{
                "key":  vinf[0],
                "expr": vinf[3]
            }]

    def get_json(self):
        js = None
        with open(self.file, 'r') as f1:
            js = json.loads(f1.read())
        for var in self.para['json']:
            if var['expr'] in js:
                self.vars[var['key']] = js[var['expr']]
            else:
                self.error = True
                self.etype = "ValueError"
                self.logger.error("\033[91m {}: Variable {} cannot be found in the output file {}.\033[0m".format(
                    self.etype, var['expr'], self.file))
                break

    def set_oup_file(self, vinf):
        if not (len(vinf) == 3 or len(vinf) == 4):
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "File" method, 3 or 4 items ( Name ID, Input file, ID, Method, "save" ) need to be provided.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "file",
                "save": False
            }
            if len(vinf) == 4 and vinf[3].lower() == "save":
                res['save'] = True
            self.vars.append(res)
            self.set_file_with_pkgname(os.path.basename(self.file))
    
    def set_file_with_pkgname(self, basename):
            basename = basename.split(".")
            if len(basename) >= 2:
                basename[-2] += f"@{self.pkg}"
            else:
                basename[-1] += f"@{self.pkg}"
            basename = ".".join(basename)
            return basename


    def set_par_file(self, vinf):
        par = {
            "expr": vinf[0],
            "save": len(vinf) == 4 and "save" in vinf
        }
        self.para['file'] = par

    def get_file(self):
        from shutil import move
        if self.para['file']['save']:
            move(self.file, os.path.join(
                self.samppath, self.set_file_with_pkgname(os.path.basename(self.file))))
            self.file = os.path.join(
                self.samppath, self.set_file_with_pkgname(os.path.basename(self.file)))
            self.vars[self.para['file']['expr']] = self.file
        else:
            move(self.file, os.path.join(self.samppath,
                 "temp", self.set_file_with_pkgname(os.path.basename(self.file))))
            self.file = os.path.join(
                self.samppath, "temp", self.set_file_with_pkgname(os.path.basename(self.file)))
            self.vars[self.para['file']['expr']] = self.file

    def set_oup_slha(self, vinf, model="mssm"):
        if len(vinf) < 5:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "SLHA" method, 4 items ( Name ID, Input file ID, Method, SLHA information ) need to be provived.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        else:
            from slha import get_variable_info
            inf = get_variable_info(vinf)
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "slha",
            }
            if vinf:
                for kk, vv in inf.items():
                    res[kk] = vv
            else:
                self.logger.error(
                    '\033[91m Variable "{}" in program "{}" with "SLHA" method can not be decode, please check your variable setting\033[0m'.format(vinf[0], self.pkg))
                sys.exit(1)
            self.vars.append(res)

    def set_par_slha(self, vinf, model="mssm"):
        res = {
            "expr": vinf[0]
        }
        from slha import get_variable_info
        inf = get_variable_info(vinf)
        res.update(inf)
        if "slha" in self.para:
            self.para['slha'].append(res)
        else:
            self.para['slha'] = [res]

    def get_slha(self):
        try:
            spectr = xslha.read(self.file)
        except:
            self.error = True
            self.etype = "FileError"
            self.logger.error(
                "\033[91m ParseError: File {} cannot be read by SLHA reader \033[0m".format(self.file))
        from slha import read_variable_value_by_xSLHA
        for var in self.para['slha']:
            try:
                self.vars[var['expr']] = read_variable_value_by_xSLHA(
                    spectr, var)
            except:
                self.error = True
                self.etype = "ValueError"
                self.logger.error("\033[91m {}: Variable {} cannot be found in the output file {}.\033[0m".format(
                    self.etype, var['expr'], self.file))
                break

    def set_oup_xml(self, vinf):
        if len(vinf) < 4:
            self.logger.error(
                '\033[91m Variable "{}" in program "{}" with "XML" method, 4 items ( Name ID, Input file ID, Method, Match ) need to be provided.\033[0m'.format(vinf[0], self.pkg))
            sys.exit(1)
        elif len(vinf):
            res = {
                "expr": vinf[0],
                "file": self.file,
                "meth": "xml",
                "patt": vinf[3],
                "posi": int(vinf[4]),
                "code": vinf[5:]
            }
            self.vars.append(res)

    def set_par_xml(self, vinf):
        res = {
            "expr": vinf[0],
            "patt": vinf[3],
            "posi": int(vinf[4]),
            "code": vinf[5:]
        }
        if "xml" in self.para:
            self.para['xml'].append(res)
        else:
            self.para['xml'] = [res]

    def get_xml(self):
        try:
            with open(self.file, 'r') as f1:
                cont = xmltodict.parse(f1.read())
        except:
            self.error = True
            self.etype = "FileError"
            self.logger.error(
                "\033[91m ParseError: File {} cannot be read by XML reader \033[0m".format(self.file))
        from xmlreader import read_variable_value
        for var in self.para['xml']:
            try:
                self.vars[var['expr']] = read_variable_value(cont, var)
            except:
                self.error = True
                self.etype = "ValueError"
                self.logger.error("\033[91m {}: Variable {} cannot be found in the output file {}.\033[0m".format(
                    self.etype, var['expr'], self.file))

    def get_variables(self):
        if "position" in self.para:
            self.get_position()
        if "json" in self.para:
            self.get_json()
        if "file" in self.para:
            self.get_file()
        if "slha" in self.para:
            self.get_slha()
        if "xml" in self.para:
            self.get_xml()
        if "yoda" in self.para:
            self.get_yoda()


class Variables():
    def __init__(self) -> None:
        self.info = None
        self.sinfo = None

    def set_info(self, info):
        self.info = info


class YodaFile():
    def __init__(self) -> None:
        self.file = None
        self.vinf = None
        self.vars = {}
        self.hist1d = {}
        self.hist2d = {}
        self.type = None
        self.hist1ddata = {}

    def decode_histo1d(self, item):
        from numpy import loadtxt
        if self.type == "dat":
            hname = item[0].split()[-1].strip()
            culine = ["xlow","xhigh","val","errminus","errplus"]
            cuid = 0
            self.hist1ddata[hname] = {}
            for line in item:
                if "=" in line and line[0] != "#":
                    temlist = line.split("=")
                    if temlist[1].isdigit():
                        self.hist1ddata[hname][temlist[0]] = eval(temlist[1])
                    else:
                        self.hist1ddata[hname][temlist[0]] = temlist[1]
                if list(map(str.strip, line[2:].split())) == culine:
                    cuid = item.index(line)
                    break
                # if "ScaledBy" == line[0:8]:
            self.hist1d[hname] = loadtxt(item[cuid+1:], dtype=str)
            self.hist1ddata[hname]["data"] = pd.DataFrame(loadtxt(item[cuid+1:]), columns = culine)
            self.hist1ddata[hname]['sumw'] = self.hist1ddata[hname]['data'].val.sum()
            
            
    def decode_histo2d(self, item):
        if self.type == "dat":
            from numpy import loadtxt
            self.hist2d[item[0].split()[-1].strip()] = loadtxt(item[item.index("# xlow	 xhigh	 ylow	 yhigh	 val	 errminus	 errplus") + 1:], dtype=str)
        

    def read(self):
        import re
        if self.file.split('.')[-1] == "yoda":
            self.type = "yoda"
        elif self.file.split('.')[-1] == "dat":
            self.type = "dat"
        with open(self.file, 'r') as f1:
            yd = f1.read()
            yd = yd.replace("\n", " endl ")
            p_rec = re.compile(r"BEGIN(.*?)END", re.M)
            a = re.findall(p_rec, yd)
            for item in a:
                item = item.replace(" endl ", "\n").split('\n')
                if "HISTO1D" in item[0] and "RAW" not in item[0] and "[Weight_MERGING=0.000]" not in item[0]:
                    self.decode_histo1d(item)
                elif "HISTO2D" in item[0]:
                    self.decode_histo2d(item)
        if self.vinf:
            for var in self.vinf:
                if var['type'] == "hist1d":
                    data = self.hist1d[var['name']]
                    self.vars[var['expr']] = float(data[numpy.where(numpy.all(data[:, 0:2] == numpy.array(var['code']), axis=1))][0, 2])
                elif var['type'] == "hist2d":
                    data = self.hist2d[var['name']]
                    self.vars[var['expr']] = float(data[numpy.where(numpy.all(data[:, 0:4] == numpy.array(var['code']), axis=1))][0, 4])
            
class YodaPlotFile():
    def __init__(self):
        self.file = None
        self.cs   = None
        self.hist1d = {}
        self.hist2d = {}
        self.histlist = None 
        
    def decode_histo1d_inf(self, item):
        hname = item[0].split()[-1]
        item.pop(0)
        item.remove("")
        for line in item:
            if line[0] == "#":
                item.remove(line)
        res = {}
        for line in item:
            line = line.split("=")
            res[line[0]] = line[1].strip()
        self.hist1d[hname] = deepcopy(self.cs)
        self.hist1d[hname].update(res)
                
    def read(self):
        import re 
        with open(self.file, 'r') as f1:
            yd = f1.read()
            yd = yd.replace("\n", " endl ")
            p_rec = re.compile(r"BEGIN(.*?)END", re.M)
            a = re.findall(p_rec, yd)
            for item in a:
                item = item.replace(" endl ", '\n').split('\n')
                hname = item[0].split()[-1]
                if hname in self.histlist:
                    self.decode_histo1d_inf(item)
        

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj)
        elif isinstance(obj, numpy.floating):
            return float(obj)
        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def to_file_woblk(data, pth, method="to_csv"):
    """
    Write pandas.DataFrame table into csv file without blocking the main stream,
    or write python dict into json without blocking the main stream,
    this function is using the asyncio library

    Parameters
    -------------------
    data: A pandas.DataFrame or pandas.Series object for to_csv method;
          A python dict object for to_json method.
    pth:  save file path.
    method: default is "to_csv", or "to_json"

    Return
    -------------------
    Nothing 
    """
    import fcntl, asyncio, aiofiles, time, threading
    async def to_csv(data, pth):
        async with aiofiles.open(pth, 'w') as f1:
            fcntl.flock(f1, fcntl.LOCK_EX)
            data.to_csv(f1, index=False)
            fcntl.flock(f1, fcntl.LOCK_UN)

    async def to_json(data, pth):
        async with aiofiles.open(pth, 'w') as f1:
            await f1.write(json.dumps(data, indent=4))

    def run_coroutine(mtd):
        lp = asyncio.new_event_loop()
        asyncio.set_event_loop(lp)

        if mtd == "to_csv":
            lp.run_until_complete(to_csv(data=data, pth=pth))
        elif mtd == "to_json":
            lp.run_until_complete(to_json(data=data, pth=pth))
        lp.close()

    t = threading.Thread(target=run_coroutine, args=(method,))
    t.start()
