#!/usr/bin/env python3

import  os, sys 
sys.path.append("{}/NTools/".format(os.path.abspath(os.path.dirname(__file__))))
from NMSSMTools import NToolsOutput
import pandas as pd 
import numpy as np 
import configparser 
import xslha, pyslha
import time 



def convert_str_2_float(num):
    def is_number(num):
        try:
            float(num)
            return True
        except:
            return False

    if type(num) == str:
        if is_number(num):
            return float(num)
        bb = num.replace("d", 'e', 1)
        if is_number(bb):
            return float(bb)
        else:
            return num
    else:
        return num

class reader():
    def __init__(self):
        pass 

    def get_inf(self):
        self.cf = configparser.ConfigParser()
        self.cf.read(self.parser)
        self.logfile.write("\tInfo -> Using input parser file:\n\t\t{}\n".format(self.parser))

    def set_parser(self, inf):
        self.pwd = os.path.abspath(os.path.dirname(__file__))
        self.start = time.time()
        self.parser = inf 
        self.logfile = open('.'.join(self.parser.split('.')[:-1]) + '.info', 'w')
        with open(os.path.join(self.pwd, "Info/logo"), 'r') as f1:
            self.logfile.write(f1.read())
            self.logfile.write("\n\tInfo -> Recording time: \n\t\t{}\n".format(time.asctime()))
            self.logfile.flush()
        self.get_inf()
        self.check_parser()
        self.set_value()
        self.logfile.close()
        self.record_mass_info()
        # self.printInfo()

    def check_parser(self):
        if not self.cf.has_section("Info"):
            print("Error -> No Section\tInfo founded!!!\n\tPlease check your input parser file -> {}".format(self.parser))
            sys.exit(0)
        if not os.path.exists("{}/Info/parser.json".format(self.pwd)):
            print("Error -> Can not load the parser checker ...")
            sys.exit(0)
        import json
        with open("{}/Info/parser.json".format(self.pwd), 'r') as f1:
            self.parser_checker = json.load(f1)
        with open("{}/Info/PDG.json".format(self.pwd), 'r') as f1:
            self.pdgset = json.load(f1)
        for key, value in self.parser_checker.items():
            if not value == "False":
                if not self.cf.has_option("Info", key):
                    print("Error -> No Option\t'{}'\t found in parser file\n\tPlease check your parser file -> {}".format(key, self.parser))
                    sys.exit(0)
        self.task = {}
        if self.cf.get("Info", "format") in self.parser_checker['format']:
            print("Info -> Using SLHA format : {}".format(self.cf.get("Info", "format")))
            self.logfile.write("\tInfo -> Using SLHA format : {}\n".format(self.cf.get("Info", "format")))
            self.task["format"] =  self.cf.get("Info", "format")
        if self.cf.get("Info", "method") == 'Spectrum':
            self.task["method"] = "Spectrum" 
        elif self.cf.get("Info", "method") == "Folder":
            self.task['method'] = "Folder"
        if not (self.cf.get("Info", "model") in self.parser_checker["model"][self.task['format']]):
            print("Error -> Model :\t{}\t not match to '{}' format\n\tPlease check your parser file -> {}".format(self.cf.get("Info", "model"), self.task['format'], self.parser))
            sys.exit(0)
        self.task['path'] = self.cf.get("Info", "path")
        if "$pwd$" in self.task['path']:
            self.task['path'] = self.task['path'].replace("$pwd$", os.path.dirname(self.parser))
        self.getcsvpath()
        self.load_pdg()
        self.spectrpath()
        self.get_variables()
        self.record()

    def getcsvpath(self):
        self.datacsv = self.cf.get("Info", "saveto")
        if "$pwd$" in self.datacsv:
            self.datacsv = self.datacsv.replace("$pwd$", os.path.dirname(self.parser))
        if not os.path.exists(os.path.dirname(self.datacsv)):
            os.makedirs(os.path.dirname(self.datacsv))
        print("Info -> result data will be saved in the file : \n\t{}\n".format(self.datacsv))
        self.logfile.write("\tInfo -> result data will be saved in the file : \n\t\t{}\n".format(self.datacsv))

    def spectrpath(self):
        def matchpattern(file, dfpatt):
            k, index = "", ""
            for key, value in dfpatt.items():
                value = value.split("$ID$")
                index = file
                for v in value:
                    index = index.replace(v, "")
                temp = file.split(index)
                if temp == value:
                    k = key 
                    break 
            # print("SpectrMatch", file, k, index)
            return k, index
        if self.task['format'] == "NMSSMTools" and os.path.isdir(self.task['path']) and self.task['method'] == "Spectrum":
            pattern = self.cf.get("Info", "pattern").split("\n")
            self.task['pattern'] ={}
            self.task['documents'] = {}
            for p in pattern:
                patt = p.split(",")[0].strip()
                if self.parser_checker["pattern"][self.task['format']]:
                    if patt not in self.parser_checker["pattern"]['NMSSMTools']:
                        print("Error -> Pattern\t{}\tfor NMSSMTools method not available!!!".format(patt))
                        sys.exit(0)
                self.task['pattern'][p.split(',')[0].strip()] = p.split(',')[1].strip()
            for p in pattern:
                patt, index = matchpattern(p.split(',')[2].strip(), self.task['pattern'])
                slha = os.path.join(self.task['path'], p.split(",")[2].strip())
                if os.path.exists(slha):
                    if index not in self.task['documents'].keys():
                        self.task['documents'][index] = {}
                    self.task['documents'][index][patt] = slha 
        elif self.task['format'] == "NMSSMTools" and os.path.isdir(self.task['path']) and self.task['method'] == "Folder":
            pattern = self.cf.get("Info", "pattern").split("\n")
            self.task['pattern'] = {}
            for ii in range(len(pattern)):
                self.task['pattern'][pattern[ii].split(",")[0].strip()] = pattern[ii].split(",")[1].strip()
            self.task['documents'] = {}
            for file in os.listdir(self.task['path']):
                patt, index = matchpattern(file, self.task['pattern'])
                if patt not in self.parser_checker["pattern"]['NMSSMTools']:
                    print("Error -> Pattern\t{}\tfor NMSSMTools method not available!!!".format(patt))
                if self.task['path'] == '/':
                    slha = "{}{}".format(self.task['path'], file)
                else:
                    slha = "{}/{}".format(self.task['path'], file)
                if os.path.exists(slha):
                    # print(slha)
                    if index not in self.task['documents'].keys():
                        self.task['documents'][index] = {}
                    self.task['documents'][index][patt]=slha
            # print(self.task['documents'])
        elif self.task['format'] == "Standard SLHA" and os.path.isdir(self.task['path']) and self.task['method'] == "Spectrum":
            pattern = self.cf.get("Info", "pattern").split("\n")
            self.task['pattern'] =  {}
            self.task['documents'] = {}
            for p in pattern:
                patt = p.split(",")[0].strip()
                self.task['pattern'][p.split(',')[0].strip()] = p.split(',')[1].strip()
            for p in pattern:
                patt, index = matchpattern(p.split(',')[2].strip(), self.task['pattern'])
                slha = os.path.join(self.task['path'], p.split(',')[2].strip())
                if os.path.exists(slha):
                    if index not in self.task['documents'].keys():
                        self.task['documents'][index] = {}
                    self.task['documents'][index][patt] = slha
                else:
                    print("Error -> Pattern\t{}\tSLHA file can not be found! \n\t{}".format(patt, slha))
                    sys.exit(0)
        elif self.task['format'] == "Standard SLHA" and os.path.isdir(self.task['path']) and self.task['method'] == "Folder":
            pattern = self.cf.get("Info", "pattern").split("\n")
            self.task['pattern'] = {}
            for p in pattern:
                self.task['pattern'][p.split(',')[0].strip()] = p.split(',')[1].strip()
            self.task['documents'] = {}
            for file in os.listdir(self.task['path']):
                patt, index = matchpattern(file, self.task['pattern'])
                slha = os.path.join(self.task['path'], file)
                if index not in self.task['documents'].keys():
                    self.task['documents'][index] = {}
                self.task['documents'][index][patt] = slha
            # print(self.task['documents'])

    def printInfo(self):
        # print(self.var)
        if self.isBP:
            print(self.particle_Mass)
        # for item in self.var['var']:
            # print(item)
        # with open("{}/")

    def load_path(self, path):
        if "$pwd$" in path:
            path = os.path.abspath(os.path.join(os.path.dirname(self.parser), path.replace("$pwd$/", "")))
        return path
        

    def record_mass_info(self):
        if self.isBP:
            mjpath = self.load_path(self.cf.get("Info", "savemassinfo").split(",")[1].strip())
            with open(mjpath, 'w') as f1:
                import json 
                f1.write(json.dumps(self.particle_Mass))
            

    def record(self):
        def record_info_line(info, value):
            self.logfile.write("\tInfo -> {} :\n\t\t{}\n".format(info, value))
        def record_var_info(variable):
            # print(variable)
            if variable['meth'] == "SLHA":
                if variable['info'] == "WIDTH":
                    self.logfile.write("{},\t{},\tSLHA,\tDECAY,\t{}\n\t\t\t".format(variable['patt'], variable['name'], variable['code']))
                elif variable['info'] == "DECAY":
                    for fstate in variable['finalstates']:
                        self.logfile.write("{},\t{},\tSLHA,\tDECAY,\t{},\t{},\t{}\n\t\t\t".format(variable['patt'], variable['name'], variable['pdg'], len(fstate), ",\t".join(map(str, fstate))))
                else:
                    if type(variable['code']) == tuple:
                        self.logfile.write("{},\t{},\tSLHA,\tBLOCK,\t{},\t{}\n\t\t\t".format(variable['patt'], variable['name'], variable['info'], ",\t".join(map(str, variable['code']))))
                    else:
                        self.logfile.write("{},\t{},\tSLHA,\tBLOCK,\t{},\t{}\n\t\t\t".format(variable['patt'], variable['name'], variable['info'], variable['code']))
            if variable['meth'] == "FindBR":
                self.logfile.write("{},\t{},\tFindBR,\t{},\t{}\n\t\t\t".format(variable['patt'], variable['name'], variable['info'], variable['code']))

        record_info_line("Using reader method", self.cf.get("Info", "method"))
        record_info_line("Using input SLHA file", self.cf.get("Info", "path"))

        self.logfile.write("\nRecording the variables information, which can be used directly in the parser file:\n\n[Info]\nvariable =\t")
        for var in self.var['var']:
            record_var_info(var)
        self.logfile.flush()

    def set_csv_info(self):
        self.data['comn'] = {
            "Index":    "ID"
        }
        self.data['columns'] = ["Index"]
        for item in self.var['var']:
            self.data['comn'][item['name']] = np.NaN
            self.data['columns'].append(item['name'])

    def set_value(self):
        self.isBP = False
        if self.cf.has_option("Info", "savemassinfo"):
            if eval(self.cf.get("Info", "savemassinfo").split(',')[0]):
                self.particle_Mass = []
                self.isBP = True
        self.data = {}
        self.data['counting'] = 0
        self.set_csv_info()
        self.data['data'] = pd.DataFrame(columns=self.data['columns'])
        for index, files in self.task['documents'].items():
            data = pd.Series(self.data['comn'])
            data['Index'] = index
            data = self.get_data_from_files(data, files)
            self.data['data'] = self.data['data'].append(data, ignore_index=True)
            self.data['counting'] += 1
            print("Finish \t{} :Time consuming -> {:.4f} second".format(self.data['counting'], time.time()-self.start), end='\r')
        self.data['data'].to_csv(self.datacsv, index=False)

    def get_data_from_files(self, data, files):
        def getvalue(spectr, var, readlib):
            if var['meth'] == "SLHA":
                if var['info'] == "MASS" and self.isBP and self.task['method'] == "Spectrum":
                    self.particle_Mass.append({
                        "particle": var['name'],
                        "PDGID":    var['code']
                    })
                if readlib == "PySLHA":
                    if var['info'] == "WIDTH":
                        width = 0
                        try:
                            width = spectr.decays[var["code"]].totalwidth
                        except:
                            pass 
                        return width
                    elif var['info'] == "DECAY":
                        br = 0. 
                        # print(var)
                        from itertools import permutations
                        try:
                            for item in spectr.decays[var['pdg']].decays:
                                for fstate in var['finalstates']:
                                    for ids in permutations(fstate):
                                        if item.ids == list(ids):
                                            br += item.br
                        except:
                            pass 
                        return br 
                    else:
                        try:
                            value = spectr.blocks[var['info']][var['code']]
                        except:
                            # print(spectr.blocks)
                            value = np.NaN
                        value = convert_str_2_float(value)
                        return value 
                elif readlib == "xSLHA":
                    if var['info'] == 'DECAY':
                        br = 0. 
                        from itertools import permutations
                        try:
                            for fstate in var['finalstates']:
                                br += spectr.Value("BR", [var['pdg'], fstate])
                        except:
                            br += 0. 
                        return br 
                    elif var['info'] == "WIDTH":
                        width = 0. 
                        try:
                            width = spectr.Value("WIDTH", var['code'])
                        except:
                            pass 
                        return width
                    else:
                        if type(var['code']) == int:
                            value = spectr.Value(var['info'], [var['code']])
                        else:
                            value = spectr.Value(var['info'], var['code'])

                        return float(value)
            elif var['meth'] == "FindSLHA":
                # print(var)
                if readlib == "PySLHA":
                    try:                      
                        for bb in self.pdgset[self.modelid]['MixingMatrix']:
                            if bb['block'] == var['mixingblock'].upper():
                                dim = bb['Dimension']
                        templist = []
                        for ii in range(dim[0]):
                            templist.append(spectr.blocks[var['mixingblock']][tuple([ii+1, var['code']])]**2 )
                        index = templist.index(max(templist))
                        value = spectr.blocks[var['info']][var['target'][index]]
                        if var['info'] == "MASS" and self.isBP and self.task['method'] == "Spectrum":
                            self.particle_Mass.append({
                                "particle": var['name'],
                                "PDGID":    var['target'][index]
                            })
                        return value 
                    except:
                        return np.NaN
                elif readlib == "xSLHA":
                    try:
                        for bb in self.pdgset[self.modelid]['MixingMatrix']:
                            if bb['block'] == var['mixingblock'].upper():
                                dim = bb['Dimension']
                        templist = []
                        for ii in range(dim[0]):
                            templist.append( spectr.Value(var['mixingblock'], tuple([ii+1, var['code']]))**2 )
                        index = templist.index(max(templist))

                        if type(var['target'][index]) == int:
                            value = spectr.Value(var['info'], [var['target'][index]])
                        else:
                            value = spectr.Value(var['info'], var['target'][index])
                        if var['info'] == "MASS" and self.isBP and self.task['method'] == "Spectrum":
                            self.particle_Mass.append({
                                "particle": var['name'],
                                "PDGID":    var['target'][index]
                            })
                        return value
                    except:
                        return np.NaN
            elif var['meth'] == "FindBR":
                if readlib == "PySLHA":
                    br = 0.
                    try:
                        for item in spectr.decays[var["info"]].decays:
                            if var['code'] in tuple(map(abs, item.ids)):
                                br += item.br
                            # print(item.br, list(map(abs, item.ids)))
                    except:
                        pass 
                    # if br > 0.1:
                        # print(spectr.decays[var["info"]].decays  )
                        # print("Info: particle {} -> {} : {}".format(var['info'], var['code'], br))
                    return br 
                elif readlib == "xSLHA":
                    # print(var)
                    br = 0. 
                    try:
                        for fs, vv in spectr.br[var["info"]].items():
                            if var['code'] in tuple(map(abs, fs)):
                                br += vv
                    except:
                        pass 
                    return br

        def getvalue_from_NTools(spectr, var):
            if var['meth'] == "SLHA":
                # print(var)
                if var['info'] == "DECAY":
                    value = 0. 
                    from itertools import permutations
                    for fstate in var['finalstates']:
                        for ids in permutations(fstate):
                            br = 0.
                            try:
                                br = spectr(var['info'], var['pdg'], ids)
                                value += br
                            except:
                                pass
                    return value
                elif var['info'] == "WIDTH":
                    try: 
                        value = spectr(var['info'], var['code'])
                    except:
                        value = 0. 
                    return value
                else:
                    # print(var)
                    value = spectr(var['info'], var['code'])
                    return value
                    # print("value -> ", value, "for var -> ", var)
            elif var['meth'] == "FindSLHA":
                try:
                    for bb in self.pdgset[self.modelid]['MixingMatrix']:
                        if bb['block'] == var['mixingblock'].upper():
                            dim = bb['Dimension']
                    templist = []
                    for ii in range(dim[0]):
                        templist.append( spectr(var['mixingblock'], tuple([ii+1, var['code']]))**2 )
                    index = templist.index(max(templist))
                    value = spectr(var['info'], var['target'][index])
                    return value 
                except:
                    return np.NaN


        if self.task['format'] == "Standard SLHA":
            for patt, slha in files.items():
                readlib = "xSLHA"
                # readlib = "PySLHA"
                try:
                    spectr = xslha.read(slha)
                    # spectr = pyslha.read(slha, ignorenomass=True)
                    # print("Read method -> PySLHA for file {}".format(slha))
                except:
                    # spectr = xslha.read(slha)
                    spectr = pyslha.read(slha, ignorenomass=True)
                    # print("Read method -> xSLHA for file {}".format(slha))
                    readlib = "PySLHA"
                    # readlib = "xSLHA"
                for var in self.var['var']:
                    if var['patt'] == patt:
                        data[var['name']] = getvalue(spectr, var, readlib)
        elif self.task['format'] == "NMSSMTools":
            if "spectr" in files.keys():
                if "omega" in files.keys():
                    spectr = NToolsOutput(files['spectr'], files['omega'])
                else:
                    spectr = NToolsOutput(files['spectr'])

                for var in self.var['var']:
                    data[var['name']] = getvalue_from_NTools(spectr, var)
        return data 

    def load_pdg(self):
        self.modelid = 0
        for ii in range(len(self.pdgset)):
            if self.pdgset[ii]['Model'] == self.cf.get("Info", "model").strip():
                self.modelid = ii
        self.particle = self.pdgset[self.modelid]['particles']
        if self.pdgset[self.modelid]['Include']:
            for model in self.pdgset[self.modelid]['Include']:
                for ii in range(len(self.pdgset)):
                    if self.pdgset[ii]['Model'] == model:
                        for p in self.pdgset[ii]['particles']:
                            tag = True
                            for pp in self.particle:
                                if p["PDG_code"] == pp['PDG_code']:
                                    tag = False
                            if tag:
                                self.particle.append(p)

    def get_variable_name(self, varinfo):
        if varinfo[2].strip().upper() == "SLHA":
            return varinfo[1].strip()
        elif varinfo[1].strip().upper() == "SLHA":
            name = ""
            if varinfo[2].strip().upper() == "BLOCK":
                for key, ls in self.pdgset[self.modelid]["Namespace"].items():
                    if varinfo[3].strip().upper() in ls:
                        name += key
                        break
                if varinfo[3].strip().upper() in self.pdgset[self.modelid]['Block']["BLOCK"].keys():
                    if len(varinfo[4:]) == 1:
                        if varinfo[4].strip() in self.pdgset[self.modelid]['Block']["BLOCK"][varinfo[3].strip().upper()].keys():
                            name += varinfo[3].strip().upper()[0:2]
                            name +=  self.pdgset[self.modelid]['Block']["BLOCK"][varinfo[3].strip().upper()][varinfo[4].strip()]
                        # if len(varinfo[4:]) == 1:
                elif (varinfo[3].strip().upper() == "MASS"):
                    tag = True
                    for pp in self.particle:
                        if varinfo[4].strip() == pp['PDG_code']:
                            tag = False
                            name += pp['name']
                            break
                    if tag:
                        name += "pdg"
                        name += varinfo[4].strip()
            elif varinfo[2].strip().upper() == "DECAY":
                if (len(varinfo) == 4):
                    name += "W"
                    tag = True
                    for pp in self.particle:
                        if varinfo[3].strip() == pp['PDG_code']:
                            tag = False
                            name += pp['name']
                            break
                    if tag:
                        name += varinfo[3].strip()
                if (len(varinfo) == 7):
                    name += "DK"
                    mp = varinfo[3].strip()
                    sp1 = varinfo[5].strip()
                    sp2 = varinfo[6].strip()
                    for pp in self.particle:
                        if mp == pp['PDG_code']:
                            mp = pp['name']
                        if sp1[0] == '-' and sp1[1:] == pp['PDG_code']:
                            sp1 = pp['antiname']
                        elif sp1 == pp['PDG_code']:
                            sp1 = pp['name']
                        if sp2[0] == '-' and sp2[1:] == pp['PDG_code']:
                            sp2 = pp['antiname']
                        elif sp2 == pp['PDG_code']:
                            sp2 = pp['name']
                    name += mp + "2" + sp1 + sp2
                if (len(varinfo) == 8):
                    name += "DK"
                    mp = varinfo[3].strip()
                    sp1 = varinfo[5].strip()
                    sp2 = varinfo[6].strip()
                    sp3 = varinfo[7].strip()
                    for pp in self.particle:
                        if mp == pp['PDG_code']:
                            mp = pp['name']
                        if sp1[0] == '-' and sp1[1:] == pp['PDG_code']:
                            sp1 = pp['antiname']
                        elif sp1 == pp['PDG_code']:
                            sp1 = pp['name']           
                        if sp2[0] == '-' and sp2[1:] == pp['PDG_code']:
                            sp2 = pp['antiname']
                        elif sp2 == pp['PDG_code']:
                            sp2 = pp['name']                            
                        if sp3[0] == '-' and sp3[1:] == pp['PDG_code']:
                            sp3 = pp['antiname']
                        elif sp3 == pp['PDG_code']:
                            sp3 = pp['name']    
                    name += mp + "2" + sp1 + sp2 + sp3                        
        elif varinfo[2].strip().upper() == "FINDSLHA":
            return varinfo[1].strip()
        elif varinfo[2].strip().upper() == "FINDBR":
            return varinfo[1].strip()
        elif varinfo[1].strip().upper():
            name = "obs_{}".format(self.var['order'])
            self.var['order'] += 1                    
        return name

    def get_variables(self):
        variable = self.cf.get("Info", "variable").split("\n")
        # print(variable)
        if self.task['format'] == "NMSSMTools":
            self.var = {
                "patternused":  [],
                "var":  [],
                "order":    1
            }
            for var in variable:
                dfname = 0
                varinfo = var.split(',')
                if varinfo[0].strip() not in self.task['pattern'].keys():
                    print("Error -> pattern '{}' not founded for variable".format(varinfo[0], var))
                    sys.exit(0)
                elif varinfo[0].strip() not in self.var['patternused']:
                    self.var['patternused'].append(varinfo[0].strip())
                if varinfo[1].strip().upper() == "SLHA" or varinfo[1].strip() == "FindSLHA" or varinfo[1].strip() == "FindBR":
                    # print("Info -> Using default name for variable:\n\t{}".format(var))
                    pass 
                    # print()
                elif varinfo[2].strip().upper() == "SLHA" or varinfo[2].strip() == "FindSLHA" or varinfo[2].strip() == "FindBR":
                    dfname = 1
                else:
                    print("Error -> Only 'SLHA' method support in this version")
                    sys.exit(0)
                if varinfo[1 + dfname].strip() == "SLHA":
                    if varinfo[2 + dfname].strip().upper() == "BLOCK":
                        if len(varinfo[3+dfname:]) == 2:
                            self.var['var'].append({
                                "name": self.get_variable_name(varinfo),
                                "info": varinfo[3+dfname].strip().upper(),
                                "patt": varinfo[0].strip(),
                                "meth": varinfo[1+dfname].strip(),
                                "code": int(varinfo[4+dfname].strip())
                            })
                        elif len(varinfo[3+dfname:]) > 2:
                            self.var['var'].append({
                                "name": self.get_variable_name(varinfo),
                                "info": varinfo[3+dfname].strip().upper(),
                                "patt": varinfo[0].strip(),
                                "meth": varinfo[1+dfname].strip(),
                                "code": tuple(list(map(int, varinfo[4+dfname:])))
                            })
                        elif len(varinfo[3+dfname:]) == 1:
                            # print(self.pdgset[self.modelid]['MixingMatrix'])
                            for bb in self.pdgset[self.modelid]['MixingMatrix']:
                                if bb['block'] == varinfo[3+dfname].strip().upper():
                                    for ii in range(bb['Dimension'][0]):
                                        for jj in range(bb['Dimension'][1]):
                                            if dfname:
                                                self.var['var'].append({
                                                    "name": "{}{}{}".format(varinfo[1].strip(), ii+1, jj+1),
                                                    "meth": varinfo[1+dfname].strip(),
                                                    "patt": varinfo[0].strip(),
                                                    "info": varinfo[3+dfname].strip().upper(),
                                                    "code": (ii+1, jj+1)
                                                })
                                            else:
                                                self.var['var'].append({
                                                    "name": "{}{}{}".format(bb['name'], ii+1, jj+1),
                                                    "info": varinfo[3+dfname].strip().upper(),
                                                    "patt": varinfo[0].strip(),
                                                    "meth": varinfo[1+dfname].strip(),
                                                    "code": (ii+1, jj+1)
                                                })
                    elif varinfo[2 + dfname].strip().upper() == "DECAY":
                        if len(varinfo[3+dfname:]) == 1:
                            self.var['var'].append({
                                "name": self.get_variable_name(varinfo),
                                "info": "WIDTH",
                                "patt": varinfo[0].strip(),
                                "meth": varinfo[1+dfname].strip(),
                                "code": (int(varinfo[3+dfname].strip()))
                            })
                        elif len(varinfo[5+dfname:]) >= 2:
                            dktag = True
                            for ii in range(len(self.var['var'])):
                                if (self.var['var'][ii]['info'] == "DECAY") and (self.get_variable_name(varinfo) == self.var['var'][ii]['name']) and self.var['var'][ii]['name'] != '':
                                    self.var['var'][ii]['finalstates'].append(tuple(list(map(int, varinfo[5+dfname:]))))
                                    dktag = False
                                    break
                            if dktag:
                                self.var['var'].append({
                                    "name": self.get_variable_name(varinfo),
                                    "patt": varinfo[0].strip(),
                                    "info": "DECAY",
                                    "meth": varinfo[1+dfname].strip(),
                                    "pdg":  int(varinfo[3+dfname]),
                                    "finalstates":  [tuple(list(map(int, varinfo[5+dfname:])))]
                                })
                        else:
                            print("Error -> Variable information Unvaliable:\n\t{}".format(var))
                            sys.exit(0)
                elif varinfo[1+dfname].strip() == "FindSLHA":
                    self.var['var'].append({
                        "name":         self.get_variable_name(varinfo),
                        "mixingblock":  varinfo[2+dfname].strip().upper(),
                        "info":         varinfo[4+dfname].strip(),
                        "patt":         varinfo[0].strip(),
                        "meth":         varinfo[1+dfname].strip(),
                        "code":         int(varinfo[3+dfname].strip()),
                        "target":       eval(','.join(varinfo[5+dfname:]))
                    })
                elif varinfo[1+dfname].strip() == "FindBR":
                    self.var['var'].append({
                        "name":     self.get_variable_name(varinfo),
                        "meth":     varinfo[1+dfname].strip(),
                        "info":     int(varinfo[2+dfname].strip()),
                        "code":     int(varinfo[3+dfname].strip()),
                        "patt":     varinfo[0].strip()
                    })
        if self.task['format'] == "Standard SLHA":
            self.var = {
                "var":  [],
                "order":    1
            }
            for var in variable:
                dfname = 0
                varinfo = var.split(',')
                if varinfo[0].strip() not in self.task['pattern'].keys():
                    print("Error -> pattern '{}' not founded for variable".format(varinfo[0], var))
                    sys.exit(0)
                if varinfo[1].strip() == "SLHA" or varinfo[1].strip() == "FindSLHA" or varinfo[1].strip() == "FindBR":
                    # print("Info -> Using default name for variable:\n\t{}".format(var))
                    pass 
                elif varinfo[2].strip() == "SLHA" or varinfo[2].strip() == "FindSLHA" or varinfo[2].strip() == "FindBR":
                    dfname = 1
                else:
                    print("Error -> Only 'SLHA' and 'FindSLHA' method support in this version")
                    sys.exit(0)
                if varinfo[1+dfname].strip() == "SLHA":
                    if varinfo[2+dfname].strip().upper() == "BLOCK":
                        if len(varinfo[3+dfname:]) == 2:
                            self.var["var"].append({
                                "name": self.get_variable_name(varinfo),
                                "info": varinfo[3+dfname].strip().upper(),
                                "meth": varinfo[1+dfname].strip(),
                                "code": (int(varinfo[4+dfname].strip())),
                                "patt": varinfo[0].strip()
                            })
                        elif len(varinfo[3+dfname:]) > 2:
                            self.var['var'].append({
                                "name": self.get_variable_name(varinfo),
                                "info": varinfo[3+dfname].strip().upper(),
                                "meth": varinfo[1+dfname].strip(),
                                "code": tuple(list(map(int, varinfo[4+dfname:]))),
                                "patt": varinfo[0].strip()

                            })
                        elif len(varinfo[3+dfname:]) == 1:
                            for bb in self.pdgset[self.modelid]['MixingMatrix']:
                                if bb['block'] == varinfo[3+dfname].strip().upper():
                                    for ii in range(bb['Dimension'][0]):
                                        for jj in range(bb['Dimension'][1]):
                                            if dfname:
                                                tempname = varinfo[1].strip()   
                                            else:
                                                tempname = bb['name']
                                            self.var['var'].append(
                                                {
                                                    "name": "{}{}{}".format(tempname, ii+1, jj+1),
                                                    "info": varinfo[3+dfname].strip().upper(),
                                                    "meth": varinfo[1+dfname].strip(),
                                                    "code": (ii+1, jj+1),
                                                    "patt": varinfo[0].strip()
                                                }
                                            )

                    elif varinfo[2+dfname].strip().upper() == "DECAY":
                        # print("->", var)
                        if len(varinfo[3+dfname:]) == 1:
                            self.var['var'].append({
                                "name": self.get_variable_name(varinfo),
                                "info": "WIDTH",
                                "meth": varinfo[1+dfname].strip(),
                                "code": (int(varinfo[3+dfname].strip())),
                                "patt": varinfo[0].strip()
                            })
                        elif len(varinfo) >= 7:
                            dktag = True
                            for ii in range(len(self.var['var'])):
                                if (self.var['var'][ii]['info'] == "DECAY") and (self.get_variable_name(varinfo) == self.var['var'][ii]['name']) and self.var['var'][ii]['name'] != '':
                                    self.var['var'][ii]['finalstates'].append(tuple(list(map(int, varinfo[5+dfname:]))))
                                    dktag = False
                                    break
                            if dktag:
                                self.var['var'].append({
                                    "name": self.get_variable_name(varinfo),
                                    "info": "DECAY",
                                    "meth": varinfo[1+dfname].strip(),
                                    "pdg":  int(varinfo[3+dfname]),
                                    "patt": varinfo[0].strip(),
                                    "finalstates":  [tuple(list(map(int, varinfo[5+dfname:])))]
                                })
                        else:
                            print("Error -> Variable infomation Unvaliable:\n\t{}".format(var))
                            sys.exit(0)
                elif varinfo[1+dfname].strip() == "FindSLHA":
                    # print(varinfo)
                    self.var['var'].append({
                        "name":         self.get_variable_name(varinfo),
                        "mixingblock":  varinfo[2+dfname].strip().upper(), 
                        "info":         varinfo[4+dfname].strip(),
                        "meth":         varinfo[1+dfname].strip(),
                        "patt":         varinfo[0].strip(),
                        "code":         int(varinfo[3+dfname].strip()),
                        "target":       eval(','.join(varinfo[5+dfname:]))
                    })
                elif varinfo[1+dfname].strip() == "FindBR":
                    self.var['var'].append({
                        "name":     self.get_variable_name(varinfo),
                        "meth":     varinfo[1+dfname].strip(),
                        "info":     int(varinfo[2+dfname].strip()),
                        "code":     int(varinfo[3+dfname].strip()),
                        "patt":     varinfo[0].strip()
                    })
        # print(self.var['var'])
 