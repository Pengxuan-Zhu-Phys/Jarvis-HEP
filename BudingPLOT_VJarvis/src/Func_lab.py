#!/usr/bin/env python3

import os
import sys
import subprocess
import json
import time
pwd = os.path.abspath(os.path.dirname(__file__))


def init(argv):
    package = {'matplotlib': None, 'numpy': None, 'pandas': None, 'scipy': None, 'sympy': None,
               'python-ternary': None, 'opencv-python': None, 'seaborn': None, 'emoji': None, "adjustText": None, "pyslha": None, "xslha": None, "shapely": None}
    out = subprocess.getoutput("pip3 freeze").split('\n')
    tag = False
    for ii in package.keys():
        de = False
        for line in out:
            if line.split('==')[0] == ii:
                package[ii] = line.split('==')[1]
                de = True
        if not de:
            pipcmd = "pip3 install {} --user".format(ii)
            if len(argv) > 1:
                if argv[1] == '-s':
                    pipcmd += " -i https://pypi.tuna.tsinghua.edu.cn/simple/"
            so = subprocess.getstatusoutput(pipcmd)
            print("Installing Python Package {}".format(ii))
            if so[0]:
                tag = True
       
    if not tag:
        with open("{}/config.json".format(os.path.abspath(os.path.dirname(__file__))), 'w') as f1:
            json.dump(package, f1)
        import cv2
        output = subprocess.getoutput("xrandr | grep '*'").split('\n')
        sw, sh = 0, 0
        for line in output:
            if line.strip()[-1] == '*':
                sw, sh = line.split()[0].split('x')
        cv2.namedWindow("Welcome to BudingPLOT")
        img = cv2.imread('{}/image-src/BPicon.png'.format(pwd))
        width, hight, ch = img.shape
        cv2.moveWindow("Welcome to BudingPLOT", int(
            0.5*(int(sw) - hight)), int(0.38*(int(sh) - width)))
        img = cv2.resize(img, (hight, width))
        cv2.imshow("Welcome to BudingPLOT", img)
        cv2.waitKey(1)
        time.sleep(5)
        cv2.destroyAllWindows()
    alplt = 'alias Plot="{}"'.format(os.path.abspath(os.path.join(pwd, "..", "Plot")))
    print("\nPlease add line this line in ~/.bashrc file -> \n\n\t{} \n\nand run command\n\n\t{} ".format(alplt, "source ~/.bashrc"))

def helptool(argv):
    if len(argv) > 1:
        if argv[1] == "BP":
            print("BudingPLOT base path: {}".format(os.path.abspath(os.path.join(pwd, "..", "Plot"))))
        # print(argv)
    elif len(argv) == 1:
        if argv[0] == '-V':
            with open(os.path.join(pwd, "image-src/neofetch.dat"), 'r') as f1:
                print(f1.read())       
        else:
            print("Help (Beta Version) ")


def make_slhaReader_mass_input(fig):
    def set_mass_variable():
        variable = "variable    =    "
        for sect in fig['particle']:
            for pp in fig['particle'][sect]:
                # print("Func_lab", pp, fig['particle'][sect][pp]['PDG'])
                if fig['particle'][sect][pp]["PDG"] == "FindSLHA":
                    variable += "{},\t{},\tFindSLHA,\t{}\n                 ".format(
                        fig['slha']['patt'], pp, fig['particle'][sect][pp]['Mass'])
                else:
                    variable += "{},\t{},\tSLHA,\tBLOCK,\tMASS,\t{}\n                 ".format(
                        fig['slha']['patt'], pp, fig['particle'][sect][pp]["PDG"])
        return variable

    with open(fig['slha']['rcfg'], 'w') as f1:
        with open(fig['model']['slhaReadertemplet']) as f2:
            templet = f2.read()
            templet = templet.replace(
                ">>>SPECTRUM_DIR<<<", os.path.dirname(fig['slha']['file']))
            templet = templet.replace(">>>PATTERN<<<", ", ".join(
                [fig['slha']['patt'], fig['slha']['norm'], os.path.basename(fig['slha']['file'])]))
            templet = templet.replace(">>>SAVE_DIR<<<", fig['slha']['svdr'])
            templet = templet.replace(">>>SPECTRUM_ID<<<", fig['slha']['ID'])
            variable = set_mass_variable()
            f1.write(templet)
            f1.write(variable)


def check_slhaReader_pattern(patt):
    return os.path.exists(patt['file'])


def reset_slhaReader_mass_input(fig):
    def reset_particle_mass_pdg(pp, massinfo):
        for item in fig['slha']['pdginfo']:
            if item['particle'] == pp[0]:
                pp[1]['PDG'] = item['PDGID']
                break
        pp[1]['Mass'] = massinfo[pp[0]]

    def assign_particle_decay_list(particle, decay_list):
        for ss, sect in fig['particle'].items():
            for pp, inf in sect.items():
                if pp == particle:
                    fig['particle'][ss][pp]['Decay_list'] = decay_list

    def make_particle_decay_settting(massinfo):
        for kk in range(len(massinfo.keys())):
            # print(key, massinfo[key])
            dkl = []
            for ii in range(kk):
                dkl.append({
                    "fstate":       massinfo.keys()[ii],
                    "fspdg":        get_particle_pdg(massinfo.keys()[ii]),
                    "name":         "DK{}@{}".format(massinfo.keys()[kk], ii)
                })
            # print(kk, dkl)
            assign_particle_decay_list(massinfo.keys()[kk], dkl)

    def get_particle_pdg(pp):
        for item in fig['slha']['pdginfo']:
            if item['particle'] == pp:
                return item['PDGID']

    with open(os.path.join(fig['slha']['svdr'], "particles{}.json".format(fig['slha']['ID'])), "r") as f1:
        fig['slha']['pdginfo'] = json.loads(f1.read())
    import pandas as pd
    massinfo = pd.read_csv(os.path.join(
        fig['slha']['svdr'], "slhaReaderOutPut{}.csv".format(fig['slha']['ID'])))
    massinfo = pd.Series(massinfo.iloc[0])
    massinfo = massinfo.drop(labels=["Index"])
    massinfo = massinfo.apply(abs)
    massinfo = massinfo.sort_values(ascending=True)
    for ss, sect in fig['particle'].items():
        for pp in sect.items():
            reset_particle_mass_pdg(pp, massinfo)
    import configparser
    SR_config = configparser.ConfigParser()
    SR_config.read(fig['slha']['rcfg'])
    variable = ""
    for ss, sect in fig['particle'].items():
        for pp, inf in sect.items():
            variable += "{},\t{},\tSLHA,\tBLOCK,\tMASS,\t{}\n       ".format(
                fig['slha']['patt'], pp, inf["PDG"])
    make_particle_decay_settting(massinfo)
    for ss, sect in fig['particle'].items():
        for pp, inf in sect.items():
            for fs in inf['Decay_list']:
                variable += "{},\t{},\tFindBR,\t{},\t{}\n       ".format(
                    fig['slha']['patt'], fs['name'], inf["PDG"], fs['fspdg'])
    SR_config.set("Info", "variable", variable)
    SR_config.write(open(fig['slha']['rcfg'], 'w'))
    # for item in fig:
    # print(fig[item])


def block_screen_print():
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')


def enable_screen_print():
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


def load_particle_info(fig):
    import pandas as pd
    info = pd.read_csv(os.path.join(
        fig['slha']['svdr'], "slhaReaderOutPut{}.csv".format(fig['slha']['ID'])))
    info = pd.Series(info.iloc[0])
    info = info.drop(labels=['Index'])
    for ss, sect in fig['particle'].items():
        for pp, part in sect.items():
            for fs in part['Decay_list']:
                fs['br'] = info[fs['name']]

    
def load_rivet_plot_info(ps):
    import re
    ps = ps.replace("\n", " endl " )            
    # print(fig['ps']['setting'])
    p_rec = re.compile(r"BEGIN PLOT(.*?)END PLOT", re.M)
    a = re.findall(p_rec, ps)
    b = {}
    for item in a:
        res = {}
        item = item.replace(" endl ", "\n").split("\n")
        title = item[0].strip()
        item = item[1:]
        for ll in item:
            if ll:
                if ll[0] != "#":
                    ll = ll.split("=")
                    res[ll[0].strip()] = ll[1].strip()
        # print(res)
        b[title] = res
    return b

def load_rivet_data_info(ds):
    ds = ds.split('\n')
    cum = []
    sep = []
    for lin in ds:
        lin = lin.split(',')
        res = {
            "meth": lin[0].strip(),
            "wigh": float(lin[1].strip()),
            "yoda": lin[2].strip(),
            "titl": lin[3].strip()
        }
        if "Cum" in res['meth']:
            cum.append(res)
        elif "Sep" in res['meth']:
            sep.append(res)
    return cum, sep

def load_yoda_data(yoda):
    import re
    plots = []
    from numpy import loadtxt
    from pandas import DataFrame
    with open(yoda, 'r') as f1:
        yd = f1.read()
        yd = yd.replace("\n", " endl " )
        p_rec = re.compile(r"BEGIN(.*?)END", re.M)
        a = re.findall(p_rec, yd)
        for item in a:
            item = item.replace(" endl ", "\n").split('\n')
            if ("Weight_MERGING" in item[0]) or ("RAW" in item[0]):
                continue
            elif "YODA_HISTO1D_V2" in item[0]:
                res = {
                    "name": item[0].split()[1].strip(),
                    "type": item[0].split()[0].strip(),
                    "sumw": float(item[item.index("---") + 4].split()[2]),
                    "data": DataFrame(loadtxt(item[item.index("---") + 8:]), columns=item[item.index("---") + 7][1:].split())
                }
                plots.append(res)
            elif "YODA_HISTO2D_V2" in item[0]:
                res = {
                    "name": item[0].split()[1].strip(),
                    "type": item[0].split()[0].strip(),
                    "sumw": float(item[item.index("---") + 4].split()[2]),
                    "data": DataFrame(loadtxt(item[item.index("---") + 7:]), columns=item[item.index("---") + 7][1:].split())
                }
                plots.append(res)
    return plots
                
