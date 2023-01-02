#!/usr/bin/env python3

import os
import sys
import pandas as pd
import time
import configparser
import math
import sympy
import json
import emoji
import numpy as np 

pwd = os.path.abspath(os.path.dirname(__file__))
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.abspath(os.path.join(pwd, "BudingPLOT")))


class Plot():
    def __init__(self) -> None:
        self.ax = None
        self.fig = None
        self.cs = None
        self.data = None
        self.cf = None
        self.scf = None
        self.type = None
        self.path = {}
        self.inf = None
        self.figs = {}

    def read_config(self, cf):
        self.cf = configparser.ConfigParser()
        self.cf.read(cf)
        if "&PWD" in self.cf.get('PLOT_CONFI', 'path').strip().upper():
            self.path['path'] = self.cf.get('PLOT_CONFI', 'path').replace("&pwd", os.path.abspath(os.path.dirname(cf)))
            # self.path['path'] = os.path.abspath(os.path.dirname(cf))
        else:
            self.path['path'] = os.path.abspath(
                self.cf.get('PLOT_CONFI', 'path'))
        self.path['jarvis'] = jpath
        self.path['pwd'] = os.path.abspath(os.path.join(pwd, "BudingPLOT"))

    def decode_path(self, path):
        if "~" in path:
            path = path.replace("~", os.path.expanduser("~"))
        if "&BP" in path:
            path = path.replace("&BP", self.path['pwd'])
        if "&J" in path:
            path = path.replace("&J", self.path['jarvis'])
        return path

    def figures_inf(self):
        if self.cf.has_option("PLOT_CONFI", "save dir"):
            self.path['save dir'] = self.decode_path(os.path.join(
                self.path['path'], self.cf.get('PLOT_CONFI', "save dir")))
        elif self.cf.has_option("PLOT_CONFI", "scan config"):
            self.scf = configparser.ConfigParser()
            self.scf.read(os.path.join(
                self.path['path'], self.cf.get("PLOT_CONFI", "scan config")))
            self.path['save dir'] = os.path.abspath(os.path.join(
                self.decode_path(self.scf.get("Scan", 'save dir')), "PLOTs"))

        if not os.path.exists(self.path['save dir']):
            os.makedirs(self.path['save dir'])

        if not self.cf.has_option("PLOT_CONFI", "colorsetting"):
            cspath = os.path.abspath(os.path.join(
                self.path['pwd'], "Cards/preference.json"))
        else:
            cspath = self.decode_path(
                self.cf.get("PLOT_CONFI", "colorsetting"))
        with open(cspath, 'r') as f1:
            self.cs = json.loads(f1.read())

        if self.cf.has_option("PLOT_CONFI", "plot"):
            ppt = self.cf.get("PLOT_CONFI", "plot").split("\n")
            for line in ppt:
                line = list(map(str.strip, line.split(",")))
                self.figs[line[0]] = []
                if len(line) == 1 or (len(line) == 2 and "all" in line):
                    for sec in self.cf.sections():
                        if line[0] in sec:
                            self.figs[line[0]].append(sec)
                elif len(line) > 1:
                    for sec in line[1:]:
                        if "PLOT_{}_{}".format(line[0], sec) in self.cf.sections():
                            self.figs[line[0]].append(
                                "PLOT_{}_{}".format(line[0], sec))

    def set_fig_info(self, fig):
        fig.cf = self.cf
        fig.path = self.path
        fig.scf = self.scf

    def plot(self):
        self.figures_inf()
        for ftype, figs in self.figs.items():
            for plot in figs:
                if ftype == "Yoda_histo1d":
                    from RivetYODA import PlotYoda1D
                    fig = PlotYoda1D()
                    fig.type = "histo1d"
                    fig.cs = self.decode_path(self.cs['Yoda']['Setting'])
                    fig.pconfig = dict(dict(self.cf.items())["PLOT_CONFI"].items())
                    fig.inf = dict(dict(self.cf.items())[plot].items())
                    fig.inf['sect'] = plot
                    self.set_fig_info(fig)
                    if "name" not in fig.inf:
                        fig.inf['name'] = plot
                    fig.plot()
                elif ftype == "Jarvis":
                    from JarvisFlow import JarvisFlow
                    fig = JarvisFlow()
                    fig.type = "JarvisFlow"
                    fig.cs = self.decode_path(self.cs['JarvisFlow']['Setting'])
                    fig.pconfig = dict(dict(self.cf.items())[
                                       "PLOT_CONFI"].items())
                    fig.inf = dict(dict(self.cf.items())[plot].items())
                    fig.inf['sect'] = plot
                    self.set_fig_info(fig)
                    if "name" not in fig.inf:
                        fig.inf['name'] = plot
                    fig.plot()
                elif ftype == "Voronoi2D":
                    from Voronoi2D import Voronoi2D
                    fig = Voronoi2D()
                    fig.type = "Voronoi2D"
                    fig.pconfig = dict(dict(self.cf.items())[
                                       "PLOT_CONFI"].items())
                    fig.cs = self.decode_path(self.cs['Voronoi2D']['Setting'])
                    fig.inf = dict(dict(self.cf.items())[plot].items())
                    fig.inf['sect'] = plot
                    self.set_fig_info(fig)
                    if "name" not in fig.inf:
                        fig.inf['name'] = plot
                    fig.plot()

    def print_figure(self, plt):
        if "print mode" in self.inf.keys():
            if "show" not in self.inf['print mode'].lower() and "save" not in self.inf['print mode'].lower():
                self.inf['print mode'] = self.cs['default']['print mode']
            elif "show" in self.inf['print mode'].lower() and "save" not in self.inf['print mode'].lower():
                self.inf['print mode'] = ["show"]
            elif "show" not in self.inf['print mode'].lower() and "save" in self.inf['print mode'].lower():
                self.inf['print mode'] = ["save"]
            else:
                self.inf['print mode'] = ['show', 'save']
        else:
            self.inf['print mode'] = self.cs['default']['print mode']
        if "save" in self.inf['print mode']:
            ffg = plt
            if "/" == self.inf['name'][0]:
                self.inf['name'] = self.inf['name'][1:]
            ffp = os.path.join(self.path['save dir'], self.inf['name'])
            if self.format is None:
                self.format = self.cs['default']['save format']
                print(emoji.emojize(
                    '\t:ghost::ghost::ghost: figure save format is not specified by the user, using the default -> {} !!'.format(self.format), language="alias"))
            support_fmt_list = ['eps', 'pdf', 'pgf', 'png', 'raw',
                                'rgba', 'svg', 'svgz', 'jpg', 'jpeg', 'tif', 'tiff']
            self.format = list(set(self.format) & set(support_fmt_list))
            if not os.path.exists(os.path.dirname(ffp)):
                os.makedirs(os.path.dirname(ffp))
            print(emoji.emojize("\t:clock2: {:.2f} Sec;  Figure {} saved in the path\n\t\t-> {} \n\t\t>> {}.{}".format(
                time.time()-self.time, self.inf['name'], os.path.dirname(ffp), os.path.basename(ffp), ", >> {}.".format(os.path.basename(ffp)).join(self.format)), language="alias"))
            for fmt in self.format:
                if (fmt == 'pdf'):
                    ffg.savefig("{}.pdf".format(ffp, format='pdf'))
                else:
                    ffg.savefig("{}.{}".format(ffp, fmt), dpi=300)
        if "show" in self.inf['print mode']:
            plt.show(block=False)
            plt.pause(1)
            input(emoji.emojize("\n\t:clock2: {:.2f} Sec;  Press 'Enter' to continue ...\n".format(
                time.time()-self.time), language="alias"))
            plt.close()


class Figure(Plot):
    def __init__(self):
        self.type = None
        self.cf = None
        self.scf = None
        self.cs = None
        self.name = None
        self.path = None
        self.time = time.time()
        self.format = None
        self.vars = {}

    def load_variable(self):
        def load_var_data(xx):
            # print(self.data)
            # print("Loading variables info", self.vars[xx]['expr'])
            self.vars[xx]['data'] = self.data.eval(self.vars[xx]['expr'])
                        
        def load_var_info(xx):
            if "{}_label".format(xx) in self.inf:
                self.vars[xx]['label'] = self.inf['{}_label'.format(xx)]
            else:
                self.vars[xx]['label'] = self.vars[xx]['expr']
            if "{}_lim".format(xx) in self.inf:
                tlim = [self.vars[xx]['data'].min(), self.vars[xx]['data'].max()]
                lim = self.inf['{}_lim'.format(xx)]
                if "AUTO" in lim:
                    count = float(lim[5:])
                    tlim[0] = (tlim[0] // count) * count
                    tlim[1] = (tlim[1] // count + 1) * count
                    self.vars[xx]['lim'] = tlim
                else:
                    self.vars[xx]['lim'] = list(map(float, map(str.strip, lim.split(","))))  
            if "{}_ticks".format(xx) in self.inf:
                tick = self.inf['{}_ticks'.format(xx)]
                if "AUTO" == tick.upper()[0:4]:
                    a = float(tick.split("_")[-1])
                    low = self.vars[xx]['lim'][0] // a 
                    upp = self.vars[xx]['lim'][1] // a + 1
                    self.vars[xx]['ticks'] = np.linspace(low*a, upp*a, int(upp-low+1))
            if "{}_scale".format(xx) in self.inf:
                if self.inf['{}_scale'.format(xx)].strip().lower() in ['log', 'linear']:
                    self.vars[xx]['scale'] = self.inf['{}_scale'.format(xx)].strip().lower()
                else:
                    print("Illegal '{}_scale' in config file, using 'flat' scale as default".format(xx))
                    self.vars[xx]['scale'] = 'linear'
            if "{}_bin".format(xx) in self.inf:
                self.vars[xx]['bin'] = float(self.inf["{}_bin".format(xx)])
            
        for xx in ['x', 'y']:
            if "{}_variable".format(xx) in self.inf:
                self.vars[xx] = {"expr": self.inf['{}_variable'.format(xx)]}
                load_var_data(xx)
                load_var_info(xx)
        
        if "c_variable" in self.inf:
            self.vars["c"] = {"expr": self.inf['c_variable']}
            load_var_data('c')
            load_var_info('c')
            self.cbar = True

                    
            
            # print(self.inf['x_variable'])

    def plot(self):
        print("=== Ploting Figure : {} ===".format(self.inf['name']))
        self.drawpicture()

    def drawpicture(self):
        pass
