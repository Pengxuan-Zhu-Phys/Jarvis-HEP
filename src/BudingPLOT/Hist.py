#!/usr/bin/env python3
import os
import sys
import matplotlib
import pandas as pd
import numpy as np
import time
import configparser
import math
from pyparsing import originalTextFor
import sympy
import scipy 
from plot import Figure
import json
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import emoji
from matplotlib.ticker import LinearLocator, FixedLocator, AutoMinorLocator, MaxNLocator
from sklearn.neighbors import KernelDensity


config = {
    "font.family": ["serif", "Times New Roman"],
    "mathtext.fontset": 'stix',
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)


class Hist1d(Figure):
    def __init__(self):
        super().__init__()
        self.ternary = {}
        self.mode = None

    def load(self):
        if self.cf.has_option("PLOT_CONFI", "save format"):
            self.format = list(map(str.strip, self.cf.get(
                "PLOT_CONFI", 'save format').split(",")))
        if self.cf.has_option("PLOT_CONFI", "result file"):
            resf = self.cf.get("PLOT_CONFI", "result file").split("\n")
            self.data = {}
            for line in resf:
                if len(line.split(",")) != 2:
                    print("Illegal format of result file")
                    sys.exit(1)
                line = list(map(str.strip, line.split(",")))
                # line[1] = os.path.join(self.path['path'], line[1])
                line[1] = self.decode_path(line[1])
                self.data[line[0]] = line[1]


    def load_json(self):
        if "config file" in self.inf:
            self.inf['config file'] = self.decode_path(self.inf['config file'])
            if not os.path.exists(self.inf['config file']):
                print("No config file found -> {}".format(self.inf['config file']))
                sys.exit(1)
            with open(self.inf['config file'], "r") as f1:
                self.inf.update(json.load(f1))
                # print(self.inf.keys())
                if "name" not in self.inf:
                    self.inf['name'] = self.inf['sect']
        if os.path.exists(self.cs):
            with open(self.cs, 'r') as f1:
                self.cs = json.loads(f1.read())

    def load_variable(self):
        if self.mode == "1D_CSV":
            if not set(['x_variable', 'x_bins']) <= set(self.inf.keys()):
                print("Illegal configure input, please check and retry!\n\t-> No 'x_variable' or 'x_bins' defined !")
                sys.exit(1)
            self.inf['x'] = {
                    "label":    self.inf['x_label'],
                    "lim":      [self.inf['x_bins'][0], self.inf['x_bins'][1]]
            }
            if "sep" in self.inf:
                self.vars['sep'] = {}
                # print(self.vars)
                self.vars['norm'] = None
                nn = 0
                for xx, vv in self.inf['sep'].items():
                    nn += 1
                    res = {
                        "label":    vv['label'],
                        "data":   self.load_var_data(vv),
                        'expr':     self.inf['x_variable'],
                        'sumW':     vv['sumW'],
                        "scnm":     "SEP{}".format(nn % 10)
                    }
                    res['data']['vals'] = res["data"].eval(res['expr'])
                    res['sumW_MC'] = res['data']['weight'].sum()
                    res['N_MC']  = res['data'].shape[0]
                    res['XSect']  = res['N_MC'] / res['sumW_MC']
                    res['hdat']  = {}
                    res['hdat']['h_MC'], res['hdat']['bin_edge'] = np.histogram(
                        res['data']['vals'],
                        range=self.inf['x_bins'][0:2],
                        bins=self.inf['x_bins'][-1],
                        weights=res['data']['weight']
                    )
                    res['hdat']['hist']  = res['hdat']['h_MC'] / res['sumW_MC'] * res['sumW'] * res['XSect']
                    res['hdat']['hN_MC'] = res['hdat']['h_MC'] / res['sumW_MC'] * res['N_MC']
                    res['hdat']['hErr_MC'] = np.sqrt(res['hdat']['hN_MC']) / res['hdat']['hN_MC']
                    if self.cs['default']['kde']:
                        res['kde'] = {
                            "vals": res['hdat']['hist'],
                            "bins": (res['hdat']['bin_edge'][0:-1] + res['hdat']['bin_edge'][1:]) / 2.0
                        }
    
                        res['kde']['expr'] = scipy.stats.gaussian_kde(
                            res['data']['vals'], 
                            weights=res['data']['weight'],
                            bw_method=0.005 
                        )

                    self.vars['sep'][xx] = res
                    if self.vars['norm'] is None:
                        self.vars['norm'] = xx 


    def load_var_data(self, vv):
        if set(['set', 'Norm']) <= set(vv.keys()):
            data = []
            with open(self.decode_path(vv['Norm']), 'r') as f1:
                norm = json.loads(f1.read())
            for lb in vv['set']:
                df = pd.read_csv(self.data[lb])
                df['weight'] = norm[lb] / df.shape[0]
                data.append(df)
                # print(df.shape, norm[lb], "\n", df.iloc[0])
            if len(data) > 1:
                data = pd.concat(data)
            else:
                data = data[0]
            return data




    def make_canvas(self):
        self.fig = plt.figure(**self.cs['CS']["figsize"])
        self.ax  = self.fig.add_axes(**self.cs['CS']["axsize"])
        self.axr = self.fig.add_axes(**self.cs['CS']['axrsize'])



    def draw_plot(self):
        # print(self.vars['sep'].keys())
        # self.axr.plot(
        #     np.array([self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][0:-1], self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
        #     np.array([
        #         self.vars['sep'][self.vars['norm']]['hdat']['hist'] / self.vars['sep'][self.vars['norm']]['hdat']['hist'], 
        #         self.vars['sep'][self.vars['norm']]['hdat']['hist'] / self.vars['sep'][self.vars['norm']]['hdat']['hist']
        #     ]).ravel(order='F'), 
        # )
        # self.axr.fill_between(
        #     np.array([self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][0:-1], self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
        #     np.array([
        #         1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'], 
        #         1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']]).ravel(order='F'), 
        #     np.array([
        #         1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'], 
        #         1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']]).ravel(order='F'), 
        #     alpha=0.2
        # )
        self.ax.plot(
            np.array([self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][0:-1], self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
            np.array([self.vars['sep'][self.vars['norm']]['hdat']['hist'], self.vars['sep'][self.vars['norm']]['hdat']['hist']]).ravel(order='F'), 
            "-",
            label=self.vars['sep'][self.vars['norm']]['label'],
            **self.cs['sep'][self.vars['sep'][self.vars['norm']]['scnm']]['cline']
        )
        if self.cs['default']['kde']:
            xx = self.vars['sep'][self.vars['norm']]['kde']['bins']
            yy = self.vars['sep'][self.vars['norm']]['kde']['expr'](xx)
            ff = scipy.interpolate.interp1d(xx, yy, kind='cubic')
            yy = yy / sum(yy) * self.vars['sep'][self.vars['norm']]['sumW'] * self.vars['sep'][self.vars['norm']]['XSect']
            self.ax.plot(xx, yy, "-",  **self.cs['css']['NormLine']['edgecss'])

        # yyp = self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'])
        # yym = self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'])
        self.axr.plot(
            self.inf['x']['lim'], [1.0, 1.0], "-", **self.cs['sep'][self.vars['sep'][self.vars['norm']]['scnm']]['cline']
        )
        self.axr.fill_between(
            np.array([self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][0:-1], self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
            np.array([
                1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'], 
                1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']]).ravel("F"), 
            np.array([
                1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'], 
                1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']]).ravel("F"), 
            **self.cs['sep'][self.vars['sep'][self.vars['norm']]['scnm']]['err']
        )
        self.ax.fill_between(
            np.array([self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][0:-1], self.vars['sep'][self.vars['norm']]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
            np.array([
                self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']), 
                self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 + self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'])]).ravel(order='F'), 
            np.array([
                self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC']), 
                self.vars['sep'][self.vars['norm']]['hdat']['hist'] * (1.0 - self.vars['sep'][self.vars['norm']]['hdat']['hErr_MC'])]).ravel(order='F'), 
            **self.cs['sep'][self.vars['sep'][self.vars['norm']]['scnm']]['err']
        )


        for var in self.vars['sep']:
            if var != self.vars['norm']:

                self.ax.plot(
                    np.array([self.vars['sep'][var]['hdat']['bin_edge'][0:-1], self.vars['sep'][var]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
                    np.array([self.vars['sep'][var]['hdat']['hist'], self.vars['sep'][var]['hdat']['hist']]).ravel(order='F'), 
                    "-", 
                    label=self.vars['sep'][var]['label'], 
                    **self.cs['sep'][self.vars['sep'][var]['scnm']]['cline']
                )
                self.ax.fill_between(
                    np.array([self.vars['sep'][var]['hdat']['bin_edge'][0:-1], self.vars['sep'][var]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
                    np.array([
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 + self.vars['sep'][var]['hdat']['hErr_MC']), 
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 + self.vars['sep'][var]['hdat']['hErr_MC'])]).ravel(order='F'), 
                    np.array([
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 - self.vars['sep'][var]['hdat']['hErr_MC']), 
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 - self.vars['sep'][var]['hdat']['hErr_MC'])]).ravel(order='F'), 
                    **self.cs['sep'][self.vars['sep'][var]['scnm']]['err']
                )
                self.axr.plot(
                    np.array([
                        self.vars['sep'][var]['hdat']['bin_edge'][0:-1], 
                        self.vars['sep'][var]['hdat']['bin_edge'][1:]
                    ]).ravel(order="F"), 
                    np.array([
                        self.vars['sep'][var]['hdat']['hist'] / self.vars['sep'][self.vars['norm']]['hdat']['hist'], 
                        self.vars['sep'][var]['hdat']['hist'] / self.vars['sep'][self.vars['norm']]['hdat']['hist']
                    ]).ravel(order='F'), 
                    "-",
                    **self.cs['sep'][self.vars['sep'][var]['scnm']]['cline']
                )
                self.axr.fill_between(
                    np.array([self.vars['sep'][var]['hdat']['bin_edge'][0:-1], self.vars['sep'][var]['hdat']['bin_edge'][1:]]).ravel(order="F"), 
                    np.array([
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 + self.vars['sep'][var]['hdat']['hErr_MC']) / self.vars['sep'][self.vars['norm']]['hdat']['hist'], 
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 + self.vars['sep'][var]['hdat']['hErr_MC']) / self.vars['sep'][self.vars['norm']]['hdat']['hist']
                    ]).ravel(order='F'), 
                    np.array([
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 - self.vars['sep'][var]['hdat']['hErr_MC']) / self.vars['sep'][self.vars['norm']]['hdat']['hist'], 
                        self.vars['sep'][var]['hdat']['hist'] * (1.0 - self.vars['sep'][var]['hdat']['hErr_MC']) / self.vars['sep'][self.vars['norm']]['hdat']['hist']
                    ]).ravel(order='F'), 
                    **self.cs['sep'][self.vars['sep'][var]['scnm']]['err']
                )                
        
        self.ax.set_xticklabels([])
        self.draw_axis()
        self.axr.set_ylim(1.0 - self.cs['default']['PRatio_dY'], 1.0 + self.cs['default']['PRatio_dY'])
        self.ax.tick_params(**self.cs['css']['axtick']['major'])
        self.ax.tick_params(**self.cs['css']['axtick']['minor'])
        self.ax.tick_params(**self.cs['css']['axtick']['both'])
        self.axr.tick_params(**self.cs['css']['axtick']['major'])
        self.axr.tick_params(**self.cs['css']['axtick']['minor'])
        self.axr.tick_params(**self.cs['css']['axtick']['both'])
        self.ax.legend(fontsize=16)

    def draw_axis(self):
        self.ax.set_xlim(self.inf['x']['lim'])
        self.axr.set_xlim(self.inf['x']['lim'])
        if "y_lim" in self.inf.keys():
            self.ax.set_ylim(self.inf['y_lim'])
        if "y_label" in self.inf.keys():
            self.ax.set_xlabel(r"{}".format(self.inf['y_label']), **self.cs['css']['YLabel'])
        from matplotlib.ticker import AutoMinorLocator, MultipleLocator, MaxNLocator, AutoLocator
        if "x_tick" in self.inf.keys():
            self.ax.xaxis.set_minor_locator(AutoMinorLocator())
            self.axr.xaxis.set_minor_locator(AutoMinorLocator())
        else: 
            self.ax.xaxis.set_minor_locator(AutoMinorLocator())
            self.axr.xaxis.set_minor_locator(AutoMinorLocator())
        
        self.ax.yaxis.set_major_locator(MaxNLocator(6))
        self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.axr.yaxis.set_minor_locator(AutoMinorLocator())

        self.axr.set_xlabel(r"{}".format(self.inf['x_label']), **self.cs['css']['XLabel'])
        self.axr.xaxis.set_label_coords(**self.cs['css']['xlabelcoords'])


    def drawpicture(self):
        print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
            time.time()-self.time, self.inf['name']), language="alias"))
        self.load()
        self.load_variable()
        self.make_canvas()
        self.draw_plot()
        self.print_figure(plt)
