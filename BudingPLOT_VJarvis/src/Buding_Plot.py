#!/usr/bin/env python3

import pandas as pd
import re
import emoji
import subprocess
import os
import sys
import json

pwd = os.path.abspath(os.path.dirname(__file__))

import time
import configparser
import math
from scipy.interpolate import Rbf
import copy
from matplotlib.patches import FancyBboxPatch
from matplotlib.patches import Rectangle
from scipy.interpolate import interp1d
import sympy

from matplotlib import rc, rcParams
import matplotlib.pyplot as plt
from matplotlib.image import NonUniformImage
from matplotlib import cm, ticker
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib
# fm = matplotlib.font_manager.json_load("{}/font/fontlist.json".format(pwd))
# fm.findfont("serif", rebuild_if_missing=False)
# print(rcParams)
# rc('text', usetex=True)
config = {
    "font.family":["serif", "Times New Roman"],
    "font.size": 20,
    "mathtext.fontset":'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

class Figure():
    def __init__(self):
        pass

    def load_path(self, path, df="path"):
        df = self.cf.get("PLOT_CONFI", df)
        if "~" in df:
            from os.path import expanduser
            df = df.replace("~", expanduser("~"))
        if "&BP" in path:
            path = path.replace("&BP", pwd)
            return path
        else:
            path = os.path.abspath(os.path.join(
                df, path))
            # print(os.path.join(self.cf.get("PLOT_CONFI", df), path))
            return path

    def get_inf(self, inf):
        self.cf = configparser.ConfigParser()
        self.cf.read(inf)
        if self.cf.get("PLOT_CONFI", "path").strip().upper() == "&PWD":
            self.cf.set("PLOT_CONFI", "path", os.path.dirname(inf))

    def load_data(self):
        self.data = []
        for opt in self.cf.get("PLOT_CONFI", 'result_file').split('\n'):
            dat = pd.read_csv(self.load_path(opt.split()[1]))
            # dat = dat.drop(["Index"], axis=1)
            self.data.append(dat.dropna(axis=0, how='any'))
        self.data = pd.concat(self.data, axis=0, join='outer',
                              ignore_index=True, sort=False)

    def load_colormap(self, figset, cmname):
        # print(figset['colormappath'])
        cmapjs = self.load_path(figset['colormappath'])
        with open(cmapjs, 'r') as f1:
            clist = json.loads(f1.read())
        if cmname in clist.keys():
            cmap = clist[cmname]
            if cmap['type'] == "self":
                from matplotlib.colors import LinearSegmentedColormap
                cmap['cmap'] = LinearSegmentedColormap.from_list(
                    cmname, tuple(cmap['code']), cmap['nbin'])
            elif cmap['type'] == "inner":
                cmap['cmap'] = plt.get_cmap(cmap['code'])
        else:
            cmap = {}
            cmap['name'] = cmname
            cmap['cmap'] = plt.get_cmap(cmname)
        return cmap

    def plot(self):
        self.figures_inf()
        self.Analytic_funcs()
        for fig in self.figs:
            if fig['type'] == "SLHA":
                self.plotSLHA(fig)
            elif fig['type'] == "Jarvis":
                self.plotJarvis(fig)
            elif fig['type'] == "RivetYoda":
                self.plotRivet(fig)
            else:
                self.load_fig_setting(fig)
                self.load_data()
                self.drawpicture(fig)

    def plotRivet(self, fig):
        print("=== Ploting Figure : {} ===".format(fig['name']))
        fig['start'] = time.time()
        self.load_fig_setting(fig)
        self.loadRivetInfo(fig)
        from copy import deepcopy
        fig['axis'] = {}
        for pname in fig['ax'].keys():
            print(fig['ps']['setting'][pname])
            if fig['ax'][pname]['type'] == "YODA_HISTO1D_V2":    
                fig['fig'] = plt.figure(
                    figsize=(10, 6)
                )
                ax = fig['fig'].add_axes([0.13, 0.16, 0.81, 0.8])
                cud = []
                cmd = None
                fig['axis']['xlim'] = []
                xcud = []
                xcue = []
                for sub in fig['ax'][pname]['Cum']:
                    if cmd is not None:
                        tmd = deepcopy(sub['data'])
                        tmd['cum'] = cmd
                        # tmd['sumw'] = cmd + tmd['sumw']
                        cmd += tmd['sumw']
                        cud.append(tmd)
                    else:
                        tmd = deepcopy(sub['data'])
                        if not fig['axis']['xlim']:
                            fig['axis']['xlim'] = [float(tmd['xlow'].min()), float(tmd['xhigh'].max())]
                        tmd['cum'] = 0.
                        cmd = tmd['sumw']
                        xcud = tmd['xlow']
                        xcue = tmd['xhigh']
                        cud.append(tmd)
                    # break
                n = 1

                for sub in fig['ax'][pname]['Sep']:
                    tmd = deepcopy(sub['data'])

                    # ax.plot([tmd['xlow'].min(), tmd['xlow'].min()], [0., yy[0]], c="red", linewidth=3, solid_joinstyle="miter")
                    ax.step(tmd['xlow'], tmd['sumw'], c="red", linewidth=3, where='post', solid_joinstyle="miter")
                    ax.step(tmd['xhigh'], tmd['sumw'], c="red", linewidth=3, where='pre', solid_joinstyle="miter")

                xcud = list(xcud)
                cmd  = list(cmd)

                ax.step(xcud, cmd, color='black', where='post')
                ax.step(xcue, cmd, color='black', where='pre')

                prop_cycle = plt.rcParams['axes.prop_cycle']
                color = prop_cycle.by_key()['color']

                xbin = None
                for ds in cud:
                    n += 1
                    for row in ds.iterrows():
                        row = row[1]
                        if xbin == None:
                            xbin = row['xhigh'] - row['xlow']
                        ax.fill(
                            [row['xlow'], row['xlow'], row['xhigh'], row['xhigh'], row['xlow']], 
                            [row['cum'], row['cum'] + row['sumw'],  row['cum'] + row['sumw'], row['cum'], row['cum']], 
                            facecolor=color[n],
                            edgecolor=None
                        ) 
                if "LogY" in fig['ps']['setting'][pname].keys() and eval(fig['ps']['setting'][pname]['LogY']):
                    ax.set_yscale("log")
                elif fig['colorset']['HISTO1D']['LogY']:
                    ax.set_yscale("log")
                ylim = [fig['colorset']['HISTO1D']['YMin'], fig['colorset']['HISTO1D']['YMax']] 
                if "YMin" in fig['ps']['setting'][pname].keys():
                    ylim[0] = eval(fig['ps']['setting'][pname]['YMin'])
                if "YMax" in fig['ps']['setting'][pname].keys():
                    ylim[1] = eval(fig['ps']['setting'][pname]['YMax'])
                ax.set_ylim(ylim[0], ylim[1])
                ax.set_xlim(fig['axis']['xlim'][0], fig['axis']['xlim'][1])
                
                from matplotlib.ticker import AutoMinorLocator, FixedLocator, StrMethodFormatter, MaxNLocator, MultipleLocator
                if "XPi" in fig['ps']['setting'][pname].keys():
                    print(">>>>Tagng", fig['ps']['setting'][pname]['XPi'])
                    if eval(fig['ps']['setting'][pname]['XPi']):
                        xbin = eval(fig['ps']['setting'][pname]['XPi'])
                        ax.xaxis.set_major_locator(MultipleLocator(np.pi * xbin * 5))
                        ax.xaxis.set_minor_locator(MultipleLocator(np.pi * xbin))
                        x_tick = np.arange(fig['axis']['xlim'][0], fig['axis']['xlim'][1], np.pi*xbin*5)
                        x_ratio = np.arange(fig['axis']['xlim'][0]/np.pi, fig['axis']['xlim'][1]/np.pi, xbin*5 )
                        print(x_ratio)
                        print((xbin*5).as_integer_ratio())
                        x_label = []
                        for xx in x_ratio:
                            if xx != 0.:
                                c = list(xx.as_integer_ratio())
                                if c[0] == 1:
                                    if c[1] != 1:
                                        x_label.append(r"$\frac{\pi}{" + str(c[1])+ r"}$")
                                    else:
                                        x_label.append(r"$\pi$")
                                elif c[1] == 1:    
                                    x_label.append(r"$" + c[0] + "\pi$")
                                else:    
                                    x_label.append(r"$\frac{" + str(c[0]) + r"\pi}{" + str(c[1])+ r"}$")
                            else:
                                x_label.append(r'0')
                        print(x_label)
                        ax.set_xticks(x_tick)
                        ax.set_xticklabels(x_label)
                        # ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
                    else:
                        ax.xaxis.set_major_locator(MultipleLocator(5*xbin))
                        ax.xaxis.set_minor_locator(MultipleLocator(xbin))
                else:
                    ax.xaxis.set_major_locator(MultipleLocator(5*xbin))
                    ax.xaxis.set_minor_locator(MultipleLocator(xbin))   
                    
                yminor = matplotlib.ticker.LogLocator(base=10.0, subs = np.arange(1.0, 10.0)*0.1, numticks=10)
                ax.yaxis.set_minor_locator(yminor)
                ymajor = matplotlib.ticker.LogLocator(base=10.0, numticks=(math.log10(ylim[1])-math.log10(ylim[0])))
                ax.yaxis.set_major_locator(ymajor)
                
                ax.tick_params(labelsize=fig['colorset']['HISTO1D']['ticks']['labelsize'],  direction=fig['colorset']['HISTO1D']['ticks']['direction'],  bottom=fig['colorset']['HISTO1D']['ticks']
                               ['bottom'],  left=fig['colorset']['HISTO1D']['ticks']['left'],  top=fig['colorset']['HISTO1D']['ticks']['top'],  right=fig['colorset']['HISTO1D']['ticks']['right'], which='both')
                ax.tick_params(which='major', length=fig['colorset']['HISTO1D']['ticks']
                               ['majorlength'], color=fig['colorset']['HISTO1D']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['HISTO1D']['ticks']
                               ['minorlength'], color=fig['colorset']['HISTO1D']['ticks']['minorcolor'])
                
                if "XLabel" in fig['ps']['setting'][pname].keys():
                    ax.set_xlabel(
                        r"{}".format(fig['ps']['setting'][pname]['XLabel']),
                        loc='right', 
                        fontsize=fig['colorset']['canvas']['x_label']['size'],
                        color=fig['colorset']['canvas']['x_label']['color']                   
                    )
                if "YLabel" in fig['ps']['setting'][pname].keys():
                    ax.set_ylabel(
                        r"{}".format(fig['ps']['setting'][pname]['YLabel']),
                        loc='top', 
                        fontsize=fig['colorset']['canvas']['yl_label']['size'],
                        color=fig['colorset']['canvas']['yl_label']['color'],               
                    )
                # ax.set_xlim(-1, fig['ax']['xlim'][1])
                if 'save' in self.cf.get(fig['section'], 'print_mode'):
                    fig['fig'] = plt
                    if "/" == pname[0]:
                        plname = pname[1:]
                    fig['file'] = os.path.join(self.figpath, plname)
                    fig['name'] = os.path.basename(plname)
                    if not os.path.exists(os.path.dirname(fig['file'])):
                        print(os.path.dirname(fig['file']))
                        os.makedirs(os.path.dirname(fig['file']))
                    # print(self.figpath, fig["file"], os.path.basename(pname), os.path.dirname(pname))
                    self.savefig(fig, plt)
                    plt.close()
                if 'show' in self.cf.get(fig['section'], 'print_mode'):
                    plt.show(block=False)
                    plt.pause(1)
                    input("\n\tPress 'Enter' to continue ...")
                    plt.close()
                # break
                
            
            
        
        
    def loadRivetInfo(self, fig):
        fig['ps'] = {}
        if self.cf.has_option(fig['section'], "plot_setting"):
            fig['ps']['cont'] = self.cf.get(fig['section'], "plot_setting")
            from Func_lab import load_rivet_plot_info
            with open(self.load_path(fig['ps']['cont'], "result_path"), 'r') as f1:
                fig['ps']['setting'] = load_rivet_plot_info(f1.read())
        else:
            fig['ps'] = None 
        from Func_lab import load_rivet_data_info
        from Func_lab import load_yoda_data
        cum, sep = load_rivet_data_info(self.cf.get(fig['section'], "data"))
        for yd in cum:
            yd['yoda'] = self.load_path(yd['yoda'], "result_path")
            yd['data'] = load_yoda_data(yd['yoda'])
        for yd in sep:
            yd['yoda'] = self.load_path(yd['yoda'], "result_path")
            yd['data'] = load_yoda_data(yd['yoda'])
        fig['data'] = {
            "cum":  cum,
            "sep":  sep
        }

        fig['ax'] = {}

        nd = len(fig['data']['cum'])
        np = len(fig['data']['cum'][0]['data'])
        
        for jj in range(np):
            for ii in range(nd):
                pname = fig['data']['cum'][ii]['data'][jj]['name']
                if pname not in fig['ax'].keys():
                    fig['ax'][pname] = {
                        "Cum": [fig['data']['cum'][ii]['data'][jj]],
                        "type": fig['data']['cum'][ii]['data'][jj]['type']
                        }
                else:
                    fig['ax'][pname]['Cum'].append(fig['data']['cum'][ii]['data'][jj])
        
        nd = len(fig['data']['sep'])
        for yd in fig['data']['sep']:
            np = len(yd['data'])
            for jj in range(np):
                pname = yd['data'][jj]['name']
                if pname not in fig['ax'].keys():
                    fig['ax'][pname] = {
                        "Sep": [yd['data'][jj]],
                        "type":     yd['data'][jj]['type']
                        }
                else:
                    if "Sep" not in fig['ax'][pname].keys():
                        fig['ax'][pname]['Sep']=[yd['data'][jj]]
                    else:
                        fig['ax'][pname]['Sep'].append(yd['data'][jj])


    def load_fig_setting(self, fig):
        with open(self.load_path(self.cf.get("COLORMAP", "colorsetting")), 'r') as settfile:
            self.colors = json.loads(settfile.read())
        for item in self.colors:
            if item['name'] in self.cf.get(fig['section'], "colorset").replace("&", ""):
                fig['colorset'] = item
        with open(self.load_path(fig['colorset']['figureSetting']), 'r') as f1:
            fig['colorset'].update(json.loads(f1.read()))
        # print(fig['colorset'])
        for key in fig['colorset'].keys():
            if "path" in key:
                fig['colorset'][key] = self.load_path(fig['colorset'][key])
        if "colormap" in fig['colorset'].keys():
            fig['colorset']['cmap'] = self.load_colormap(
                fig['colorset'], fig['colorset']['colormap'])

    def plotJarvis(self, fig):
        print("=== Ploting Figure : {} ===".format(fig['name']))
        fig['start'] = time.time()
        self.load_fig_setting(fig)
        self.load_Jarvis_flow(fig)
        self.plot_Jarvis_flow(fig)
        
        
    def load_Jarvis_flow(self, fig):
        fig['layers'] = {}
        with open(self.cf.get(fig['section'], 'info_path'), 'r') as f1:
            dat = json.loads(f1.read())
        # print(dat.keys())
        fs = fig['colorset']['figsize']
        n = 1
        for ii in range(len(dat['plot']['layer'])):
            nn_currentlayer = len(dat['plot']['layer']) - ii - 1
            layer = dat['plot']['layer'][nn_currentlayer]
            # print("info", ii, nn_currentlayer, layer)
            fig['layers'][nn_currentlayer] = {
                "nodes":    layer,
                "var_nextlayer":    {}
            }
            
            for node in layer:
                hfinp   = len(dat['plot']['nodes'][node]['input file']) * fs['h3']
                hsinp   = (len(dat['plot']['nodes'][node]['input file']) - 1) * fs['h6']
                hvinp   = 0. 
                for ff in dat['plot']['nodes'][node]['input file']:                    
                    hvinp += len(ff['vars']) * fs['hv']

                hfoup   = len(dat['plot']['nodes'][node]['output file']) * fs['h3']
                hsoup   = (len(dat['plot']['nodes'][node]['input file']) - 1) * fs['h6']
                hvoup   = 0. 
                for ff in dat['plot']['nodes'][node]['output file']:
                    for var in ff['vars']:
                        n_nextlayer = self.find_var_nextlayer_number(var, nn_currentlayer, dat)
                        # print(var, n_nextlayer)
                        nnl = n_nextlayer['nNL'] + n_nextlayer['nNNL']
                        if nnl == 0:
                            hvoup += fs['hv']
                        else:
                            hvoup += nnl * fs['hv']
        
                hinp    = hfinp + hsinp + hvinp
                houp    = hfoup + hsoup + hvoup
                htxt    = fs['wt'] * (len(node)+1.5)

                fig['layers'][nn_currentlayer][node] = {
                    "hinp":    hinp,
                    "houp":    houp
                }
                
                locs    = {}
                hmax = max(hinp, houp, htxt)
                y0 = -(hmax - hinp)/2.0 
                y1 = -(hmax - houp)/2.0 
                    
                for ff in dat['plot']['nodes'][node]['input file']:                    
                    locs[ff['name']] = {
                        "type": "inpfile",
                        "x":    0., 
                        "y":    y0 - fs['h3'],
                        "w":    fs['wv'] + fs['wns'], 
                        "h":    fs['h3']
                    }
                    y0 -= fs['h3']
                    endvv = ""
                    for vv, pp in ff['vars'].items():
                        locs[vv] = {
                            "type": "inpvar",
                            "posi": "in",
                            "meth": pp,
                            "x":    fs['wsfv'],       
                            "y":    y0 - fs['hv'],
                            'w':    fs['wv'] - fs['wsfv'],
                            'h':    fs['hv']
                        }
                        endvv = vv 
                        y0 -= fs['hv']
                    y0 -= fs['h6']
                    if endvv != "":
                        locs[endvv]['posi'] = "side"
                
                for ff in dat['plot']['nodes'][node]['output file']:
                    locs[ff['name']] = {
                        "type": "outfile",
                        "x":    fs['wv'] + fs['wns'] + fs['wn'],
                        "y":    y1 - fs['h3'],
                        "w":    fs['wv'] + fs['wns'],
                        "h":    fs['h3']
                    }
                    y1 -= fs['h3']
                    endvv = ""
                    for vv, pp in ff['vars'].items():
                        n_nextlayer = self.find_var_nextlayer_number(vv, nn_currentlayer, dat)
                        nnl = n_nextlayer['nNL'] + n_nextlayer['nNNL']
                        if nnl == 0:
                            nnl = 1 
                        locs[vv] = {
                            "type": "oupvar",
                            "posi": "in",
                            "meth": pp,  
                            'x':    fs['wv'] + 2*fs['wns'] + fs['wn'],
                            'y':    y1 - nnl * fs['hv'],
                            'w':    fs['wv'] - fs['wsfv'],
                            "h":    nnl * fs['hv']
                        }
                        endvv = vv
                        y1 -= nnl * fs['hv']
                    y1 -= fs['h6']
                    locs[endvv]['posi'] = "side"
                            
                locs[node] = {
                    "type": "node",
                    "x":    fs['wv'] + fs['wns'],
                    "y":    -hmax, 
                    "w":    fs['wn'],
                    "h":    hmax
                }
                            
                            
                fig['layers'][nn_currentlayer][node]['locs'] = locs 
                fig['layers'][nn_currentlayer][node]['height'] = hmax 
                                
        nn_currentlayer = -1 
        fig['layers'][nn_currentlayer] = {
            "nodes":    [],
            "var_nextlayer":    {},
            "scan variables":   {}
        }
        y0 = 0. 
        y0 -= fs['h2']
        
        locs    = {
            "ID":   {
                "type": "outvar",
                'meth': "UUID(4)",
                "x":    0.,
                "y":    y0 - fs['hv'],
                'w':    fs['wv'],
                'h':    fs['hv']
            },
            "Status":   {
                "type": "outvar",
                "meth": "free",
                "x":    0.,
                "y":    y0 - 2 * fs['hv'],
                "w":    fs['wv'],
                "h":    fs['hv']
            }
        }
        y0 -= 2* fs['hv']
        for var in dat['scan variables']:
            n_nextlayer = self.find_var_nextlayer_number(var['name'], nn_currentlayer, dat)
            nnl = n_nextlayer['nNL'] + n_nextlayer['nNNL']
            locs[var['name']] = {
                "type": "outvar",
                "meth": var['prior'],
                "x":    0., 
                "y":    y0 - nnl * fs['hv'],
                'w':    fs['wv'],
                'h':    nnl * fs['hv'],
            }
            y0 -= nnl * fs['hv']
            if n_nextlayer['nNNL'] > 0:
                fig['layers'][nn_currentlayer + 1]['var_nextlayer'][var['name']] = {"nhv":    n_nextlayer['nNNL']}
        
        fig['layers'][nn_currentlayer]['scan variables']['locs'] = locs 
        fig['layers'][nn_currentlayer]['scan variables']['height'] = -y0  
        fig['layers'][nn_currentlayer]['height'] = -y0 
        
        for nn_currentlayer in range(len(dat['plot']['layer'])):
            layer = dat['plot']['layer'][nn_currentlayer]
            for node in layer: 
                for ff in dat['plot']['nodes'][node]['output file']:
                    for var in ff['vars']:
                        n_nextlayer = self.find_var_nextlayer_number(var, nn_currentlayer, dat)
                        if n_nextlayer['nNNL'] > 0:
                            fig['layers'][nn_currentlayer + 1]['var_nextlayer'][var] = {"nhv":  n_nextlayer['nNNL']}

        nn = -1
        for ii in range(len(fig['layers'])):
            kk = nn 
            nn += 1
            vv = fig['layers'][kk]
            if kk != -1:
                for var in vv['var_nextlayer']:
                    n_nextlayer = self.find_var_nextlayer_number(var, kk, dat)
                    if n_nextlayer['nNNL'] > 0:
                        fig['layers'][kk+1]['var_nextlayer'][var] = {"nhv": n_nextlayer['nNNL']}

        for kk, vv in fig['layers'].items():
            if kk != -1:
                y0      = 0. 
                for var in vv['var_nextlayer']:
                    vv['var_nextlayer'][var].update({
                        "x":    fs['wv'] + (kk + 1)*fs['ws'] + kk*(2*fs['wv'] + 2*fs['wns'] + fs['wn']),
                        'y':    y0 - vv['var_nextlayer'][var]['nhv'] * fs['hv'],
                        "w":    fs['wn'],
                        'h':    vv['var_nextlayer'][var]['nhv'] * fs['hv'] 
                    })
                    y0 -= vv['var_nextlayer'][var]['nhv'] * fs['hv']
                    y0 -= fs['h5']
                
                if len(vv['var_nextlayer']) > 0:
                    y0 += fs['h5']
                    y0 -= fs['h4']
                
                # fig['layers']['height'] = 
                for node in vv['nodes']:
                    vv[node]['base'] = {
                        "x0":   fs['wv'] + (kk + 1)*fs['ws'] + kk*(2*fs['wv'] + 2*fs['wns'] + fs['wn']),
                        "y0":   y0
                    }
                    y0 -= vv[node]['height'] 
                    y0 -= fs['h5']
                y0 += fs['h5']
                fig['layers'][kk]['height'] = -y0 
        fig['size'] = {
            "h":    0.,
            "w":    (len(fig['layers']) - 1) * (fs['ws'] + 2*fs['wv'] + 2*fs['wns'] + fs['wn'] ) + fs['wv']
        }
        for kk, vv in fig['layers'].items():
            if fig['layers'][kk]['height'] > fig['size']['h']:
                fig['size']['h'] = fig['layers'][kk]['height']
            
    def plot_Jarvis_flow(self, fig):
        fs = fig['colorset']['figsize']

        fig['fig'] = plt.figure(
            figsize=(
                fig['size']['w'] + 2*fs['w0'], 
                fig['size']['h'] + 2*fs['h0'])
        )
        fig['fig'].text(0., 0., "Test", color="None")
        fig['fig'].text(1., 1., "Test", color="None")
        fig['ax'] = fig['fig'].add_axes([0., 0., 1., 1.])
        fig['ax'].axis("off")
        fig['ax'].set_xlim(-fs['w0'], fig['size']['w'] + fs['w0'])
        fig['ax'].set_ylim(-fig['size']['h'] - fs['h0'], fs['h0'])
        
        self.draw_scan_variable(fig)
        nly = len(fig['layers']) - 1 
        for ii in range(nly):
            self.draw_layer(fig, fig['layers'][ii], ii)
        
        if 'save' in self.cf.get(fig['section'], 'print_mode'):
            fig['fig'] = plt
            fig['file'] = os.path.join(self.figpath, fig['name'])
            self.savefig(fig, plt)
        if 'show' in self.cf.get(fig['section'], 'print_mode'):
            plt.show(block=False)
            plt.pause(1)
            input("\n\tPress 'Enter' to continue ...")
            plt.close()        
          
    def get_jarvis_var_source(self, fig, vname, kk, h):
        fs = fig['colorset']['figsize']
        if kk == 0:
            out = fig['layers'][kk-1]['scan variables']['locs']
            res = {
                "x":    out[vname]['x'] + fs['wvc'],
                "y":    out[vname]['y'] + out[vname]['h'] - h,
                "c":    out[vname]['c'] 
            }
            out[vname].update({
                "h":    out[vname]['h'] - h
            })
            return res 
        else:
            out = None
            if vname in fig['layers'][kk-1]['var_nextlayer'].keys():
                outv = fig['layers'][kk-1]['var_nextlayer'][vname]
                res = {
                    "x":    outv['x'] + fs['wn'],
                    "y":    outv['y'] + outv['h'] - h,
                    'c':    outv['c']
                }
                outv.update({
                    "h":    outv['h'] - h 
                })
                return res 
            else:
                for node in fig['layers'][kk-1]['nodes']:
                    if vname in fig['layers'][kk-1][node]['locs'].keys():
                        outv = fig['layers'][kk-1][node]['locs'][vname]
                        res = {
                            "x":    outv['x'] + fs['wvc'],
                            "y":    outv['y'] + outv['h'] - h,
                            "c":    outv['c']
                        }
                        outv.update({
                            "h":    outv['h'] - h
                        })
                        break 
                return res 
            
    def draw_one_flow(self, fig, s, e):
        x = e['x'] - s['x']
        y = e['y'] - s['y']
        omega = np.pi / x 
        mag = - y/2. 
        tt = np.linspace(0., 1., 100)
        xx = tt * x + s['x']
        yy0 = mag * np.cos(tt*np.pi) + s['y'] + 0.5*y 
        yy1 = mag * np.cos(tt*np.pi) + s['y'] + 0.5*y + e['h']
        
        fig['ax'].fill_between(
            xx, yy0, yy1, 
            edgecolor=None, facecolor=e['c'],
            alpha=0.3
        )
        
    def draw_layer(self, fig, vv, kk):
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']

        y0 = - 0.5 * (fig['size']['h'] - vv['height'])
        # print(kk, vv['var_nextlayer'].keys())
        for nn, va in vv['var_nextlayer'].items():
            va.update({"x": va['x'] + fs['wv'] + fs['wns'], "y": va['y'] + y0})
            sinf = self.get_jarvis_var_source(fig, nn, kk, va['h'])
            va.update({"c": sinf['c']})
            self.draw_one_flow(fig, sinf, va)
            fig['ax'].fill(
                [va['x'], va['x'] + va['w'], va['x'] + va['w'], va['x'], va['x']],
                [va['y'], va['y'], va['y'] + va['h'], va['y'] + va['h'], va['y']],
                alpha=0.6,
                facecolor=va['c'],
                edgecolor=None
            ) 
            fl = fig['colorset']['nodelabel']
            txt = fig['ax'].text(
                va['x'] + 0.5*va['w'],
                va['y'] + 0.5*va['h'],
                nn, 
                fontsize=fl['fontsize'],
                color='w',
                rotation=90,
                fontname=fl['fontname'],
                fontweight=fl['fontweight'],
                fontstyle=fl['fontstyle'],
                va = 'center',
                ha = 'center',
                wrap = True
            )
            import matplotlib.patheffects as PathEffects
            txt.set_path_effects([PathEffects.withStroke(linewidth=1.2, foreground=va['c'])])
        for nn in vv['nodes']:
            va = vv[nn]
            width = 2*fs['wv'] + 2*fs['wns'] + fs['wn']
            # print(va)
            va['base'].update({"y0": va['base']['y0'] + y0})
            # fig['ax'].plot(
            #     [va['base']['x0'], va['base']['x0'] + width, va['base']['x0'] + width, va['base']['x0'], va['base']['x0']],
            #     [va['base']['y0'], va['base']['y0'], va['base']['y0'] - va['height'], va['base']['y0'] - va['height'], va['base']['y0']]
            # )

            for nn, var in va['locs'].items():
                if var['type'] == "inpvar":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    sinf = self.get_jarvis_var_source(fig, nn, kk, var['h'])
                    var.update({"c": sinf['c']})
                    if var['posi'] == "in":
                        self.draw_inpvar_label(fig, var)
                    elif var['posi'] == "side":
                        self.draw_endinpvar_label(fig, var)

                    p = fig['ax'].fill(
                        [var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                         var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                         var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                         var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                         var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb']],
                        [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']],
                        facecolor=var['c']
                    )
                    var['x'] += var['w'] - fs['wvc']
                    self.draw_one_flow(fig, sinf, var)
                    fig['ax'].text(
                        var['x'] + vs['offset'], 
                        var['y'] + 0.5*var['h'], 
                        nn, 
                        fontsize= vs['fontsize'],
                        color=vs['fontcolor'],
                        fontname=vs['namefont'],
                        fontstyle=vs['namestyle'],
                        fontweight=vs['nameweight'],
                        horizontalalignment='left', 
                        verticalalignment='center'            
                    )            
                    fig['ax'].text(
                        var['x'] - vs['offset'], 
                        var['y'] + 0.5*var['h'], 
                        var['meth'], 
                        fontsize= vs['fontsize'],
                        color=vs['fontcolor'],
                        fontname=vs['methfont'],
                        fontstyle=vs['methstyle'],
                        fontweight=vs['methweight'],
                        horizontalalignment='right', 
                        verticalalignment='center'            
                    )
                elif var['type'] == "inpfile":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    self.draw_inpfile(fig, nn, var)
                elif var['type'] == "node":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    p = fig['ax'].fill(
                        [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                        [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                        alpha=0.5
                    )
                    var.update({"c": p[0].get_facecolor()[0:3]})
                    fig['ax'].plot(
                        [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                        [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                        '-',
                        linewidth=0.5,
                        c=var['c']
                    )

                    fl = fig['colorset']['nodelabel']
                    txt = fig['ax'].text(
                        var['x'] + 0.5*var['w'],
                        var['y'] + 0.5*var['h'],
                        nn, 
                        fontsize=fl['fontsize'],
                        color='w',
                        rotation=90,
                        fontname=fl['fontname'],
                        fontweight=fl['fontweight'],
                        fontstyle=fl['fontstyle'],
                        va = 'center',
                        ha = 'center',
                        wrap = True
                    )
                    import matplotlib.patheffects as PathEffects
                    txt.set_path_effects([PathEffects.withStroke(linewidth=1.2, foreground=var['c'])])
                elif var['type'] == "outfile":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    # fig['ax'].plot(
                    #     [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                    #     [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                    #     "-",
                    # )
                    self.draw_oupfile(fig, nn, var)
                elif var['type'] == "oupvar":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    p = fig['ax'].fill(
                        [var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x']  + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb']],
                        [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']]
                    )
                    var.update({"c": p[0].get_facecolor()[0:3]})
                    if var['posi'] == "in":
                        self.draw_outvar_label(fig, var)
                    elif var['posi'] == "side":
                        self.draw_endoutvar_label(fig, var)
                    
                    fig['ax'].text(
                        var['x'] + fs['wvc'] - vs['offset'], 
                        var['y'] + 0.5*var['h'], 
                        nn, 
                        fontsize= vs['fontsize'],
                        color=vs['fontcolor'],
                        fontname=vs['namefont'],
                        fontstyle=vs['namestyle'],
                        fontweight=vs['nameweight'],
                        horizontalalignment='right', 
                        verticalalignment='center'            
                    )            
                    fig['ax'].text(
                        var['x'] + fs['wvc'] + vs['offset'], 
                        var['y'] + 0.5*var['h'], 
                        var['meth'], 
                        fontsize= vs['fontsize'],
                        color=vs['fontcolor'],
                        fontname=vs['methfont'],
                        fontstyle=vs['methstyle'],
                        fontweight=vs['methweight'],
                        horizontalalignment='left', 
                        verticalalignment='center'            
                    )
                                
    def draw_inpfile(self, fig, fname, var):
        fs = fig['colorset']['figsize']
        fl = fig['colorset']['filelabel']
        def draw_logo():
            po = np.array([var['x'], var['y']])
            pa1 = po + np.array([fs['wvcb'], 0.])
            pb1 = pa1 + np.array([0., var['h']])
            pa2 = pa1 + np.array([var['h'], 0. ])
            pb2 = pb1 + np.array([var['h'], 0. ])
            pd1 = po  + np.array([var['w'], 0. ])
            pd2 = pd1 + np.array([0., var['h']])
            pc1 = po + np.array([fs['wvcb'], fs['wvcb']])
            pc2 = pc1 + np.array([0., var['h'] - 2* fs['wvcb']])
            pc4 = po + np.array([0., fs['wvcb']])
            pc3 = pc4 + np.array([0., var['h'] - 2* fs['wvcb']])
            from shapely.geometry import Polygon, Point 
            rt = Polygon([pc1, pc2, pc3, pc4, pc1])
            c1 = Point(pc1).buffer(fs['wvcb'])
            c2 = Point(pc2).buffer(fs['wvcb'])
            from shapely.ops import unary_union
            shp = unary_union([rt, c1, c2])
            rt0 = Polygon([pa1, pa2, pb2, pb1, pa1])
            shp = shp.difference(rt0)
            
            xx, yy = shp.exterior.xy 
            p = fig['ax'].fill(xx, yy, alpha=0.7)
            var.update({"c": p[0].get_facecolor()[0:3]})
            
            rt1 = Polygon([pa1, pd1, pd2, pb1, pa1])
            xx, yy = rt1.exterior.xy 
            fig['ax'].plot(xx, yy, "-", linewidth=0.5, c=var['c'])
            
            pe1 = pd1 + np.array([0., 0.5*var['h']])
            pf1 = pe1 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg1 = pe1 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe2 = pe1 + np.array([-fs['wvcb'], 0.])
            pf2 = pe2 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg2 = pe2 + np.array([-fs['wvcb'], -0.5*var['h']])            
            pe3 = pe2 + np.array([-fs['wvcb'], 0.])
            pf3 = pe3 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg3 = pe3 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe4 = pe3 + np.array([-fs['wvcb'], 0.])
            pf4 = pe4 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg4 = pe4 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe5 = pe4 + np.array([-fs['wvcb'], 0.])
            pf5 = pe5 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg5 = pe5 + np.array([-fs['wvcb'], -0.5*var['h']])
            
            s1 = Polygon([pe1, pd1, pg1, pe1])
            s2 = Polygon([pe1, pd2, pf1, pe1])
            s3 = Polygon([pe2, pf2, pf3, pe3, pg3, pg2, pe2])
            s4 = Polygon([pe4, pf4, pf5, pe5, pg5, pg4, pe4])
            # shp = unary_union([s1, s2, s3])
            xx, yy = s1.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
            xx, yy = s2.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
            xx, yy = s3.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)
            xx, yy = s4.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
                                
        draw_logo()
        ch, cs = len(fname), 16
        fname = [fname[i: i+cs] for i in range(0, ch, cs)]
        fname = "\n".join(fname)        
        
        fig['ax'].text(
            var['x'] + fs['wvcb'] + fl['offset'],
            var['y'] + 0.5*var['h'],
            fname, 
            fontsize=fl['fontsize'],
            color=var['c'],
            fontname=fl['fontname'],
            fontweight=fl['fontweight'],
            fontstyle=fl['fontstyle'],
            va = 'center',
            ha = 'left'
        )
        
    def draw_oupfile(self, fig, fname, var):
        fs = fig['colorset']['figsize']
        fl = fig['colorset']['filelabel']
        def draw_logo():
            po = np.array([var['x']+var['w'], var['y']])
            pa1 = po + np.array([-fs['wvcb'], 0.])
            pb1 = pa1 + np.array([0., var['h']])
            pa2 = pa1 - np.array([var['h'], 0. ])
            pb2 = pb1 - np.array([var['h'], 0. ])
            pd1 = po  - np.array([var['w'], 0. ])
            pd2 = pd1 + np.array([0., var['h']])
            pc1 = po + np.array([-fs['wvcb'], fs['wvcb']])
            pc2 = pc1 + np.array([0., var['h'] - 2* fs['wvcb']])
            pc4 = po + np.array([0., fs['wvcb']])
            pc3 = pc4 + np.array([0., var['h'] - 2* fs['wvcb']])
            from shapely.geometry import Polygon, Point 
            rt = Polygon([pc1, pc2, pc3, pc4, pc1])
            c1 = Point(pc1).buffer(fs['wvcb'])
            c2 = Point(pc2).buffer(fs['wvcb'])
            from shapely.ops import unary_union
            shp = unary_union([rt, c1, c2])
            rt0 = Polygon([pa1, pa2, pb2, pb1, pa1])
            shp = shp.difference(rt0)
            
            xx, yy = shp.exterior.xy 
            p = fig['ax'].fill(xx, yy, alpha=0.7)
            var.update({"c": p[0].get_facecolor()[0:3]})
            
            rt1 = Polygon([pa1, pd1, pd2, pb1, pa1])
            xx, yy = rt1.exterior.xy 
            fig['ax'].plot(xx, yy, "-", linewidth=0.5, c=var['c'])
            
            pe1 = pd1 + np.array([fs['wvcb'], 0.5*var['h']])
            pf1 = pe1 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg1 = pe1 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe2 = pe1 + np.array([fs['wvcb'], 0.])
            pf2 = pe2 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg2 = pe2 + np.array([-fs['wvcb'], -0.5*var['h']])            
            pe3 = pe2 + np.array([fs['wvcb'], 0.])
            pf3 = pe3 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg3 = pe3 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe4 = pe3 + np.array([fs['wvcb'], 0.])
            pf4 = pe4 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg4 = pe4 + np.array([-fs['wvcb'], -0.5*var['h']])
            pe5 = pe4 + np.array([fs['wvcb'], 0.])
            pf5 = pe5 + np.array([-fs['wvcb'], 0.5*var['h']])
            pg5 = pe5 + np.array([-fs['wvcb'], -0.5*var['h']])
            
            s1 = Polygon([pe1, pd1, pd2, pe1])
            s3 = Polygon([pe2, pf2, pf3, pe3, pg3, pg2, pe2])
            s4 = Polygon([pe4, pf4, pf5, pe5, pg5, pg4, pe4])
            # shp = unary_union([s1, s2, s3])
            xx, yy = s1.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)                    
            xx, yy = s3.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)
            xx, yy = s4.exterior.xy 
            fig['ax'].fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
                                
        draw_logo()     

        ch, cs = len(fname), 16
        fname = [fname[i: i+cs] for i in range(0, ch, cs)]
        fname = "\n".join(fname)

        fig['ax'].text(
            var['x'] + 5*fs['wvcb'] + fl['offset'],
            var['y'] + 0.5*var['h'],
            fname, 
            fontsize=fl['fontsize'],
            color=var['c'],
            fontname=fl['fontname'],
            fontweight=fl['fontweight'],
            fontstyle=fl['fontstyle'],
            va = 'center',
            ha = 'left',
            wrap = True
        )
        
    def draw_scan_variable(self, fig):
        x0 = 0. 
        y0 = -(fig['size']['h'] - fig['layers'][-1]['height']) / 2.0 
        ts = fig['colorset']['scanvarlabel']
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']
        lsx, lsy = self.get_label_shape([x0, x0 + fs['wvc']], [y0-fs['h2'], y0-fs['h2']], 2*fs['wvcb'], fs['h2'])
        fig['ax'].fill(lsx, lsy, facecolor=ts['labelcolor'])
        fig['ax'].plot(lsx, lsy, '-', linewidth=0.5, color=ts['labelcolor'])
        fig['ax'].text(
            0.5 * fs['wvc'], 
            y0 - 0.5 * ts['offset'] - 0.5 * fs['h2'], 
            "SCAN PARAMETERS", 
            fontsize=ts['fontsize'],
            color=ts['color'],
            fontweight=ts['fontweight'],
            fontstyle=ts['fontstyle'],
            fontname=ts['fontname'],
            horizontalalignment='center', 
            verticalalignment='center'
        )

        for kk, var in fig['layers'][-1]['scan variables']['locs'].items():
            var.update({
                "x":    var['x'] + x0, 
                "y":    var['y'] + y0, 
            })
            p = fig['ax'].fill(
                [var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x']  + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb']],
                [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']]
            )
            var.update({"c": p[0].get_facecolor()[0:3]})
          
            fig['ax'].text(
                var['x'] + fs['wvc'] - vs['offset'], 
                var['y'] + 0.5*var['h'], 
                kk, 
                fontsize= vs['fontsize'],
                color=vs['fontcolor'],
                fontname=vs['namefont'],
                fontstyle=vs['namestyle'],
                fontweight=vs['nameweight'],
                horizontalalignment='right', 
                verticalalignment='center'            
            )            
            fig['ax'].text(
                var['x'] + fs['wvc'] + vs['offset'], 
                var['y'] + 0.5*var['h'], 
                var['meth'], 
                fontsize= vs['fontsize'],
                color=vs['fontcolor'],
                fontname=vs['methfont'],
                fontstyle=vs['methstyle'],
                fontweight=vs['methweight'],
                horizontalalignment='left', 
                verticalalignment='center'            
            )
        nscvar = len(fig['layers'][-1]['scan variables']['locs'])
        ii = 0
        for kk, var in fig['layers'][-1]['scan variables']['locs'].items():
            if ii != nscvar - 1:
                self.draw_outvar_label(fig, var)
            else:
                self.draw_endoutvar_label(fig, var)
            ii += 1 
        # print(fig['layers'][-1]['scan variables'], fig['size'])
          
    def draw_outvar_label(self, fig, loc):
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']
        fig['ax'].plot(
            [loc['x'], loc['x'], loc['x'] + fs['wvc'] + 0.4*fs['wvcb']],
            [loc['y']+loc['h'], loc['y'], loc['y']],
            '-',
            linewidth=0.5, 
            color=vs['labelcolor'],
            zorder=0
        )
        
    def draw_inpvar_label(self, fig, loc):
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']
        fig['ax'].plot(
            [loc['x'] + loc['w'], loc['x'] + loc['w'], loc['x'] + loc['w'] - fs['wvc'] - 0.4*fs['wvcb']],
            [loc['y']+loc['h'], loc['y'], loc['y']],
            '-',
            linewidth=0.5, 
            color=vs['labelcolor'],
            zorder=0
        )
        
    def draw_endoutvar_label(self, fig, loc):
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']
        from shapely.geometry import Polygon, Point 
        po  = np.array([loc['x'] + fs['wvcb'], loc['y'] + fs['wvcb']])
        pc1 = po + np.array([0., -fs['wvcb']])
        pc3 = po + np.array([fs['wvc']-fs['wvcb'], 0.])
        pc2 = pc1 + np.array([fs['wvc']-fs['wvcb'], 0.])
        pc7 = po + np.array([-fs['wvcb'], 0.])
        pc6 = pc7 + np.array([0., loc['h']-fs['wvcb']])
        pc5 = po + np.array([0., loc['h']-fs['wvcb']])
        pc4 = pc3 + np.array([0., loc['h']-fs['wvcb']])
        
        rt1 = Polygon([pc7, po, pc5, pc6, pc7])
        rt2 = Polygon([po, pc3, pc4, pc5, po])
        rt3 = Polygon([pc1, pc2, pc3, po, pc1])
        c1  = Point(po).buffer(fs['wvcb'])
        
        from shapely.ops import unary_union
        usp = unary_union([rt1, rt2, rt3, c1])
        xx, yy = usp.exterior.xy 
        fig['ax'].plot(xx, yy, "-", linewidth=0.5, color=vs['labelcolor'], zorder=0)
               
    def draw_endinpvar_label(self, fig, loc):
        fs = fig['colorset']['figsize']
        vs = fig['colorset']['varlabel']
        from shapely.geometry import Polygon, Point 
        po  = np.array([loc['x'] + loc['w'] - fs['wvcb'], loc['y'] + fs['wvcb']])
        pc1 = po + np.array([0., -fs['wvcb']])
        pc3 = po - np.array([fs['wvc']-fs['wvcb'], 0.])
        pc2 = pc1 - np.array([fs['wvc']-fs['wvcb'], 0.])
        pc7 = po + np.array([fs['wvcb'], 0.])
        pc6 = pc7 + np.array([0., loc['h']-fs['wvcb']])
        pc5 = po + np.array([0., loc['h']-fs['wvcb']])
        pc4 = pc3 + np.array([0., loc['h']-fs['wvcb']])
        
        rt1 = Polygon([pc7, po, pc5, pc6, pc7])
        rt2 = Polygon([po, pc3, pc4, pc5, po])
        rt3 = Polygon([pc1, pc2, pc3, po, pc1])
        c1  = Point(po).buffer(fs['wvcb'])
        
        from shapely.ops import unary_union
        usp = unary_union([rt1, rt2, rt3, c1])
        xx, yy = usp.exterior.xy 
        fig['ax'].plot(xx, yy, "-", linewidth=0.5, color=vs['labelcolor'], zorder=0)
                  
    def get_label_shape(self, x, y, r, h):
        from shapely.geometry import Polygon, Point 
        pa  = np.array([x[0], y[0]])
        pb  = np.array([x[1], y[1]])
        v_ab= np.array([x[1]-x[0], y[1] - y[0]])
        l   = np.linalg.norm(v_ab)
        e_ab= v_ab / l
        nv  = np.array([-e_ab[1], e_ab[0]])
        # l_ab = np.linalg.mag(v_ab)

        pd1 = pb + (h-r)*nv 
        pd2 = pa + (h-r)*nv 
        pc2 = pd1 - r*e_ab 
        pc1 = pd2 + r*e_ab 
        pc4 = pc1 + r*nv 
        pc3 = pc2 + r*nv 
        
        c1  = Point(pc1).buffer(r)
        c2  = Point(pc2).buffer(r)
        rt1 = Polygon([pc1, pc2, pc3, pc4, pc1])
        rt2 = Polygon([pa, pb, pd1, pd2, pa])
        
        from shapely.ops import unary_union
        shp = unary_union([c1, c2, rt1, rt2])
        return shp.exterior.xy 
            
    def find_var_nextlayer_number(self, var, nc, dat):
        n_NL = 0 
        n_NNL = 0  
        for ii in range(len(dat['plot']['layer'])):
            if ii == nc + 1:
                for node in dat['plot']['layer'][ii]:
                    node = dat['plot']['nodes'][node]['input file']
                    for jj in range(len(node)):
                        for vname in node[jj]['vars'].keys():
                            if vname == var:
                                n_NL += 1
            elif ii > nc + 1:
                for node in dat['plot']['layer'][ii]:
                    node = dat['plot']['nodes'][node]['input file']
                    for jj in range(len(node)):
                        for vname in node[jj]['vars'].keys():
                            if vname == var:
                                n_NNL += 1
        return {"nNL":  n_NL, "nNNL":   n_NNL}
        
    def plotSLHA(self, fig):
        print("=== Ploting Figure : {} ===".format(fig['name']))
        fig['start'] = time.time()
        self.slha_check_config(fig)
        self.slha_load_colorset(fig)
        if fig['standard']:
            self.read_slha_info(fig)
            self.make_slha_canvas(fig)
            self.draw_slha_particle(fig)
            if self.cf.has_option(fig['section'], "Text"):
                self.drawtext(fig, fig['ax'])
            if 'save' in self.cf.get(fig['section'], 'print_mode'):
                fig['fig'] = plt
                fig['file'] = os.path.join(self.figpath, fig['name'])
                self.savefig(fig, plt)
            if 'show' in self.cf.get(fig['section'], 'print_mode'):
                plt.show(block=False)
                plt.pause(1)
                input("\n\tPress 'Enter' to continue ...")
                plt.close()

    def refine_particle_list(self, fig, sect, od):
        res = {}
        res['order'] = od
        res['charged'] = []
        res['neutral'] = []
        res['showlist'] = []
        for pp, info in sect.items():
            if info['Mass'] > fig['yaxis']['lim'][0] and info['Mass'] < fig['yaxis']['lim'][1]:
                if info['ElectricCharge']:
                    res['charged'].append(pp)
                else:
                    res['neutral'].append(pp)
                res['showlist'].append(pp)
        res['dfpos'] = True
        if len(res['charged'])*len(res['neutral']):
            res['dfpos'] = False
        for pp, info in sect.items():
            if pp in res['showlist']:
                fig['show'][pp] = info
                fig['show'][pp]['labelva'] = "top"
                fig['show'][pp]['color'] = fig['colors'][info['color']]
                if res['dfpos']:
                    if pp in res['charged']:
                        fig['show'][pp]['labelha'] = "right"
                        fig['show'][pp]['labelpos'] = fig['fs']['mass']['defaultpos'][0] + \
                            res['order'] - fig['fs']['mass']['labeloffline']
                    elif pp in res['neutral']:
                        fig['show'][pp]['labelha'] = 'left'
                        fig['show'][pp]['labelpos'] = fig['fs']['mass']['defaultpos'][1] + \
                            res['order'] + fig['fs']['mass']['labeloffline']
                    # fig['show'][pp]['color']        = fig['fs']['mass']['defaultcolor']
                    fig['show'][pp]['masspos'] = [fig['fs']['mass']['defaultpos'][0] +
                                                  res['order'], fig['fs']['mass']['defaultpos'][1] + res['order']]
                    fig['show'][pp]['centerpos'] = 0.5 * \
                        sum(fig['show'][pp]['masspos'])
                else:
                    if pp in res['charged']:
                        fig['show'][pp]['labelha'] = "right"
                        fig['show'][pp]['labelpos'] = fig['fs']['mass']['chargedpos'][0] + \
                            res['order'] - fig['fs']['mass']['labeloffline']
                        fig['show'][pp]['masspos'] = [fig['fs']['mass']['chargedpos'][0] +
                                                      res['order'], fig['fs']['mass']['chargedpos'][1] + res['order']]
                        fig['show'][pp]['centerpos'] = 0.5 * \
                            sum(fig['show'][pp]['masspos'])
                        # fig['show'][pp]['color']        = fig['fs']['mass']['chargedcolor']
                    elif pp in res['neutral']:
                        fig['show'][pp]['labelha'] = "left"
                        fig['show'][pp]['labelpos'] = fig['fs']['mass']['neutralpos'][1] + \
                            res['order'] + fig['fs']['mass']['labeloffline']
                        fig['show'][pp]['masspos'] = [fig['fs']['mass']['neutralpos'][0] +
                                                      res['order'], fig['fs']['mass']['neutralpos'][1] + res['order']]
                        fig['show'][pp]['centerpos'] = 0.5 * \
                            sum(fig['show'][pp]['masspos'])
                        # fig['show'][pp]['color']        = fig['fs']['mass']['neutralcolor']

    def draw_slha_particle(self, fig):
        from adjustText import adjust_text
        if fig['slha']['PDGConfirm']:
            step = 0
            fig['show'] = {}
            tboxc = []
            tboxn = []
            for ss, sect in fig['particle'].items():
                self.refine_particle_list(fig, sect, step)
                step += 1
            for pp, info in fig['show'].items():
                fig['ax'].plot(
                    info['masspos'],
                    [info['Mass'], info['Mass']],
                    '-',
                    linewidth=fig['fs']['mass']['linewidth'],
                    color=info['color'],
                    zorder=99
                )
                if info['labelha'] == 'left':
                    tboxn.append(fig['ax'].text(
                        info['labelpos'], info['Mass'],
                        r"${}$".format(info['LaTeX']),
                        fontsize=fig['fs']['mass']['labelsize'],
                        horizontalalignment=info['labelha'],
                        verticalalignment=info['labelva'],
                        color=fig['fs']["mass"]['labelcolor'],
                        zorder=99
                    ))
                if info['labelha'] == 'right':
                    tboxc.append(fig['ax'].text(
                        info['labelpos'], info['Mass'],
                        r"${}$".format(info['LaTeX']),
                        fontsize=fig['fs']['mass']['labelsize'],
                        horizontalalignment=info['labelha'],
                        verticalalignment=info['labelva'],
                        color=fig['fs']["mass"]['labelcolor'],
                        zorder=99
                    ))
            adjust_text(tboxc, autoalign='y', only_move={
                        "points": "x", 'text': "y", 'objects': 'y'}, ha="right")
            adjust_text(tboxn, autoalign='y', only_move={
                        "points": "x", 'text': "y", 'objects': 'y'}, ha="left")
        if fig['slha']['DecayConfirm']:
            fig['decay'] = fig['fs']['decay']
            fig['decay']['cmap'] = self.load_colormap(
                fig['fs'], fig['fs']['decay']['colormap'])

            for pp, info in fig['show'].items():
                for fs in info["Decay_list"]:
                    if fs['br'] > fig['decay']['minbr'] and fs['fstate'] in fig['show'].keys():
                        fs['width'] = (fs['br'] - fig['decay']['minbr']) * (
                            fig['decay']['maxwidth'] - fig['decay']['minwidth']) + fig['decay']['minwidth']
                        self.arrow(
                            fig['ax'], [info['centerpos'], info['Mass']],
                            [fig['show'][fs['fstate']]['centerpos'], fig['show'][fs['fstate']]['Mass']],
                            fig['decay'],
                            info, fs, fig['fs']
                        )

    def ax_transform_axes(self, ax, point):
        trans = ax.transData.transform(tuple(point))
        trans = ax.transAxes.inverted().transform(trans)
        return trans

    def arrow(self, ax, sp, ep, cmap, info, fs, ffs):
        sp = self.ax_transform_axes(ax, sp)
        ep = self.ax_transform_axes(ax, ep)
        ax.plot(
            [sp[0], ep[0]],
            [sp[1], ep[1]],
            linestyle = cmap['linestyle'],
            color = cmap['cmap']['cmap'](fs['br']),
            alpha = fs['br'],
            linewidth = fs['width'],
            transform=ax.transAxes
        )
        lens = math.sqrt((ep - sp)[0]**2 + (ep-sp)[1]**2)
        dirc = (ep - sp) / lens 
        vert = (dirc[1], -dirc[0]* ffs['canvas']['axratio'] * ffs['canvas']['axratio'])

        if lens > ffs['decay']['minlens']:
            arcood = [
                [ep[0], ep[0]-(ffs['decay']['head_length'] + ffs['decay']['head_outer'])*dirc[0]-ffs['decay']['head_width']*vert[0], ep[0]-ffs['decay']['head_length']*dirc[0], ep[0]-(ffs['decay']['head_length'] + ffs['decay']['head_outer'])*dirc[0] + ffs['decay']['head_width']*vert[0], ep[0]],
                [ep[1], ep[1]-(ffs['decay']['head_length'] + ffs['decay']['head_outer'])*dirc[1]-ffs['decay']['head_width']*vert[1], ep[1]-ffs['decay']['head_length']*dirc[1], ep[1]-(ffs['decay']['head_length'] + ffs['decay']['head_outer'])*dirc[1] + ffs['decay']['head_width']*vert[1], ep[1]]
            ]
            ax.fill(
                arcood[0],
                arcood[1],
                color =     cmap['cmap']['cmap'](fs['br']),
                alpha =     fs['br'],
                edgecolor = None,
                transform = ax.transAxes
            )
        # if lens >= 0.01:
            # dir 
        # print(sp, ep, ep-sp, math.sqrt((ep - sp)[0]**2 + (ep-sp)[1]**2))
        
    def make_slha_canvas(self, fig):
        with open(self.load_path(fig['colorset']['figureSetting']), "r") as f1:
            fig['fs'] = json.loads(f1.read())
        fig['fs']['canvas']['width'] = len(
            fig["particle"])*fig['fs']['canvas']['block'] + fig['fs']['canvas']['left'] + fig['fs']['canvas']['right']
        fig['fig'] = plt.figure(
            figsize=(fig['fs']['canvas']['width'], fig['fs']['canvas']['height']))
        fig['fs']['canvas']['axsize'] = {
            "xmin": fig['fs']['canvas']['left']/fig['fs']['canvas']['width'],
            "ymin": fig['fs']['canvas']['bottom']/fig['fs']['canvas']['height'],
            "xwidth": len(fig["particle"])*fig['fs']['canvas']['block']/fig['fs']['canvas']['width'],
            "ywidth": 1.0 - (fig['fs']['canvas']['bottom'] + fig['fs']['canvas']['top'])/fig['fs']['canvas']['height']
        }
        fig['fig'].text(0., 0., "Test", color="None")
        fig['fig'].text(1., 1., "Test", color="None")
        fig['cv'] = fig['fig'].add_axes([0., 0., 1., 1.])
        fig['cv'].axis("off")
        fig['cv'].set_xlim(0, 1)
        fig['cv'].set_ylim(0, 1)
        fig['ax'] = fig['fig'].add_axes([
            fig['fs']['canvas']['axsize']['xmin'],
            fig['fs']['canvas']['axsize']['ymin'],
            fig['fs']['canvas']['axsize']['xwidth'],
            fig['fs']['canvas']['axsize']['ywidth']
        ])
        fig['fs']['canvas']['axratio'] = (fig['fs']['canvas']['axsize']['xwidth']*fig['fs']['canvas']['width']) / (fig['fs']['canvas']['axsize']['ywidth']*fig['fs']['canvas']['height'] )
        fig['ax'].text(
            -0.11, 0.5,
            r"{}".format(fig['fs']['yaxis']['label']),
            fontsize=fig['fs']['yaxis']['fontsize'],
            horizontalalignment='right',
            verticalalignment='center',
            rotation=fig['fs']['yaxis']['rotation'],
            color='black',
            transform=fig['ax'].transAxes
        )
        fig['ax'].set_xlim(0, len(fig['particle']))
        fig['ax'].spines['bottom'].set_color(None)
        fig['ax'].spines['top'].set_color(None)
        fig['ax'].spines['right'].set_color(None)
        fig['ax'].spines['left'].set_color(fig['fs']['yaxis']['fc'])
        from matplotlib.ticker import AutoMinorLocator, FixedLocator
        if self.cf.get(fig['section'], 'm_scale').strip().lower() == 'flat':
            if not self.cf.get(fig['section'], 'm_ticks')[0:4] == "Manu":
                fig['ax'].set_yticks(fig['yaxis']['ticks'])
                fig['ax'].yaxis.set_minor_locator(AutoMinorLocator())
        elif self.cf.get(fig['section'], 'm_scale').strip().lower() == 'log':
            fig['ax'].set_yscale('log')
        if self.cf.get(fig['section'], 'm_ticks')[0:4] == "Manu":
            fig['ax'].yaxis.set_major_locator(
                ticker.FixedLocator(fig['yaxis']['ticks'][0]))
            fig['ax'].set_yticklabels(fig['yaxis']['ticks'][1])
            if self.cf.get(fig['section'], 'm_scale').strip().lower() == 'flat':
                fig['ax'].yaxis.set_minor_locator(AutoMinorLocator())
        fig['ax'].set_ylim(fig['yaxis']['lim'][0], fig['yaxis']['lim'][1])
        xaxis = fig['ax'].get_xaxis()
        xaxis.set_visible(False)
        fig['cv'].arrow(
            fig['fs']['canvas']['axsize']['xmin'],
            fig['fs']['canvas']['axsize']['ymin'] -
            fig['fs']['yaxis']['extend'],
            0.,
            fig['fs']['canvas']['axsize']['ywidth'] +
            2.0*fig['fs']['yaxis']['extend'],
            head_width=fig['fs']['yaxis']['head_width'],
            head_length=fig['fs']['yaxis']['head_length'],
            fc=fig['fs']['yaxis']['fc'],
            ec=fig['fs']['yaxis']['ec'],
            shape=fig['fs']['yaxis']['shape'],
            joinstyle=fig['fs']['yaxis']['joinstyle']
        )
        fig['ax'].tick_params(
            which="major",
            direction=fig['fs']['yticks']["direction"],
            length=fig['fs']['yticks']["majorlength"],
            color=fig['fs']['yticks']["color"],
            width=fig['fs']['yticks']['width'],
            labelsize=fig['fs']['yticks']['labelsize']
        )
        fig['ax'].tick_params(
            which="minor",
            direction=fig['fs']['yticks']["direction"],
            length=fig['fs']['yticks']["minorlength"],
            color=fig['fs']['yticks']["color"],
            width=fig['fs']['yticks']['width']
        )

    def load_slha_model(self, fig):
        fig['particle'] = {}
        with open(fig['model']['code'], 'r') as f1:
            import json
            fig['model']['setup'] = json.loads(f1.read())
        for pp in self.cf.get(fig['section'], "particle").split('+'):
            pp = pp.strip()
            if pp in fig['model']['setup']['Sectors']:
                fig['particle'][pp] = {}
            else:
                print(
                    "\tWarning: No particle list ->  '{}' found in model '{}'".format(pp, fig['model']['name']))

    def load_particle_info(self, fig):
        print("\tTimer: {:.2f} Second;  Message from {} : particle loading from file ->\n\t\t{} ".format(
            time.time()-fig['start'], fig['section'], fig['model']['PDGid']))
        with open(fig['model']['PDGid'], 'r') as f1:
            PDGids = json.loads(f1.read())
            for pp in fig['particle']:
                for item in fig['model']['setup']['Sections'][pp]['particles']:
                    fig['particle'][pp][item] = PDGids[item]
                # print(fig['model']['setup']['Sections'][pp])
                if fig['model']['setup']['Sections'][pp]['subsection']:
                    # print(fig['model']['setup']['Sections'][pp]['subsection'])
                    for sect in fig['model']['setup']['Sections'][pp]['subsection']:
                        for item in fig['model']['setup']['Sections'][sect]['particles']:
                            if item in PDGids:
                                # print(PDGids[item])
                                fig['particle'][pp][item] = PDGids[item]
                            else:
                                print(emoji.emojize('\t:ghost::ghost::ghost: Multilist Particle {} not found in Model particle list !!\n\tPlease check your model file \n'.format(
                                    item), use_aliases=True))
                                sys.exit(0)

    def slha_load_colorset(self, fig):
        if fig['standard']:
            mtag = True
            self.load_fig_setting(fig)
            for model in fig['colorset']['models']:
                if self.cf.get(fig['section'], "model").strip() == model['name']:
                    model['code'] = self.load_path(model['code'])
                    model['PDGid'] = self.load_path(model['PDGid'])
                    model['slhaReadertemplet'] = self.load_path(
                        model['slhaReadertemplet'])
                    fig['model'] = model
                    with open(fig['colorset']['colorcardpath'], 'r') as f1:
                        fig['colors'] = json.loads(f1.read())
                    mtag = False
                    break
            if mtag:
                fig['standard'] = False
                print(emoji.emojize('\t:ghost::ghost::ghost: Model not support in the current version -> {} !! \n\t\tPlease check your configure file'.format(
                    self.cf.get(fig['section'], "model").strip()), use_aliases=True))

    def read_slha_info(self, fig):
        self.load_slha_model(fig)
        self.load_particle_info(fig)
        from Func_lab import make_slhaReader_mass_input
        make_slhaReader_mass_input(fig)
        sys.path.append(fig['colorset']['slhareaderpath'])
        from Func_lab import block_screen_print
        block_screen_print()
        from reader import reader
        slhareader = reader()
        slhareader.set_parser(fig['slha']['rcfg'])
        from Func_lab import reset_slhaReader_mass_input
        reset_slhaReader_mass_input(fig)
        slhareader.set_parser(fig['slha']['rcfg'])
        from Func_lab import enable_screen_print
        enable_screen_print()
        from Func_lab import load_particle_info
        load_particle_info(fig)
        print("\tTimer: {:.2f} Second;  Message from SLHAReader : Reading particle Complete !!!".format(
            time.time()-fig['start'], fig['model']['PDGid']))

    def slha_check_config(self, fig):
        fig['standard'] = True
        if not self.cf.has_option(fig['section'], 'colorset'):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "colorset" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        if not self.cf.has_option("COLORMAP", "colorsetting"):
            fig['standard'] = False
            print(emoji.emojize(
                '\t:ghost::ghost::ghost: No option "colorsetting" found in Section [COLORMAP] !!\n\tPlease check your configure file \n', use_aliases=True))
        if not self.cf.has_option(fig['section'], 'model'):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "model" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        if not self.cf.has_option(fig['section'], "particle"):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "particle" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        if not self.cf.has_option(fig['section'], 'spectra'):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "spectra" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        fig['slha'] = {}
        fig['slha']['patt'] = self.cf.get(
            fig['section'], "spectra").split(",")[0].strip()
        fig['slha']['norm'] = self.cf.get(
            fig['section'], "spectra").split(",")[1].strip()
        fig['slha']['file'] = self.load_path(self.cf.get(
            fig['section'], "spectra").split(",")[2].strip())
        fig['slha']['svdr'] = os.path.join(self.figpath, fig['name'])
        if not os.path.exists(fig['slha']['svdr']):
            os.makedirs(fig['slha']['svdr'])
        fig['slha']['rcfg'] = os.path.join(
            fig['slha']['svdr'], "slhaReader.ini")
        file = fig['slha']['norm'].replace("$ID$", "")
        fig['slha']['ID'] = os.path.basename(
            fig['slha']['file']).replace(file, "")
        fig['slha']['PDGConfirm'] = True
        fig['slha']['DecayConfirm'] = True
        from Func_lab import check_slhaReader_pattern
        if not check_slhaReader_pattern(fig['slha']):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: Wrong grammer "spectra" detected !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        if not self.cf.has_option(fig['section'], "mode"):
            fig['slha']['DecayConfirm'] = False
        else:
            mode = self.cf.get(fig['section'], "mode").lower()
            if "mass" not in mode:
                fig['slha']['PDGConfirm'] = False
            if "decay" not in mode:
                fig['slha']['DecayConfirm'] = False
        fig['yaxis'] = {}
        if not self.cf.has_option(fig['section'], "m_range"):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "m_range" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        else:
            lim = self.cf.get(fig['section'], "m_range").split(',')
            fig['yaxis']['lim'] = [
                float(lim[0]),
                float(lim[1])
            ]
        if not self.cf.has_option(fig['section'], "m_ticks"):
            fig['standard'] = False
            print(emoji.emojize('\t:ghost::ghost::ghost: No option "m_ticks" found in Section [{}] !!\n\tPlease check your configure file \n'.format(
                fig['section']), use_aliases=True))
        else:
            fig['yaxis']['ticks'] = {}
            tick = self.cf.get(fig['section'], '{}_ticks'.format("m"))
            if 'AUTO' in tick:
                a = float(tick.split('_')[-1])
                low = fig['yaxis']['lim'][0] // a
                upp = fig['yaxis']['lim'][1] // a + 1
                fig['yaxis']['ticks'] = np.linspace(
                    low*a, upp*a, int(upp-low+1))
            elif tick[0:4] == 'Manu':
                p_rec = re.compile(r'[[].*?[]]', re.S)
                a = re.findall(p_rec, tick[5:].strip())
                tk = []
                label = []
                for it in a:
                    it = it.strip().strip('[').strip(']')
                    tk.append(float(it.split(',')[0]))
                    label.append(r"{}".format(it.split(',')[1].strip()))
                fig['yaxis']['ticks'] = tuple([tk, label])
            else:
                tick = tick.split(',')
                for ii in range(len(tick)):
                    tick[ii] = tick[ii].strip()
                fig['yaxis']['ticks'] = np.linspace(
                    float(tick[0]), float(tick[1]), int(tick[2]))

    # Plot Method is write in the drawpicture :
    # Add " elif fig['type'] == $FIGURE_TYPE$: "
    def drawpicture(self, fig):
        fig['start'] = time.time()
        fig['fig'] = plt.figure(
            figsize=(
                fig['colorset']['figureSize']['width'], 
                fig['colorset']["figureSize"]['height'])
        )
        fig['fig'].text(0., 0., "Test", color="None")
        fig['fig'].text(1., 1., "Test", color="None")
        print("\n=== Ploting Figure : {} ===".format(fig['name']))
        self.basic_selection(fig)
        if fig['data'].shape[0] < 1:
            print(emoji.emojize(
                '\t:ghost::ghost::ghost: No data selected!!\n\tPlease check your Data or the selection condition \n', use_aliases=True))
        else:
            print(emoji.emojize("\t:space_invader::space_invader::space_invader:\tSelected Data is    -> {}\trows".format(
                fig['data'].shape[0]), use_aliases=True))
            from matplotlib.ticker import AutoMinorLocator, FixedLocator
            if   fig['type'] == "2D_Stat_Profile":
                ax = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axbox']['x0'], 
                        fig['colorset']['figureSize']['axbox']['y0'], 
                        fig['colorset']['figureSize']['axbox']['width'], 
                        fig['colorset']['figureSize']['axbox']['height'] 
                ])                
                axc = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axcbox']['x0'], 
                        fig['colorset']['figureSize']['axcbox']['y0'], 
                        fig['colorset']['figureSize']['axcbox']['width'], 
                        fig['colorset']['figureSize']['axcbox']['height'] 
                ])
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                self.GetStatData(fig)
                fig['var']['BestPoint'] = {
                    'x':    fig['var']['data'].loc[fig['var']['data'].Stat.idxmin()].x,
                    'y':    fig['var']['data'].loc[fig['var']['data'].Stat.idxmin()].y,
                    'Stat': fig['var']['data'].loc[fig['var']['data'].Stat.idxmin()].Stat
                }
                # Two Method to achieve Profile Likelihood Data in Pandas! Lambda expr is more compact in code
                # fig['var']['data'] = fig['var']['data'].assign( PL=np.exp( -0.5 * fig['var']['data'].Stat )/np.exp(-0.5 * fig['var']['BestPoint']['Stat'])  )
                fig['var']['data'] = fig['var']['data'].assign(
                    PL=lambda x: np.exp(-0.5 * x.Stat)/np.exp(-0.5 * fig['var']['BestPoint']['Stat']))
                fig['var']['lim'] = {
                    'x': [fig['var']['data'].x.min(), fig['var']['data'].x.max()],
                    'y': [fig['var']['data'].y.min(), fig['var']['data'].y.max()],
                }
                fig['ax'] = {}
                self.ax_setlim(fig, 'xyc')
                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat" and self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    XI, YI = np.meshgrid(
                        np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], int(
                            self.cf.get(fig['section'], 'x_nbin'))+1),
                        np.linspace(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1], int(
                            self.cf.get(fig['section'], 'y_nbin'))+1)
                    )
                    fig['var']['AxesBP'] = {
                        'x':    (fig['var']['BestPoint']['x'] - fig['ax']['lim']['x'][0])/(fig['ax']['lim']['x'][1] - fig['ax']['lim']['x'][0]),
                        'y':    (fig['var']['BestPoint']['y'] - fig['ax']['lim']['y'][0])/(fig['ax']['lim']['y'][1] - fig['ax']['lim']['y'][0])
                    }
                    fig['ax']['var'] = {
                        'x':    XI,
                        'y':    YI,
                        'dx':   (fig['ax']['lim']['x'][1] - fig['ax']['lim']['x'][0])/int(self.cf.get(fig['section'], 'x_nbin')),
                        'dy':   (fig['ax']['lim']['y'][1] - fig['ax']['lim']['y'][0])/int(self.cf.get(fig['section'], 'y_nbin'))
                    }
                    fig['ax']['grid'] = pd.DataFrame(
                        index=np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], int(
                            self.cf.get(fig['section'], 'x_nbin'))+1),
                        columns=np.linspace(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1], int(
                            self.cf.get(fig['section'], 'y_nbin'))+1)
                    ).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
                    fig['ax']['grid'] = pd.DataFrame({
                        "xi":   fig['ax']['grid']['xx'],
                        'yi':   fig['ax']['grid']['yy'],
                        'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x > tt['xx'] - 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].x < tt['xx'] + 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].y > tt['yy'] - 0.5*fig['ax']['var']['dy']) & (fig['var']['data'].y < tt['yy'] + 0.5*fig['ax']['var']['dy'])].PL.max(axis=0, skipna=True), axis=1)
                    }).fillna({'PL': -0.1})
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log" and self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    XI, YI = np.meshgrid(
                        np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), int(self.cf.get(fig['section'], 'x_nbin'))+1),
                        np.linspace(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1], int(
                            self.cf.get(fig['section'], 'y_nbin'))+1)
                    )
                    fig['var']['AxesBP'] = {
                        'x':    (math.log10(fig['var']['BestPoint']['x']) - math.log10(fig['ax']['lim']['x'][0]))/(math.log10(fig['ax']['lim']['x'][1]) - math.log10(fig['ax']['lim']['x'][0])),
                        'y':    (fig['var']['BestPoint']['y'] - fig['ax']['lim']['y'][0])/(fig['ax']['lim']['y'][1] - fig['ax']['lim']['y'][0])
                    }
                    fig['ax']['var'] = {
                        'x':    XI,
                        'y':    YI,
                        'dx':   (math.log10(fig['ax']['lim']['x'][1]) - math.log10(fig['ax']['lim']['x'][0]))/int(self.cf.get(fig['section'], 'x_nbin')),
                        'dy':   (fig['ax']['lim']['y'][1] - fig['ax']['lim']['y'][0])/int(self.cf.get(fig['section'], 'y_nbin'))
                    }
                    fig['ax']['grid'] = pd.DataFrame(
                        index=np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), int(self.cf.get(fig['section'], 'x_nbin'))+1),
                        columns=np.linspace(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1], int(
                            self.cf.get(fig['section'], 'y_nbin'))+1)
                    ).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
                    fig['ax']['grid'] = pd.DataFrame({
                        "xi":   fig['ax']['grid']['xx'],
                        'yi':   fig['ax']['grid']['yy'],
                        'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x > 10**(math.log10(tt['xx']) - 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].x < 10**(math.log10(tt['xx']) + 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].y > tt['yy'] - 0.5*fig['ax']['var']['dy']) & (fig['var']['data'].y < tt['yy'] + 0.5*fig['ax']['var']['dy'])].PL.max(axis=0, skipna=True), axis=1)
                    }).fillna({'PL': -0.1})
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat" and self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    XI, YI = np.meshgrid(
                        np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], int(
                            self.cf.get(fig['section'], 'x_nbin'))+1),
                        np.logspace(math.log10(fig['ax']['lim']['y'][0]), math.log(
                            fig['ax']['lim']['y'][1]), int(self.cf.get(fig['section'], 'y_nbin'))+1)
                    )
                    fig['var']['AxesBP'] = {
                        'x':    (fig['var']['BestPoint']['x'] - fig['ax']['lim']['x'][0])/(fig['ax']['lim']['x'][1] - fig['ax']['lim']['x'][0]),
                        'y':    (math.log10(fig['var']['BestPoint']['y']) - math.log10(fig['ax']['lim']['y'][0]))/(math.log10(fig['ax']['lim']['y'][1]) - math.log10(fig['ax']['lim']['y'][0]))
                    }
                    fig['ax']['var'] = {
                        'x':    XI,
                        'y':    YI,
                        'dx':   (fig['ax']['lim']['x'][1] - fig['ax']['lim']['x'][0])/int(self.cf.get(fig['section'], 'x_nbin')),
                        'dy':   (math.log10(fig['ax']['lim']['y'][1]) - math.log10(fig['ax']['lim']['y'][0]))/int(self.cf.get(fig['section'], 'y_nbin'))
                    }
                    fig['ax']['grid'] = pd.DataFrame(
                        index=np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], int(
                            self.cf.get(fig['section'], 'x_nbin'))+1),
                        columns=np.logspace(math.log10(fig['ax']['lim']['y'][0]), math.log10(
                            fig['ax']['lim']['y'][1]), int(self.cf.get(fig['section'], 'y_nbin'))+1)
                    ).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
                    fig['ax']['grid'] = pd.DataFrame({
                        "xi":   fig['ax']['grid']['xx'],
                        'yi':   fig['ax']['grid']['yy'],
                        'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x > tt['xx'] - 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].x < tt['xx'] + 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].y > 10**(math.log10(tt['yy']) - 0.5*fig['ax']['var']['dy'])) & (fig['var']['data'].y < 10**(math.log10(tt['yy']) + 0.5*fig['ax']['var']['dy']))].PL.max(axis=0, skipna=True), axis=1)
                    }).fillna({'PL': -0.1})
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log" and self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    XI, YI = np.meshgrid(
                        np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), int(self.cf.get(fig['section'], 'x_nbin'))+1),
                        np.logspace(math.log10(fig['ax']['lim']['y'][0]), math.log(
                            fig['ax']['lim']['y'][1]), int(self.cf.get(fig['section'], 'y_nbin'))+1)
                    )
                    fig['var']['AxesBP'] = {
                        'x':    (math.log10(fig['var']['BestPoint']['x']) - math.log10(fig['ax']['lim']['x'][0]))/(math.log10(fig['ax']['lim']['x'][1]) - math.log10(fig['ax']['lim']['x'][0])),
                        'y':    (math.log10(fig['var']['BestPoint']['y']) - math.log10(fig['ax']['lim']['y'][0]))/(math.log10(fig['ax']['lim']['y'][1]) - math.log10(fig['ax']['lim']['y'][0]))
                    }
                    fig['ax']['var'] = {
                        'x':    XI,
                        'y':    YI,
                        'dx':   (math.log10(fig['ax']['lim']['x'][1]) - math.log10(fig['ax']['lim']['x'][0]))/int(self.cf.get(fig['section'], 'x_nbin')),
                        'dy':   (math.log10(fig['ax']['lim']['y'][1]) - math.log10(fig['ax']['lim']['y'][0]))/int(self.cf.get(fig['section'], 'y_nbin'))
                    }
                    fig['ax']['grid'] = pd.DataFrame(
                        index=np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), int(self.cf.get(fig['section'], 'x_nbin'))+1),
                        columns=np.logspace(math.log10(fig['ax']['lim']['y'][0]), math.log10(
                            fig['ax']['lim']['y'][1]), int(self.cf.get(fig['section'], 'y_nbin'))+1)
                    ).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
                    fig['ax']['grid'] = pd.DataFrame({
                        "xi":   fig['ax']['grid']['xx'],
                        'yi':   fig['ax']['grid']['yy'],
                        'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x > 10**(math.log10(tt['xx']) - 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].x < 10**(math.log10(tt['xx']) + 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].y > 10**(math.log10(tt['yy']) - 0.5*fig['ax']['var']['dy'])) & (fig['var']['data'].y < 10**(math.log10(tt['yy']) + 0.5*fig['ax']['var']['dy']))].PL.max(axis=0, skipna=True), axis=1)
                    }).fillna({'PL': -0.1})
                else:
                    print("No such Mode For X,Y scale -> {},{}".format(self.cf.get(
                        fig['section'], 'x_scale'), self.cf.get(fig['section'], 'y_scale')))
                    sys.exit(0)
                XI, YI = np.meshgrid(
                    np.linspace(0., 1., int(
                        self.cf.get(fig['section'], 'x_nbin'))+1),
                    np.linspace(0., 1., int(
                        self.cf.get(fig['section'], 'y_nbin'))+1)
                )
                fig['ax']['mashgrid'] = pd.DataFrame(
                    index=np.linspace(0., 1., int(
                        self.cf.get(fig['section'], 'x_nbin'))+1),
                    columns=np.linspace(0., 1., int(
                        self.cf.get(fig['section'], 'y_nbin'))+1)
                ).unstack().reset_index().rename(columns={'level_0': 'yy', 'level_1': 'xx', 0: 'z'})
                fig['ax']['mashgrid'] = pd.DataFrame({
                    'x':  fig['ax']['mashgrid']['xx'],
                    'y':  fig['ax']['mashgrid']['yy'],
                    'PL':   fig['ax']['grid']['PL']
                })
                self.ax_setticks(fig, 'xyc')

                from matplotlib.ticker import MaxNLocator
                levels = MaxNLocator(nbins=100).tick_values(
                    fig['ax']['lim']['c'][0], fig['ax']['lim']['c'][1])
                self.ax_setcmap(fig)
                from matplotlib.tri import Triangulation, TriAnalyzer, UniformTriRefiner
                fig['ax']['tri'] = Triangulation(
                    fig['ax']['mashgrid']['x'], fig['ax']['mashgrid']['y'])
                fig['ax']['refiner'] = UniformTriRefiner(fig['ax']['tri'])
                fig['ax']['tri_refine_PL'], fig['ax']['PL_refine'] = fig['ax']['refiner'].refine_field(
                    fig['ax']['mashgrid']['PL'], subdiv=3)
                fig['ax']['PL_refine'] = (fig['ax']['PL_refine'] > 0.) * fig['ax']['PL_refine'] / (
                    np.max(fig['ax']['PL_refine'])) + 0. * (fig['ax']['PL_refine'] < 0.)

                print("\tTimer: {:.2f} Second;  Message from '{}' -> Data analysis completed".format(
                    time.time()-fig['start'], fig['section']))

                a1 = ax.tricontourf(fig['ax']['tri_refine_PL'], fig['ax']['PL_refine'], cmap=fig['colorset']
                                    ['cmap']['cmap'], levels=levels, zorder=1, transform=ax.transAxes)
                plt.colorbar(
                    a1, axc, ticks=fig['ax']['ticks']['c'], orientation='vertical')
                if 'curve' in fig['colorset'].keys():
                    for curve in fig['colorset']['curve']:
                        ct = ax.tricontour(fig['ax']['tri_refine_PL'], fig['ax']['PL_refine'], [math.exp(-0.5 * curve['value'])],
                                           colors=curve['linecolor'], linewidths=curve['linewidth'], zorder=(10-curve['tag'])*4, transform=ax.transAxes)

                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                    if not self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                        ax.set_xticks(fig['ax']['ticks']['x'])
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                    ax.set_xscale('log')
                if self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                    ax.xaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['x'][0]))
                    ax.set_xticklabels(fig['ax']['ticks']['x'][1])
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.set_xlim(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1])

                if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        ax.set_yticks(fig['ax']['ticks']['y'])
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    ax.set_yscale('log')
                if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                    ax.yaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                    ax.set_yticklabels(fig['ax']['ticks']['y'][1])
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.set_ylim(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])

                if self.cf.has_option(fig['section'], 'BestPoint'):
                    self.drawBestPoint(fig, ax)
                    # if eval(self.cf.get(fig['section'], 'BestPoint')) and 'bestpoint' in fig['colorset'].keys():
                    # ax.scatter(fig['var']['BestPoint']['x'], fig['var']['BestPoint']['y'], 300, marker='*', color=fig['colorset']['bestpoint'][0], zorder=2000)
                    # ax.scatter(fig['var']['BestPoint']['x'], fig['var']['BestPoint']['y'], 50, marker='*', color=fig['colorset']['bestpoint'][1], zorder=2100)
                    # print(fig['data'][fig['data']["chi2_h1"] == fig['var']['BestPoint']['Stat']])

                ax.tick_params(
                    labelsize=fig['colorset']['ticks']['labelsize'],
                    direction=fig['colorset']['ticks']['direction'],
                    bottom=fig['colorset']['ticks']['bottom'],
                    left=fig['colorset']['ticks']['left'],
                    top=fig['colorset']['ticks']['top'],
                    right=fig['colorset']['ticks']['right'],
                    which='both'
                )
                ax.tick_params(which='major', length=fig['colorset']['ticks']
                               ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['ticks']
                               ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                axc.tick_params(
                    labelsize=fig['colorset']['colorticks']['labelsize'],
                    direction=fig['colorset']['colorticks']['direction'],
                    bottom=fig['colorset']['colorticks']['bottom'],
                    left=fig['colorset']['colorticks']['left'],
                    top=fig['colorset']['colorticks']['top'],
                    right=fig['colorset']['colorticks']['right'],
                    color=fig['colorset']['colorticks']['color']
                )
                ax.set_xlabel(r"{}".format(self.cf.get(
                    fig['section'], 'x_label')), fontsize=30)
                ax.set_ylabel(r"{}".format(self.cf.get(
                    fig['section'], 'y_label')), fontsize=30)
                axc.set_ylabel(r"{}".format(self.cf.get(
                    fig['section'], 'c_label')), fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.068)
                if self.cf.has_option(fig['section'], 'Line_draw'):
                    self.drawline(fig, ax)
                if self.cf.has_option(fig['section'], "Text"):
                    self.drawtext(fig, ax)
            
            elif fig['type'] == "2D_Scatter":
                ax = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axxbox']['x0'], 
                        fig['colorset']['figureSize']['axxbox']['y0'], 
                        fig['colorset']['figureSize']['axxbox']['width'], 
                        fig['colorset']['figureSize']['axxbox']['height'] 
                ]) 
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                self.Get2DData(fig)
                # print(fig['var']['data'])
                fig['var']['lim'] = {
                    'x': [fig['var']['data'].x.min(), fig['var']['data'].x.max()],
                    'y': [fig['var']['data'].y.min(), fig['var']['data'].y.max()]
                }
                fig['ax'] = {}
                self.ax_setlim(fig, 'xy')
                self.ax_setcmap(fig)
                # print(fig['colorset'])
                if self.cf.has_option(fig['section'], 'marker'):
                    with open(self.load_path(fig['colorset']['markercodepath']), 'r') as f1:
                        fig['colorset']['scattermarker'] = json.loads(
                            f1.read())
                    with open(self.load_path(fig['colorset']["colorpath"]), 'r') as f1:
                        fig['colorset']['scattercolor'] = json.loads(f1.read())
                    lines = self.cf.get(fig['section'], 'marker').split('\n')
                    if len(lines) == 1:
                        marker = lines[0].split(',')
                        if len(marker) == 2:
                            fig['ax']['markercolor'] = marker[0].strip()
                            fig['ax']['markertype'] = marker[1].strip()
                            ax.scatter(
                                fig['var']['data'].x,
                                fig['var']['data'].y,
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                c=fig['colorset']['scattercolor'][fig['ax']
                                                                  ['markercolor']],
                                edgecolor=fig['colorset']['marker']['edgecolor'],
                                s=fig['colorset']['marker']['size']**2,
                                alpha=fig['colorset']['marker']['alpha']
                            )
                        elif len(marker) > 2:
                            if marker[2].strip()[0:3] == "&Bo":
                                fig['ax']['markercolor'] = marker[0].strip()
                                fig['ax']['markertype'] = marker[1].strip()
                                self.scatter_classify_data(
                                    fig, marker[2].strip()[4:])
                                ax.scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    c=fig['colorset']['scattercolor'][fig['ax']
                                                                      ['markercolor']],
                                    edgecolor=fig['colorset']['marker']['edgecolor'],
                                    s=fig['colorset']['marker']['size']**2,
                                    alpha=fig['colorset']['marker']['alpha']
                                )
                    else:
                        for line in lines:
                            marker = line.split(',')
                            if len(marker) > 2:
                                fig['ax']['markercolor'] = marker[0].strip()
                                fig['ax']['markertype'] = marker[1].strip()
                                self.scatter_classify_data(
                                    fig, marker[2].strip()[4:])
                                ax.scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    c=fig['colorset']['scattercolor'][fig['ax']
                                                                      ['markercolor']],
                                    edgecolor=fig['colorset']['marker']['edgecolor'],
                                    s=fig['colorset']['marker']['size']**2,
                                    alpha=fig['colorset']['marker']['alpha']
                                )
                else:
                    fig['ax']['markercolor'] = 'Blue'
                    fig['ax']['markertype'] = 'round'
                    ax.scatter(
                        fig['var']['data'].x,
                        fig['var']['data'].y,
                        marker=fig['colorset']['scattermarker'][fig['ax']
                                                                ['markertype']],
                        c=fig['colorset']['scattercolor'][fig['ax']
                                                          ['markercolor']],
                        edgecolor=fig['colorset']['marker']['edgecolor'],
                        s=fig['colorset']['marker']['size']**2,
                        alpha=fig['colorset']['marker']['alpha']
                    )
                self.ax_setticks(fig, 'xy')

                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                    if not self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                        ax.set_xticks(fig['ax']['ticks']['x'])
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                    ax.set_xscale('log')
                if self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                    ax.xaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['x'][0]))
                    ax.set_xticklabels(fig['ax']['ticks']['x'][1])
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.set_xlim(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1])

                if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        ax.set_yticks(fig['ax']['ticks']['y'])
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    ax.set_yscale('log')
                if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                    ax.yaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                    ax.set_yticklabels(fig['ax']['ticks']['y'][1])
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.set_ylim(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])

                ax.tick_params(
                    labelsize=fig['colorset']['ticks']['labelsize'],
                    direction=fig['colorset']['ticks']['direction'],
                    bottom=fig['colorset']['ticks']['bottom'],
                    left=fig['colorset']['ticks']['left'],
                    top=fig['colorset']['ticks']['top'],
                    right=fig['colorset']['ticks']['right'],
                    which='both'
                )
                ax.tick_params(which='major', length=fig['colorset']['ticks']
                               ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['ticks']
                               ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                ax.set_xlabel(r"{}".format(self.cf.get(
                    fig['section'], 'x_label')), fontsize=30)
                ax.set_ylabel(r"{}".format(self.cf.get(
                    fig['section'], 'y_label')), fontsize=30)
                ax.xaxis.set_label_coords(0.5, -0.068)

                if self.cf.has_option(fig['section'], 'Line_draw'):
                    self.drawline(fig, ax)

                if self.cf.has_option(fig['section'], "Text"):
                    self.drawtext(fig, ax)

                if self.cf.has_option(fig['section'], "fill"):
                    self.drawfillarea(fig, ax)
            
                if self.cf.has_option(fig['section'], "draw_legend"):
                    self.draw_scatter_legend(fig, ax)

            elif fig['type'] == "2DC_Scatter":
                ax = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axbox']['x0'], 
                        fig['colorset']['figureSize']['axbox']['y0'], 
                        fig['colorset']['figureSize']['axbox']['width'], 
                        fig['colorset']['figureSize']['axbox']['height'] 
                ])                
                axc = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axcbox']['x0'], 
                        fig['colorset']['figureSize']['axcbox']['y0'], 
                        fig['colorset']['figureSize']['axcbox']['width'], 
                        fig['colorset']['figureSize']['axcbox']['height'] 
                ])
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                self.Get3DData(fig)
                fig['var']['lim'] = {
                    'x': [fig['var']['data'].x.min(), fig['var']['data'].x.max()],
                    'y': [fig['var']['data'].y.min(), fig['var']['data'].y.max()],
                    'c': [fig['var']['data'].c.min(), fig['var']['data'].c.max()]
                }
                fig['ax'] = {}
                self.ax_setlim(fig, 'xyc')
                self.ax_setcmap(fig)

                fig['ax']['a1'] = []
                if self.cf.has_option(fig['section'], 'marker'):
                    lines = self.cf.get(fig['section'], 'marker').split('\n')
                    with open(self.load_path(fig['colorset']['markercodepath']), 'r') as f1:
                        fig['colorset']['scattermarker'] = json.loads(
                            f1.read())
                    with open(self.load_path(fig['colorset']["colorpath"]), 'r') as f1:
                        fig['colorset']['scattercolor'] = json.loads(f1.read())
                for line in lines:
                    marker = line.split(',')
                    if marker[0].strip() == "&Color":
                        fig['ax']['markertype'] = marker[1].strip()
                        if len(marker) == 2:
                            fig['ax']['a1'].append(ax.scatter(
                                fig['var']['data'].x,
                                fig['var']['data'].y,
                                c=fig['var']['data'].c,
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                edgecolor=fig['colorset']['marker']['edgecolor'],
                                s=fig['colorset']['marker']['size']**2,
                                alpha=fig['colorset']['marker']['alpha'],
                                cmap=fig['colorset']['cmap']['cmap'],
                                vmin=fig['ax']['lim']['c'][0],
                                vmax=fig['ax']['lim']['c'][1]
                            ))
                        elif len(marker) > 2:
                            if marker[2].strip()[0:3] == '&Bo':
                                self.scatter_classify_data(
                                    fig, marker[2].strip()[4:])
                                fig['ax']['a1'].append(ax.scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    c=fig['classify']['c'],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    edgecolor=fig['colorset']['marker']['edgecolor'],
                                    s=fig['colorset']['marker']['size']**2,
                                    alpha=fig['colorset']['marker']['alpha'],
                                    cmap=fig['colorset']['cmap']['cmap'],
                                    vmin=fig['ax']['lim']['c'][0],
                                    vmax=fig['ax']['lim']['c'][1]
                                ))
                    else:
                        fig['ax']['markercolor'] = marker[0].strip()
                        fig['ax']['markertype'] = marker[1].strip()
                        if len(marker) == 2:
                            ax.scatter(
                                fig['var']['data'].x,
                                fig['var']['data'].y,
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                c=fig['colorset']['scattercolor'][fig['ax']
                                                                  ['markercolor']],
                                edgecolor=fig['colorset']['marker']['edgecolor'],
                                s=fig['colorset']['marker']['size']**2,
                                alpha=fig['colorset']['marker']['alpha']
                            )
                        elif len(marker) > 2:
                            if marker[2].strip()[0:3] == '&Bo':
                                self.scatter_classify_data(
                                    fig, marker[2].strip()[4:])
                                ax.scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    c=fig['colorset']['scattercolor'][fig['ax']
                                                                      ['markercolor']],
                                    edgecolor=fig['colorset']['marker']['edgecolor'],
                                    s=fig['colorset']['marker']['size']**2,
                                    alpha=fig['colorset']['marker']['alpha']
                                )

                self.ax_setticks(fig, 'xyc')
                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                    if not self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                        ax.set_xticks(fig['ax']['ticks']['x'])
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                    ax.set_xscale('log')
                if self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                    ax.xaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['x'][0]))
                    ax.set_xticklabels(fig['ax']['ticks']['x'][1])
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.set_xlim(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1])

                if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        ax.set_yticks(fig['ax']['ticks']['y'])
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    ax.set_yscale('log')
                if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                    ax.yaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                    ax.set_yticklabels(fig['ax']['ticks']['y'][1])
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.set_ylim(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])

                if self.cf.has_option(fig['section'], 'c_scale'):
                    if self.cf.get(fig['section'], 'c_scale').strip().lower() == 'log':
                        from matplotlib.colors import LogNorm
                        plt.colorbar(fig['ax']['a1'][0], axc, norm=LogNorm(
                            fig['var']['lim']['c'][0], fig['var']['lim']['c'][1]), orientation='vertical', extend='neither')
                        axc.set_yscale("log")
                    elif self.cf.get(fig['section'], 'c_scale').strip().lower() == "flat":
                        plt.colorbar(fig['ax']['a1'][0], axc, ticks=fig['ax']
                                     ['ticks']['c'], orientation='vertical', extend='neither')
                # axc.set_ylim(fig['ax']['lim']['c'][0], fig['ax']['lim']['c'][1])

                ax.tick_params(
                    labelsize=fig['colorset']['ticks']['labelsize'],
                    direction=fig['colorset']['ticks']['direction'],
                    bottom=fig['colorset']['ticks']['bottom'],
                    left=fig['colorset']['ticks']['left'],
                    top=fig['colorset']['ticks']['top'],
                    right=fig['colorset']['ticks']['right'],
                    which='both'
                )
                ax.tick_params(which='major', length=fig['colorset']['ticks']
                               ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['ticks']
                               ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                axc.tick_params(
                    which="both",
                    labelsize=fig['colorset']['colorticks']['labelsize'],
                    direction=fig['colorset']['colorticks']['direction'],
                    bottom=fig['colorset']['colorticks']['bottom'],
                    left=fig['colorset']['colorticks']['left'],
                    top=fig['colorset']['colorticks']['top'],
                    right=fig['colorset']['colorticks']['right'],
                    color=fig['colorset']['colorticks']['color']
                )
                # axc.tick_params(which='major', length=fig['colorset']['ticks']['majorlength'], color=fig['colorset']['ticks']['majorcolor'])

                ax.set_xlabel(
                    r"{}".format(self.cf.get(fig['section'], 'x_label')), 
                    fontsize=fig['colorset']['canvas']['x_label']['size'],
                    color=fig['colorset']['canvas']['x_label']['color']                   
                )
                ax.xaxis.set_label_coords(0.5, fig['colorset']['canvas']['x_label']['offline'])
                ax.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], 'y_label')), 
                    fontsize=fig['colorset']['canvas']['yl_label']['size'],
                    color=fig['colorset']['canvas']['yl_label']['color'],
                )
                ax.yaxis.set_label_coords(fig['colorset']['canvas']['yl_label']['offline'], 0.5)
                axc.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], 'c_label')), 
                    fontsize=fig['colorset']['canvas']['yr_label']['size'],
                    color=fig['colorset']['canvas']['yr_label']['color'],
                )
                axc.yaxis.set_label_coords(fig['colorset']['canvas']['yr_label']['offline'], 0.5)

                if self.cf.has_option(fig['section'], 'Line_draw'):
                    self.drawline(fig, ax)

                if self.cf.has_option(fig['section'], "Text"):
                    self.drawtext(fig, ax)

                if self.cf.has_option(fig['section'], "fill"):
                    self.drawfillarea(fig, ax)

            elif fig['type'] == "1D_Stat":
                self.Get1DStatData(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                if 'x' in fig['var']['type'] and 'CHI2' in fig['var']['type'] and 'PDF' in fig['var']['type']:
                    ax = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axbox']['x0'], 
                        fig['colorset']['figureSize']['axbox']['y0'], 
                        fig['colorset']['figureSize']['axbox']['width'], 
                        fig['colorset']['figureSize']['axbox']['height'] 
                    ])
                    fig['var']['BestPoint'] = {
                        'x':    fig['var']['data'].loc[fig['var']['data'].CHI2.idxmin()].x,
                        'CHI2': fig['var']['data'].loc[fig['var']['data'].CHI2.idxmin()].CHI2,
                        'PDF':  fig['var']['data'].loc[fig['var']['data'].CHI2.idxmin()].PDF
                    }
                    fig['var']['data'] = fig['var']['data'].assign(
                        PL=lambda x: np.exp(-0.5 * x.CHI2)/np.exp(-0.5 * fig['var']['BestPoint']['CHI2']))
                    fig['var']['lim'] = {
                        'x':    [fig['var']['data'].x.min(), fig['var']['data'].x.max()]
                    }
                    fig['ax'] = {}
                    self.ax_setlim(fig, 'x')
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        XXGrid = np.linspace(
                            0, 1, int(self.cf.get(fig['section'], 'x_nbin'))+1)
                        XI = np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], int(
                            self.cf.get(fig['section'], 'x_nbin'))+1)
                        fig['ax']['var'] = {
                            'x':    XI,
                            'dx':   (fig['ax']['lim']['x'][1] - fig['ax']['lim']['x'][0])/int(self.cf.get(fig['section'], 'x_nbin'))
                        }
                        fig['ax']['grid'] = pd.DataFrame({
                            "xx":   XI
                        })
                        fig['ax']['grid'] = pd.DataFrame({
                            "xi":   fig['ax']['grid']['xx'],
                            "xxgrid":   XXGrid,
                            'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x >= tt['xx'] - 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].x < tt['xx'] + 0.5*fig['ax']['var']['dx'])].PL.max(axis=0, skipna=True), axis=1),
                            'PDF':  fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x >= tt['xx'] - 0.5*fig['ax']['var']['dx']) & (fig['var']['data'].x < tt['xx'] + 0.5*fig['ax']['var']['dx'])].PDF.sum(), axis=1)
                        }).fillna({"PL": 0.0, "PDF": 0.0})
                        from scipy.stats import gaussian_kde
                        # pdfkde = gaussian_kde(fig['ax']['grid'].xi, bw_method=0.03*fig['ax']['var']['dx'], weights=fig['ax']['grid'].PDF)
                        if self.cf.has_option(fig['section'], 'pdf_kde_bw'):
                            fig['ax']['pdf_kde_bw'] = float(
                                self.cf.get(fig['section'], 'pdf_kde_bw'))
                        else:
                            fig['ax']['pdf_kde_bw'] = 'silverman'
                        fig["ax"]["pdfkde"] = gaussian_kde(
                            fig['ax']['grid'].xxgrid, bw_method=fig['ax']['pdf_kde_bw'], weights=fig['ax']['grid'].PDF)
                        xgrid = np.linspace(0, 1, 10000)
                        fig['ax']['pdfkdedata'] = pd.DataFrame({
                            'xx':   np.linspace(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], 10000),
                            'pdf':  fig["ax"]["pdfkde"].evaluate(xgrid),
                        })
                        fig['var']['BestPoint']['PDF'] = fig['ax']['pdfkde'].evaluate(
                            fig['var']['BestPoint']['x'])/fig['ax']['pdfkdedata']['pdf'].max()
                        fig['ax']['pdfkdedata']['pdf'] = fig['ax']['pdfkdedata']['pdf'] / \
                            fig['ax']['pdfkdedata']['pdf'].max()
                        fig['ax']['pdfpara'] = {
                            'norm':                     sum(fig['ax']['pdfkdedata'].pdf),
                            '1sigma_critical_prob':     self.find_critical_prob(fig['ax']['pdfkdedata'], 0.675),
                            '2sigma_critical_prob':     self.find_critical_prob(fig['ax']['pdfkdedata'], 0.95),
                            'mode':                     fig['ax']['pdfkdedata'].iloc[fig['ax']['pdfkdedata'].pdf.idxmax()].xx
                        }
                        axpl = interp1d(
                            fig['ax']['grid']['xi'], fig['ax']['grid']['PL'], kind='linear')
                        plgrid = axpl(np.linspace(
                            fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1], 10000))
                    elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                        XXGrid = np.linspace(
                            0, 1, int(self.cf.get(fig['section'], 'x_nbin'))+1)
                        XI = np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), int(self.cf.get(fig['section'], 'x_nbin'))+1)
                        fig['ax']['var'] = {
                            'x':    XI,
                            'dx':   (math.log10(fig['ax']['lim']['x'][1]) - math.log10(fig['ax']['lim']['x'][0]))/int(self.cf.get(fig['section'], 'x_nbin'))
                        }
                        fig['ax']['grid'] = pd.DataFrame({
                            "xx":   XI
                        })
                        fig['ax']['grid'] = pd.DataFrame({
                            "xi":   fig['ax']['grid']['xx'],
                            "xxgrid":   XXGrid,
                            'PL':   fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x >= 10**(math.log10(tt['xx']) - 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].x < 10**(math.log10(tt['xx']) + 0.5*fig['ax']['var']['dx']))].PL.max(axis=0, skipna=True), axis=1),
                            'PDF':  fig['ax']['grid'].apply(lambda tt: fig['var']['data'][(fig['var']['data'].x >= 10**(math.log10(tt['xx']) - 0.5*fig['ax']['var']['dx'])) & (fig['var']['data'].x < 10**(math.log10(tt['xx']) + 0.5*fig['ax']['var']['dx']))].PDF.sum(), axis=1)
                        }).fillna({'PL': 0.0, 'PDF': 0.0})
                        # print(fig['ax']['grid'])
                        from scipy.stats import gaussian_kde
                        # pdfkde = gaussian_kde(fig['ax']['grid'].xi, bw_method=0.03*fig['ax']['var']['dx'], weights=fig['ax']['grid'].PDF)
                        if self.cf.has_option(fig['section'], 'pdf_kde_bw'):
                            fig['ax']['pdf_kde_bw'] = float(
                                self.cf.get(fig['section'], 'pdf_kde_bw'))
                        else:
                            fig['ax']['pdf_kde_bw'] = 'silverman'
                        fig["ax"]["pdfkde"] = gaussian_kde(
                            fig['ax']['grid'].xxgrid, bw_method=fig['ax']['pdf_kde_bw'], weights=fig['ax']['grid'].PDF)
                        xgrid = np.linspace(0, 1, 10000)
                        fig['ax']['pdfkdedata'] = pd.DataFrame({
                            'xx':   np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(fig['ax']['lim']['x'][1]), 10000),
                            'pdf':  fig["ax"]["pdfkde"].evaluate(xgrid),
                        })

                        fig['var']['BestPoint']['PDF'] = fig['ax']['pdfkde'].evaluate(
                            fig['var']['BestPoint']['x'])/fig['ax']['pdfkdedata']['pdf'].max()
                        fig['ax']['pdfkdedata']['pdf'] = fig['ax']['pdfkdedata']['pdf'] / \
                            fig['ax']['pdfkdedata']['pdf'].max()
                        fig['ax']['pdfpara'] = {
                            'norm':                     sum(fig['ax']['pdfkdedata'].pdf),
                            '1sigma_critical_prob':     self.find_critical_prob(fig['ax']['pdfkdedata'], 0.675),
                            '2sigma_critical_prob':     self.find_critical_prob(fig['ax']['pdfkdedata'], 0.95),
                            'mode':                     fig['ax']['pdfkdedata'].iloc[fig['ax']['pdfkdedata'].pdf.idxmax()].xx
                        }
                        axpl = interp1d(
                            fig['ax']['grid']['xi'], fig['ax']['grid']['PL'], kind='linear')
                        plgrid = axpl(np.logspace(math.log10(fig['ax']['lim']['x'][0]), math.log10(
                            fig['ax']['lim']['x'][1]), 10000))
                        ax.set_xscale("log")
                    self.ax_setcmap(fig)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], -0.09, -0.12, where=fig['ax']['pdfkdedata']['pdf'] > fig['ax']['pdfpara']['2sigma_critical_prob'], edgecolor=fig['colorset']
                                    ['1dpdf']['2sigma']['edge'], facecolor=fig['colorset']['1dpdf']['2sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['2sigma']['alpha'], zorder=5)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], -0.09, -0.12, where=fig['ax']['pdfkdedata']['pdf'] > fig['ax']['pdfpara']['1sigma_critical_prob'], edgecolor=fig['colorset']
                                    ['1dpdf']['1sigma']['edge'], facecolor=fig['colorset']['1dpdf']['1sigma']['facecolor'], alpha=fig['colorset']['1dpdf']['1sigma']['alpha'], zorder=6)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], 0, fig['ax']['pdfkdedata']['pdf'],  where=fig['ax']['pdfkdedata']['pdf'] > fig['ax']['pdfpara']['2sigma_critical_prob'],
                                    edgecolor=fig['colorset']['1dpdf']['2sigma']['edge'], facecolor=fig['colorset']['1dpdf']['2sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['2sigma']['alpha'], zorder=5)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], 0, fig['ax']['pdfkdedata']['pdf'],  where=fig['ax']['pdfkdedata']['pdf'] > fig['ax']['pdfpara']['1sigma_critical_prob'],
                                    edgecolor=fig['colorset']['1dpdf']['1sigma']['edge'], facecolor=fig['colorset']['1dpdf']['1sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['1sigma']['alpha'], zorder=6)
                    ax.set_xlim(fig['ax']['lim']['x'][0],
                                fig['ax']['lim']['x'][1])
                    ax.plot([fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1]], [
                            0, 0], '-', color='grey', linewidth=0.8, zorder=0)
                    ax.plot([fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1]], [
                            0.1354, 0.1354], '--', color='#ea4702', linewidth=1.2, zorder=0)
                    ax.text(
                        0.025, 0.237, r"$95.4\%~{\rm C.L.}$", fontsize=12, transform=ax.transAxes)
                    ax.plot([fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1]], [
                            0.60583, 0.60583], '--', color='#ea4702', linewidth=1.2, zorder=0)
                    ax.text(
                        0.025, 0.595, r"$68.3\%~{\rm C.L.}$", fontsize=12, transform=ax.transAxes)
                    ax.set_ylim(-0.16, 1.15)
                    if self.cf.has_option(fig['section'], "BestPoint"):
                        for item in fig['colorset']['1dpdf']['bestpoint']:
                            ax.plot([fig['var']['BestPoint']['x'], fig['var']['BestPoint']['x']], [-0.069, -0.041],
                                    '-', linewidth=item['width'], color=item['color'], alpha=item['alpha'], zorder=7)
                        for item in fig['colorset']['1dpdf']['pdfmode']:
                            ax.plot([fig['ax']['pdfpara']['mode'], fig['ax']['pdfpara']['mode']], [
                                    0.001, 0.999], '-', linewidth=item['width'], color=item['color'], alpha=item['alpha'], zorder=7)
                            ax.plot([fig['ax']['pdfpara']['mode'], fig['ax']['pdfpara']['mode']], [-0.09, -0.12],
                                    '-', linewidth=item['width'], color=item['color'], alpha=item['alpha'], zorder=7)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], 0, fig['ax']['pdfkdedata']['pdf'],
                                    color="w", alpha=fig['colorset']['1dpdf']['line']['alpha'], zorder=1)
                    ax.plot(fig['ax']['pdfkdedata']['xx'], fig['ax']['pdfkdedata']['pdf'], '-', linewidth=fig['colorset']['1dpdf']['line']
                            ['width'], color=fig['colorset']['1dpdf']['line']['color'], alpha=fig['colorset']['1dpdf']['line']['alpha'], zorder=10)
                    ax.fill_between(fig['ax']['grid']['xi'], 0, fig['ax']['grid']
                                    ['PL'], color='red', step='mid', alpha=0.1, zorder=0)
                    ax.step(fig['ax']['grid']['xi'], fig['ax']['grid']
                            ['PL'], color='red', where='mid', zorder=9)
                    # axpl = interp1d(fig['ax']['grid']['xi'], fig['ax']['grid']['PL'], kind='linear')
                    # plgrid = axpl(xgrid)
                    # ax.fill_between( fig['ax']['grid']['xi'], -0.07, -0.04, where=fig['ax']['grid']['PL']>0.1354, color='#f8a501', step='mid', alpha=0.8 )
                    # ax.fill_between( fig['ax']['grid']['xi'], -0.07, -0.04, where=fig['ax']['grid']['PL']>0.60583, color='#10a6e7', step='mid', alpha=0.7 )
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], -0.07, -0.04,
                                    where=plgrid > 0.1354, color='#f8a501', step='mid', alpha=0.8)
                    ax.fill_between(fig['ax']['pdfkdedata']['xx'], -0.07, -0.04,
                                    where=plgrid > 0.60583, color='#10a6e7', step='mid', alpha=0.7)
                self.ax_setticks(fig, 'x')

                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                    if not self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                        ax.set_xticks(fig['ax']['ticks']['x'])
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                    ax.set_xscale('log')
                if self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                    ax.xaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['x'][0]))
                    ax.set_xticklabels(fig['ax']['ticks']['x'][1])
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.set_xlim(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1])
                ax.yaxis.set_minor_locator(FixedLocator(np.linspace(0, 1, 26)))

                ax.tick_params(labelsize=fig['colorset']['ticks']['labelsize'],  direction=fig['colorset']['ticks']['direction'],  bottom=fig['colorset']['ticks']
                               ['bottom'],  left=fig['colorset']['ticks']['left'],  top=fig['colorset']['ticks']['top'],  right=fig['colorset']['ticks']['right'], which='both')
                ax.tick_params(which='major', length=fig['colorset']['ticks']
                               ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['ticks']
                               ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                ax.set_xlabel(
                    r"{}".format(self.cf.get(fig['section'], 'x_label')), 
                    fontsize=fig['colorset']['canvas']['x_label']['size'],
                    color=fig['colorset']['canvas']['x_label']['color']                   
                )
                ax.xaxis.set_label_coords(0.5, fig['colorset']['canvas']['x_label']['offline'])
                ax.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], 'chi2_label')), 
                    fontsize=fig['colorset']['canvas']['yl_label']['size'],
                    color=fig['colorset']['canvas']['yl_label']['color'],
                )
                ax.yaxis.set_label_coords(fig['colorset']['canvas']['yl_label']['offline'], 0.5)
                plax = ax.secondary_yaxis(1-fig['colorset']['canvas']['yr_label']['offline'], color="None")
                plax.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], "pdf_label")), 
                    fontsize=fig['colorset']['canvas']['yr_label']['size'], 
                    color=fig['colorset']['canvas']['yr_label']['color']
                )
                # plax.set_label_coords(fig['colorset']['canvas']['yr_label']['offline'], 0.5)
                plax.yaxis.set_major_locator(FixedLocator([]))
                plax.tick_params(labelsize=3, labelcolor="None",
                                 direction=fig['colorset']['ticks']['direction'])

                if self.cf.has_option(fig['section'], 'Line_draw'):
                    self.drawline(fig, ax)

                if self.cf.has_option(fig['section'], "Text"):
                    self.drawtext(fig, ax)

                if self.cf.has_option(fig['section'], "drawlegend"):
                    if eval(self.cf.get(fig['section'], "drawlegend")):
                        print("\tTimer: {:.2f} Second;  Drawing default format legend ...".format(
                            time.time()-fig['start']))
                        axlegend = fig['fig'].add_axes(
                            [0.14, 0.875, 0.72, 0.07])
                        axlegend.spines['top'].set_visible(False)
                        axlegend.spines['bottom'].set_visible(False)
                        axlegend.spines['left'].set_visible(False)
                        axlegend.spines['right'].set_visible(False)
                        from matplotlib.ticker import NullLocator
                        axlegend.xaxis.set_major_locator(NullLocator())
                        axlegend.yaxis.set_major_locator(NullLocator())
                        axlegend.set_xlim(0, 100)
                        axlegend.set_ylim(-2.5, 2)
                        axlegend.step([2, 3, 5, 7, 9, 10], [
                                      0.5, 0.5, 1.3, 1.5, 0.2, 0.2], color='red', where='mid', zorder=0)
                        axlegend.fill_between([2, 3, 5, 7, 9, 10], 0, [
                                              0.5, 0.5, 1.3, 1.5, 0.2, 0.2], color='red', step='mid', alpha=0.1)
                        axlegend.text(
                            12, 0.2, "Profile Likelihood", fontsize=11)
                        axlegend.fill_between(
                            [2, 10], -0.7, -2, color='#f8a501', step='mid', alpha=0.8)
                        axlegend.fill_between(
                            [4, 8], -0.7, -2,  color='#10a6e7', step='mid', alpha=0.7)
                        axlegend.text(
                            12, -1.9, "Best-fit point, $1\sigma$ & $2\sigma$ confidence interval", fontsize=11)
                        legendxi = np.linspace(52, 60, 100)
                        legendyy = 1.5 * \
                            np.exp(-(legendxi-55.0)**2/2) + 0.6 * \
                            np.exp(-(legendxi-58)**2 / 1.5)
                        axlegend.plot(legendxi, legendyy, '-', linewidth=fig['colorset']['1dpdf']['line']['width'],
                                      color=fig['colorset']['1dpdf']['line']['color'], alpha=fig['colorset']['1dpdf']['line']['alpha'], zorder=10)
                        axlegend.fill_between(
                            legendxi, 0, legendyy, color='w', zorder=1)
                        axlegend.fill_between(legendxi, 0, legendyy, where=legendyy > 0.23, edgecolor=fig['colorset']['1dpdf']['2sigma']['edge'], facecolor=fig[
                                              'colorset']['1dpdf']['2sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['2sigma']['alpha'], zorder=5)
                        axlegend.fill_between(legendxi, 0, legendyy, where=legendyy > 0.72, edgecolor=fig['colorset']['1dpdf']['1sigma']['edge'], facecolor=fig[
                                              'colorset']['1dpdf']['1sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['1sigma']['alpha'], zorder=6)
                        axlegend.fill_between(legendxi, -0.7, -2, where=legendyy > 0.23, edgecolor=fig['colorset']['1dpdf']['2sigma']['edge'],
                                              facecolor=fig['colorset']['1dpdf']['2sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['2sigma']['alpha'], zorder=5)
                        axlegend.fill_between(legendxi, -0.7, -2, where=legendyy > 0.72, edgecolor=fig['colorset']['1dpdf']['1sigma']['edge'],
                                              facecolor=fig['colorset']['1dpdf']['1sigma']['facecolor'],  alpha=fig['colorset']['1dpdf']['1sigma']['alpha'], zorder=6)
                        if self.cf.has_option(fig['section'], "BestPoint"):
                            for item in fig['colorset']['1dpdf']['pdfmode']:
                                axlegend.plot([55.0, 55.0], [
                                              0.05, 1.45], '-', linewidth=item['width'], color=item['color'], alpha=item['alpha'], zorder=7)
                                axlegend.plot([55.0, 55.0], [-0.8, -1.95], '-', linewidth=item['width'],
                                              color=item['color'], alpha=item['alpha'], zorder=7)
                            for item in fig['colorset']['1dpdf']['bestpoint']:
                                axlegend.plot([7.0, 7.0], [-0.8, -1.95], '-', linewidth=item['width'],
                                              color=item['color'], alpha=item['alpha'], zorder=7)
                        axlegend.text(
                            62, 0.2, "Posterior PDF", fontsize=11)
                        axlegend.text(
                            62, -1.9, "PDF mode, $1\sigma$ & $2\sigma$ credible region", fontsize=11)

            elif fig['type'] == "TernaryRGB_Scatter":
                self.load_ternary_configure(fig)
                self.makecanvas(fig)
                self.init_ternary(fig, 'axt')
                self.init_ternary(fig, 'axr')
                self.init_ternary(fig, 'axg')
                self.init_ternary(fig, 'axb')
                self.init_ternary(fig, 'axl')
                self.draw_ternary_legend(fig, 'axl')
                self.mark_ternary_channel(fig, 'axr', 'r')
                self.mark_ternary_channel(fig, 'axg', 'g')
                self.mark_ternary_channel(fig, 'axb', 'b')
                self.getTernaryRGBData(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                fig['ax']['axt'].scatter(fig['var']['axdata']['x'], fig['var']['axdata']['y'],
                                         marker='^', color=fig['var']['axdata']['c'], s=6, zorder=1, alpha=1)

                fig['ax']['axr'].scatter(fig['var']['axdata']['x'], fig['var']['axdata']['y'],
                                         marker='^', color=fig['var']['axdata']['r'], s=1.5, zorder=1, alpha=0.8)
                fig['ax']['axg'].scatter(fig['var']['axdata']['x'], fig['var']['axdata']['y'],
                                         marker='^', color=fig['var']['axdata']['g'], s=1.5, zorder=1, alpha=0.8)
                fig['ax']['axb'].scatter(fig['var']['axdata']['x'], fig['var']['axdata']['y'],
                                         marker='^', color=fig['var']['axdata']['b'], s=1.5, zorder=1, alpha=1)

                self.set_Ternary_label(fig, 'axt')
            
            elif fig['type'] == "Ternary_Scatter":
                self.load_ternary_configure(fig)
                self.makecanvas(fig)
                self.init_ternary(fig, 'axt')
                self.ax_setcmap(fig)
                self.getTernaryData(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                lines = self.cf.get(fig['section'], "marker").split("\n")
                with open(self.load_path(fig['colorset']['markercodepath']), 'r') as f1:
                    fig['colorset']['scattermarker'] = json.loads(f1.read())
                with open(self.load_path(fig['colorset']["colorpath"]), 'r') as f1:
                    fig['colorset']['scattercolor'] = json.loads(f1.read())
                for line in lines:
                    marker = line.split(',')
                    if len(marker) == 2:
                        fig['ax']['markercolor'] = marker[0].strip()
                        fig['ax']['markertype'] = marker[1].strip()
                        fig['ax']['axt'].scatter(
                            fig['var']['axdata']['x'],
                            fig['var']['axdata']['y'],
                            marker=fig['colorset']['scattermarker'][fig['ax']
                                                                    ['markertype']],
                            color=fig['colorset']['scattercolor'][fig['ax']
                                                                  ['markercolor']],
                            s=fig['colorset']['marker']['size'],
                            zorder=1, alpha=0.8
                        )
                    elif len(marker) >= 3:
                        if marker[2].strip()[0:3] == "&Bo":
                            marker_selection = ','.join(marker[2:]).strip()[4:]
                            fig['ax']['markercolor'] = marker[0].strip()
                            fig['ax']['markertype'] = marker[1].strip()
                            self.Ternary_classify_data(fig, marker_selection)
                            fig['ax']['axt'].scatter(
                                fig['classify']['x'],
                                fig['classify']['y'],
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                color=fig['colorset']['scattercolor'][fig['ax']
                                                                      ['markercolor']],
                                s=fig['colorset']['marker']['size'],
                                zorder=1, alpha=0.8
                            )
                self.set_Ternary_label(fig, 'axt')
            
            elif fig['type'] == "TernaryC_Scatter":
                self.load_ternary_configure(fig)
                axc = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axcbox']['x0'], 
                        fig['colorset']['figureSize']['axcbox']['y0'], 
                        fig['colorset']['figureSize']['axcbox']['width'], 
                        fig['colorset']['figureSize']['axcbox']['height'] 
                    ])
                self.makecanvas(fig)
                self.init_ternary(fig, "axt")
                self.ax_setcmap(fig)
                self.getTernaryCData(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                fig['var']['lim'] = {
                    'c':    [fig['var']['color'].min(), fig['var']['color'].max()]
                }
                self.ax_setlim(fig, "c")
                lines = self.cf.get(fig['section'], "marker").split('\n')
                fig['ax']['a1'] = []
                with open(self.load_path(fig['colorset']['markercodepath']), 'r') as f1:
                    fig['colorset']['scattermarker'] = json.loads(f1.read())
                # with open(self.load_path(fig['colorset']["colorpath"]), 'r') as f1:
                    # fig['colorset']['scattercolor']  = json.loads(f1.read())
                for line in lines:
                    marker = line.split(",")
                    if marker[0].strip() == "&Color":
                        fig['ax']['markertype'] = marker[1].strip()
                        if len(marker) == 2:
                            fig['ax']['a1'].append(fig['ax']['axt'].scatter(
                                fig['var']['axdata']['x'],
                                fig['var']['axdata']['y'],
                                c=fig['var']['axdata'].c,
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                s=fig['colorset']['marker']['size'],
                                alpha=fig['colorset']['marker']['alpha'],
                                cmap=fig['colorset']['cmap']['cmap'],
                                vmin=fig['ax']['lim']['c'][0],
                                vmax=fig['ax']['lim']['c'][1]
                            ))
                        elif len(marker) >= 3:
                            if marker[2].strip()[0:3] == "&Bo":
                                self.Ternary_classify_data(
                                    fig, marker[2].strip()[4:])
                                fig['ax']['a1'].append(fig['ax']['axt'].scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    c=fig['classify']['c'],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    s=fig['colorset']['marker']['size'],
                                    alpha=fig['colorset']['marker']['alpha'],
                                    cmap=fig['colorset']['cmap']['cmap'],
                                    vmin=fig['ax']['lim']['c'][0],
                                    vmax=fig['ax']['lim']['c'][1]
                                ))
                    else:
                        fig['ax']['markertype'] = marker[1].strip()
                        fig['ax']['markercolor'] = marker[0].strip()
                        if len(marker) == 2:
                            fig['ax']['axt'].scatter(
                                fig['var']['axdata']['x'],
                                fig['var']['axdata']['y'],
                                c=fig['colorset']['scattercolor'][fig['ax']
                                                                  ['markercolor']],
                                marker=fig['colorset']['scattermarker'][fig['ax']
                                                                        ['markertype']],
                                s=fig['colorset']['marker']['size'],
                                alpha=fig['colorset']['marker']['alpha'],
                            )
                        elif len(marker) >= 3:
                            if marker[2].strip()[0:3] == "&Bo":
                                self.Ternary_classify_data(
                                    fig, marker[2].strip()[4:])
                                fig['ax']['axt'].scatter(
                                    fig['classify']['x'],
                                    fig['classify']['y'],
                                    c=fig['colorset']['scattercolor'][fig['ax']
                                                                      ['markercolor']],
                                    marker=fig['colorset']['scattermarker'][fig['ax']
                                                                            ['markertype']],
                                    s=fig['colorset']['marker']['size'],
                                    alpha=fig['colorset']['marker']['alpha'],
                                )

                self.ax_setticks(fig, 'c')
                if self.cf.has_option(fig['section'], 'c_scale'):
                    if self.cf.get(fig['section'], 'c_scale').strip().lower() == 'log':
                        from matplotlib.colors import LogNorm
                        plt.colorbar(fig['ax']['a1'][0], axc, norm=LogNorm(
                            fig['ax']['lim']['c'][0], fig['ax']['lim']['c'][1]), orientation='vertical', extend='neither')
                        axc.set_yscale('log')
                    elif self.cf.get(fig['section'], 'c_scale').strip().lower() == "flat":
                        plt.colorbar(fig['ax']['a1'][0], axc, ticks=fig['ax']
                                     ['ticks']['c'], orientation='vertical', extend='neither')

                axc.set_ylabel(r"{}".format(self.cf.get(
                    fig['section'], 'c_label')), fontsize=24)
                axc.yaxis.set_label_coords(2, 0.5)
                axc.tick_params(
                    axis='y',
                    labelleft=True,
                    which="both",
                    labelsize=fig['colorset']['colorticks']['labelsize'],
                    direction=fig['colorset']['colorticks']['direction'],
                    bottom=fig['colorset']['colorticks']['bottom'],
                    left=fig['colorset']['colorticks']['left'],
                    right=fig['colorset']['colorticks']['right'],
                    top=fig['colorset']['colorticks']['top'],
                    color=fig['colorset']['colorticks']['color'],
                    labelright=False
                )
                self.set_Ternary_label(fig, 'axt')
            
            elif fig['type'] == "Violin":
                self.get_violin_data(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                self.ax_setcmap(fig)
                self.set_violin_canvas(fig)
                self.ax_setlim(fig, 'y')
                fig['ax']['data'] = []
                fig['ax']['label'] = []
                for var in fig['varlist']:
                    fig['ax']['data'].append(fig['var'][var['name']])
                    fig['ax']['label'].append(var['label'])
                parts = fig['ax']['ax'].violinplot(
                    fig['ax']['data'],
                    widths=fig['colorset']['violin']['width'],
                    showmeans=False,
                    showmedians=False,
                    showextrema=False
                )
                for pc in parts['bodies']:
                    pc.set_facecolor(fig['colorset']['violin']['facecolor'])
                    pc.set_edgecolor(fig['colorset']['violin']['edgecolor'])
                    pc.set_alpha(fig['colorset']['violin']['alpha'])
                    pc.set_linewidth(fig['colorset']['violin']['linewidth'])
                if fig['colorset']['colorfulSchema']:
                    ccd = 0.
                    for pc in parts['bodies']:
                        pc.set_facecolor(fig['colorset']['cmap']['cmap'](ccd)[0:3])
                        ccd += 1./(len(parts['bodies'])-1)

                inds = np.arange(1, len(fig['varlist']) + 1)
                quartile1, quartile2, medians, quartile3, quartile4 = np.percentile(
                    fig['ax']['data'], [2.5, 16, 50, 84, 97.5], axis=1)
                fig['ax']['ax'].vlines(
                    inds, quartile2, quartile3,
                    color=fig['colorset']['violin']['1slcolor'],
                    lw=fig['colorset']['violin']['1slwidth'],
                    linestyle='-',
                    zorder=19
                )
                fig['ax']['ax'].vlines(
                    inds, quartile1, quartile4,
                    color=fig['colorset']['violin']['2slcolor'],
                    lw=fig['colorset']['violin']['2slwidth'],
                    linestyle='-',
                    zorder=18
                )
                fig['ax']['ax'].scatter(
                    inds, medians,
                    marker=fig['colorset']['violin']['mdmarker'],
                    color=fig['colorset']['violin']['mdcolor'],
                    s=fig['colorset']['violin']['mdsize'],
                    zorder=20
                )

                fig['ax']['ax'].set_ylim(
                    fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])
                fig['ax']['ax'].set_xlim(0.5 + fig['colorset']['violin']['xmin'],
                                         (len(fig["varlist"])+0.5-fig['colorset']['violin']['xmin']))
                fig['ax']['ax'].set_xticks(
                    np.arange(1, len(fig['ax']['label']) + 1))
                fig['ax']['ax'].set_xticklabels(fig['ax']['label'])

                fig['ax']['ax'].spines['bottom'].set_color(
                    fig['colorset']['violin']['spinesc'])
                fig['ax']['ax'].spines['left'].set_color(
                    fig['colorset']['violin']['spinesc'])
                fig['ax']['ax'].spines['top'].set_color(
                    fig['colorset']['violin']['spinesc'])
                fig['ax']['ax'].spines['right'].set_color(
                    fig['colorset']['violin']['spinesc'])

                fig['ax']['ax'].tick_params(
                    axis="y",
                    labelsize=fig['colorset']['ticks']['labelsize'],
                    direction=fig['colorset']['ticks']['direction'],
                    bottom=fig['colorset']['ticks']['bottom'],
                    left=fig['colorset']['ticks']['left'],
                    top=fig['colorset']['ticks']['top'],
                    right=fig['colorset']['ticks']['right'],
                    which='both'
                )
                fig['ax']['ax'].tick_params(which='major', length=fig['colorset']['ticks']
                                            ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                fig['ax']['ax'].tick_params(which='minor', length=fig['colorset']['ticks']
                                            ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                fig['ax']['ax'].tick_params(
                    axis='x',
                    labelsize=24,
                    bottom=False,
                    pad=12
                )

                self.ax_setticks(fig, 'y')
                if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        fig['ax']['ax'].set_yticks(fig['ax']['ticks']['y'])
                        fig['ax']['ax'].yaxis.set_minor_locator(
                            AutoMinorLocator())
                elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    fig['ax']['ax'].set_yscale('log')
                if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                    fig['ax']['ax'].yaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                    fig['ax']['ax'].set_yticklabels(fig['ax']['ticks']['y'][1])
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        fig['ax']['ax'].yaxis.set_minor_locator(
                            AutoMinorLocator())
                fig['ax']['ax'].set_ylim(
                    fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])
                fig['ax']['ax'].set_ylabel(r"{}".format(
                    self.cf.get(fig['section'], 'y_label')), fontsize=30)

            elif fig['type'] == "BiViolin":
                # print(self.data)
                if "right_selection" in self.cf.options(fig['section']) and "left_selection" in self.cf.options(fig['section']):
                    fig['Biviolin'] = True
                else:
                    fig['Biviolin'] = False
                if fig['Biviolin']:
                    self.get_violin_data(fig)
                    self.set_violin_canvas(fig)                    
                    self.ax_setcmap(fig)
                    print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                        time.time()-fig['start'], fig['section']), use_aliases=True))
                    self.ax_setlim(fig, 'y')
                    for var in fig['varlist']:
                        fig['ax']['data']['sum'].append(fig['var'][var['name']])
                        fig['ax']['data']['label'].append(var['label'])
                    parts1 = fig['ax']['ax'].violinplot(
                        fig['ax']['data']['left'],
                        widths=fig['colorset']['violin']['width'],
                        showmeans=False,
                        showmedians=False,
                        showextrema=False
                    )
                    parts2 = fig['ax']['ax'].violinplot(
                        fig['ax']['data']['right'],
                        widths=fig['colorset']['violin']['width'],
                        showmeans=False,
                        showmedians=False,
                        showextrema=False
                    )
                    for pc in parts1['bodies']:
                        m = np.mean(pc.get_paths()[0].vertices[:,0])
                        pc.get_paths()[0].vertices[:, 0] = np.clip(pc.get_paths()[0].vertices[:, 0], -np.inf, m)
                        pc.set_facecolor(fig['colorset']['violin']['BiSidefacecolor'])
                        pc.set_edgecolor(fig['colorset']['violin']['BiSideedgecolor'])
                        pc.set_alpha(fig['colorset']['violin']['alpha'])
                        pc.set_linewidth(fig['colorset']['violin']['linewidth'])
                    for pc in parts2['bodies']:
                        m = np.mean(pc.get_paths()[0].vertices[:,0])
                        pc.get_paths()[0].vertices[:, 0] = np.clip(pc.get_paths()[0].vertices[:, 0], m, np.inf)
                        pc.set_facecolor(fig['colorset']['violin']['BiSidefacecolor'])
                        pc.set_edgecolor(fig['colorset']['violin']['BiSideedgecolor'])
                        pc.set_alpha(fig['colorset']['violin']['alpha'])
                        pc.set_linewidth(fig['colorset']['violin']['linewidth'])

                    if fig['colorset']['violin']['BiSidecolorfulface'] == "left":
                        if fig['colorset']['colorfulSchema']:
                            ccd = 0.
                            for pc in parts1['bodies']:
                                pc.set_facecolor(fig['colorset']['cmap']['cmap'](ccd)[0:3])
                                pc.set_edgecolor(fig['colorset']['violin']['edgecolor'])
                                ccd += 1./(len(parts1['bodies'])-1)
                                pc.set_zorder(17)
                    elif fig['colorset']['violin']['BiSidecolorfulface'] == "right":
                        if fig['colorset']['colorfulSchema']:
                            ccd = 0.
                            for pc in parts2['bodies']:
                                pc.set_facecolor(fig['colorset']['cmap']['cmap'](ccd)[0:3])
                                pc.set_edgecolor(fig['colorset']['violin']['edgecolor'])
                                ccd += 1./(len(parts2['bodies'])-1)
                                pc.set_zorder(17)

                    inds = np.arange(1, len(fig['varlist']) + 1)
                    quartile1, quartile2, medians, quartile3, quartile4 = np.percentile(
                        fig['ax']['data']['sum'], [2.5, 16, 50, 84, 97.5], axis=1)
                        
                    fig['ax']['ax'].vlines(
                        inds, quartile2, quartile3,
                        color=fig['colorset']['violin']['1slcolor'],
                        lw=fig['colorset']['violin']['1slwidth'],
                        linestyle='-',
                        zorder=19
                    )
                    fig['ax']['ax'].vlines(
                        inds, quartile1, quartile4,
                        color=fig['colorset']['violin']['2slcolor'],
                        lw=fig['colorset']['violin']['2slwidth'],
                        linestyle='-',
                        zorder=18
                    )
                    fig['ax']['ax'].scatter(
                        inds, medians,
                        marker=fig['colorset']['violin']['mdmarker'],
                        color=fig['colorset']['violin']['mdcolor'],
                        s=fig['colorset']['violin']['mdsize'],
                        zorder=20
                    )

                    fig['ax']['ax'].set_ylim(
                        fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])
                    fig['ax']['ax'].set_xlim(0.5 + fig['colorset']['violin']['xmin'],
                                             (len(fig["varlist"])+0.5-fig['colorset']['violin']['xmin']))
                    fig['ax']['ax'].set_xticks(
                        np.arange(1, len(fig['ax']['data']['label']) + 1))
                    fig['ax']['ax'].set_xticklabels(fig['ax']['data']['label'])

                    fig['ax']['ax'].spines['bottom'].set_color(
                        fig['colorset']['violin']['spinesc'])
                    fig['ax']['ax'].spines['left'].set_color(
                        fig['colorset']['violin']['spinesc'])
                    fig['ax']['ax'].spines['top'].set_color(
                        fig['colorset']['violin']['spinesc'])
                    fig['ax']['ax'].spines['right'].set_color(
                        fig['colorset']['violin']['spinesc'])

                    fig['ax']['ax'].tick_params(
                        axis="y",
                        labelsize=fig['colorset']['ticks']['labelsize'],
                        direction=fig['colorset']['ticks']['direction'],
                        bottom=fig['colorset']['ticks']['bottom'],
                        left=fig['colorset']['ticks']['left'],
                        top=fig['colorset']['ticks']['top'],
                        right=fig['colorset']['ticks']['right'],
                        which='both'
                    )
                    fig['ax']['ax'].tick_params(which='major', length=fig['colorset']['ticks']
                                                ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                    fig['ax']['ax'].tick_params(which='minor', length=fig['colorset']['ticks']
                                                ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                    fig['ax']['ax'].tick_params(
                        axis='x',
                        labelsize=24,
                        bottom=False,
                        pad=12
                    )

                    self.ax_setticks(fig, 'y')
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                            fig['ax']['ax'].set_yticks(fig['ax']['ticks']['y'])
                            fig['ax']['ax'].yaxis.set_minor_locator(
                                AutoMinorLocator())
                    elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                        fig['ax']['ax'].set_yscale('log')
                    if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        fig['ax']['ax'].yaxis.set_major_locator(
                            ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                        fig['ax']['ax'].set_yticklabels(fig['ax']['ticks']['y'][1])
                        if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                            fig['ax']['ax'].yaxis.set_minor_locator(
                                AutoMinorLocator())
                    fig['ax']['ax'].set_ylim(
                        fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])
                    fig['ax']['ax'].set_ylabel(r"{}".format(
                        self.cf.get(fig['section'], 'y_label')), fontsize=30)

            elif fig['type'] == "Decay":
                self.get_DK_data(fig)
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                self.make_decays_canvas(fig)
                self.make_decays_title(fig)
                for axis in fig['ax']['axs']:
                    for ii, pp in enumerate(fig['DKs']['data'][axis['name']]):
                        axis['1dd'].text(
                            fig['colorset']['label']['x'], axis['ymax']-ii-0.95, pp['label'], 
                            fontsize=fig['colorset']['label']['fontsize'],
                            color=fig['colorset']['label']['color'],
                            horizontalalignment="left",
                            verticalalignment="bottom"
                        )
                        barwid = pp['data'][pp['data'] < fig['colorset']['canvas']['xlim'][0]].count() / fig['data'].shape[0]
                        axis['bar'].text(
                            1-barwid/1.5 -0.01, axis['ymax']-ii-0.95, "{:.0%}".format(barwid),
                            fontsize=fig['colorset']['bar']['labelsize'],
                            fontfamily = "monospace",
                            color=fig['colorset']['bar']['labelcolor'],
                            horizontalalignment="right",
                            verticalalignment="bottom"
                        )
                        onedd = fig['colorset']['1dd']
                        grid, edge = np.histogram(pp['data'], bins=fig['colorset']['1dd']['bin']-1, range=(fig['colorset']['canvas']['xlim'][0], fig['colorset']['canvas']['xlim'][1]), density=False)
                        dx = (fig['colorset']['canvas']['xlim'][1] - fig['colorset']['canvas']['xlim'][0])/(fig['colorset']['1dd']['bin'] - 1)
                        xxgrid = np.linspace(fig['colorset']['canvas']['xlim'][0]+0.5*dx, fig['colorset']['canvas']['xlim'][1]-0.5*dx, fig['colorset']['1dd']['bin'] - 1)
                        from scipy.stats import gaussian_kde
                        from scipy.interpolate import Rbf
                        bw = fig['colorset']['1dd']['kde_bw']
                        if max(grid) > 0:
                            axis['bar'].fill(
                                [1-barwid/1.5, 1, 1, 1-barwid/1.5], [axis['ymax']-ii-1, axis['ymax']-ii-1, axis['ymax']-ii-0.1, axis['ymax']-ii-0.1], 
                                facecolor=fig['colorset']['cmap']['cmap'](1 - ii/axis['ymax'])[0:3]
                            )
                            grid = grid / max(grid) 
                            xx = np.linspace(fig['colorset']['canvas']['xlim'][0], fig['colorset']['canvas']['xlim'][1], 1000)
                            try:
                                from Func_lab import block_screen_print
                                block_screen_print()
                                grdkde = gaussian_kde(xxgrid, bw_method=bw, weights=grid)
                                pp['pdf'] = grdkde.evaluate(xx)
                                pp['pdf'] = pp['pdf'] / pp['pdf'].max() * fig['colorset']['figureSize']['1dheight']
                                axis['1dd'].fill_between(
                                    xx, pp['pdf'] + axis['ymax']-ii-1, axis['ymax'] - ii-1,
                                    where=pp['pdf']>0,
                                    facecolor=fig['colorset']['cmap']['cmap'](1 - ii/axis['ymax'])[0:3],
                                    alpha=0.95,
                                    zorder = 100+ii
                                )
                                from Func_lab import enable_screen_print
                                enable_screen_print()
                                axis['1dd'].plot(
                                    xx, pp['pdf'] + axis['ymax']-ii-1, 
                                    '-', 
                                    color=fig['colorset']['1dd']['color'],
                                    linewidth=fig['colorset']['1dd']['linewidth'],
                                    zorder = 100 + ii
                                )
                            except:
                                grdkde = Rbf(
                                    xxgrid, grid
                                )
                                pp['pdf']=grdkde(xx)
                                pp['pdf'] = pp['pdf'] / pp['pdf'].max() * fig['colorset']['figureSize']['1dheight']
                                pp['pdf'] = np.where(pp['pdf'] > 0, pp['pdf'], 0.0*pp['pdf'])
                                axis['1dd'].fill_between(
                                    xx, pp['pdf'] + axis['ymax']-ii-1, axis['ymax'] - ii-1,
                                    where=pp['pdf']>0,
                                    facecolor=fig['colorset']['cmap']['cmap'](1 - ii/axis['ymax'])[0:3],
                                    alpha=0.95,
                                    zorder = 100+ii
                                )
                                axis['1dd'].plot(
                                    xx, pp['pdf'] + axis['ymax']-ii-1, 
                                    '-', 
                                    color=fig['colorset']['1dd']['color'],
                                    linewidth=fig['colorset']['1dd']['linewidth'],
                                    zorder = 100+ii
                                )
                        else:
                            axis['bar'].fill(
                                [1-barwid/1.5, 1, 1, 1-barwid/1.5], [axis['ymax']-ii-1, axis['ymax']-ii-1, axis['ymax']-ii-0.1, axis['ymax']-ii-0.1], 
                                facecolor=fig['colorset']['bar']['fullcolor']
                            )
                            axis['1dd'].plot(
                                fig['colorset']['canvas']['xlim'], [axis['ymax']-ii-1, axis['ymax']-ii-1], 
                                '-', 
                                color=fig['colorset']['1dd']['color'],
                                linewidth=fig['colorset']['1dd']['linewidth'],
                                zorder = 100+ii
                            )
                    axis['1dd'].tick_params(
                        labelsize=fig['colorset']['ticks']['labelsize'],  
                        direction=fig['colorset']['ticks']['direction'],  
                        which='both',
                        zorder=500
                    )
                        # for ii in range(len(grid)):
                        #     print(xxgrid[ii], grid[ii], edge[ii], edge[ii+1])

            elif fig['type'] == "Grid":
                ax = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axbox']['x0'], 
                        fig['colorset']['figureSize']['axbox']['y0'], 
                        fig['colorset']['figureSize']['axbox']['width'], 
                        fig['colorset']['figureSize']['axbox']['height'] 
                ])                
                axc = fig['fig'].add_axes([
                        fig['colorset']['figureSize']['axcbox']['x0'], 
                        fig['colorset']['figureSize']['axcbox']['y0'], 
                        fig['colorset']['figureSize']['axcbox']['width'], 
                        fig['colorset']['figureSize']['axcbox']['height'] 
                ])
                print(emoji.emojize("    :beginner:\tTimer: {:.2f} Second;  Message from '{}' \n\t\t-> Data loading completed".format(
                    time.time()-fig['start'], fig['section']), use_aliases=True))
                sect = dict(self.cf)[fig['section']]
                self.Get3DData(fig)
                fig['var']['lim'] = {
                    'x': [fig['var']['data'].x.min(), fig['var']['data'].x.max()],
                    'y': [fig['var']['data'].y.min(), fig['var']['data'].y.max()],
                    'c': [fig['var']['data'].c.min(), fig['var']['data'].c.max()]
                }
                fig['ax'] = {}
                self.ax_setlim(fig, 'xyc')
                self.ax_setcmap(fig)
                fig['ax']['a1'] = []
                if self.cf.has_option(fig['section'], "x_bin") and self.cf.has_option(fig['section'], "y_bin"):
                    fig['var']['Wbin'] = {
                        "x":    float(self.cf.get(fig['section'], "x_bin"))/2.,
                        "y":    float(self.cf.get(fig['section'], "y_bin"))/2.
                    }
                    for index, row in fig['var']['data'].iterrows():
                        row = pd.Series(row)
                        if self.cf.get(fig['section'], 'c_scale').strip().lower() == 'log':
                            if row['c'] <= fig['ax']['lim']['c'][0]:
                                row['c'] = fig['ax']['lim']['c'][0]
                            cc = (math.log10(row['c']) - math.log10(fig['ax']['lim']['c'][0])) / (math.log10(fig['ax']['lim']['c'][1]) - math.log10(fig['ax']['lim']['c'][0]))
                            fc = fig['colorset']['cmap']['cmap'](cc)[0:3]
                        elif self.cf.get(fig['section'], 'c_scale').strip().lower() == 'flat':
                            cc = (row['c'] - fig['ax']['lim']['c'][0]) / (fig['ax']['lim']['c'][1] - fig['ax']['lim']['c'][0])
                            fc = fig['colorset']['cmap']['cmap'](cc)[0:3]

                        fig['ax']['a1'].append(
                            ax.fill(
                                [row['x'] - fig['var']['Wbin']['x'], row['x'] + fig['var']['Wbin']['x'], row['x'] + fig['var']['Wbin']['x'], row['x'] - fig['var']['Wbin']['x'], row['x'] - fig['var']['Wbin']['x']],
                                [row['y'] - fig['var']['Wbin']['y'], row['y'] - fig['var']['Wbin']['y'], row['y'] + fig['var']['Wbin']['y'], row['y'] + fig['var']['Wbin']['y'], row['y'] - fig['var']['Wbin']['y']],
                                facecolor=fc,
                                edgecolor=None,
                            )
                        )           
                if self.cf.get(fig['section'], 'c_scale').strip().lower() == "flat":          
                    fig['ax']['a2'] = ax.scatter(
                            fig['var']['data'].x,
                            fig['var']['data'].y,
                            c=fig['var']['data'].c,
                            marker='s',
                            edgecolor=None,
                            s=1,
                            cmap=fig['colorset']['cmap']['cmap'],
                            vmin=fig['ax']['lim']['c'][0],
                            vmax=fig['ax']['lim']['c'][1],
                            zorder=0
                        )
                elif self.cf.get(fig['section'], 'c_scale').strip().lower() == "log":
                    fig['ax']['a2'] = ax.scatter(
                            fig['var']['data'].x,
                            fig['var']['data'].y,
                            c=fig['var']['data'].c,
                            marker='s',
                            edgecolor=None,
                            s=1,
                            cmap=fig['colorset']['cmap']['cmap'],
                            norm=matplotlib.colors.LogNorm(vmin=fig['ax']['lim']['c'][0],vmax=fig['ax']['lim']['c'][1]),
                            zorder=0
                        )
                
                self.ax_setticks(fig, 'xyc')
                if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                    if not self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                        ax.set_xticks(fig['ax']['ticks']['x'])
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'x_scale').strip().lower() == "log":
                    ax.set_xscale('log')
                if self.cf.get(fig['section'], 'x_ticks')[0:4] == 'Manu':
                    ax.xaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['x'][0]))
                    ax.set_xticklabels(fig['ax']['ticks']['x'][1])
                    if self.cf.get(fig['section'], 'x_scale').strip().lower() == "flat":
                        ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.set_xlim(fig['ax']['lim']['x'][0], fig['ax']['lim']['x'][1])

                if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                    if not self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                        ax.set_yticks(fig['ax']['ticks']['y'])
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                elif self.cf.get(fig['section'], 'y_scale').strip().lower() == 'log':
                    ax.set_yscale('log')
                if self.cf.get(fig['section'], 'y_ticks')[0:4] == "Manu":
                    ax.yaxis.set_major_locator(
                        ticker.FixedLocator(fig['ax']['ticks']['y'][0]))
                    ax.set_yticklabels(fig['ax']['ticks']['y'][1])
                    if self.cf.get(fig['section'], 'y_scale').strip().lower() == 'flat':
                        ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.set_ylim(fig['ax']['lim']['y'][0], fig['ax']['lim']['y'][1])

                if self.cf.has_option(fig['section'], 'c_scale'):
                    if self.cf.get(fig['section'], 'c_scale').strip().lower() == 'log':
                        from matplotlib.colors import LogNorm
                        plt.colorbar(fig['ax']['a2'], axc, norm=LogNorm(vmin=fig['ax']['lim']['c'][0], vmax=fig['ax']['lim']['c'][1]), orientation='vertical', extend='neither')
                        axc.set_yscale("log")
                    elif self.cf.get(fig['section'], 'c_scale').strip().lower() == "flat":
                        plt.colorbar(fig['ax']['a2'], cax=axc, ticks=fig['ax']['ticks']['c'], orientation='vertical', extend='neither')
                # axc.set_ylim(fig['ax']['lim']['c'][0], fig['ax']['lim']['c'][1])

                ax.tick_params(
                    labelsize=fig['colorset']['ticks']['labelsize'],
                    direction=fig['colorset']['ticks']['direction'],
                    bottom=fig['colorset']['ticks']['bottom'],
                    left=fig['colorset']['ticks']['left'],
                    top=fig['colorset']['ticks']['top'],
                    right=fig['colorset']['ticks']['right'],
                    which='both'
                )
                ax.tick_params(which='major', length=fig['colorset']['ticks']
                               ['majorlength'], color=fig['colorset']['ticks']['majorcolor'])
                ax.tick_params(which='minor', length=fig['colorset']['ticks']
                               ['minorlength'], color=fig['colorset']['ticks']['minorcolor'])
                axc.tick_params(
                    which="both",
                    labelsize=fig['colorset']['colorticks']['labelsize'],
                    direction=fig['colorset']['colorticks']['direction'],
                    bottom=fig['colorset']['colorticks']['bottom'],
                    left=fig['colorset']['colorticks']['left'],
                    top=fig['colorset']['colorticks']['top'],
                    right=fig['colorset']['colorticks']['right'],
                    color=fig['colorset']['colorticks']['color']
                )
                # axc.tick_params(which='major', length=fig['colorset']['ticks']['majorlength'], color=fig['colorset']['ticks']['majorcolor'])

                ax.set_xlabel(
                    r"{}".format(self.cf.get(fig['section'], 'x_label')), 
                    fontsize=fig['colorset']['canvas']['x_label']['size'],
                    color=fig['colorset']['canvas']['x_label']['color']                   
                )
                ax.xaxis.set_label_coords(0.5, fig['colorset']['canvas']['x_label']['offline'])
                ax.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], 'y_label')), 
                    fontsize=fig['colorset']['canvas']['yl_label']['size'],
                    color=fig['colorset']['canvas']['yl_label']['color'],
                )
                ax.yaxis.set_label_coords(fig['colorset']['canvas']['yl_label']['offline'], 0.5)
                axc.set_ylabel(
                    r"{}".format(self.cf.get(fig['section'], 'c_label')), 
                    fontsize=fig['colorset']['canvas']['yr_label']['size'],
                    color=fig['colorset']['canvas']['yr_label']['color'],
                )
                axc.yaxis.set_label_coords(fig['colorset']['canvas']['yr_label']['offline'], 0.5)

                if self.cf.has_option(fig['section'], 'Line_draw'):
                    self.drawline(fig, ax)

                if self.cf.has_option(fig['section'], "Text"):
                    self.drawtext(fig, ax)

                if self.cf.has_option(fig['section'], "fill"):
                    self.drawfillarea(fig, ax)
                
            if 'save' in self.cf.get(fig['section'], 'print_mode'):
                fig['fig'] = plt
                fig['file'] = os.path.join(self.figpath, fig['name'])
                self.savefig(fig, plt)
            if 'show' in self.cf.get(fig['section'], 'print_mode'):
                plt.show(block=False)
                plt.pause(1)
                input("\n\tPress 'Enter' to continue ...")
                plt.close()

        # print(vars(self).keys())

    def get_DK_data(self, fig):
        def get_conbine_name(it):
            name = "BCdk2"
            tex  = ""
            it['fstate'] = it['finalstates'][0]
            it['res_pp'] = list(set(it['fstate']) - set(it['cb']['item']))
            if it['res_pp']:
                for pp in it['res_pp']:
                    for itp in fig['particle']:
                        if fig['particle'][itp]['PDG'] ==  pp or fig['particle'][itp]['PDG'] == -pp:
                            tex += "{} ".format(fig['particle'][itp]['LaTeX'])
                            name += itp
            name += it['cb']['name']
            tex  += it['cb']['TeX']
            return name, "${}$".format(tex)


        def load_BRComb():
            if fig['DKs']['BRcomb']:
                if fig['DKs']['BRcomb'] in fig['colorset']['BRComb']:
                    fig['DKs']['Combcf'] = self.load_path(fig['colorset']['BRComb'][fig['DKs']['BRcomb']])
                    print("\tTimer: {:.2f} Second; Message from {}\n\t\tLoading Combine Schema -> {} : {}".format(time.time()-fig['start'], fig['section'], fig['DKs']['BRcomb'], fig['DKs']['Combcf']))
                    with open(fig['DKs']['Combcf'], 'r') as f1:
                        fig['DKs']['Combcf'] = json.loads(f1.read())
                    fig['DKs']['CS'] = {}
                    for scmic in fig['DKs']['Combcf']['Schema_include']:
                        for scm in fig['DKs']['Combcf'][scmic]:
                            fig['DKs']['CS'][scm] = fig['DKs']['Combcf'][scm]
            fig['DKs']['cbs']       = []
            fig['DKs']['cbslist']   = []
            fig['DKs']['inflist']   = []
            for ii, item in enumerate(fig['DKs']['info']):
                intag = False
                for jj, cs in fig['DKs']['CS'].items():
                    for cb in cs['Combs']:
                        if set(cb) <= set(item['finalstates'][0]):
                            intag = True
                            item['cb'] = {
                                "name":   jj,
                                "item":   cb,
                                "TeX":    cs['TeX']
                            }
                            break
                if intag:
                    name, tex = get_conbine_name(item)
                    if name not in fig['DKs']['cbslist']:
                        fig['DKs']['cbslist'].append(name)
                        dset = {
                            "name": name,
                            "label":  tex,
                            "data": {
                                item['name']:   fig['data'][item['name']]
                            },
                            "fstate":   item['finalstates']
                        }
                        fig['DKs']['cbs'].append(dset)
                    else:
                        for cb in fig['DKs']['cbs']:
                            if cb['name'] == name:
                                cb['data'][item['name']] = fig['data'][item['name']]
                                cb['fstate'] += item['finalstates']
                                break 
                else:
                    fig['DKs']['inflist'].append(item['name'])
            for cb in fig['DKs']['cbs']:
                cb['data'] = pd.DataFrame(cb['data']).sum(axis=1)

        def set_columns_data():
            if fig['colorset']['figureSize']['columns'] == 1:
                fig['DKs']['data'] = {
                    "ax0":  []
                }
            else:
                fig['DKs']['data'] = {}
                for ii in range(fig['colorset']['figureSize']['columns']):
                    fig['DKs']['data']['ax{}'.format(ii)] = []
            lencb = 0
            if self.cf.has_option(fig['section'], "BRcomb"):
                lencb += len(fig['DKs']['cbs'])
                for ii, item in enumerate(fig['DKs']['cbs']):
                    axno = ii // fig['colorset']['figureSize']['maxin1C']
                    fig['DKs']['data']['ax{}'.format(axno)].append(item)
                nn = 0
                for ii, item in enumerate(fig['DKs']['info']):
                    if item['name'] not in fig['DKs']['inflist']:
                        pass 
                    else:
                        axno = (nn + lencb) // fig['colorset']['figureSize']['maxin1C']
                        dset = {
                            "name":     item['name'],
                            "label":    get_fslabel(item['finalstates']),
                            "data":     fig['data'][item['name']],
                            "fstate":   item['finalstates']
                        }
                        fig['DKs']['data']['ax{}'.format(axno)].append(dset)
                        nn += 1
            else:
                axno = (nn + lencb) // fig['colorset']['figureSize']['maxin1C']
                dset = {
                    "name":     item['name'],
                    "label":    get_fslabel(item['finalstates']),
                    "data":     fig['data'][item['name']],
                    "fstate":   item['finalstates']
                }
                fig['DKs']['data']['ax{}'.format(axno)].append(dset)
            # print(fig['DKs']['data'])

        def get_fslabel(fs):
            tex = ""
            for pp in fs[0]:
                for it in fig['particle']:
                    if fig['particle'][it]['PDG'] ==  pp or fig['particle'][it]['PDG'] == -pp:
                        tex += "{} ".format(fig['particle'][it]['LaTeX'])
                        break
            return "${}$".format(tex)

        dkpath = self.load_path(self.cf.get(fig['section'], "table_path"))
        DKdata = pd.read_csv(dkpath)
        self.data = pd.merge(self.data, DKdata, on=['Index'], how='outer')
        self.basic_selection(fig)
        fig['DKs'] = {}
        ifpath = self.load_path(self.cf.get(fig['section'], "info_path"))
        fig['colorset']['PDGid'] = self.load_path(fig['colorset']['PDGid'])
        with open(fig['colorset']['PDGid'], 'r') as f1:
            fig['particle'] = json.loads(f1.read())
        with open(ifpath, 'r') as f1:
            fig['DKs']['info'] = json.loads(f1.read())
        fig['ax'] = {
            "xlim":  fig['colorset']['canvas']['xlim'],
            "axs":  []
        }
        if self.cf.has_option(fig['section'], "BRcomb"):
            fig['DKs']['BRcomb'] = self.cf.get(fig['section'], "BRcomb").strip()
            load_BRComb()
            fig['colorset']['figureSize']['rows'] = len(fig['DKs']['cbslist']) + len(fig['DKs']['inflist'])
        else:
            fig['colorset']['figureSize']['rows'] = len(fig['DKs']['info'])
        fig['colorset']['figureSize']['ends'] = fig['colorset']['figureSize']['rows'] % fig['colorset']['figureSize']['maxin1C'] + fig['colorset']['figureSize']['1dheight'] - 1
        fig['colorset']['figureSize']['columns'] = fig['colorset']['figureSize']['rows'] // fig['colorset']['figureSize']['maxin1C'] + 1
        if fig['colorset']['figureSize']['rows'] % fig['colorset']['figureSize']['maxin1C'] == 0:
            fig['colorset']['figureSize']['columns'] -= 1
            fig['colorset']['figureSize']['ends'] = fig['colorset']['figureSize']['maxin1C']+ fig['colorset']['figureSize']['1dheight']-1
        if fig['colorset']['figureSize']['columns'] == 1:
            fig['ax']['ylim'] = [0.0, fig['colorset']['figureSize']['rows'] + fig['colorset']['figureSize']['1dheight']]
            fig['colorset']['figureSize']['axheight'] = (fig['colorset']['figureSize']['rows'] + fig['colorset']['figureSize']['1dheight']-1) * fig['colorset']['figureSize']['1dhunit'] 
        elif fig['colorset']['figureSize']['columns'] > 1:
            fig['ax']['ylim'] = [0.0, fig['colorset']['figureSize']['maxin1C'] + fig['colorset']['figureSize']['1dheight']]
            fig['colorset']['figureSize']['axheight'] = (fig['colorset']['figureSize']['maxin1C']+fig['colorset']['figureSize']['1dheight']-1) * fig['colorset']['figureSize']['1dhunit'] 
        set_columns_data()

    def make_decays_title(self, fig):
        x0 = 0.04 * (fig['colorset']['figureSize']['text'] + fig['colorset']['figureSize']['barwidth'] + fig['colorset']['figureSize']['1dwidth']) / fig['colorset']['figureSize']['width']
        x2 = x0 + 0.92*(fig['colorset']['figureSize']['text'] + fig['colorset']['figureSize']['barwidth'] + fig['colorset']['figureSize']['1dwidth'])/ fig['colorset']['figureSize']['width']
        y0 = (fig['colorset']['figureSize']['bottom'] + fig['colorset']['figureSize']['axheight'])/ fig['colorset']['figureSize']['height']
        y1 = y0 + 0.8*fig['colorset']['figureSize']['tiheight']/fig['colorset']['figureSize']['height']
        x1 = x0 + (y1-y0) * fig['colorset']['figureSize']['height'] / fig['colorset']['figureSize']['width']
        fig['ax']['title'] = fig['fig'].add_axes([x1, y0, x2-x1, y1-y0])
        fig['ax']['logo'] = fig['fig'].add_axes([x0, y0, x1-x0, y1-y0])
        fig['title'] = {
            "xy_ratio": ((x2-x1)*fig['colorset']['figureSize']['width'] )/ ((y1-y0)*fig['colorset']['figureSize']['height'])
        }
        fig['ax']['title'].spines['top'].set_visible(False)
        fig['ax']['title'].spines['bottom'].set_visible(False)
        fig['ax']['title'].spines['left'].set_visible(False)
        fig['ax']['title'].spines['right'].set_visible(False)
        fig['ax']['title'].get_xaxis().set_visible(False)
        fig['ax']['title'].get_yaxis().set_visible(False)        
        fig['ax']['logo'].spines['top'].set_visible(False)
        fig['ax']['logo'].spines['bottom'].set_visible(False)
        fig['ax']['logo'].spines['left'].set_visible(False)
        fig['ax']['logo'].spines['right'].set_visible(False)
        fig['ax']['logo'].get_xaxis().set_visible(False)
        fig['ax']['logo'].get_yaxis().set_visible(False)
        fig['ax']['title'].set_xlim(0, fig['title']['xy_ratio'])

        fig['ax']['title'].set_ylim(0, 1)
        fig['ax']['logo'].set_xlim(0, 674)
        fig['ax']['logo'].set_ylim(0, 674)
        import matplotlib.patches as mpatches
        fancybox = mpatches.FancyBboxPatch(
            [0.275*674, 0.275*674], 0.45*674, 0.45*674,
            boxstyle=mpatches.BoxStyle("Round", pad=0.27*674),
            color="None"
        )
        fig['ax']['logo'].add_patch(fancybox)
        fig['ax']['title'].text(0.5, 0.06, "BudingPLOT", 
                            fontsize=5,
                            color='w', fontfamily = "monospace",
                            horizontalalignment="center",
                            verticalalignment="bottom"
        )
        import matplotlib.cbook as cbook
        with cbook.get_sample_data(fig['colorset']['logopath']) as image_file:
            logo = plt.imread(image_file)
        im = fig['ax']['logo'].imshow(logo)
        im.set_clip_path(fancybox)
        if self.cf.has_option(fig['section'], "title"):
            fig['title']['label'] = self.cf.get(fig['section'], "title")
        else:
            pp = int(os.path.basename(self.load_path(self.cf.get(fig['section'], "info_path"))).replace(fig['colorset']['title']['dfjs'], ""))
            for itp in fig['particle']:
                if abs(fig['particle'][itp]['PDG']) ==  abs(pp):
                    fig['title']['pp'] = "${}$".format(fig['particle'][itp]['LaTeX'])
                    break
            fig['title']['label'] = "{}{}".format(fig['title']['pp'], fig['colorset']['title']['ends'])
        fig['ax']['title'].text(
            fig['colorset']['title']['offleft'], fig['colorset']['title']['offtop'], fig['title']['label'],
            fontsize=fig['colorset']['title']['fontsize'],
            color=fig['colorset']['title']['fontcolor'],
            horizontalalignment="left",
            verticalalignment="center"
        )

    def make_decays_canvas(self, fig):
        fig['colorset']['figureSize']['height'] = fig['colorset']['figureSize']['axheight'] + fig['colorset']['figureSize']['tiheight'] + fig['colorset']['figureSize']['bottom']
        fig['colorset']['figureSize']['width'] = (fig['colorset']['figureSize']['text'] + fig['colorset']['figureSize']['barwidth'] + fig['colorset']['figureSize']['1dwidth']) * fig['colorset']['figureSize']['columns']
        fig['fig'].set_figheight(fig['colorset']['figureSize']['height'])
        fig['fig'].set_figwidth(fig['colorset']['figureSize']['width'])
        for ii in range(fig['colorset']['figureSize']['columns']):
            xtot = ii / fig['colorset']['figureSize']['columns']
            x0 = (1.0 - fig['colorset']['figureSize']['bar1dtot'])/2.0 + xtot
            x1 = fig['colorset']['figureSize']['barwidth'] / fig['colorset']['figureSize']['width'] * fig['colorset']['figureSize']['bar1dtot'] + x0 
            x2 = x1 + fig['colorset']['figureSize']['barin1d'] / fig['colorset']['figureSize']['columns']
            x3 = x2 + fig['colorset']['figureSize']['1dwidth'] / fig['colorset']['figureSize']['width'] * fig['colorset']['figureSize']['bar1dtot']
            ymin = fig['colorset']['figureSize']['bottom'] / fig['colorset']['figureSize']['height']
            ymax = (fig['colorset']['figureSize']['bottom'] + fig['colorset']['figureSize']['axheight'])/ fig['colorset']['figureSize']['height']
            yend = (fig['colorset']['figureSize']['bottom'] + (fig['colorset']['figureSize']['maxin1C'] - fig['colorset']['figureSize']['ends']+ fig['colorset']['figureSize']['1dheight']-1) * fig['colorset']['figureSize']['1dhunit'] )/ fig['colorset']['figureSize']['height']
            if ii + 1 == fig['colorset']['figureSize']['columns'] and ii != 0:
                axes = {
                    "name": "ax{}".format(ii),
                    "bar":  [x0, yend, x1-x0, ymax-yend],
                    "1dd":  [x2, yend, x3-x2, ymax-yend],
                    "yupp": fig['colorset']['figureSize']['ends'],
                    "ymax": fig['colorset']['figureSize']['ends'] - fig['colorset']['figureSize']['1dheight']+1
                }
            elif ii + 1 == fig['colorset']['figureSize']['columns'] and ii == 0:
                fig['colorset']['figureSize']['maxin1C'] = fig['colorset']['figureSize']['rows']
                yend = ymax
                axes = {
                    "name": "ax{}".format(ii),
                    "bar":  [x0, ymin, x1-x0, ymax-ymin],
                    "1dd":  [x2, ymin, x3-x2, ymax-ymin],
                    "yupp": fig['colorset']['figureSize']['maxin1C'] + fig['colorset']['figureSize']['1dheight']-1,
                    "ymax": fig['colorset']['figureSize']['maxin1C']
                }
            else:
                axes = {
                    "name": "ax{}".format(ii),
                    "bar":  [x0, ymin, x1-x0, ymax-ymin],
                    "1dd":  [x2, ymin, x3-x2, ymax-ymin],
                    "yupp": fig['colorset']['figureSize']['maxin1C'] + fig['colorset']['figureSize']['1dheight']-1,
                    "ymax":  fig['colorset']['figureSize']['maxin1C']
                }
            fig['ax']['axs'].append(axes)
        for axis in fig['ax']['axs']:
            axis['bar'] = fig['fig'].add_axes(axis['bar'])
            axis['bar'].set_xlim(0, 1)
            axis['bar'].set_ylim(0, axis['yupp'])
            axis['bar'].spines['top'].set_visible(False)
            axis['bar'].spines['bottom'].set_visible(False)
            axis['bar'].spines['left'].set_visible(False)
            axis['bar'].spines['right'].set_visible(False)
            axis['bar'].get_xaxis().set_visible(False)
            axis['bar'].get_yaxis().set_visible(False)

            axis['1dd'] = fig['fig'].add_axes(axis['1dd'])
            # axis['1dd'].set_xscale("log")
            axis['1dd'].set_xlim(fig['ax']['xlim'])
            axis['1dd'].set_ylim(0, axis['yupp'])
            axis['1dd'].spines['top'].set_visible(False)
            axis['1dd'].spines['left'].set_visible(False)
            axis['1dd'].spines['right'].set_visible(False)
            axis['1dd'].get_yaxis().set_visible(False)
            from matplotlib.ticker import MultipleLocator, PercentFormatter
            axis['1dd'].set_xticks([0.05, 0.25, 0.5, 0.75, 1])
            # axis['1dd'].xaxis.set_major_locator(MultipleLocator(0.25))
            axis['1dd'].xaxis.set_major_formatter(PercentFormatter(xmax=1))
            axis['1dd'].xaxis.set_minor_locator(MultipleLocator(0.05))

            # axis['1dd'].xaxis.set_subminor_locator(MultipleLocator(0.02))
            # axis['1dd'].set_xticklabels(["0.05", "0.2", "0.4", "0.6", "0.8", "1"])

    def draw_scatter_legend(self, fig, ax):
        def get_lg_info(lcd):
            lgset = {}
            with open(self.load_path(fig['colorset']['markercodepath']), 'r') as f1:
                fig['colorset']['scattermarker'] = json.loads(f1.read())
            with open(self.load_path(fig['colorset']["colorpath"]), 'r') as f1:
                fig['colorset']['scattercolor'] = json.loads(f1.read())
            mline = self.cf.get(fig['section'], 'marker').split('\n')
            for line in lcd:
                p_rec = re.compile(r'[\[].*?[\]]', re.S)
                a = re.findall(p_rec, line.strip())[0]
                line = line.replace(a, "").split(',')
                lname = line[0].strip()
                lst = {
                    'loc':  line[1].strip(),
                    'include':  eval(a),
                    'marker':   [],
                    'label':    line[-1].strip()
                }
                for ids in lst['include']:
                    mks     = mline[ids-1].split(',')[0:2]
                    lst['marker'].append({
                        'type':         fig['colorset']['scattermarker'][mks[1].strip()],     
                        'color':      fig['colorset']['scattercolor'][mks[0].strip()],
                        'edgecolor':   fig['colorset']['legend']['edgecolor'],
                        'size':        fig['colorset']['legend']['markersize'],
                        'alpha':        fig['colorset']['legend']['alpha']
                    })
                if lname not in lgset.keys():
                    lgset[lname] = [lst]
                else:
                    lgset[lname].append(lst)
            return lgset

        def get_legend_data():
            fig['legend']['data'] = {}
            renderer1 = fig['fig'].canvas.get_renderer()
            fig['legend']['ax'] = {
                "axbb":     ax.get_window_extent(renderer1)
            }
            for kk, lg in fig['legend']['info'].items():
                fig['legend']['data'][kk] = {
                    "txtwidth":    [],
                    "txtheight":   []
                }
                for it in lg:
                    txt = ax.text(
                        0., 0., it['label'], color="None",
                        fontsize=fig['colorset']['legend']['fontsize'],
                        horizontalalignment='left',
                        verticalalignment='top', 
                        bbox=dict(facecolor='None', alpha=0., edgecolor='None')
                    )
                    bbox = txt.get_window_extent(renderer1)
                    it['bbox'] = {
                        "width":    bbox.width/fig['legend']['ax']['axbb'].width,
                        "height":    bbox.height/fig['legend']['ax']['axbb'].height
                    }
                    fig['legend']['data'][kk]['txtwidth'].append(bbox.width/fig['legend']['ax']['axbb'].width)
                    fig['legend']['data'][kk]['txtheight'].append(bbox.height/fig['legend']['ax']['axbb'].height)
                fig['legend']['data'][kk]['x_txt'] = max(fig['legend']['data'][kk]['txtwidth'])
                fig['legend']['data'][kk]['y_txt'] = max(fig['legend']['data'][kk]['txtheight'])
                fig['legend']['data'][kk]['xx'] = 1.0 - fig['legend']['data'][kk]['x_txt'] - fig['colorset']['legend']['offset']

            ymin= 0.
            for ll in fig['legend']['info'].keys():
                for ii in range(len(fig['legend']['info'][ll])):
                    it = fig['legend']['info'][ll][ii]
                    yy = ymin + fig['colorset']['legend']['offset'] + (0.5 + ii) * (fig['legend']['data'][ll]['y_txt'] + fig['colorset']['legend']['dy'])
                    txt = ax.text(
                        fig['legend']['data'][ll]['xx'], yy, it['label'], 
                        fontsize=fig['colorset']['legend']['fontsize'],
                        color=fig['colorset']['legend']['fontcolor'],
                        horizontalalignment='left',
                        verticalalignment='center',
                        transform=ax.transAxes
                    )
                    for jj in range(len(it['marker'])):
                        xx = fig['legend']['data'][ll]['xx'] - fig['colorset']['legend']['markerbttext'] - (0.5 + jj) * fig['colorset']['legend']['markersplit']
                        ax.plot(
                            xx, yy, it['marker'][jj]['type'],
                            color=it['marker'][jj]['color'],
                            markersize=it['marker'][jj]['size'],
                            transform=ax.transAxes
                        )
                xx = fig['legend']['data'][kk]['x_txt'] + fig['colorset']['legend']['markerbttext'] + (len(fig['legend']['info'][ll])+2.5) * fig['colorset']['legend']['markersplit']
                yy = (fig['legend']['data'][ll]['y_txt'] + fig['colorset']['legend']['dy']) * len(fig['legend']['info'][ll]) 
                import matplotlib.patches as mpatches
                fancybox = mpatches.FancyBboxPatch(
                    [
                        1.0 - xx - fig['colorset']['legend']['offset'], 
                        fig['colorset']['legend']['offset'] + ymin
                    ],
                    xx,
                    yy,
                    boxstyle=mpatches.BoxStyle(
                        "Round", 
                        pad=fig['colorset']['legend']['pad']
                    ),
                    edgecolor=fig['colorset']['legend']['padlc'],
                    facecolor=fig['colorset']['legend']['padcolor'],
                    alpha=fig['colorset']['legend']['padalpha'],
                    transform=ax.transAxes
                )
                ax.add_patch(fancybox)
                ymin += (fig['legend']['data'][ll]['y_txt'] + fig['colorset']['legend']['dy']) * len(fig['legend']['info'][ll]) + fig['colorset']['legend']['offset']




        if self.cf.has_option(fig['section'], "draw_legend"):
            legcmd = self.cf.get(fig['section'], "draw_legend")
            fig['legend'] = {}
            if legcmd.lower() == "true":
                pass
            elif self.cf.has_option(fig['section'], 'marker'):
                legcmd = legcmd.split("\n")
                fig['legend']['info'] = get_lg_info(legcmd)
                get_legend_data()
            
                
    def mark_ternary_channel(self, fig, ax, cmode):
        chkv = {'r': "Red ", 'g': "Green", "b": 'Blue'}
        axs = "{}size".format(ax)
        ptk = [fig['cf'][axs]['x']-0.02, fig['cf'][axs]
               ['y']+fig['cf'][axs]['height'] - 0.015]
        ptv = [fig['cf'][axs]['x']-0.02, fig['cf'][axs]
               ['y']+fig['cf'][axs]['height'] - 0.018]
        ktex = r"$\rm {}~Channel:$".format(chkv[cmode])
        fig['cv'].text(ptk[0], ptk[1], ktex, fontsize=8,
                       horizontalalignment="left", verticalalignment="bottom")
        fig['cv'].text(ptv[0], ptv[1], r"{}".format(self.cf.get(fig['section'], "{}_label".format(
            cmode))), fontsize=12, horizontalalignment="left", verticalalignment="top")

    def draw_ternary_legend(self, fig, ax):
        def pick_tri_color(verts):
            x = 0.9 * ((verts[0][0] + verts[1][0] +
                        verts[2][0]) / 3.0)**0.25 + 0.1
            y = 0.9 * ((verts[0][1] + verts[1][1] +
                        verts[2][1]) / 3.0)**0.25 + 0.1
            z = 0.9 * ((verts[0][2] + verts[1][2] +
                        verts[2][2]) / 3.0)**0.25 + 0.1
            return (x, y, z, 1.)

        def corrdinate(vert):
            return vert[1] + 0.5*vert[2], vert[2]

        def draw_triangle(vertex, pm):
            if pm == "+":
                order = "abc"
            elif pm == '-':
                order = "acb"
            xl, yl = [], []
            verts = []
            for arr in order:
                vertex += fig['legend']['arraw'][arr]
                x, y = corrdinate(vertex)
                xl.append(x)
                yl.append(y)
                verts.append(list(vertex))

            fig['ax'][ax].fill(xl, yl, facecolor=pick_tri_color(verts))

        def draw_color_space():
            top = np.array([0., 0., 1.])
            vertex = copy.copy(top)
            while top[2] > 1.e-4:
                if vertex[0] < 1.e-4:
                    top += fig['legend']['arraw']["c"]
                    vertex = copy.copy(top)
                while vertex[0] > 1.e-4:
                    draw_triangle(vertex, "+")
                    if vertex[2] > 1.e-4:
                        draw_triangle(vertex, '-')
                    vertex += fig['legend']['arraw']["a"]

        def draw_colorbar(axis, pos, cmode):
            # print(axis)
            sp, ep, lp = axis
            # print(sp, ep)
            yUnitArraw = (np.array(ep) - np.array(sp))/fig['legend']['scale']
            xUnitArraw = np.array(lp) - 0.5*(np.array(ep)+np.array(sp))
            # print(xUnitArraw, yUnitArraw)
            for ii in range(fig['legend']['scale']):
                pa = sp - 0.02 * xUnitArraw + ii * yUnitArraw
                pb = sp - 0.07 * xUnitArraw + ii * yUnitArraw
                pc = sp - 0.07 * xUnitArraw + (ii+1) * yUnitArraw
                pd = sp - 0.02 * xUnitArraw + (ii+1) * yUnitArraw
                cc = (ii + 0.5)/fig['legend']['scale']
                # print(pa, pb, pc, pd)
                if cmode == 'r':
                    color = tuple([cc, 0.38*cc, 0.38*cc])
                elif cmode == 'g':
                    color = tuple([0.44*cc, cc, 0.44*cc])
                elif cmode == 'b':
                    color = tuple([0.42*cc, 0.42*cc, cc])
                fig['cv'].fill(
                    [pa[0], pb[0], pc[0], pd[0]],
                    [pa[1], pb[1], pc[1], pd[1]],
                    facecolor=color
                )
            for ii in np.linspace(0., 1., 6):
                pa = sp - 0.07 * xUnitArraw + ii * \
                    fig['legend']['scale'] * yUnitArraw
                pb = sp - 0.04 * xUnitArraw + ii * \
                    fig['legend']['scale'] * yUnitArraw
                pt = sp - 0.12 * xUnitArraw + ii * \
                    fig['legend']['scale'] * yUnitArraw
                fig['cv'].plot(
                    [pa[0], pb[0]],
                    [pa[1], pb[1]],
                    '-',
                    linewidth=fig['legend']['sty']['ticks']['width'],
                    color=fig['legend']['sty']['ticks']['color'],
                    solid_capstyle=fig['legend']['sty']['ticks']['end'],
                    zorder=100
                )
                fig['cv'].text(
                    pt[0],
                    pt[1],
                    r"${:.1f}$".format(ii),
                    fontsize=fig['legend']['sty']['ticks']['fontsize'],
                    rotation=fig['legend']['sty']['ticks']['fontangle'][pos],
                    horizontalalignment=fig['legend']['sty']['ticks']["labelshift"][pos][2],
                    verticalalignment=fig['legend']['sty']['ticks']["labelshift"][pos][3]
                )
            pa = sp - 0.07 * xUnitArraw
            pb = ep - 0.07 * xUnitArraw
            fig['cv'].plot(
                [pa[0], pb[0]],
                [pa[1], pb[1]],
                '-',
                linewidth=fig['legend']['sty']['ticks']['width'],
                color=fig['legend']['sty']['ticks']['color'],
                solid_capstyle=fig['legend']['sty']['ticks']['end'],
                zorder=100
            )

            pl = 0.5 * np.array(sp) + 0.5 * np.array(ep) - 0.26 * xUnitArraw
            fig['cv'].text(
                pl[0],
                pl[1],
                r"{}".format(self.cf.get(
                    fig['section'], "{}_label".format(cmode))),
                fontsize=fig['legend']['sty']['label']['fontsize'],
                color=fig['legend']['sty']['label']['color'],
                rotation=fig['legend']['sty']['label']['rotation'][pos],
                horizontalalignment="center",
                verticalalignment="center"
            )

        fig['legend'] = {}
        axs = "{}size".format(ax)
        fig['legend']['scale'] = fig['cf'][axs]['legendgrids']
        fig['legend']['unit'] = 1.0 / fig['legend']['scale']
        fig['legend']['sty'] = fig['cf'][fig['cf'][axs]['axisstyle']]
        fig['legend']['arraw'] = {
            "a":    np.array([-fig['legend']['unit'], fig['legend']['unit'], 0.0]),
            "b":    np.array([0.0, -fig['legend']['unit'], fig['legend']['unit']]),
            "c":    np.array([fig['legend']['unit'], 0.0, -fig['legend']['unit']])
        }
        draw_color_space()

        pointA = [fig['cf'][axs]['x'], fig['cf'][axs]['y']]
        pointB = [fig['cf'][axs]['x'] + fig['cf']
                  [axs]['width'], fig['cf'][axs]['y']]
        pointC = [fig['cf'][axs]['x']+0.5*fig['cf'][axs]['width'],
                  fig['cf'][axs]['y']+fig['cf'][axs]['height']]
        # print(pointA, pointB, pointC)
        draw_colorbar((pointA, pointB, pointC), "bottom", 'g')
        draw_colorbar((pointB, pointC, pointA), "right", 'b')
        draw_colorbar((pointC, pointA, pointB), "left", 'r')

    def set_violin_canvas(self, fig):
        fig['fig'].set_figheight(fig['colorset']['figureSize']['height'])
        totalwidth = (fig['colorset']['figureSize']['left'] + fig['colorset']['figureSize']
                      ['right'] + len(fig['varlist']))*fig['colorset']['figureSize']['Unit']
        fig['fig'].set_figwidth(totalwidth)
        fig['ax']['ax'] = fig['fig'].add_axes([fig['colorset']['figureSize']['left']/totalwidth, fig['colorset']['figureSize']['bottom']/fig['colorset']['figureSize']['height'], 1-(
            fig['colorset']['figureSize']['left']+fig['colorset']['figureSize']['right'])/totalwidth, 1-(fig['colorset']['figureSize']['bottom']+fig['colorset']['figureSize']['top'])/fig['colorset']['figureSize']['height']])

    def get_violin_data(self, fig):
        def get_variable_list():
            lvar = []
            for opt in self.cf.options(fig['section']):
                if opt[0:8] == 'variable':
                    no = opt[8:]
                    if "label{}".format(no) in self.cf.options(fig['section']):
                        lvar.append({
                            "name":     opt,
                            "label":    self.cf.get(fig['section'], "label{}".format(no))
                        })
            return lvar

        def get_selection_data(bo):
            if bo[0:3] == "&Bo":
                if "&V" in bo:
                    bo = self.get_defined_var(fig, bo, 'data')
                bo = bo[4:]
                x_sel = self.var_symbol(bo)
                for x in x_sel:
                    bo = bo.replace("_{}".format(
                        x), "self.data['{}']".format(x))
                if "&FC_" in bo:
                    for ii in range(len(self.funcs)):
                        bo = bo.replace(
                            self.funcs[ii]['name'], "self.funcs[{}]['expr']".format(ii))
                if "*Bool*" not in self.data.columns:
                    bool_list = np.ones(self.data.shape[0], dtype=np.bool)
                    self.data['*Bool*'] = bool_list
                    bo = bo + "& self.data['*Bool*']"
                else:
                    bool_list = np.ones(self.data.shape[0], dtype=np.bool)
                    self.data['abc**cba**bool**loob**alphabeta'] = bool_list
                    bo = bo + "& self.data['abc**cba**bool**loob**alphabeta']"
                res = fig['data'].loc[eval(bo)]
                return res

        def get_axvardata(fig, name, varinf, dataset):
            if varinf[0:3] == '&Eq':
                varinf = varinf[4:]
                if "&V" in varinf:
                    varinf = self.get_defined_var(fig, varinf, 'data')
                x_sel = self.var_symbol(varinf)
                for x in x_sel:
                    varinf = varinf.replace("_{}".format(
                        x), "fig['axdata']['{}']['{}']".format(dataset, x))
                if "&FC_" in varinf:
                    for ii in range(len(self.funcs)):
                        varinf = varinf.replace(
                            self.funcs[ii]['name'], "self.funcs[{}]['expr']".format(ii))
                res = eval(varinf)
            elif varinf in fig['data'].columns.values:
                res = fig['axdata'][dataset][varinf]
            else:
                print("No Variable {} found in Data!".format(varinf))
                sys.exit(0)
            return res

        fig['ax'] = {}
        fig['var'] = {}
        fig['varlist'] = get_variable_list()
        if fig['type'] == "BiViolin":
            if fig['Biviolin']:
                fig['axdata'] = {
                    "left": get_selection_data(self.cf.get(fig['section'], "left_selection")),
                    "right":get_selection_data(self.cf.get(fig['section'], "right_selection"))
                }
                print("\t\tLeft side contains  -> {} \trows\n\t\tRight side contains -> {} \trows".format(fig['axdata']['left'].shape[0], fig['axdata']['right'].shape[0]))
                fig['ax']['data'] = {
                    "left":     [],
                    "right":    [],
                    "label":    [],
                    "sum":      []
                }
                for var in fig['varlist']:
                    fig['ax']['data']['left'].append(get_axvardata(fig, var['name'], self.cf.get(fig['section'], var['name']), 'left'))
                    fig['ax']['data']['right'].append(get_axvardata(fig, var['name'], self.cf.get(fig['section'], var['name']), 'right'))
                    self.get_variable_data(
                        fig, var['name'], self.cf.get(fig['section'], var['name'])
                    )
        elif fig['type'] == "Violin":
            fig['ax']['data'] = []
            for var in fig['varlist']:
                self.get_variable_data(
                    fig, var['name'], self.cf.get(fig['section'], var['name']))
            fig['limvar'] = {
                "min":  [],
                "max":  []
            }
            for kk in fig['var'].keys():
                fig['limvar']['min'].append(fig['var'][kk].min())
                fig['limvar']['max'].append(fig['var'][kk].max())
            fig['var']['lim'] = {}
            fig['var']['lim']['y'] = [
                min(fig['limvar']['min']), max(fig['limvar']['max'])]

    def savefig(self, fig, plt):
        from matplotlib.backends.backend_pdf import PdfPages
        support_fmt_list = ['ps', 'eps', 'pdf', 'pgf', 'png', 'raw',
                            'rgba', 'svg', 'svgz', 'jpg', 'jpeg', 'tif', 'tiff']
        if self.cf.has_option("PLOT_CONFI", 'save_format'):
            save_format = self.cf.get("PLOT_CONFI", 'save_format').split(",")
        else:
            save_format = ['pdf']
        unsupport = []
        support = []
        file_list = []
        for fmt in save_format:
            if fmt.strip().lower() in support_fmt_list:
                support.append(fmt.strip().lower())
                file_list.append("'{}.{}'".format(
                    fig['name'], fmt.strip().lower()))
            else:
                unsupport.append("'*.{}'".format(fmt.strip().lower()))
        if len(support) == 0:
            support = ['pdf']
            file_list = ["'{}.{}'".format(fig['name'], 'pdf')]
            print("\tTimer: {:.2f} Second;  Message from '{}' -> No support picture format found in configure file! Default format '*.pdf' used in this plot".format(
                time.time()-fig['start'], fig['section']))
        for fmt in support:
            if (fmt == 'ps'):
                pp = plt
                pp.savefig("{}.pdf".format(fig['file']), format='pdf')
                self.compress_figure_to_PS(fig['file'])
            elif fmt == "pdf" and ('ps' not in save_format):
                pp = plt
                pp.savefig("{}.{}".format(fig['file'], fmt))
            else:
                pp = plt
                pp.savefig("{}.{}".format(fig['file'], fmt), dpi=300)

        if ('pdf' not in support) and ('ps' in support):
            os.remove("{}.pdf".format(fig['file']))

        print("\tTimer: {:.2f} Second;  Figure {} saved in the path\n\t\t-> {} \n\t\t>> {}.".format(
            time.time()-fig['start'], fig['name'], os.path.dirname(fig['file']), ", >> ".join(file_list)))
        if unsupport:
            print(emoji.emojize('\t:ghost::ghost::ghost: Figure format unsupport -> {}. '.format(
                ", ".join(unsupport)), use_aliases=True))

    def drawBestPoint(self, fig, ax):
        import matplotlib.patches as patches
        theta = np.linspace(0, math.pi, 50)
        xx1 = fig['colorset']['bestpoint']["size"] * \
            (0.25*np.sin(2*theta)) + fig['var']['AxesBP']['x']
        yy1 = fig['colorset']['bestpoint']["size"] * \
            (np.sin(theta)) + fig['var']['AxesBP']['y']
        path1 = self.get_path(xx1, yy1)
        patch1 = patches.PathPatch(
            path1, facecolor=fig['colorset']['bestpoint']['color'][0], lw=0, transform=ax.transAxes, zorder=2000)
        xx2 = fig['colorset']['bestpoint']["size"] * \
            (0.13*np.cos(2*theta)) + fig['var']['AxesBP']['x']
        yy2 = fig['colorset']['bestpoint']["size"] * \
            (0.157*np.sin(2*theta) + 0.7) + fig['var']['AxesBP']['y']
        path2 = self.get_path(xx2, yy2)
        patch2 = patches.PathPatch(
            path2, facecolor=fig['colorset']['bestpoint']['color'][1], lw=0, transform=ax.transAxes, zorder=2001)
        ax.add_patch(patch1)
        ax.add_patch(patch2)

    def get_path(self, xx, yy):
        from matplotlib.path import Path
        verts, codes = [], []
        for ii in range(len(xx)):
            verts.append((xx[ii], yy[ii]))
            codes.append(Path.LINETO)
        codes[0] = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY

        path = Path(verts, codes)
        return path

    def find_critical_prob(self, data, prob):
        norm = data['pdf'].sum()
        temp_pdf = data.sort_values(by=['pdf'], ascending=True)
        temp_pdf = temp_pdf.reset_index()
        temp_pdf.pop("index")
        idx = {
            'min':  0,
            'upp':  temp_pdf.shape[0]-1,
            'id':   temp_pdf.shape[0]//2,
            'flag': True
        }
        while idx['flag']:
            if idx['upp'] - idx['min'] < 2:
                idx['flag'] = False
            elif temp_pdf['pdf'][idx['id']:].sum() / norm > prob:
                idx['min'] = idx['id']
                idx['id'] = (idx['upp']+idx['min'])//2
            else:
                idx['upp'] = idx['id']
                idx['id'] = (idx['upp']+idx['min'])//2
        return temp_pdf['pdf'][idx['id']]

    def set_Ternary_label(self, fig, ax):
        def setlabel(pa, pb, pc, label, sty):
            posx = 0.5*(pa[0]+pb[0])
            posy = 0.5*(pa[1]+pb[1])
            posx = posx + sty['label']['offline']*(posx - pc[0])
            posy = posy + sty['label']['offline']*(posy - pc[1])
            if fig['type'] == "TernaryRGB_Scatter" and label == "bottom":
                lb = "{} $+$ {} $+$ {}".format(self.cf.get(fig['section'], 'r_label'), self.cf.get(
                    fig['section'], 'g_label'), self.cf.get(fig['section'], 'b_label'))
                fig['cv'].text(
                    posx,
                    posy,
                    r"{}".format(lb),
                    fontsize=sty['label']['fontsize'],
                    color=sty['label']['color'],
                    rotation=-1.0 * sty['label']['rotation'][label],
                    horizontalalignment="center",
                    verticalalignment="center"
                )
            else:
                fig['cv'].text(
                    posx,
                    posy,
                    r"{}".format(self.cf.get(
                        fig['section'], "{}_label".format(label))),
                    fontsize=sty['label']['fontsize'],
                    color=sty['label']['color'],
                    rotation=-1.0 * sty['label']['rotation'][label],
                    horizontalalignment="center",
                    verticalalignment="center"
                )

        axs = "{}size".format(ax)
        pointA = [fig['cf'][axs]['x'], fig['cf'][axs]['y']]
        pointB = [fig['cf'][axs]['x'] + fig['cf']
                  [axs]['width'], fig['cf'][axs]['y']]
        pointC = [fig['cf'][axs]['x']+0.5*fig['cf'][axs]['width'],
                  fig['cf'][axs]['y']+fig['cf'][axs]['height']]
        axsty = fig['cf'][fig['cf'][axs]["axisstyle"]]
        setlabel(pointA, pointB, pointC, "bottom", axsty)
        setlabel(pointB, pointC, pointA, "right", axsty)
        setlabel(pointC, pointA, pointB, "left", axsty)

    def load_ternary_configure(self, fig):
        fig['cf'] = fig['colorset']
        fig['fig'].set_figheight(fig['cf']['figureSize']['height'])
        fig['fig'].set_figwidth(fig['cf']['figureSize']['width'])

    def makecanvas(self, fig):
        fig['fig'].set_figheight(fig['cf']['figureSize']['height'])
        fig['fig'].set_figwidth(fig['cf']['figureSize']['width'])
        fig['cv'] = fig['fig'].add_axes([0., 0., 1., 1.])
        fig['cv'].axis("off")
        fig['cv'].set_xlim(0, 1)
        fig['cv'].set_ylim(0, 1)
        fig['ax'] = {}

    def init_ternary(self, fig, ax):
        def ticks(axs, spl, epl, info):
            axsty = fig['cf'][fig['cf'][axs]['axisstyle']]['ticks']
            if eval(axsty['switch']):
                xl = np.linspace(spl[0], epl[0], fig['cf'][axs]["multiple"]+1)
                yl = np.linspace(spl[1], epl[1], fig['cf'][axs]["multiple"]+1)
                from math import sin, cos, radians
                for ii in range(fig['cf'][axs]["multiple"]+1):
                    fig['cv'].plot(
                        [xl[ii], xl[ii]+axsty['length'] *
                            cos(radians(axsty['angle'][info]))],
                        [yl[ii], yl[ii]+axsty['length'] *
                            sin(radians(axsty['angle'][info]))],
                        '-',
                        linewidth=axsty['width'],
                        color=axsty['color'],
                        solid_capstyle=axsty['end'],
                        zorder=10000
                    )
                    fig['cv'].text(
                        xl[ii]+axsty['length'] *
                        cos(radians(axsty['angle'][info])) +
                        axsty["labelshift"][info][0],
                        yl[ii]+axsty['length'] *
                        sin(radians(axsty['angle'][info])) +
                        axsty["labelshift"][info][1],
                        r"${:.1f}$".format(
                            fig['cf'][axs]['scale'][0] + ii/fig['cf'][axs]["multiple"] * fig['cf'][axs]['scale'][1]),
                        fontsize=axsty['fontsize'],
                        rotation=axsty['fontangle'][info],
                        horizontalalignment=axsty["labelshift"][info][2],
                        verticalalignment=axsty["labelshift"][info][3]
                    )

        def gridline(axs, pa, pb, pc, info):
            axsty = fig['cf'][fig['cf'][axs]['axisstyle']]['gridline']
            if eval(axsty['switch']):
                xd = np.linspace(pa[0], pb[0], fig['cf'][axs]['multiple']+1)
                yd = np.linspace(pa[1], pb[1], fig['cf'][axs]['multiple']+1)
                xe = np.linspace(pc[0], pb[0], fig['cf'][axs]['multiple']+1)
                ye = np.linspace(pc[1], pb[1], fig['cf'][axs]['multiple']+1)
            for ii in range(fig['cf'][axs]['multiple']-1):
                plt.plot(
                    [xd[ii+1], xe[ii+1]],
                    [yd[ii+1], ye[ii+1]],
                    axsty[info]['linestyle'],
                    linewidth=axsty[info]['linewidth'],
                    color=axsty[info]['color'],
                    zorder=9999,
                    transform=fig['cv'].transAxes,
                    alpha=axsty['alpha']
                )

        axs = "{}size".format(ax)
        fig['ax'][ax] = fig['fig'].add_axes([
            fig['cf'][axs]['x'],
            fig['cf'][axs]['y'],
            fig['cf'][axs]['width'],
            fig['cf'][axs]['height']
        ])
        fig['ax'][ax].axis('off')

        fig['cv'].fill(
            [fig['cf'][axs]['x'], fig['cf'][axs]['x']+fig['cf'][axs]['width'], fig['cf'][axs]['x']+0.5 *
                fig['cf'][axs]['width'], fig['cf'][axs]['x'], fig['cf'][axs]['x']+0.5*fig['cf'][axs]['width']],
            [fig['cf'][axs]['y'], fig['cf'][axs]['y'], fig['cf'][axs]['y'] +
                fig['cf'][axs]['height'], fig['cf'][axs]['y'], fig['cf'][axs]['y']],
            facecolor=fig['cf'][fig['cf'][axs]
                                ["axisstyle"]]['bkground']['bkgcolor'],
            zorder=9999
        )
        fig['cv'].plot(
            [fig['cf'][axs]['x'], fig['cf'][axs]['x']+fig['cf'][axs]['width'], fig['cf'][axs]['x']+0.5 *
                fig['cf'][axs]['width'], fig['cf'][axs]['x'], fig['cf'][axs]['x']+0.5*fig['cf'][axs]['width']],
            [fig['cf'][axs]['y'], fig['cf'][axs]['y'], fig['cf'][axs]['y'] +
                fig['cf'][axs]['height'], fig['cf'][axs]['y'], fig['cf'][axs]['y']],
            '-',
            linewidth=fig['cf'][fig['cf'][axs]
                                ["axisstyle"]]['axisline']['width'],
            color=fig['cf'][fig['cf'][axs]["axisstyle"]]['axisline']['color'],
            solid_joinstyle='miter',
            transform=fig['cv'].transAxes,
            zorder=10000
        )
        pointA = [fig['cf'][axs]['x'], fig['cf'][axs]['y']]
        pointB = [fig['cf'][axs]['x'] + fig['cf']
                  [axs]['width'], fig['cf'][axs]['y']]
        pointC = [fig['cf'][axs]['x']+0.5*fig['cf'][axs]['width'],
                  fig['cf'][axs]['y']+fig['cf'][axs]['height']]
        ticks(axs, pointC, pointA, "left")
        ticks(axs, pointB, pointC, "right")
        ticks(axs, pointA, pointB, "bottom")
        gridline(axs, pointC, pointA, pointB, "left")
        gridline(axs, pointB, pointC, pointA, "right")
        gridline(axs, pointA, pointB, pointC, "bottom")
        fig['ax'][ax].set_xlim(fig['cf'][axs]['scale'][0],
                               fig['cf'][axs]['scale'][1])
        fig['ax'][ax].set_ylim(fig['cf'][axs]['scale'][0],
                               fig['cf'][axs]['scale'][1])

    def scatter_classify_data(self, fig, bo):
        x_sel = self.var_symbol(bo)
        for x in x_sel:
            bo = bo.replace("_{}".format(x), "fig['data']['{}']".format(x))
        if "&FC_" in bo:
            for ii, func in enumerate(self.funcs):
                bo = bo.replace(
                    func['name'], "self.funcs[{}]['expr']".format(ii))
        fig['classify'] = {}
        fig['classify']['data'] = fig['data'][eval(bo)].reset_index()
        x_info = self.cf.get(fig['section'], 'x_variable')
        if x_info[0: 3] == '&Eq':
            x_info = x_info[4:]
            if "&V" in x_info:
                x_info = self.get_defined_var(fig, x_info, "classify']['data")
            x_sel = self.var_symbol(x_info)
            for x in x_sel:
                x_info = x_info.replace("_{}".format(
                    x), "fig['classify']['data']['{}']".format(x))
            if "&FC_" in x_info:
                for ii, func in enumerate(self.funcs):
                    x_info = x_info.replace(
                        func['name'], "self.funcs[{}]['expr']".format(ii))
            fig['classify']['x'] = eval(x_info)
        elif x_info in fig['classify']['data'].columns.values:
            fig['classify']['x'] = fig['classify']['data'][x_info]
        y_info = self.cf.get(fig['section'], 'y_variable')
        if y_info[0: 3] == '&Eq':
            y_info = y_info[4:]
            if "&V" in y_info:
                y_info = self.get_defined_var(fig, y_info, "classify']['data")
            x_sel = self.var_symbol(y_info)
            for x in x_sel:
                y_info = y_info.replace("_{}".format(
                    x), "fig['classify']['data']['{}']".format(x))
            if "&FC_" in y_info:
                for ii, func in enumerate(self.funcs):
                    y_info = y_info.replace(
                        func['name'], "self.funcs[{}]['expr']".format(ii))
            fig['classify']['y'] = eval(y_info)
        elif y_info in fig['classify']['data'].columns.values:
            fig['classify']['y'] = fig['classify']['data'][y_info]
        if self.cf.has_option(fig['section'], 'c_variable'):
            c_info = self.cf.get(fig['section'], 'c_variable')
            if c_info[0: 3] == '&Eq':
                c_info = c_info[4:]
                if "&V" in c_info:
                    c_info = self.get_defined_var(
                        fig, c_info, "classify']['data")
                x_sel = self.var_symbol(c_info)
                for x in x_sel:
                    c_info = c_info.replace("_{}".format(
                        x), "fig['classify']['data']['{}']".format(x))
                if "&FC_" in c_info:
                    for ii, func in enumerate(self.funcs):
                        c_info = c_info.replace(
                            func['name'], "self.funcs[{}]['expr']".format(ii))
                fig['classify']['c'] = eval(c_info)
            elif c_info in fig['classify']['data'].columns.values:
                fig['classify']['c'] = fig['classify']['data'][c_info]

    def Ternary_classify_data(self, fig, bo):
        def var_normalize(l, r, b):
            res = pd.DataFrame({
                "left":     l,
                "right":    r,
                "bottom":   b
            })
            res['sum'] = res['left'] + res['right'] + res['bottom']
            res['left'] = res['left']/res['sum']
            res['right'] = res['right']/res['sum']
            res['bottom'] = res['bottom']/res['sum']
            return res

        x_sel = self.var_symbol(bo)
        for x in x_sel:
            bo = bo.replace("_{}".format(x), "fig['data']['{}']".format(x))
        if "&FC_" in bo:
            for ii, func in enumerate(self.funcs):
                bo = bo.replace(
                    func['name'], "self.funcs[{}]['expr']".format(ii))
        fig['classify'] = {}
        fig['classify']['data'] = fig['data'][eval(bo)].reset_index()

        if fig['type'] == "Ternary_Scatter":
            variable_list = ['left', 'right', 'bottom']
        elif fig['type'] == "TernaryC_Scatter":
            variable_list = ['left', "right", "bottom", "c"]
        for varname in variable_list:
            var_info = self.cf.get(
                fig['section'], '{}_variable'.format(varname)).strip()
            if var_info[0: 3] == "&Eq":
                var_info = var_info[4:]
                x_sel = self.var_symbol(var_info)
                for x in x_sel:
                    var_info = var_info.replace("_{}".format(
                        x), "fig['classify']['data']['{}']".format(x))
                if "&FC_" in var_info:
                    for ii, func in enumerate(self.funcs):
                        var_info = var_info.replace(
                            func['name'], "self.funcs[{}]['expr']".format(ii))
                fig['classify'][varname] = eval(var_info)
            elif var_info in fig['classify']['data'].columns.values:
                fig['classify'][varname] = fig['classify']['data'][varname]

        fig['classify']['data'] = var_normalize(
            fig['classify']['left'], fig['classify']['right'], fig['classify']['bottom'])
        fig['classify']['x'] = fig['classify']['data']['bottom'] + \
            0.5 * fig['classify']['data']['right']
        fig['classify']['y'] = fig['classify']['data']['right']

    def get_Linestyle(self, fig, style):
        # style_file = self.load_path(self.cf.get("COLORMAP", "StyleSetting"))
        style_file = self.load_path(fig['colorset']['linestylepath'])
        style = style.strip()
        if (not style[0] == '&') or (not os.path.exists(style_file)):
            print(
                "\tLine info Error: Unvaliable line format {} \n\t Default Line Format used".format(style))
            return {'width':    2., 'color':    '#b71c1c', 'style':    '-', 'alpha':    1.0, 'marker':   None, 'markersize':   5}
        else:
            style = style.split('_')
            style = {
                'name':     style[0].replace("&", ''),
                'label':    style[1].strip()
            }
            with open(style_file, 'r') as f1:
                default_style = json.loads(f1.read())
            style_tag = False
            for n, v in default_style.items():
                if style['name'] == n:
                    for item in v:
                        if style['label'] == item['label']:
                            style.update(item['LineStyle'])
                            style_tag = True
                            break
            if style_tag:
                return style
            else:
                print(
                    "\tLine info Error: Unvaliable line format {} \n\t Default Line Format used".format(style))
                return {'width':    2., 'color':    '#b71c1c', 'style':    '-', 'alpha':    1.0, 'marker':   None, 'markersize':   5}

    def drawline(self, fig, ax):
        def get_variable(info):
            info = info.replace('{', '').replace('}', '').strip()
            info = info.split('|')
            info = {
                'varname':  info[0].strip(),
                'vardata':  np.linspace(float(info[1].strip()), float(info[2].strip()), int(info[3].strip()))
            }
            return info

        def get_data(expr, var):
            for ii, fc in enumerate(self.funcs):
                expr = expr.replace(
                    fc['name'], "self.funcs[{}]['expr']".format(ii))
            y = []
            for ii in var['vardata']:
                y.append(
                    eval(expr.replace('_{}'.format(var['varname']), str(ii))))
            return np.array(y)

        def get_function(var, func):
            var = get_variable(var)
            decode = re.compile(r'[{](.*?)[}]', re.S)
            func = re.findall(decode, func)
            for ii, fc in enumerate(func):
                line = fc.split(':')
                func[ii] = {
                    'varname':  line[0].strip(),
                    'vardata':  get_data(line[1], var)
                }
            func.append(var)
            return func

        def get_line_info(line):
            info = line.split(',')
            if info[0].strip() == 'parametric':
                res = {}
                res['method'] = info.pop(0)
                res['zorder'] = float(info.pop(0).strip())
                res['var'] = info.pop(0)
                res['style'] = self.get_Linestyle(fig, info.pop(-1).strip())
                res['Func'] = ','.join(info)
                res['var'] = get_function(res['var'], res['Func'])
                xtag, ytag = False, False
                for ii, var in enumerate(res['var']):
                    if var['varname'] == 'x':
                        xx = var['vardata']
                        xtag = True
                    if var['varname'] == 'y':
                        yy = var['vardata']
                        ytag = True
                if xtag and ytag:
                    res['data'] = pd.DataFrame({
                        'x':  xx,
                        'y':  yy
                    })
                elif not xtag:
                    print(
                        "\tLine Info Error: No x coordinate founded in Line Setting\n\t{}".format(line))
                elif not ytag:
                    print(
                        "\tLine Info Error: No y coordinate founded in Line Setting\n\t{}".format(line))
                return res
            elif info[0].strip() == "Equation":
                res = {}
                res['method'] = info.pop(0)
                res['zorder'] = float(info.pop(0).strip())
                res['var'] = info.pop(0)
                res['style'] = self.get_Linestyle(fig, info.pop(-1))
                res['Func'] = ','.join(info)
                res['var'] = get_function(res['var'], res['Func'])
                xtag, ytag = False, False
                for ii, var in enumerate(res['var']):
                    if var['varname'] == 'x':
                        xx = var['vardata']
                        xtag = True
                    if var['varname'] == 'y':
                        yy = var['vardata']
                        ytag = True
                if xtag and ytag:
                    res['data'] = pd.DataFrame({
                        'x':  xx,
                        'y':  yy
                    })
                elif not xtag:
                    print(
                        "\tLine Info Error: No x coordinate founded in Line Setting\n\t{}".format(line))
                elif not ytag:
                    print(
                        "\tLine Info Error: No y coordinate founded in Line Setting\n\t{}".format(line))
                return res
            else:
                print(
                    "Line Drawing Error: No such line function methed {}\n\t-> {}".format(info[0], line))
                sys.exit(1)

        def draw(lineinfo):
            if lineinfo['style']['marker'] == 'None':
                ax.plot(lineinfo['data']['x'], lineinfo['data']['y'], linewidth=lineinfo['style']['width'],
                        color=lineinfo['style']['color'], linestyle=lineinfo['style']['style'], zorder=lineinfo['zorder'])

        fig['lineinfo'] = self.cf.get(fig['section'], 'Line_draw').split('\n')
        for ii, line in enumerate(fig['lineinfo']):
            draw(get_line_info(line))

    def get_TextStyle(self, style):
        style_file = self.load_path(self.cf.get("COLORMAP", "StyleSetting"))
        style = style.strip()
        if style[0] != '&' and type(eval(style)) == dict:
            style = eval(style)
            kys = ['FontSize', 'color', 'alpha']
            style_tag = True
            for k in kys:
                if 'FontSize' not in style.keys():
                    style_tag = False
            if style_tag:
                return style
            else:
                print(
                    "\tLine info Error: Unvaliable Text format {} \n\t Default Line Text Format used".format(style))
                return {"FontSize": 20, "color": "#b71c1c", "alpha": 1.0}
        elif (not style[0] == '&') or (not os.path.exists(style_file)):
            print(
                "\tLine info Error: Unvaliable Text format {} \n\t Default Line Text Format used".format(style))
            return {"FontSize": 20, "color": "#b71c1c", "alpha": 1.0}
        else:
            style = style.split('_')
            style = {
                'name':     style[0].replace("&", '').strip(),
                'label':    style[1].strip()
            }
            with open(style_file, 'r') as f1:
                default_style = json.load(f1)
            style_tag = False
            for n, v in default_style.items():
                if style['name'] == n:
                    for item in v:
                        if style['label'] == item['label']:
                            style = item['TextStyle']
                            style_tag = True
                            break
            if style_tag:
                return style
            else:
                print(
                    "\tLine info Error: Unvaliable Text format {} \n\t Default Line Text Format used".format(style))
                return {"FontSize": 20, "color": "#b71c1c", "alpha": 1.0}

    def drawtext(self, fig, ax):
        def get_text_info(line):
            decode = re.compile(r'[(](.*?)[)]', re.S)
            pos = re.findall(decode, line)[0]
            line = line.replace('({})'.format(pos), '').strip()
            pos = pos.split(',')
            line = line.lstrip('|')
            line = line.split('|')
            res = {
                'pos':      [float(pos[0].strip()), float(pos[1].strip())],
                'rotation': float(line.pop(0).strip()),
                'style':    self.get_TextStyle(line.pop(-1).strip()),
                'text':     '|'.join(line).strip()
            }
            return res

        def draw(info):
            ax.text(info['pos'][0], info['pos'][1], r"{}".format(info['text']), fontsize=info['style']['FontSize'], color=info['style']['color'],
                    alpha=info['style']['alpha'], rotation=info['rotation'], horizontalalignment='left', verticalalignment='bottom', zorder=1999)

        fig['textinfo'] = self.cf.get(fig['section'], 'Text').split('\n')
        for ii, line in enumerate(fig['textinfo']):
            draw(get_text_info(line))

    def compress_figure_to_PS(self, figpath):
        os.system('pdf2ps {}.pdf {}.ps'.format(figpath, figpath))

    def ax_setcmap(self, fig):
        if self.cf.has_option(fig['section'], 'colorset'):
            cname = self.cf.get(fig['section'], 'colorset')
            if cname[0] == '&':
                for ii in self.colors:
                    if cname[1:] == ii['name']:
                        fig['colorset'] = ii
                        break

    def ax_setlim(self, fig, label):
        fig['ax']['lim'] = {}
        for aa in label:
            tem = self.cf.get(fig['section'], '{}_lim'.format(aa)).split(',')
            fig['ax']['lim'][aa] = tem
            for it in tem:
                if ('AUTO' in it) and (aa in fig['var']['lim'].keys()):
                    fig['ax']['lim'][aa][tem.index(it)] = float(it.split('_')[1].strip(
                    )) * (fig['var']['lim'][aa][1] // float(it.split('_')[1].strip()) + 1)
                else:
                    fig['ax']['lim'][aa][tem.index(it)] = float(it.strip())

    def ax_setticks(self, fig, axislabel):
        fig['ax']['ticks'] = {}
        for aa in axislabel:
            tick = self.cf.get(fig['section'], '{}_ticks'.format(aa))
            if 'AUTO' in tick:
                a = float(tick.split('_')[-1])
                low = fig['ax']['lim'][aa][0] // a
                upp = fig['ax']['lim'][aa][1] // a + 1
                fig['ax']['ticks'][aa] = np.linspace(
                    low*a, upp*a, int(upp-low+1))
            elif tick[0:4] == 'Manu':
                p_rec = re.compile(r'[\[].*?[\]]', re.S)
                a = re.findall(p_rec, tick[5:].strip())
                tk = []
                label = []
                for it in a:
                    it = it.strip().strip('[').strip(']')
                    tk.append(float(it.split(',')[0]))
                    label.append("{}".format(it.split(',')[1].strip()))
                fig['ax']['ticks'][aa] = tuple([tk, label])
            else:
                tick = tick.split(',')
                for ii in range(len(tick)):
                    tick[ii] = tick[ii].strip()
                fig['ax']['ticks'][aa] = np.linspace(
                    float(tick[0]), float(tick[1]), int(tick[2]))
            if aa == 'y':
                if fig['ax']['lim']['y'][0] in fig['ax']['ticks']['y']:
                    fig['ax']['ticks']['y'] = fig['ax']['ticks']['y'][np.where(
                        fig['ax']['ticks']['y'] != fig['ax']['lim']['y'][0])]

    def GetStatData(self, fig):
        if self.cf.has_option(fig['section'], 'stat_variable'):
            fig['var'] = {}
            if self.cf.get(fig['section'], 'stat_variable').split(',')[0] == 'CHI2':
                self.get_variable_data(fig, 'x', self.cf.get(
                    fig['section'], 'x_variable'))
                self.get_variable_data(fig, 'y', self.cf.get(
                    fig['section'], 'y_variable'))
                self.get_variable_data(fig, 'Stat', ",".join(self.cf.get(
                    fig['section'], 'stat_variable').split(',')[1:]).strip())
                fig['var']['data'] = pd.DataFrame({
                    'x':    fig['var']['x'],
                    'y':    fig['var']['y'],
                    'Stat': fig['var']['Stat']})
                fig['var'].pop('x')
                fig['var'].pop('y')
                fig['var'].pop('Stat')

    def Get1DStatData(self, fig):
        if self.cf.has_option(fig['section'], "stat_variable"):
            fig['var'] = {}
            self.get_variable_data(fig, 'x', self.cf.get(
                fig['section'], 'x_variable'))
            for line in self.cf.get(fig['section'], 'stat_variable').split('\n'):
                if line.split(",")[0].strip() == 'CHI2':
                    self.get_variable_data(
                        fig, 'CHI2', ','.join(line.split(",")[1:]).strip())
                if line.split(',')[0].strip() == 'PDF':
                    self.get_variable_data(
                        fig, 'PDF', ','.join(line.split(',')[1:]).strip())
            if 'CHI2' in fig['var'].keys() and 'PDF' in fig['var'].keys():
                fig['var']['data'] = pd.DataFrame({
                    "x":    fig['var']['x'],
                    "CHI2": fig['var']['CHI2'],
                    "PDF":  fig['var']['PDF']
                })
                fig['var'].pop('x')
                fig['var'].pop("PDF")
                fig['var'].pop('CHI2')
                fig['var']['type'] = ["x", "PDF", "CHI2"]
            elif 'CHI2' in fig['var'].keys():
                fig['var']['data'] = pd.DataFrame({
                    "x":    fig['var']['x'],
                    "CHI2": fig['var']['CHI2']
                })
                fig['var'].pop('x')
                fig['var'].pop('CHI2')
                fig['var']['type'] = ['x', 'CHI2']
            elif 'PDF' in fig['var'].keys():
                fig['var']['data'] = pd.DataFrame({
                    "x":    fig['var']['x'],
                    "PDF":  fig['var']['PDF']
                })
                fig['var'].pop('x')
                fig['var'].pop("PDF")
                fig['var']['type'] = ['x', 'PDF']

    def getTernaryRGBData(self, fig):
        def var_normalize(a, b, c, d, e):
            res = pd.DataFrame({
                'left':   a,
                'right':   b,
                'r':   c,
                'g':   d,
                'b':   e
            })
            res['sum'] = res['left'] + res['right'] + \
                res['r'] + res['g'] + res['b']
            res['left'] = res['left']/res['sum']
            res['right'] = res['right']/res['sum']
            res['r'] = res['r']/res['sum']
            res['g'] = res['g']/res['sum']
            res['b'] = res['b']/res['sum']
            res['bottom'] = res['r'] + res['g'] + res['b']
            maxrgb = max(max(res['r']), max(res['g']))
            res['r'] = np.power(res['r'], 0.25)
            res['g'] = np.power(res['g'], 0.25)
            res['b'] = np.power(res['b'], 0.25)

            return res

        if self.cf.has_option(fig['section'], 'left_variable') & self.cf.has_option(fig['section'], 'right_variable') & self.cf.has_option(fig['section'], 'r_variable') & self.cf.has_option(fig['section'], 'g_variable') & self.cf.has_option(fig['section'], 'b_variable'):
            fig['var'] = {}
            self.get_variable_data(fig, 'left', self.cf.get(
                fig['section'], 'left_variable'))
            self.get_variable_data(fig, 'right', self.cf.get(
                fig['section'], 'right_variable'))
            self.get_variable_data(fig, 'r', self.cf.get(
                fig['section'], 'r_variable'))
            self.get_variable_data(fig, 'g', self.cf.get(
                fig['section'], 'g_variable'))
            self.get_variable_data(fig, 'b', self.cf.get(
                fig['section'], 'b_variable'))

            fig['var']['oridata'] = var_normalize(
                fig['var']['left'], fig['var']['right'], fig['var']['r'], fig['var']['g'], fig['var']['b'])
            fig['var']['axdata'] = pd.DataFrame({
                'x':    fig['var']['oridata']['bottom'] + 0.5*fig['var']['oridata']['right'],
                'y':    fig['var']['oridata']['right']
            })
            fig['var']['axdata']['c'] = fig['var']['oridata'].apply(
                lambda x: tuple([0.9*x['r']+0.1, 0.1+0.9*x['g'], 0.1+0.9*x['b']]), axis=1)
            fig['var']['axdata']['r'] = fig['var']['oridata'].apply(
                lambda x: tuple([x['r'], 0.38*x['r'], 0.38*x['r']]), axis=1)
            fig['var']['axdata']['g'] = fig['var']['oridata'].apply(
                lambda x: tuple([0.44*x['g'], x['g'], 0.44*x['g']]), axis=1)
            fig['var']['axdata']['b'] = fig['var']['oridata'].apply(
                lambda x: tuple([0.42*x['b'], 0.42*x['b'], x['b']]), axis=1)

    def getTernaryData(self, fig):
        def var_normalize(l, r, b):
            res = pd.DataFrame({
                "left":     l,
                "right":    r,
                "bottom":   b
            })
            res['sum'] = res['left'] + res['right'] + res['bottom']
            res['left'] = res['left']/res['sum']
            res['right'] = res['right']/res['sum']
            res['bottom'] = res['bottom']/res['sum']
            return res

        if self.cf.has_option(fig['section'], "left_variable") and self.cf.has_option(fig['section'], "right_variable") and self.cf.has_option(fig['section'], "bottom_variable"):
            fig['var'] = {}
            self.get_variable_data(fig, 'left', self.cf.get(
                fig['section'], "left_variable"))
            self.get_variable_data(fig, 'right', self.cf.get(
                fig['section'], "right_variable"))
            self.get_variable_data(fig, 'bottom', self.cf.get(
                fig['section'], "bottom_variable"))

            fig['var']['oridata'] = var_normalize(
                fig['var']['left'], fig['var']['right'], fig['var']['bottom'])
            fig['var']['axdata'] = pd.DataFrame({
                "x":    fig['var']['oridata']['bottom'] + 0.5*fig['var']['oridata']['right'],
                "y":    fig['var']['oridata']['right']
            })

    def get_defined_var(self, fig, vars, ty):
        def set_defined_var_data(var):
            vainfo = self.cf.get("VARIBLES", var).strip()
            if vainfo.split(',')[0].strip() == "Sort":
                vainfo = ','.join(vainfo.split(',')[1:]).strip()
                decode = re.compile(r'[\[](.*?)[\]]', re.S)
                dset = re.findall(decode, vainfo)
                vdd = {
                    "name": var,
                    "sourceinfo":   dset[0],
                    "sourcedata":   {},
                    "targetdata":   {}
                }
                if len(dset) == 2:
                    vdd['target'] = dset[1]
                    vdd['type'] = "scalar"
                    vdd['code'] = int(vainfo.split(',')[-1].strip())-1
                elif len(dset) == 1:
                    vdd['target'] = dset[0]
                    vdd['type'] = "vector"
                else:
                    print("\tInvaliade variable for {}!!\n".format(var))
                    sys.exit(0)
                if vdd['type'] == "vector":
                    for ii in range(len(vdd['sourceinfo'].split(','))):
                        item = vdd['sourceinfo'].split(",")[ii].strip()
                        vdd['sourcedata']["a{}".format(
                            ii)] = self.get_listData(fig, item)
                    vdd['sourcedata'] = pd.DataFrame(vdd['sourcedata'])
                    for ii in range(len(vdd['sourceinfo'].split(','))):
                        fig['data']["@V{}{}".format(
                            var, ii+1)] = vdd["sourcedata"].apply(lambda x: np.sort(x)[ii], axis=1)
                elif vdd['type'] == "scalar":
                    lopt = len(vdd['sourceinfo'].split(','))
                    for ii in range(lopt):
                        sitem = vdd['sourceinfo'].split(',')[ii].strip()
                        titem = vdd['target'].split(',')[ii].strip()
                        vdd['sourcedata'][ii] = self.get_listData(fig, sitem)
                        vdd['sourcedata'][ii +
                                          lopt] = self.get_listData(fig, titem)
                    vdd['sourcedata'] = pd.DataFrame(vdd['sourcedata'])
                    fig['data']["@V{}".format(var)] = vdd['sourcedata'].apply(
                        lambda x: x[list(x).index(np.sort(x[0:lopt])[vdd['code']])+lopt], axis=1)

        def match_variable(fig, vars, ty):
            vlist = vars.split()
            for vv in vlist:
                if vv[0:2] == "&V":
                    va = "@V"+vv[2:]
                    if not va in fig['data'].columns:
                        for opt in self.cf.options("VARIBLES"):
                            if "&V{}".format(opt) in vv:
                                set_defined_var_data(opt)
                    vars = vars.replace(vv, "fig['{}']['{}']".format(ty, va))
            return vars

        if not self.cf.has_section("VARIBLES"):
            print(emoji.emojize(
                '\t:ghost::ghost::ghost: Section [VARIBLES] not found!!\n', use_aliases=True))
        vars = match_variable(fig, vars, ty)
        return vars

    def getTernaryCData(self, fig):
        def var_normalize(l, r, b):
            res = pd.DataFrame({
                "left":     l,
                "right":    r,
                "bottom":   b
            })
            res['sum'] = res['left'] + res['right'] + res['bottom']
            res['left'] = res['left']/res['sum']
            res['right'] = res['right']/res['sum']
            res['bottom'] = res['bottom']/res['sum']
            return res

        if self.cf.has_option(fig['section'], "left_variable") and self.cf.has_option(fig['section'], "right_variable") and self.cf.has_option(fig['section'], "bottom_variable"):
            fig['var'] = {}
            self.get_variable_data(fig, 'left', self.cf.get(
                fig['section'], "left_variable"))
            self.get_variable_data(fig, 'right', self.cf.get(
                fig['section'], "right_variable"))
            self.get_variable_data(fig, 'bottom', self.cf.get(
                fig['section'], "bottom_variable"))
            self.get_variable_data(fig, 'color', self.cf.get(
                fig['section'], 'c_variable'))

            fig['var']['oridata'] = var_normalize(
                fig['var']['left'], fig['var']['right'], fig['var']['bottom'])
            fig['var']['axdata'] = pd.DataFrame({
                "x":    fig['var']['oridata']['bottom'] + 0.5*fig['var']['oridata']['right'],
                "y":    fig['var']['oridata']['right'],
                "c":    fig['var']['color']
            })

    def Get2DData(self, fig):
        fig['var'] = {}
        self.get_variable_data(fig, 'x', self.cf.get(
            fig['section'], 'x_variable'))
        self.get_variable_data(fig, 'y', self.cf.get(
            fig['section'], 'y_variable'))
        fig['var']['data'] = pd.DataFrame({
            'x':    fig['var']['x'],
            'y':    fig['var']['y']
        })

    def Get3DData(self, fig):
        fig['var'] = {}
        self.get_variable_data(fig, 'x', self.cf.get(
            fig['section'], 'x_variable'))
        self.get_variable_data(fig, 'y', self.cf.get(
            fig['section'], 'y_variable'))
        self.get_variable_data(fig, 'c', self.cf.get(
            fig['section'], 'c_variable'))
        fig['var']['data'] = pd.DataFrame({
            'x':    fig['var']['x'],
            'y':    fig['var']['y'],
            'c':    fig['var']['c']
        })

    def get_listData(self, fig, varinf):
        if varinf in fig['data'].columns.values:
            return fig['data'][varinf]
        else:
            try:
                x_sel = self.var_symbol(varinf)
                for x in x_sel:
                    varinf = varinf.replace("_{}".format(
                        x), "fig['data']['{}']".format(x))
                if "&FC_" in varinf:
                    for ii in range(len(self.funcs)):
                        varinf = varinf.replace(
                            self.funcs[ii]['name'], "self.funcs[{}]['expr']".format(ii))
                return eval(varinf)
            except:
                print("No Variable {} found in Data!".format(varinf))
                sys.exit(0)

    def get_variable_data(self, fig, name, varinf):
        if varinf[0:3] == '&Eq':
            varinf = varinf[4:]
            if "&V" in varinf:
                varinf = self.get_defined_var(fig, varinf, 'data')
            x_sel = self.var_symbol(varinf)
            for x in x_sel:
                varinf = varinf.replace("_{}".format(
                    x), "fig['data']['{}']".format(x))
            if "&FC_" in varinf:
                for ii in range(len(self.funcs)):
                    varinf = varinf.replace(
                        self.funcs[ii]['name'], "self.funcs[{}]['expr']".format(ii))
            fig['var'][name] = eval(varinf)
        elif varinf in fig['data'].columns.values:
            fig['var'][name] = fig['data'][varinf]
        else:
            print("No Variable {} found in Data!".format(varinf))
            sys.exit(0)

    def figures_inf(self):
        self.figpath = self.load_path(self.cf.get('PLOT_CONFI', 'save_dir'))
        if not os.path.exists(self.figpath):
            os.makedirs(self.figpath)
        if self.cf.has_option('PLOT_CONFI', 'plot'):
            ppt = self.cf.get('PLOT_CONFI', 'plot').split('\n')
            self.figs = []
            for line in ppt:
                pic = {}
                pic['type'] = line.split(',')[0]
                if "ALL" in ','.join(line.split(',')[1:]).upper():
                    for item in self.cf.sections():
                        if pic['type'] in item:
                            self.figs.append({
                                'name':     self.cf.get(item, 'plot_name'),
                                'section':  item,
                                'type':     pic['type']
                            })
                else:
                    for item in line.split(',')[1:]:
                        if self.cf.has_section('PLOT_{}_{}'.format(pic['type'], item.strip())):
                            sec = 'PLOT_{}_{}'.format(
                                pic['type'], item.strip())
                            self.figs.append({
                                'section':  sec,
                                'name':     self.cf.get(sec, 'plot_name'),
                                'type':     pic['type']
                            })

    def Analytic_funcs(self):
        self.funcs = []
        for item in self.cf.sections():
            if "FUNCTION1D" in item:
                fun = {}
                fun['name'] = "&FC_{}".format(self.cf.get(item, 'name'))
                data_dir = self.load_path(self.cf.get(item, 'file'))
                fun['data'] = pd.read_csv(data_dir)
                if self.cf.has_option(item, "method"):
                    method = self.cf.get(item, "method").strip().lower()
                    if method not in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', 'next']:
                        method = "linear"
                else:
                    method = 'linear'
                fill_value = ""
                if self.cf.has_option(item, "fill_value"):
                    fill_value = self.cf.get(
                        item, "fill_value").strip().lower()
                if fill_value == "extrapolate":
                    fun['expr'] = interp1d(
                        fun['data']['x'], fun['data']['y'], kind=method, fill_value=fill_value)
                else:
                    fun['expr'] = interp1d(
                        fun['data']['x'], fun['data']['y'], kind=method)
                self.funcs.append(fun)
        self.funcs = tuple(self.funcs)

    def var_symbol(self, bo):
        x = []
        bo = bo.split()
        for it in bo:
            if it[0] == '_' and it not in x:
                x.append(it[1:])
        return x

    def basic_selection(self, fig):
        if self.cf.has_option(fig['section'], 'selection'):
            fig['data'] = self.data
            bo = self.cf.get(fig['section'], 'selection')
            if bo[0:3] == '&Bo':
                if "&V" in bo:
                    bo = self.get_defined_var(fig, bo, 'data')
                bo = bo[4:]
                x_sel = self.var_symbol(bo)
                for x in x_sel:
                    bo = bo.replace("_{}".format(
                        x), "self.data['{}']".format(x))
            if "&FC_" in bo:
                for ii in range(len(self.funcs)):
                    bo = bo.replace(
                        self.funcs[ii]['name'], "self.funcs[{}]['expr']".format(ii))
            # print("Total Data is -> {} rows".format(self.data.shape[0]))
            if "*Bool*" not in self.data.columns:
                bool_list = np.ones(self.data.shape[0], dtype=np.bool)
                self.data['*Bool*'] = bool_list
                bo = bo + "& self.data['*Bool*']"
            else:
                bool_list = np.ones(self.data.shape[0], dtype=np.bool)
                self.data['abc**cba**bool**loob**alphabeta'] = bool_list
                bo = bo + "& self.data['abc**cba**bool**loob**alphabeta']"
            fig['data'] = self.data[eval(bo)].reset_index()

    def drawfillarea(self, fig, ax):
        def decode(line):
            info = line.split(',')
            if info[0].strip() == 'between':
                res = {}
                res['method'] = info.pop(0).strip()
                res['zorder'] = float(info.pop(0).strip())
                res['var'] = info.pop(0)
                # res['style'] = self.get_Linestyle(info.pop(-1))
                res['style'] = get_fillstyle(info.pop(-1))
                res['Func'] = ','.join(info)
                res['var'] = get_function(res['var'], res['Func'])
                xtag, ytag = [], []
                res["xvars"], res["yvars"] = [], []
                res['tag'] = ""
                for ii, var in enumerate(res['var']):
                    # print(ii, var['varname'])
                    if 'x' in var['varname']:
                        xtag.append(var['varname'])
                        res["xvars"].append(var)
                    if "y" in var['varname']:
                        ytag.append(var['varname'])
                        res["yvars"].append(var)
                if len(xtag) == 2 and len(ytag) == 1:
                    res['tag'] = "fillX"
                elif len(ytag) == 2 and len(xtag) == 1:
                    res['tag'] = "fillY"
                else:
                    print(
                        "\tFill Info Error: No coordinate founded in Fill Setting\n\t{}".format(line))
                # elif not ytag:
                    # print("\tLine Info Error: No y coordinate founded in Line Setting\n\t{}".format(line))

                return res
            else:
                print(
                    "Line Drawing Error: No such line function methed {}\n\t-> {}".format(info[0], line))
                sys.exit(1)

        def get_fillstyle(info):
            with open(fig['colorset']['fillstylepath'], 'r') as f1:
                fig['colorset']['fillst'] = json.loads(f1.read())
            return fig['colorset']['fillst'][info.strip()]

        def get_variable(info):
            info = info.replace('{', '').replace('}', '').strip()
            info = info.split('|')
            info = {
                'varname':  info[0].strip(),
                'vardata':  np.linspace(float(info[1].strip()), float(info[2].strip()), int(info[3].strip()))
            }
            return info

        def get_data(expr, var):
            for ii, fc in enumerate(self.funcs):
                expr = expr.replace(
                    fc['name'], "self.funcs[{}]['expr']".format(ii))
            y = []
            for ii in var['vardata']:
                y.append(
                    eval(expr.replace('_{}'.format(var['varname']), str(ii))))
            return np.array(y)

        def get_function(var, func):
            var = get_variable(var)
            decode = re.compile(r'[{](.*?)[}]', re.S)
            func = re.findall(decode, func)
            for ii, fc in enumerate(func):
                line = fc.split(':')
                func[ii] = {
                    'varname':  line[0].strip(),
                    'vardata':  get_data(line[1], var)
                }
            func.append(var)
            return func

        def draw(info):
            if info['method'] == "between":
                if info['tag'] == "fillX":
                    x0 = info['xvars'][0]
                    x1 = info['xvars'][1]
                    if int(x1['varname'].replace('x', "")) < int(x0['varname'].replace("x", "")):
                        x0, x1 = x1, x0
                    ax.fill_betweenx(info['yvars'][0]['vardata'], x0['vardata'], x1['vardata'],
                                     color=info['style']['color'], alpha=info['style']['alpha'], zorder=info['zorder'])
                elif info['tag'] == "fillY":
                    y0 = info['yvars'][0]
                    y1 = info['yvars'][1]
                    if int(y1['varname'].replace('y', "")) < int(y0['varname'].replace("y", "")):
                        y0, y1 = y1, y0
                    ax.fill_between(info['xvars'][0]['vardata'], y0['vardata'], y1['vardata'],
                                    color=info['style']['color'], alpha=info['style']['alpha'], zorder=info['zorder'])

                    # for var in info['xvars']:

        fillcf = self.cf.get(fig['section'], 'fill').split("\n")
        for ff in fillcf:
            fc = decode(ff)
            if fc['tag']:
                # if fillcf['tag'] == "fillX":
                draw(fc)


