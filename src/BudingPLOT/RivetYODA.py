#!/usr/bin/env python3
import os, sys 
import pandas as pd 
import numpy as np 
import time 
import configparser
import math 
import sympy 
from plot import Figure
import json
import matplotlib.pyplot as plt 
from matplotlib import rc, rcParams
import emoji

config = {
    "font.family":["serif", "Times New Roman"],
    "mathtext.fontset":'stix',
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

class PlotYoda1D(Figure):
    def __init__(self) -> None:
        super().__init__()
        
    def load(self):
        if self.cf.has_option("PLOT_CONFI", "save format"):
            self.format = list(map(str.strip, self.cf.get("PLOT_CONFI", 'save format').split(",")))
            
    def load_plot(self, plot):
        if "legend" in self.plots[plot]:
            self.plots[plot]['legend'] = eval(self.plots[plot]['legend'])
        else:
            self.plots[plot]['legend'] = self.cs['default']['legend'] 
        if type(self.plots[plot]['sep_scale']) == str:
            self.plots[plot]['sep_scale'] = eval(self.plots[plot]['sep_scale'])

    def decode_plot_setting(self, hname):
        hset = dict(dict(self.cf.items())[hname].items())
        return hset

    def load_yodas(self):
        with open(self.cs, 'r') as f1:
            self.cs = json.loads(f1.read())
        if "data path" in self.inf and "accumulate data" in self.inf:
            self.inf['data path'] = self.decode_path(self.inf['data path'])
            self.cumyds = {}
            self.plots = {}
            nn = 1
            from IOs import YodaFile
            for cumd in self.inf['accumulate data'].split("\n"):
                cumd = list(map(str.strip, cumd.split(",")))
                cname = "CUM{}".format(nn)
                nn += 1
                self.cumyds[cname] = {
                    "path":      self.decode_path(os.path.join(self.inf['data path'], cumd[0])),
                    "title":     cumd[1],
                    "file":      YodaFile()
                }
                self.cumyds[cname]['file'].file = self.cumyds[cname]['path']
                self.cumyds[cname]['file'].read()
                self.cumyds[cname]["hist1d"] = self.cumyds[cname]['file'].hist1ddata
                if not self.plots:
                    for hname in self.cumyds[cname]['file'].hist1ddata.keys():
                        from copy import deepcopy
                        self.plots[hname] = deepcopy(self.cs['plot_temp']['histo1d'])
        if "data path" in self.inf and "separate data" in self.inf:
            self.inf['data path'] = self.decode_path(self.inf['data path'])
            self.sepyds = {}
            from IOs import YodaFile
            nn = 1
            for sepd in self.inf['separate data'].split("\n"):
                sped = list(map(str.strip, sepd.split(",")))
                sname = "SEP{}".format(nn)
                nn += 1
                self.sepyds[sname] = {
                    "path":      self.decode_path(os.path.join(self.inf['data path'], sped[0])),
                    "title":     ", ".join(sped[1:]),
                    "file":      YodaFile()
                }
                self.sepyds[sname]['file'].file = self.sepyds[sname]['path']
                self.sepyds[sname]['file'].read()
                self.sepyds[sname]['hist1d'] = self.sepyds[sname]['file'].hist1ddata
                if not self.plots:
                    for hname in self.sepyds[sname]['file'].hist1ddata.keys():
                        from copy import deepcopy
                        self.plots[hname] = deepcopy(self.cs['plot_temp']['histo1d'])
        if "yoda setting" in self.inf:
            from IOs import YodaPlotFile
            self.ydinf = YodaPlotFile()
            self.ydinf.cs = self.cs['plot_temp']['histo1d']
            self.ydinf.file = self.decode_path(self.inf['yoda setting'])
            self.ydinf.histlist = self.plots.keys()
            self.ydinf.read()
            self.plots = deepcopy(self.ydinf.hist1d)
        for hname in self.plots:
            if self.cf.has_section(hname):
                self.plots[hname].update(self.decode_plot_setting(hname))
                self.load_plot(hname)
            
    def draw_plot(self, plot):
        self.fig = plt.figure(**self.cs["axis"]['figsize'])
        self.ax  = self.fig.add_axes(**self.cs['axis']['axsize'])
        for ii in range(len(self.cumyds)):
            oname = "CUM{}".format(ii)
            cname = "CUM{}".format(ii+1)
            ccnm  = "CUM{}".format((ii + 1) % 10)
            self.cumyds[cname]['ccnm'] = ccnm
            if not ii:
                self.cumyds[cname]['hist1d'][plot]['data']['cum'] = 0.
            else:
                self.cumyds[cname]['hist1d'][plot]['data']['cum'] = self.cumyds[oname]['hist1d'][plot]['data']['cum'] + self.cumyds[oname]['hist1d'][plot]['data']['val']
            for jj, row in self.cumyds[cname]['hist1d'][plot]['data'].iterrows():                
                xx = [row['xlow'], row['xlow'], row['xhigh'], row['xhigh'], row['xlow']]
                yy = [row['cum'], row['cum'] + row['val'], row['val'] + row['cum'], row['cum'], row['cum']]
                self.ax.fill(xx, yy, **self.cs['cum'][ccnm])
        cumtotal = self.cumyds[cname]['hist1d'][plot]['data']['cum'] + self.cumyds[cname]['hist1d'][plot]['data']['val']
        xx = np.array([np.array(self.cumyds[cname]['hist1d'][plot]['data']['xlow']), np.array(self.cumyds[cname]['hist1d'][plot]['data']['xhigh'])]).ravel(order="F")
        yy = np.array([np.array(cumtotal), np.array(cumtotal)]).ravel(order="F")
        self.ax.plot(xx, yy, **self.cs['default']['cumLine'])
        self.plots[plot]['XMin'] = self.cumyds[cname]['hist1d'][plot]['data']['xlow'].min()
        self.plots[plot]['XMax'] = self.cumyds[cname]['hist1d'][plot]['data']['xhigh'].max()
        self.xbin = self.cumyds[cname]['hist1d'][plot]['data']['xhigh'].min() - self.cumyds[cname]['hist1d'][plot]['data']['xlow'].min()
        
        for ii in range(len(self.sepyds)):
            sname = "SEP{}".format(ii + 1)
            scnm  = "SEP{}".format((ii + 1) % 10)
            self.sepyds[sname]['scnm'] = scnm 
            xx = np.array([
                np.array(self.sepyds[sname]['hist1d'][plot]['data']['xlow']), 
                np.array(self.sepyds[sname]['hist1d'][plot]['data']['xhigh'])
            ]).ravel(order="F")
            yy = np.array([
                self.plots[plot]['sep_scale'] * np.array(self.sepyds[sname]['hist1d'][plot]['data']['val']), 
                self.plots[plot]['sep_scale'] * np.array(self.sepyds[sname]['hist1d'][plot]['data']['val'])
            ]).ravel(order="F")
            if self.cs['default']['sepLine']['edge']:
                import matplotlib.patheffects as PathEffects
                self.ax.plot(xx, yy, **self.cs['sep'][scnm], path_effects=[PathEffects.withStroke( **self.cs['default']['sepLine']['edgecss'])])       
            else:
                self.ax.plot(xx, yy, **self.cs['sep'][scnm])       
        
        self.draw_xaxis(plot)
        self.draw_yaxis(plot)
        
        self.ax.set_ylabel(self.plots[plot]['YLabel'], **self.cs['default']['YLabel'])
        self.ax.tick_params(**self.cs['default']['axtick']['major'])
        self.ax.tick_params(**self.cs['default']['axtick']['minor'])
        self.ax.tick_params(**self.cs['default']['axtick']['both'])
        
        self.draw_text(plot)
        self.draw_legend(plot)
        self.print_figure(plt)
        
    def draw_yaxis(self, plot):
        from matplotlib.ticker import AutoMinorLocator, MaxNLocator

        if "ylim" in self.plots[plot]:
            self.plots[plot]['YMin'] = eval(self.plots[plot]['ylim'])[0]
            self.plots[plot]['YMax'] = eval(self.plots[plot]['ylim'])[1]
        else:
            if type(self.plots[plot]['YMax']) is not str:
                if self.plots[plot]['YMax'] < 0:
                    self.plots[plot]['YMax'] = self.ax.get_ylim()[1]
            else:
                self.plots[plot]['YMax'] = eval(self.plots[plot]['YMax'])
            if type(self.plots[plot]['YMin']) is str:
                self.plots[plot]['YMin'] = eval(self.plots[plot]['YMin'])
            if self.plots[plot]['YMin'] < 0:
                self.plots[plot]['YMin'] = max(0., self.ax.get_ylim()[0])
                            
        if "yscale" in self.plots[plot].keys():
            self.ax.set_yscale(self.plots[plot]["yscale"].lower())
        elif type(self.plots[plot]['LogY']) is str:
            if eval(self.plots[plot]['LogY']):
                self.ax.set_yscale("log")
        elif self.plots[plot]['LogY']:
            self.ax.set_yscale("log")
        
        self.ax.set_ylim(self.plots[plot]['YMin'], self.plots[plot]['YMax'])
        if self.ax.get_yscale() == "linear":
            self.ax.yaxis.set_minor_locator(AutoMinorLocator())
            self.ax.ticklabel_format(**self.cs['default']['ylabel_format'])
            self.ax.yaxis.offsetText.set_fontsize(self.cs['default']['offset_text_fontsize'])
            self.ax.yaxis.set_major_locator(MaxNLocator(6))
        self.ax.yaxis.set_label_coords(**self.cs['default']['ylabelcoords'])

    def draw_xaxis(self, plot):
        if "xlim" in self.plots[plot]:
            self.plots[plot]['XMin'] = eval(self.plots[plot]['xlim'])[0]
            self.plots[plot]['XMax'] = eval(self.plots[plot]['xlim'])[1]        
        self.ax.set_xlim(self.plots[plot]['XMin'], self.plots[plot]['XMax'])
        # if "XTick" in self.plots[plot]:
        #     self.ax.set_xticks()
        from matplotlib.ticker import AutoMinorLocator, MultipleLocator, MaxNLocator, AutoLocator
        if self.plots[plot]['XPi']:
            xbin = eval(self.plots[plot]['XPi'])
            self.ax.xaxis.set_major_locator(MultipleLocator(np.pi * xbin * 5))
            self.ax.xaxis.set_minor_locator(MultipleLocator(np.pi * xbin))
            x_tick = np.arange(self.plots[plot]['XMin'], self.plots[plot]['XMax'], np.pi*xbin*5)
            x_ratio = np.arange(self.plots[plot]['XMin']/np.pi, self.plots[plot]['XMax']/np.pi, xbin*5 )
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
            self.ax.set_xticks(x_tick)
            self.ax.set_xticklabels(x_label)
        elif "xtick" in self.plots[plot]:
            if "AUTO" == self.plots[plot]['xtick'][0:4]:
                xbin = eval(self.plots[plot]['xtick'][5:])
                self.ax.xaxis.set_major_locator(MultipleLocator(xbin))
                self.ax.xaxis.set_minor_locator(MultipleLocator(self.xbin))
        else: 
            self.ax.xaxis.set_major_locator(MaxNLocator(9))
            self.ax.xaxis.set_minor_locator(MultipleLocator(self.xbin))
            nbin = (list(self.ax.xaxis.get_ticklocs())[1] - list(self.ax.xaxis.get_ticklocs())[0]) // (10 * self.xbin)
            if nbin:
                self.ax.xaxis.set_minor_locator(MultipleLocator((nbin+1) * self.xbin))
                self.ax.xaxis.set_major_locator(MultipleLocator(10 * (nbin+1) * self.xbin))
        self.ax.set_xlabel(self.plots[plot]['XLabel'], **self.cs['default']['XLabel'])
        self.ax.xaxis.set_label_coords(**self.cs['default']['xlabelcoords'])

    def draw_text(self, plot):
        if "title" in self.plots[plot]:
            self.ax.text(self.cs['text']['title']['pos'][0], self.cs['text']['title']['pos'][1], s=self.plots[plot]['title'], transform=self.ax.transAxes, **self.cs['text']['title']['css'])
        if "subtitle" in self.plots[plot]:
            self.ax.text(self.cs['text']['subtitle']['pos'][0], self.cs['text']['subtitle']['pos'][1], s=self.plots[plot]['subtitle'], transform=self.ax.transAxes, **self.cs['text']['subtitle']['css'])
            
    def draw_legend(self, plot):
        from copy import deepcopy
        if self.plots[plot]['legend']:
            ii, jj = 0, 0
            lgs = deepcopy(self.cs['legend'])
            self.ax.plot(
                [lgs['pos'][0] - lgs['block'][0] + ii * lgs['columnwidth'],  lgs['pos'][0] + lgs['block'][0] + ii * lgs['columnwidth']],
                [lgs['pos'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),  lgs['pos'][1]  - jj * (2 * lgs['block'][1] + lgs['hsep'] )],
                **self.cs['default']['cumLine'],
                transform=self.ax.transAxes
            ) 
            if self.plots[plot]['cumulative legend']:
                # print(self.plots[plot]['cumulative legend'], lgs['pos'][0] + 2*lgs['block'][0] + self.cs['text']['legend']['pos'][0] + ii * lgs['columnwidth'],
                    # lgs['pos'][1] + 2*lgs['block'][1] + self.cs['text']['legend']['pos'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ))
                self.ax.text(
                    lgs['pos'][0] + lgs['block'][0] + self.cs['text']['legend']['pos'][0] + ii * lgs['columnwidth'],
                    lgs['pos'][1] + self.cs['text']['legend']['pos'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                    s=self.plots[plot]['cumulative legend'],
                    transform=self.ax.transAxes,
                    **self.cs['text']['legend']['css']
                )
            ccum = list(self.cumyds.keys())
            for nn in range(len(self.cumyds)):
                jj, ii = (nn+1) // self.cs['legend']['ncolumn'], (nn+1) % self.cs['legend']['ncolumn']
                self.ax.fill(
                    [
                        lgs['pos'][0] - lgs['block'][0] + ii * lgs['columnwidth'],
                        lgs['pos'][0] + lgs['block'][0] + ii * lgs['columnwidth'],
                        lgs['pos'][0] + lgs['block'][0] + ii * lgs['columnwidth'],
                        lgs['pos'][0] - lgs['block'][0] + ii * lgs['columnwidth'],
                        lgs['pos'][0] - lgs['block'][0] + ii * lgs['columnwidth']
                    ], [
                        lgs['pos'][1] - lgs['block'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                        lgs['pos'][1] - lgs['block'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                        lgs['pos'][1] + lgs['block'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                        lgs['pos'][1] + lgs['block'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                        lgs['pos'][1] - lgs['block'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] )
                    ],
                    transform=self.ax.transAxes,
                    **self.cs['cum'][self.cumyds[ccum[nn]]['ccnm']]
                )
                self.ax.text(
                    lgs['pos'][0] + lgs['block'][0] + self.cs['text']['legend']['pos'][0] + ii * lgs['columnwidth'],
                    lgs['pos'][1] + self.cs['text']['legend']['pos'][1] - jj * (2 * lgs['block'][1] + lgs['hsep'] ),
                    s=self.cumyds[ccum[nn]]['title'],
                    transform=self.ax.transAxes,
                    **self.cs['text']['legend']['css']
                )                
            ssep = list(self.sepyds.keys())
            for nn in range(len(self.sepyds)):
                py = lgs['pos'][1]  - (nn + jj + 1) * (2 * lgs['block'][1] + lgs['hsep'] )
                xx = [lgs['pos'][0] - lgs['block'][0] ,  lgs['pos'][0] + lgs['block'][0]]
                yy = [py, py]
                if self.cs['default']['sepLine']['edge']:
                    import matplotlib.patheffects as PathEffects
                    self.ax.plot(
                        xx, yy, 
                        **self.cs['sep'][ssep[nn]], 
                        path_effects=[PathEffects.withStroke( **self.cs['default']['sepLine']['edgecss'])],
                        transform=self.ax.transAxes
                    )       
                else:
                    self.ax.plot(
                        xx, yy, 
                        **self.cs['sep'][ssep[nn]],
                        transform=self.ax.transAxes
                    )       
                self.ax.text(
                    lgs['pos'][0] + lgs['block'][0] + self.cs['text']['legend']['pos'][0],
                    py + self.cs['text']['legend']['pos'][1], 
                    s=self.sepyds[ssep[nn]]['title'],
                    transform=self.ax.transAxes,
                    **self.cs['text']['legend']['css']
                )
  

    def drawpicture(self):
        self.load_yodas()
        self.load()
        from copy import deepcopy
        self.sinf = deepcopy(self.inf)

        for plot in self.plots:
            if self.cf.has_section(plot):
                self.inf = deepcopy(self.sinf)
                self.inf['name'] = self.sinf["name"] + plot
                print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
                    time.time()-self.time, self.inf['name']), language="alias"))
                self.draw_plot(plot)
                