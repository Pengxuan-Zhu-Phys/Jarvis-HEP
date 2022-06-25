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
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

class JarvisFlow(Figure):
    def __init__(self) -> None:
        super().__init__()

    def load(self):
        if self.cf.has_option("PLOT_CONFI", "save format"):
            self.format = list(map(str.strip, self.cf.get("PLOT_CONFI", 'save format').split(",")))
            
    def load_colorsetting(self):
        with open(self.cs, 'r') as f1:
            self.cs = json.loads(f1.read())
    
    def load_layers(self):
        fs = self.cs['ruler']
        for ii in range(len(self.dat['plot']['layer'])):
            nn_currentlayer = len(self.dat['plot']['layer']) - ii - 1
            layer = self.dat['plot']['layer'][nn_currentlayer]
            self.layers[nn_currentlayer] = {
                "nodes":    layer,
                "var_nextlayer":    {}
            }
            for node in layer:
                hfinp   = len(self.dat['plot']['nodes'][node]['input file']) * fs['h3']
                hsinp   = (len(self.dat['plot']['nodes'][node]['input file']) - 1) * fs['h6']
                hvinp   = 0. 
                for ff in self.dat['plot']['nodes'][node]['input file']:                    
                    hvinp += len(ff['vars']) * fs['hv']

                hfoup   = len(self.dat['plot']['nodes'][node]['output file']) * fs['h3']
                hsoup   = (len(self.dat['plot']['nodes'][node]['input file']) - 1) * fs['h6']
                hvoup   = 0. 
                for ff in self.dat['plot']['nodes'][node]['output file']:
                    for var in ff['vars']:
                        n_nextlayer = self.find_var_nextlayer_number(var, nn_currentlayer)
                        # print(var, n_nextlayer)
                        nnl = n_nextlayer['nNL'] + n_nextlayer['nNNL']
                        if nnl == 0:
                            hvoup += fs['hv']
                        else:
                            hvoup += nnl * fs['hv']
        
                hinp    = hfinp + hsinp + hvinp
                houp    = hfoup + hsoup + hvoup
                htxt    = fs['wt'] * (len(node)+1.5)

                self.layers[nn_currentlayer][node] = {
                    "hinp":    hinp,
                    "houp":    houp
                }
                
                locs    = {}
                hmax = max(hinp, houp, htxt)
                y0 = -(hmax - hinp)/2.0 
                y1 = -(hmax - houp)/2.0    

                for ff in self.dat['plot']['nodes'][node]['input file']:                    
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

                for ff in self.dat['plot']['nodes'][node]['output file']:
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
                        n_nextlayer = self.find_var_nextlayer_number(vv, nn_currentlayer)
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
                            
                            
                self.layers[nn_currentlayer][node]['locs'] = locs 
                self.layers[nn_currentlayer][node]['height'] = hmax 

        nn_currentlayer = -1 
        self.layers[nn_currentlayer] = {
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
        for var in self.dat['scan variables']:
            n_nextlayer = self.find_var_nextlayer_number(var['name'], nn_currentlayer)
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
                self.layers[nn_currentlayer + 1]['var_nextlayer'][var['name']] = {"nhv":    n_nextlayer['nNNL']}

        self.layers[nn_currentlayer]['scan variables']['locs'] = locs 
        self.layers[nn_currentlayer]['scan variables']['height'] = -y0  
        self.layers[nn_currentlayer]['height'] = -y0 
        
        for nn_currentlayer in range(len(self.dat['plot']['layer'])):
            layer = self.dat['plot']['layer'][nn_currentlayer]
            for node in layer: 
                for ff in self.dat['plot']['nodes'][node]['output file']:
                    for var in ff['vars']:
                        n_nextlayer = self.find_var_nextlayer_number(var, nn_currentlayer)
                        if n_nextlayer['nNNL'] > 0:
                            self.layers[nn_currentlayer + 1]['var_nextlayer'][var] = {"nhv":  n_nextlayer['nNNL']}

        nn = -1
        for ii in range(len(self.layers)):
            kk = nn 
            nn += 1
            vv = self.layers[kk]
            if kk != -1:
                for var in vv['var_nextlayer']:
                    n_nextlayer = self.find_var_nextlayer_number(var, kk, dat)
                    if n_nextlayer['nNNL'] > 0:
                        self.layers[kk+1]['var_nextlayer'][var] = {"nhv": n_nextlayer['nNNL']}

        for kk, vv in self.layers.items():
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
                self.layers[kk]['height'] = -y0 
        self.size = {
            "h":    0.,
            "w":    (len(self.layers) - 1) * (fs['ws'] + 2*fs['wv'] + 2*fs['wns'] + fs['wn'] ) + fs['wv']
        }
        for kk, vv in self.layers.items():
            if self.layers[kk]['height'] > self.size['h']:
                self.size['h'] = self.layers[kk]['height']

    def find_var_nextlayer_number(self, var, nc):
        n_NL = 0 
        n_NNL = 0  
        for ii in range(len(self.dat['plot']['layer'])):
            if ii == nc + 1:
                for node in self.dat['plot']['layer'][ii]:
                    node = self.dat['plot']['nodes'][node]['input file']
                    for jj in range(len(node)):
                        for vname in node[jj]['vars'].keys():
                            if vname == var:
                                n_NL += 1
            elif ii > nc + 1:
                for node in self.dat['plot']['layer'][ii]:
                    node = self.dat['plot']['nodes'][node]['input file']
                    for jj in range(len(node)):
                        for vname in node[jj]['vars'].keys():
                            if vname == var:
                                n_NNL += 1
        return {"nNL":  n_NL, "nNNL":   n_NNL}

    def load_jarvis_flow(self):
        self.layers = {}
        with open(self.cf.get(self.inf['sect'], "info path"), 'r') as f1:
            self.dat = json.loads(f1.read())
        self.load_layers()        

    def plot_flow(self):
        fs = self.cs['ruler']

        self.fig = plt.figure(
            figsize=(
                self.size['w'] + 2*fs['w0'], 
                self.size['h'] + 2*fs['h0'])
        )
        self.fig.text(0., 0., "Test", color="None")
        self.fig.text(1., 1., "Test", color="None")
        self.ax = self.fig.add_axes([0., 0., 1., 1.])
        self.ax.axis("off")
        self.ax.set_xlim(-fs['w0'], self.size['w'] + fs['w0'])
        self.ax.set_ylim(-self.size['h'] - fs['h0'], fs['h0'])
        
        self.draw_scan_variable(fig)
        for ii in range(len(self.layers) - 1):
            self.draw_layer(fig, fig['layers'][ii], ii)

    def drawpicture(self):
        print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
                    time.time()-self.time, self.inf['name']), use_aliases=True))
        self.load()
        self.load_colorsetting()
        self.load_jarvis_flow()
        self.plot_flow()
