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
                if len(node.split("\n")) > 1:
                    htxt    = fs['wt'] * (max([len(x) for x in node.split("\n")]) + 1.5)
                else:
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
                "y":    y0 - max(nnl, 0.8) * fs['hv'],
                'w':    fs['wv'],
                'h':    max(nnl, 0.8) * fs['hv'],
            }
            y0 -= max(nnl, 0.8) * fs['hv']
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
                    n_nextlayer = self.find_var_nextlayer_number(var, kk)
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
                
                # self.layers['height'] = 
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

    def draw_scan_variable(self):
        x0 = 0. 
        y0 = -(self.size['h'] - self.layers[-1]['height']) / 2.0 
        ts = self.cs['text']['scanvar']
        fs = self.cs['ruler']
        vs = self.cs['text']['varname']
        lsx, lsy = self.get_label_shape([x0, x0 + fs['wvc']], [y0-fs['h2'], y0-fs['h2']], 2*fs['wvcb'], fs['h2'])
        self.ax.fill(lsx, lsy, **self.cs['label']['fill'])
        self.ax.plot(lsx, lsy, '-', **self.cs['label']['plot'])
        self.ax.text(
            0.5 * fs['wvc'], 
            y0 - 0.5 * self.cs['label']['scan_offset'] - 0.5 * fs['h2'], 
            "SCAN PARAMETERS", 
            **self.cs['text']['scanvar']
        )

        for kk, var in self.layers[-1]['scan variables']['locs'].items():
            var.update({
                "x":    var['x'] + x0, 
                "y":    var['y'] + y0, 
            })
            p = self.ax.fill(
                [var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x']  + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb']],
                [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']]
            )
            var.update({"c": p[0].get_facecolor()[0:3]})
          
            self.ax.text(
                var['x'] + fs['wvc'] - self.cs['label']['var_offset'], 
                var['y'] + 0.5*var['h'], 
                kk, 
                **self.cs['text']['varname'],
                horizontalalignment='right', 
                verticalalignment='center'            
            )            
            self.ax.text(
                var['x'] + fs['wvc'] + self.cs['label']['var_offset'], 
                var['y'] + 0.5*var['h'], 
                var['meth'], 
                **self.cs['text']['varmeth'],
                horizontalalignment='left', 
                verticalalignment='center'            
            )
        nscvar = len(self.layers[-1]['scan variables']['locs'])
        ii = 0
        for kk, var in self.layers[-1]['scan variables']['locs'].items():
            if ii != nscvar - 1:
                self.draw_outvar_label(var)
            else:
                self.draw_endoutvar_label(var)
            ii += 1 

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
        
        self.draw_scan_variable()
        for ii in range(len(self.layers) - 1):
            self.draw_layer(self.layers[ii], ii)
            
    def drawpicture(self):
        print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
                    time.time()-self.time, self.inf['name']), language="alias"))
        self.load()
        self.load_colorsetting()
        self.load_jarvis_flow()
        self.plot_flow()
        self.print_figure(plt)


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
    
    def draw_outvar_label(self, loc):
        fs = self.cs['ruler']
        vs = self.cs['text']['varname']
        self.ax.plot(
            [loc['x'], loc['x'], loc['x'] + fs['wvc'] + 0.4*fs['wvcb']],
            [loc['y']+loc['h'], loc['y'], loc['y']],
            '-',
            **self.cs['label']['plot'],
            zorder=0
        )
        
    def draw_endoutvar_label(self, loc):
        fs = self.cs['ruler']
        vs = self.cs['text']['varname']
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
        self.ax.plot(xx, yy, "-", **self.cs['label']['plot'], zorder=0)

    def draw_layer(self, vv, kk):
        fs = self.cs['ruler']
        # vs = fig['colorset']['varlabel']

        y0 = - 0.5 * (self.size['h'] - vv['height'])
        # print(kk, vv['var_nextlayer'].keys())
        for nn, va in vv['var_nextlayer'].items():
            va.update({"x": va['x'] + fs['wv'] + fs['wns'], "y": va['y'] + y0})
            sinf = self.get_jarvis_var_source(nn, kk, va['h'])
            va.update({"c": sinf['c']})
            self.draw_one_flow(sinf, va)
            self.ax.fill(
                [va['x'], va['x'] + va['w'], va['x'] + va['w'], va['x'], va['x']],
                [va['y'], va['y'], va['y'] + va['h'], va['y'] + va['h'], va['y']],
                alpha=0.6,
                facecolor=va['c'],
                edgecolor=None
            ) 
            txt = self.ax.text(
                va['x'] + 0.5*va['w'],
                va['y'] + 0.5*va['h'],
                nn, 
                **self.cs['text']['node'],
                rotation=90,
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
            # self.ax.plot(
            #     [va['base']['x0'], va['base']['x0'] + width, va['base']['x0'] + width, va['base']['x0'], va['base']['x0']],
            #     [va['base']['y0'], va['base']['y0'], va['base']['y0'] - va['height'], va['base']['y0'] - va['height'], va['base']['y0']]
            # )

            for nn, var in va['locs'].items():
                # if var['type'] == "inpvar" and var['meth'] != "file":
                if var['type'] == "inpvar":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    if var['meth'] != "file":
                        sinf = self.get_jarvis_var_source(nn, kk, var['h'])
                        var.update({"c": sinf['c']})
                    if var['posi'] == "in":
                        self.draw_inpvar_label(var)
                    elif var['posi'] == "side":
                        self.draw_endinpvar_label(var)

                    if "c" in var:
                        p = self.ax.fill(
                            [var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb']],
                            [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']],
                            facecolor=var['c']
                        )
                        var['x'] += var['w'] - fs['wvc']
                        self.draw_one_flow(sinf, var)
                    else:
                        p = self.ax.fill(
                            [var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] + 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb'], 
                             var['x'] + var['w'] -fs['wvc'] - 0.5*fs['wvcb']],
                            [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']]
                        )
                    self.ax.text(
                        var['x'] + self.cs['label']['var_offset'], 
                        var['y'] + 0.5*var['h'], 
                        nn, 
                        **self.cs['text']['varname'],
                        horizontalalignment='left', 
                        verticalalignment='center'            
                    )            
                    self.ax.text(
                        var['x'] - self.cs['label']['var_offset'], 
                        var['y'] + 0.5*var['h'], 
                        var['meth'], 
                        **self.cs['text']['varmeth'],
                        horizontalalignment='right', 
                        verticalalignment='center'            
                    )
                elif var['type'] == "inpfile":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    self.draw_inpfile(nn, var)
                elif var['type'] == "node":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    p = self.ax.fill(
                        [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                        [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                        alpha=0.5
                    )
                    var.update({"c": p[0].get_facecolor()[0:3]})
                    self.ax.plot(
                        [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                        [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                        '-',
                        linewidth=0.5,
                        c=var['c']
                    )

                    txt = self.ax.text(
                        var['x'] + 0.5*var['w'],
                        var['y'] + 0.5*var['h'],
                        nn, 
                        **self.cs['text']['node'],
                        rotation=90,

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
                    # self.ax.plot(
                    #     [var['x'], var['x'] + var['w'], var['x'] + var['w'], var['x'], var['x']],
                    #     [var['y'], var['y'], var['y'] + var['h'], var['y']+ var['h'], var['y']],
                    #     "-",
                    # )
                    self.draw_oupfile(nn, var)
                elif var['type'] == "oupvar":
                    var.update({
                        "x":    var['x'] + va['base']['x0'],
                        "y":    var['y'] + va['base']['y0']
                    })
                    p = self.ax.fill(
                        [var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x']  + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] + 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb'], var['x'] + fs['wvc'] - 0.5*fs['wvcb']],
                        [var['y'], var['y'], var['y']+var['h'], var['y'] + var['h'], var['y']]
                    )
                    var.update({"c": p[0].get_facecolor()[0:3]})
                    if var['posi'] == "in":
                        self.draw_outvar_label(var)
                    elif var['posi'] == "side":
                        self.draw_endoutvar_label(var)
                    
                    self.ax.text(
                        var['x'] + fs['wvc'] - self.cs['label']['var_offset'], 
                        var['y'] + 0.5*var['h'], 
                        nn, 
                        **self.cs['text']['varname'],
                        horizontalalignment='right', 
                        verticalalignment='center'            
                    )            
                    self.ax.text(
                        var['x'] + fs['wvc'] + self.cs['label']['var_offset'], 
                        var['y'] + 0.5*var['h'], 
                        var['meth'], 
                        **self.cs['text']['varmeth'],
                        horizontalalignment='left', 
                        verticalalignment='center'            
                    )

    def draw_oupfile(self, fname, var):
        fs = self.cs['ruler']
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
            p = self.ax.fill(xx, yy, alpha=0.7)
            var.update({"c": p[0].get_facecolor()[0:3]})
            
            rt1 = Polygon([pa1, pd1, pd2, pb1, pa1])
            xx, yy = rt1.exterior.xy 
            self.ax.plot(xx, yy, "-", linewidth=0.5, c=var['c'])
            
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
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)                    
            xx, yy = s3.exterior.xy 
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)
            xx, yy = s4.exterior.xy 
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
                                
        draw_logo()     

        ch, cs = len(fname), 16
        fname = [fname[i: i+cs] for i in range(0, ch, cs)]
        fname = "\n".join(fname)

        self.ax.text(
            var['x'] + 5*fs['wvcb'] + self.cs['label']['file_offset'],
            var['y'] + 0.5*var['h'],
            fname, 
            **self.cs['text']['file'],
            color=var['c'],
            va = 'center',
            ha = 'left',
            wrap = True
        )

    def draw_inpfile(self, fname, var):
        fs = self.cs['ruler']
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
            p = self.ax.fill(xx, yy, alpha=0.7)
            var.update({"c": p[0].get_facecolor()[0:3]})
            
            rt1 = Polygon([pa1, pd1, pd2, pb1, pa1])
            xx, yy = rt1.exterior.xy 
            self.ax.plot(xx, yy, "-", linewidth=0.5, c=var['c'])
            
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
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
            xx, yy = s2.exterior.xy 
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
            xx, yy = s3.exterior.xy 
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)
            xx, yy = s4.exterior.xy 
            self.ax.fill(xx, yy, facecolor=var['c'], edgecolor=None, alpha=0.7)            
                                
        draw_logo()
        ch, cs = len(fname), 16
        fname = [fname[i: i+cs] for i in range(0, ch, cs)]
        fname = "\n".join(fname)        
        
        self.ax.text(
            var['x'] + fs['wvcb'] + self.cs['label']['file_offset'],
            var['y'] + 0.5*var['h'],
            fname, 
            **self.cs['text']['file'],
            color=var['c'],
            va = 'center',
            ha = 'left'
        )
        
    def get_jarvis_var_source(self, vname, kk, h):
        fs = self.cs['ruler']
        if kk == 0:
            out = self.layers[kk-1]['scan variables']['locs']
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
            if vname in self.layers[kk-1]['var_nextlayer'].keys():
                outv = self.layers[kk-1]['var_nextlayer'][vname]
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
                for node in self.layers[kk-1]['nodes']:
                    if vname in self.layers[kk-1][node]['locs'].keys():
                        outv = self.layers[kk-1][node]['locs'][vname]
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
            
            
    def draw_inpvar_label(self, loc):
        fs = self.cs['ruler']
        self.ax.plot(
            [loc['x'] + loc['w'], loc['x'] + loc['w'], loc['x'] + loc['w'] - fs['wvc'] - 0.4*fs['wvcb']],
            [loc['y']+loc['h'], loc['y'], loc['y']],
            '-',
            **self.cs['label']['plot'],
            zorder=0
        )
 
    def draw_endinpvar_label(self, loc):
        fs = self.cs['ruler']
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
        self.ax.plot(xx, yy, "-", **self.cs['label']['plot'], zorder=0)

        
    def draw_one_flow(self, s, e):
        x = e['x'] - s['x']
        y = e['y'] - s['y']
        omega = np.pi / x 
        mag = - y/2. 
        tt = np.linspace(0., 1., 100)
        xx = tt * x + s['x']
        yy0 = mag * np.cos(tt*np.pi) + s['y'] + 0.5*y 
        yy1 = mag * np.cos(tt*np.pi) + s['y'] + 0.5*y + e['h']
        
        self.ax.fill_between(
            xx, yy0, yy1, 
            facecolor=e['c'],
            **self.cs['flow'],
        )