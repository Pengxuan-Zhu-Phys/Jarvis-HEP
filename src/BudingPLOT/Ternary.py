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
from plot import Figure
import json
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import emoji

config = {
    "font.family": ["serif", "Times New Roman"],
    "mathtext.fontset": 'stix',
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)


class Ternary(Figure):
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
            self.data = []
            for line in resf:
                if len(line.split(",")) != 2:
                    print("Illegal format of result file")
                    sys.exit(1)
                line = list(map(str.strip, line.split(",")))
                line[1] = os.path.join(self.path['path'], line[1])
                self.data.append(pd.read_csv(line[1]))
            if len(self.data) > 1:
                self.data = pd.concat(self.data)
            else:
                self.data = self.data[0]

    def load_variable(self):
        def load_var_data(xx):
            self.vars[xx]['data'] = self.data.eval(self.vars[xx]['expr'])

        def load_var_info(xx):
            if "{}_label".format(xx) in self.inf:
                self.vars[xx]['label'] = self.inf['{}_label'.format(xx)]
            else:
                self.vars[xx]['label'] = self.vars[xx]['expr']
            if "{}_lim".format(xx) in self.inf:
                tlim = [self.vars[xx]['data'].min(), self.vars[xx]
                        ['data'].max()]
                lim = self.inf['{}_lim'.format(xx)]
                if "AUTO" in lim:
                    count = float(lim[5:])
                    tlim[0] = (tlim[0] // count) * count
                    tlim[1] = (tlim[1] // count + 1) * count
                    self.vars[xx]['lim'] = tlim
                else:
                    self.vars[xx]['lim'] = list(
                        map(float, map(str.strip, lim.split(","))))
            if "{}_ticks".format(xx) in self.inf:
                tick = self.inf['{}_ticks'.format(xx)]
                if "AUTO" == tick.upper()[0:4]:
                    a = float(tick.split("_")[-1])
                    low = self.vars[xx]['lim'][0] // a
                    upp = self.vars[xx]['lim'][1] // a + 1
                    self.vars[xx]['ticks'] = np.linspace(
                        low*a, upp*a, int(upp-low+1))
            if "{}_scale".format(xx) in self.inf:
                if self.inf['{}_scale'.format(xx)].strip().lower() in ['log', 'linear']:
                    self.vars[xx]['scale'] = self.inf['{}_scale'.format(
                        xx)].strip().lower()
                else:
                    print(
                        "Illegal '{}_scale' in config file, using 'flat' scale as default".format(xx))
                    self.vars[xx]['scale'] = 'linear'
            if "{}_bin".format(xx) in self.inf:
                self.vars[xx]['bin'] = float(self.inf["{}_bin".format(xx)])

        if set(["left_variable", "right_variable", 'r_variable', 'g_variable', 'b_variable']).issubset(set(self.inf.keys())):
            self.mode = "RGB_Scatter"
            for xx in ['left', 'right', 'r', 'g', 'b']:
                if "{}_variable".format(xx) in self.inf:
                    self.vars[xx] = {
                        "expr": self.inf['{}_variable'.format(xx)]}
                    load_var_data(xx)
                    load_var_info(xx)

        elif set(["left_variable", "right_variable", 'bottom_variable', 'c_variable']).issubset(set(self.inf.keys())):
            self.mode = "C_Scatter"
            self.cmap = {}
            for xx in ['left', 'right', 'bottom', 'c']:
                if "{}_variable".format(xx) in self.inf:
                    self.vars[xx] = {
                        "expr": self.inf['{}_variable'.format(xx)]}
                    load_var_data(xx)
                    load_var_info(xx)
            self.vars["data"] = {
                "tot":  self.vars['left']['data'] + self.vars['right']['data'] + self.vars['bottom']['data']
            }
            self.vars['data']['left'] = self.vars['left']['data'] / \
                self.vars['data']['tot']
            self.vars['data']['right'] = self.vars['right']['data'] / \
                self.vars['data']['tot']
            self.vars['data']['bottom'] = self.vars['bottom']['data'] / \
                self.vars['data']['tot']
            self.vars['x'] = {"data": self.vars['data']
                              ['bottom'] + 0.5 * (self.vars['data']['right'])}
            self.vars['y'] = {"data": self.vars['data']['right']}
        elif set(["left_variable", "right_variable", 'bottom_variable']).issubset(set(self.inf.keys())):
            self.mode = "Scatter"
            for xx in ['left', 'right', 'bottom']:
                if "{}_variable".format(xx) in self.inf:
                    self.vars[xx] = {
                        "expr": self.inf['{}_variable'.format(xx)]}
                    load_var_data(xx)
                    load_var_info(xx)
            self.vars["data"] = {
                "tot":  self.vars['left']['data'] + self.vars['right']['data'] + self.vars['bottom']['data']
            }
            self.vars['data']['left'] = self.vars['left']['data'] / \
                self.vars['data']['tot']
            self.vars['data']['right'] = self.vars['right']['data'] / \
                self.vars['data']['tot']
            self.vars['data']['bottom'] = self.vars['bottom']['data'] / \
                self.vars['data']['tot']
            self.vars['x'] = {"data": self.vars['data']
                              ['bottom'] + 0.5 * (self.vars['data']['right'])}
            self.vars['y'] = {"data": self.vars['data']['right']}
        else:
            print("Illegal ternary input variables !!!")
            self.false = True
            return
        self.load_limits()

    def load_limits(self):
        if self.mode == "C_Scatter":
            self.vars['c']['limits'] = [
                self.vars['c']['data'].min(),
                self.vars['c']['data'].max()
            ]
        if self.mode == "Scatter" or self.mode == "C_Scatter":
            self.vars['left']['limits'] = [
                self.vars['data']['left'].min(),
                self.vars['data']['left'].max()
            ]
            self.vars['right']['limits'] = [
                self.vars['data']['right'].min(),
                self.vars['data']['right'].max()
            ]
            self.vars['bottom']['limits'] = [
                self.vars['data']['bottom'].min(),
                self.vars['data']['bottom'].max()
            ]
            self.vars['x']['limits'] = [
                self.vars['x']['data'].min(),
                self.vars['x']['data'].max()
            ]
            self.vars['y']['limits'] = [
                self.vars['y']['data'].min(),
                self.vars['y']['data'].max()
            ]
            self.ternary['ternary coord'] = {
                "bottom":   [
                    self.vars['bottom']['limits'][0] // 1.0,
                    self.vars['bottom']['limits'][1] // 1.0 + 1.0,
                ],
                "left":   [
                    self.vars['left']['limits'][0] // 1.0,
                    self.vars['left']['limits'][1] // 1.0 + 1.0,
                ],
                "right":   [
                    self.vars['right']['limits'][0] // 1.0,
                    self.vars['right']['limits'][1] // 1.0 + 1.0,
                ],
            }
            self.ternary['ax'] = {
                "lim": {
                    "x":    [
                        self.ternary['ternary coord']['bottom'][0] - 0.125 +
                        0.5 * self.ternary['ternary coord']['right'][0],
                        self.ternary['ternary coord']['bottom'][0] + 1.125 +
                        0.5 * self.ternary['ternary coord']['right'][0]
                    ],
                    "y":    [
                        self.ternary['ternary coord']['right'][0] -
                        0.5 * ((2.5/math.sqrt(3.0)) - 1.0),
                        self.ternary['ternary coord']['right'][1] +
                        0.5 * ((2.5/math.sqrt(3.0)) - 1.0)
                    ]
                },
                "x0": [
                    self.ternary['ternary coord']['bottom'][0] +
                    0.5 * self.ternary['ternary coord']['right'][0],
                    self.ternary['ternary coord']['right'][0]
                ],
                "tc": [
                    self.ternary['ternary coord']['bottom'][0] + 0.5 *
                    self.ternary['ternary coord']['right'][0] + 0.5,
                    self.ternary['ternary coord']['right'][0] + 1.0 / 3.0
                ]
            }
            if self.ternary['ternary coord'] != self.cs['default']['ternary coord']:
                print(emoji.emojize(
                    '\t:disguised_face::disguised_face::disguised_face: Ternary plot using the non-trival negative values, please carefully check the output plots !!', language="alias"))
                if self.vars['x']['limits'][0] > self.ternary['ternary coord']['bottom'][0]-0.5 and self.vars['x']['limits'][1] < self.ternary['ternary coord']['bottom'][1]-0.5 and self.vars['y']['limits'][0] > self.ternary['ternary coord']['right'][0] and self.vars['y']['limits'][1] < self.ternary['ternary coord']['right'][1]:
                    if self.ternary['ternary coord']['left'] != [
                        1.0 - self.ternary['ternary coord']['bottom'][0] -
                            self.ternary['ternary coord']['right'][1],
                        1.0 - self.ternary['ternary coord']['bottom'][0] -
                            self.ternary['ternary coord']['right'][0]
                    ]:
                        self.ternary['ternary coord']['left'] = [
                            1.0 - self.ternary['ternary coord']['bottom'][0] -
                            self.ternary['ternary coord']['right'][1],
                            1.0 - self.ternary['ternary coord']['bottom'][0] -
                            self.ternary['ternary coord']['right'][0]
                        ]
                    print(emoji.emojize(
                        '\t:ghost::ghost::ghost: Jarvis changes the default range \n\t\t{}\n\t\t to range\n\t\t{} !!!'.format(
                            self.cs['default']['ternary coord'], self.ternary['ternary coord']),
                        language="alias"))
                else:
                    print(emoji.emojize(
                        '\t:anguished_face::anguished_face::anguished_face: Jarvis can not deal with the ternary data, force stop plotting this figure !!!'.format(
                            self.cs['default']['ternary coord'], self.ternary['ternary coord']),
                        language="alias"))
                    self.false = True
                    return

    def load_colorsetting(self):
        with open(self.cs, "r") as f1:
            self.cs = json.loads(f1.read())
            self.cs['colormap_path'] = self.decode_path(
                self.cs['colormap_path'])
        if "colormap" in self.inf.keys():
            self.inf['colormap'] = list(
                map(str.strip, self.inf['colormap'].split(',')))

    def make_canvas(self, ax):
        if self.mode == "C_Scatter":
            self.ax.plot(
                self.ternary['ax']['x0'][0] + [0., 1., 0.5, 0.],
                self.ternary['ax']['x0'][1] + [0., 0., 1.0, 0.],
                '-',
                **self.cs['TSC']['canvas'][ax]
            )
            # self.ax.text()
            tlength = self.cs['TSC']['ticks']['length']
            for ii in range(self.cs['TSC']['ticks']['N'] + 1):
                self.ax.plot(
                    self.ternary['ax']['x0'][0] + ii /
                    (self.cs['TSC']['ticks']['N']) + [0.0, -0.5 * tlength],
                    self.ternary['ax']['x0'][1] + [0.0, -1.0 * tlength], "-",
                    **self.cs['TSC']['ticks']['tick']
                ),
                self.ax.text(
                    self.ternary['ax']['x0'][0] + ii/(
                        self.cs['TSC']['ticks']['N']) - 0.5 * self.cs['TSC']['ticks']['label_off'],
                    self.ternary['ax']['x0'][1] -
                    self.cs['TSC']['ticks']['label_off'],
                    self.cs['default']['coords_format'].format(
                        self.ternary['ternary coord']['bottom'][0] + ii/(self.cs['TSC']['ticks']['N'])),
                    **self.cs['TSC']['ticks']['bottom_tick_label']
                )
                self.ax.plot(
                    self.ternary['ax']['x0'][0] - 0.5 * ii /
                    (self.cs['TSC']['ticks']['N']) + [1.0, 1.0 + tlength],
                    self.ternary['ax']['x0'][1] + ii /
                    (self.cs['TSC']['ticks']['N']) + [0.0, 0.0], "-",
                    **self.cs['TSC']['ticks']['tick']
                )
                self.ax.text(
                    self.ternary['ax']['x0'][0] + 1.0 - 0.5 * ii /
                    (self.cs['TSC']['ticks']['N']) +
                    self.cs['TSC']['ticks']['label_off'],
                    self.ternary['ax']['x0'][1] + ii /
                    (self.cs['TSC']['ticks']['N']),
                    self.cs['default']['coords_format'].format(
                        self.ternary['ternary coord']['right'][0] + ii/(self.cs['TSC']['ticks']['N'])),
                    **self.cs['TSC']['ticks']['right_tick_label']
                )
                self.ax.plot(
                    self.ternary['ax']['x0'][0] - 0.5 * ii /
                    (self.cs['TSC']['ticks']['N']) +
                    [0.5, 0.5 - 0.5 * tlength],
                    self.ternary['ax']['x0'][1] - ii /
                    (self.cs['TSC']['ticks']['N']) + [1.0, 1.0 + tlength], "-",
                    **self.cs['TSC']['ticks']['tick']
                )
                self.ax.text(
                    self.ternary['ax']['x0'][0] + 0.5 - 0.5 * ii/(
                        self.cs['TSC']['ticks']['N']) - 0.5 * self.cs['TSC']['ticks']['label_off'],
                    self.ternary['ax']['x0'][1] + 1.0 - ii /
                    (self.cs['TSC']['ticks']['N']) +
                    self.cs['TSC']['ticks']['label_off'],
                    self.cs['default']['coords_format'].format(
                        self.ternary['ternary coord']['left'][0] + ii/(self.cs['TSC']['ticks']['N'])),
                    **self.cs['TSC']['ticks']['left_tick_label']
                )
            if self.cs['TSC']['grid']['switch']:
                for ii in range(self.cs['TSC']['ticks']['N'] - 1):
                    gg = ii + 1
                    gridlength = 1.0 - gg
                    self.ax.plot(
                        self.ternary['ax']['x0'][0] + gg/(self.cs['TSC']['ticks']['N']) + [
                            0.0, 0.5 * (1.0 - gg/(self.cs['TSC']['ticks']['N']))],
                        self.ternary['ax']['x0'][1] + [0.0,
                                                       (1.0 - gg/(self.cs['TSC']['ticks']['N']))], ":",
                        **self.cs['TSC']['grid']['style']
                    )
                    self.ax.plot(
                        self.ternary['ax']['x0'][0] + [0.5 * gg/(
                            self.cs['TSC']['ticks']['N']), 1.0 - 0.5 * gg/(self.cs['TSC']['ticks']['N'])],
                        self.ternary['ax']['x0'][1] + gg /
                        (self.cs['TSC']['ticks']['N']) + [0.0, 0.0], ":",
                        **self.cs['TSC']['grid']['style']
                    )
                    self.ax.plot(
                        self.ternary['ax']['x0'][0] + 0.5 + [- 0.5 * gg/(
                            self.cs['TSC']['ticks']['N']), 0.5 - gg/(self.cs['TSC']['ticks']['N'])],
                        self.ternary['ax']['x0'][1] + [1.0 - gg /
                                                       (self.cs['TSC']['ticks']['N']), 0.0], ":",
                        **self.cs['TSC']['grid']['style']
                    )
            self.ax.text(
                self.ternary['ax']['tc'][0],
                self.ternary['ax']['tc'][1] - 1.0 *
                self.cs['TSC']['label']['off_set'],
                r"{}".format(self.inf['bottom_label']),
                **self.cs['TSC']['label']['bottom']
            )
            self.ax.text(
                self.ternary['ax']['tc'][0] + 0.75 *
                self.cs['TSC']['label']['off_set'],
                self.ternary['ax']['tc'][1] + 0.5 *
                self.cs['TSC']['label']['off_set'],
                r"{}".format(self.inf['right_label']),
                **self.cs['TSC']['label']['right']
            )
            self.ax.text(
                self.ternary['ax']['tc'][0] - 0.75 *
                self.cs['TSC']['label']['off_set'],
                self.ternary['ax']['tc'][1] + 0.5 *
                self.cs['TSC']['label']['off_set'],
                r"{}".format(self.inf['left_label']),
                **self.cs['TSC']['label']['left']
            )

    def draw_plot(self):
        if self.mode == "Scatter":
            self.fig = plt.figure(**self.cs['TS']['figsize'])
            self.ax = self.fig.add_axes(**self.cs['TS']['axsize'])
            self.ax.scatter(self.vars['x']['data'],
                            self.vars['y']['data'], marker='^')
            self.ax.plot([-0.5, 0.5, 0.0, -0.5], [-1.0, -1.0, 0.0, -1.0])
        elif self.mode == "C_Scatter":
            self.fig = plt.figure(**self.cs['TSC']['figsize'])
            self.ax = self.fig.add_axes(**self.cs['TSC']['axsize'], zorder=1)
            self.axc = self.fig.add_axes(**self.cs['TSC']['axcsize'], zorder=2)
            self.ax.axis("off")
            self.ax.set_xlim(self.ternary['ax']['lim']['x'])
            self.ax.set_ylim(self.ternary['ax']['lim']['y'])
            self.make_canvas("ax")
            self.load_colormap()
            if self.vars['c']['scale'] == "linear":
                from matplotlib.colors import Normalize
                a1 = self.ax.scatter(self.vars['x']['data'],
                                     self.vars['y']['data'],
                                     c=self.vars['c']['data'],
                                     vmin=self.cmap['vmin'],
                                     vmax=self.cmap['vmax'],
                                     cmap=self.cmap['cmap'],
                                     **self.cs['TSC']['marker'])
                plt.colorbar(
                        a1, self.axc, ticks=self.vars['c']['ticks'], 
                        orientation='vertical', 
                        norm=Normalize(vmin=self.cmap['vmin'], vmax=self.cmap['vmax']),
                        extend='neither'
                    )
                self.axc.tick_params(**self.cs['TSC']['axc_ticks_parameter'])
                self.axc.set_ylabel(
                    self.inf['c_label'], **self.cs['TSC']['label']['c'])

            # self.ax.plot([-0.5, 0.5, 0.0, -0.5], [-1.0, -1.0, 0.0, -1.0])
            # self.ax.plot(self.ternary['ax']['x0'][0] + [0.5, 0.5],
            #              self.ternary['ax']['x0'][1] + [0., 1.])
            # self.ax.plot(self.ternary['ax']['x0'][0] + [0.0, 0.75],
            #              self.ternary['ax']['x0'][1] + [0.0, 0.5])
            # self.ax.plot(self.ternary['ax']['x0'][0] + [0.25, 1.0],
            #              self.ternary['ax']['x0'][1] + [0.5, 0.0])
            # self.ax.scatter(self.ternary['ax']['tc']
            #                 [0], self.ternary['ax']['tc'][1])

    def drawpicture(self):
        print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
            time.time()-self.time, self.inf['name']), language="alias"))
        self.load()
        self.load_colorsetting()
        self.load_variable()
        if not self.false:
            self.draw_plot()
            self.print_figure(plt)
