#!/usr/bin/env python3
import os, sys
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
    "font.family":["serif", "Times New Roman"],
    "mathtext.fontset":'stix',
    "font.serif": ['Computer Modern'],
    "text.latex.preamble": r"\usepackage{amsmath}"
}
rcParams.update(config)

class Voronoi2D(Figure):
    def __init__(self):
        super().__init__()
    
    def load(self):
        if self.cf.has_option("PLOT_CONFI", "save format"):
            self.format = list(map(str.strip, self.cf.get("PLOT_CONFI", 'save format').split(",")))
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
                
            
            
            
    def load_colorsetting(self):
        with open(self.cs, 'r') as f1:
            self.cs = json.loads(f1.read())
    
            
    def draw_plot(self):
        self.load_variable()
        if self.cbar:
            self.fig = plt.figure(**self.cs['2DC']['axis']['figsize'])
            self.ax  = self.fig.add_axes(**self.cs['2DC']['axis']["axsize"])
            self.axc = self.fig.add_axes(**self.cs['2DC']['axis']['axcsize'])
            from scipy.spatial import Voronoi, voronoi_plot_2d
            ss = pd.DataFrame({"x": self.vars['x']['data'], "y": self.vars['y']['data']}).to_numpy()
            vor = Voronoi(ss[:, :])
            self.ax.scatter(self.vars['x']['data'], self.vars['y']['data'], s=3, c='black', zorder=100)
            voronoi_plot_2d(vor, self.ax, **self.cs['default']['voronoiPlot'])

            self.draw_xaxis()
            self.draw_yaxis()

            norm = matplotlib.colors.Normalize(vmin=self.vars['c']['lim'][0], vmax=self.vars['c']['lim'][1], clip=True)
            mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap=matplotlib.cm.cool_r)
            regions, vertices = self.voronoi_finite_polygons_2d(vor)
            for ii in range(len(regions)):
                region = regions[ii]
                polygon = vertices[region]
                self.ax.fill(*zip(*polygon), color=mapper.to_rgba(self.vars['c']['data'][ii]))
            plt.colorbar(mapper, self.axc, ticks=self.vars['c']['ticks'], orientation='vertical', extend="neither")
            self.ax.tick_params(**self.cs['default']['axtick']['major'])
            self.ax.tick_params(**self.cs['default']['axtick']['minor'])
            self.ax.tick_params(**self.cs['default']['axtick']['both'])
            self.draw_caxis()
   
                
    def draw_xaxis(self):
        from matplotlib.ticker import AutoMinorLocator, MaxNLocator
        self.ax.set_xscale(self.vars['x']['scale'])
        if self.vars['x']['scale'] == "linear":
            self.ax.set_xticks(self.vars['x']['ticks'])
            self.ax.ticklabel_format(**self.cs['default']['label_format'])
            self.ax.xaxis.set_minor_locator(AutoMinorLocator())
        self.ax.set_xlim(self.vars['x']['lim'])
        self.ax.set_xlabel(r"{}".format(self.vars['x']['label']), **self.cs['default']['XLabel'])
            
    def draw_yaxis(self):
        from matplotlib.ticker import AutoMinorLocator, MaxNLocator
        self.ax.set_yscale(self.vars['y']['scale'])
        if self.vars['y']['scale'] == "linear":
            self.ax.set_yticks(self.vars['y']['ticks'])
            self.ax.ticklabel_format(**self.cs['default']['label_format'])
            self.ax.yaxis.set_minor_locator(AutoMinorLocator())
        self.ax.set_ylim(self.vars['y']['lim'])
        self.ax.set_ylabel(r"{}".format(self.vars['y']['label']), **self.cs['default']['YLabel'])        

    def draw_caxis(self):
        from matplotlib.ticker import AutoMinorLocator, MaxNLocator
        self.axc.set_yscale(self.vars['c']['scale'])
        if self.vars['c']['scale'] == "linear":
            print(self.vars['c'].items())
            self.axc.set_yticks(self.vars['c']['ticks'])
            self.axc.yaxis.set_minor_locator(AutoMinorLocator())
        self.axc.set_ylim(self.vars['c']['lim'])
        self.axc.ticklabel_format(**self.cs['default']['label_format'])
        self.axc.set_ylabel(r"{}".format(self.vars['c']['label']), **self.cs['default']['YLabel'])  
        self.axc.tick_params(**self.cs['default']['axctick']['major'])
        self.axc.tick_params(**self.cs['default']['axctick']['minor'])
        self.axc.tick_params(**self.cs['default']['axctick']['both']) 
        
    def drawpicture(self):
        print(emoji.emojize("\n\t:clock2: {:.2f} Sec;  :art::art::art: plotting {} ....".format(
                time.time()-self.time, self.inf['name']), use_aliases=True))
        self.load()
        self.load_colorsetting()
        self.draw_plot()
        self.print_figure(plt)
        
        
        
        
    def voronoi_finite_polygons_2d(self, vor, radius=None):
        """
        Reconstruct infinite voronoi regions in a 2D diagram to finite
        regions.
        Parameters
        ----------
        vor : Voronoi
            Input diagram
        radius : float, optional
            Distance to 'points at infinity'.
        Returns
        -------
        regions : list of tuples
            Indices of vertices in each revised Voronoi regions.
        vertices : list of tuples
            Coordinates for revised Voronoi vertices. Same as coordinates
            of input vertices, with 'points at infinity' appended to the
            end.
        """

        if vor.points.shape[1] != 2:
            raise ValueError("Requires 2D input")

        new_regions = []
        new_vertices = vor.vertices.tolist()

        center = vor.points.mean(axis=0)
        if radius is None:
            radius = vor.points.ptp().max()*2

        # Construct a map containing all ridges for a given point
        all_ridges = {}
        for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
            all_ridges.setdefault(p1, []).append((p2, v1, v2))
            all_ridges.setdefault(p2, []).append((p1, v1, v2))

        # Reconstruct infinite regions
        for p1, region in enumerate(vor.point_region):
            vertices = vor.regions[region]

            if all(v >= 0 for v in vertices):
                # finite region
                new_regions.append(vertices)
                continue

            # reconstruct a non-finite region
            ridges = all_ridges[p1]
            new_region = [v for v in vertices if v >= 0]

            for p2, v1, v2 in ridges:
                if v2 < 0:
                    v1, v2 = v2, v1
                if v1 >= 0:
                    # finite ridge: already in the region
                    continue

                # Compute the missing endpoint of an infinite ridge

                t = vor.points[p2] - vor.points[p1] # tangent
                t /= np.linalg.norm(t)
                n = np.array([-t[1], t[0]])  # normal

                midpoint = vor.points[[p1, p2]].mean(axis=0)
                direction = np.sign(np.dot(midpoint - center, n)) * n
                far_point = vor.vertices[v2] + direction * radius

                new_region.append(len(new_vertices))
                new_vertices.append(far_point.tolist())

            # sort region counterclockwise
            vs = np.asarray([new_vertices[v] for v in new_region])
            c = vs.mean(axis=0)
            angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
            new_region = np.array(new_region)[np.argsort(angles)]

            # finish
            new_regions.append(new_region.tolist())

        return new_regions, np.asarray(new_vertices)