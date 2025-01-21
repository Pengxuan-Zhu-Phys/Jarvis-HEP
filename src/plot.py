#!/usr/bin/env python3

import os, io 
import sys
import pandas as pd
import time
import math
import sympy
import json
import emoji
import numpy as np 
import matplotlib.pyplot as plt
import shapely as sp 
import contextlib
from base import Base
import yaml 

pwd = os.path.abspath(os.path.dirname(__file__))
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.abspath(os.path.join(pwd, "BudingPLOT")))


def draw_logo_in_square(ax):
    with io.StringIO() as buf, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):

        ax.axis("off")
        # ax.text(0.5, -0.06, "Jarvis-HEP", fontsize=12, ha="center", va='top', transform=ax.transAxes)
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-1.02, 1.02)

        # print(xx)
        rs = create_round_square((0., 0.))
        x_coords = [point[0] for point in rs.exterior.coords]
        y_coords = [point[1] for point in rs.exterior.coords]
        ax.plot(x_coords, y_coords, '-', color='grey', lw=0.3)
        ax.fill(x_coords, y_coords, fc='#34495E', ec=None)

        grid1 = [
        [False, False, False,   True, False, False, False, False],
        [False, False,  True,   True, False, False, False, False],
        [False,  True,  True,   True, False, False, False, False],
        [False, False,  True,   True, False, False, False, False],
        [False, False, False,   True, False, False, False, False],
        [True,   True,  True,   True, False, False, False, False],
        [False,  True,  True,   True, False, False, False, False],
        [False, False, False,   True, False, False, False, False]
    ]
        round1 = [
        [False, False, False, False, False, False, False, False],
        [False, False, False, False, True,  False, False, False],
        [False, False, False, False, True,  False, False, False],
        [False, False, False, False, True,  False, False, False],
        [False, False, False, False, True,  True,  False, False],
        [False, False, False, False, True,  True,  True,  True ],
        [False, False, False, False, True,  True,  False, False],
        [False, False, False, False, True,  False, False, False]
    ]

        for ii in range(8):
            for jj in range(8):
                if grid1[jj][ii]:
                    grid = create_round_square(
                        (
                            0.125 + ii * 0.25 - 1.0, 
                            1.0 - 0.125 - jj * 0.25,
                        ), 
                        rd=0.11, nn=6, frame=rs
                    )
                    x_coords = [point[0] for point in grid.exterior.coords]
                    y_coords = [point[1] for point in grid.exterior.coords]    
                    ax.fill(x_coords, y_coords, fc="#00fdff", ec=None)
                elif round1[jj][ii]:
                    grid = create_round_square(
                        (
                            0.125 + ii * 0.25 - 1.0, 
                            1.0 - 0.125 - jj * 0.25,
                        ), 
                        rd=0.11, nn=2.1, frame=rs
                    )
                    x_coords = [point[0] for point in grid.exterior.coords]
                    y_coords = [point[1] for point in grid.exterior.coords]    
                    ax.fill(x_coords, y_coords, fc="#fffc79", ec=None)    
                elif ii < 4:
                    grid = create_round_square(
                        (
                            0.125 + ii * 0.25 - 1.0, 
                            1.0 - 0.125 - jj * 0.25,
                        ), 
                        rd=0.11, nn=6, frame=rs
                    )
                    x_coords = [point[0] for point in grid.exterior.coords]
                    y_coords = [point[1] for point in grid.exterior.coords]    
                    ax.fill(x_coords, y_coords, fc="#072346", ec=None)   
                else:
                    grid = create_round_square(
                        (
                            0.125 + ii * 0.25 - 1.0, 
                            1.0 - 0.125 - jj * 0.25,
                        ), 
                        rd=0.11, nn=4, frame=rs
                    )
                    x_coords = [point[0] for point in grid.exterior.coords]
                    y_coords = [point[1] for point in grid.exterior.coords]    
                    ax.fill(x_coords, y_coords, fc="#008f00", ec=None)  

def create_round_square(ct, rd=1, nn = 3.5, frame=None):
    with io.StringIO() as buf, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        tt = np.linspace(0., 2.0 * np.pi, 314)
        xx = (np.cos(tt)) ** (2/nn)
        yy = (np.sin(tt)) ** (2/nn)
        xx[np.isnan(xx)] = 0
        yy[np.isnan(yy)] = 0
        xi = -(-np.cos(tt)) ** (2/nn)
        yi = -(-np.sin(tt)) ** (2/nn)
        xi[np.isnan(xi)] = 0
        yi[np.isnan(yi)] = 0
        xx += xi 
        yy += yi
        xx = rd * xx + ct[0] 
        yy = rd * yy + ct[1]
        cdn = np.stack((xx, yy), axis=-1)
        rs = sp.Polygon(cdn)
        if frame is not None:
            rs = rs.intersection(frame)
        return rs
   
class CustomDumper(yaml.Dumper):
    def write_line_break(self, data=None):
        super().write_line_break(data)
        if len(self.indents) == 1:
            super().write_line_break()

class BudingPLOT(Base):
    def __init__(self):
        super().__init__()
        self.info = {}
        self.sdf = None
        self.ddf = None

    def get_plot_config_from_Jarvis(self, info, scan_yaml):
        self.info.update(info)
        if scan_yaml.config['Sampling']['Method'] == "Dynesty":
            self.info['db']["nested result"] = os.path.join(self.info['sample']['task_result_dir'], "DATABASE", "dynesty_result.csv")
            self.yaml = {
                "Plot_Config":  {
                    "save_dir": self.info['plot']['save_path'],
                    "save_format":  ["pdf", "png"],
                    "dynesty":  {
                        "result":   self.info['db']['nested result']
                    },
                    "samples":  self.info['db']['out_csv'],
                    "scan_yaml":    self.info['config_file'],
                    "screen_show":  True
                },
                "Figures":  [
                    {
                        "name": "dynesty_sampling_summary",
                        "type": "dynesty_run"
                    },
                    {
                        "name": "parameter_summary",
                        "type": "dynesty_parameter",
                        "parameters":   []
                    }    
                ],
                "Variables":    []
            }
            for var in scan_yaml.config["Sampling"]["Variables"]:
                print(var)
                vardict = {
                    "name":     var['name'],
                    "label":    var['name'],
                    "scale":    var['distribution']['type'],
                    "lim":      [var['distribution']['parameters']['min'], var['distribution']['parameters']['max']]
                }
                self.yaml['Variables'].append(vardict)
                self.yaml['Figures'][1]['parameters'].append(str(var['name']))
            print(self.yaml)
            with open(self.info['plot']['config'], 'w') as file:
                yaml.dump(self.yaml, file, Dumper=CustomDumper, default_flow_style=False, allow_unicode=True)
            
    def load_config(self, filepath) -> None:
        with open(filepath, 'r') as file:
            self.yaml = yaml.safe_load(file)
    
    def plot(self) -> None: 
        if "dynesty" in self.yaml['Plot_Config']:
            self.ddf = pd.read_csv(self.yaml['Plot_Config']['dynesty']['result'])
        self.sdf = pd.read_csv(self.yaml['Plot_Config']['samples'] )
        for image in self.yaml['Figures']:
            if image['type'] == "dynesty_run":
                self.plot_dynesty_results(image['name'])
            if image['type'] == "dynesty_parameter":
                self.plot_dynesty_parameter(image)
            
                
    def plot_dynesty_parameter(self, image) -> None:
        from matplotlib.ticker import AutoMinorLocator
        from scipy.stats import gaussian_kde
        
        ndim = len(image['parameters'])    
        width = 6 + 3 * ndim
        height = 2 + 3 * ndim
        fig = plt.figure(figsize=(width, height))
        logo = fig.add_axes([1 - 1.2 / width, 1 - 1.2 / height, 1 / width, 1 / height])
        draw_logo_in_square(logo) 

        for ii in range(ndim):
            axLX = fig.add_axes([3 / width, (1 + 3*ii)/ height, 5.75 / width, 3 / height])
            a1 = axLX.scatter( - self.ddf['log_PriorVolume'], self.ddf[f"samples_u[{ii}]"], s=1, marker='.', c=self.ddf['log_Like'], cmap="plasma_r", alpha=0.7)
            axLX.set_xlim(0, max( - self.ddf['log_PriorVolume']))
            axLX.set_ylim(0, 1)
            axLX.yaxis.set_minor_locator(AutoMinorLocator())
            axLX.xaxis.set_minor_locator(AutoMinorLocator())
            axLX.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
            axLX.tick_params(which='major', length=10)
            axLX.tick_params(which='minor', length=4)
            if ii != 0:
                axLX.set_xticklabels([])
            else:
                axLX.set_xlabel(r"$-\ln(X)$", fontsize=24)
            axLX.set_yticklabels([])

            for jj in range(ndim):
                kk = ndim - jj - 1 
                if ii < jj:
                    ax2d = fig.add_axes([(8.8 + kk * 3) / width, ( 1 + ii * 3) / height, 3/width, 3/height])
                    ax2d.scatter(self.ddf[f"samples_u[{jj}]"], self.ddf[f"samples_u[{ii}]"], s=1, marker='.', c=self.ddf['log_Like'], cmap="plasma_r", alpha=0.7)
                    ax2d.set_xlim(0, 1)
                    ax2d.set_ylim(0, 1)
                    ax2d.yaxis.set_minor_locator(AutoMinorLocator())
                    ax2d.xaxis.set_minor_locator(AutoMinorLocator())
                    ax2d.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
                    ax2d.tick_params(which='major', length=10)
                    ax2d.tick_params(which='minor', length=4)
                    ax2d.set_yticklabels([])
                    if ii != 0:
                        ax2d.set_xticklabels([])
                    else:
                        ax2d.set_xlabel(f"u{jj}", fontsize=24)
                elif ii == jj:
                    ax1d = fig.add_axes([1 / width, ( 1 + ii * 3) / height, 1.95/width, 3/height])
                    wt_kde = gaussian_kde(np.array(self.ddf[f"samples_u[{ii}]"]), weights=np.exp(np.array(self.ddf['log_weight'])), bw_method="silverman")
                    yy = np.linspace(0, 5, 500)
                    wtt = wt_kde(yy)
                    wtt = wtt / max(wtt)
                    ax1d.plot(wtt, yy, '-', linewidth=1.5, color="#283593")
                    ax1d.fill_betweenx(yy, wtt, 0., edgecolor=None, facecolor="#03A9F4", alpha=0.5)
                    ax1d.yaxis.set_minor_locator(AutoMinorLocator())
                    ax1d.xaxis.set_minor_locator(AutoMinorLocator())
                    ax1d.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
                    ax1d.tick_params(which='major', length=10)
                    ax1d.tick_params(which='minor', length=4)
                    ax1d.set_ylim(0, 1)
                    ax1d.set_xlim(1.2, 0)
                    if ii != 0:
                        ax1d.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
                        ax1d.set_xticklabels([])
                    else:  
                        ax1d.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
                    ax1d.grid(visible=True, which='major', axis="x")
                    # ax1d.set_yticklabels([])
                    ax1d.set_ylabel(f"u{ii}", fontsize=24)

        ax1d.xaxis.set_label_position("top")
        ax1d.set_xlabel(r"$\text{Post-PDF}$", fontsize=20)

        axc = fig.add_axes([3 / width, (1.05 + 3*ndim )/ height, 5.75 / width, 0.2 / height])
        plt.colorbar(a1, axc, orientation="horizontal")
        axc.tick_params(labelsize=11,  direction="in", bottom=False, left=False, top=True, right=False, which='both')
        axc.tick_params(which='major', length=7)
        axc.xaxis.set_minor_locator(AutoMinorLocator())
        axc.tick_params(which='minor', length=4)
        axc.xaxis.set_label_position('top')
        axc.xaxis.tick_top()
        axc.set_xlabel(r"$\log{\mathcal{L}}$", fontsize=24)
        plt.draw()
        
        image['fig'] = plt
        image['file'] = os.path.join(self.yaml['Plot_Config']['save_dir'], image['name'])
        self.savefig(image, plt)
        # plt.close()
        if self.yaml['Plot_Config']['screen_show']:
            plt.show(block=False)
        #     plt.pause(1)
            input("\n\tPress 'Enter' to continue ...")
        plt.close()

        # Plot same plot for all samples 
        
        self.sdf.replace([np.inf, -np.inf], np.nan, inplace=True)
        self.sdf.dropna(subset=['LogL'], inplace=True)
        
        # Calculate Log(X) for samples
        from scipy.interpolate import interp1d
        lxll = interp1d(self.ddf['log_Like'], self.ddf['log_PriorVolume'], kind="linear", fill_value="extrapolate")

        lxx = lxll(self.sdf['LogL'])
        self.sdf['log_PriorVolume'] = lxx 
        self.sdf = self.sdf.sort_values(by="log_PriorVolume", ascending=False)
        self.ddf = self.ddf.sort_values(by="log_PriorVolume", ascending=False)
        
        # Calculate weight for samples 
        wt = np.exp(np.array(self.ddf["log_weight"]))
        self.ddf['CDF'] = np.cumsum(wt)

        cdf = np.insert(np.array(self.ddf['CDF']), 0, 0.)
        lnx = np.array(- self.ddf['log_PriorVolume'])
        lnx = np.insert(lnx, 0, 0.)
        cdf_func = interp1d(lnx, cdf, kind="linear")
        scdf = cdf_func(- self.sdf['log_PriorVolume'])
        aaf = np.insert(scdf, 0, 0.)
        self.sdf['weight'] = np.diff(aaf)

        self.sdf.to_csv(self.yaml['Plot_Config']['samples'], index=False)


        fig = plt.figure(figsize=(width, height))
        logo = fig.add_axes([1 - 1.2 / width, 1 - 1.2 / height, 1 / width, 1 / height])
        draw_logo_in_square(logo) 
        vardict = []
        for var in image['parameters']:
            for inp in self.yaml['Variables']:
                if inp['name'] == var:
                    vardict.append(inp)
                
        for ii in range(ndim):
            varii = vardict[ii]
            axLX = fig.add_axes([3 / width, (1 + 3*ii)/ height, 5.75 / width, 3 / height])
            a1 = axLX.scatter( - self.sdf['log_PriorVolume'], self.sdf[varii["name"]], s=1, marker='.', c=self.sdf['LogL'], cmap="plasma_r", alpha=0.7)
            axLX.set_xlim(0, max( - self.sdf['log_PriorVolume']))
            axLX.set_ylim(varii['lim'])
            axLX.yaxis.set_minor_locator(AutoMinorLocator())
            axLX.xaxis.set_minor_locator(AutoMinorLocator())
            axLX.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
            axLX.tick_params(which='major', length=10)
            axLX.tick_params(which='minor', length=4)
            if ii != 0:
                axLX.set_xticklabels([])
            else:
                axLX.set_xlabel(r"$-\ln(X)$", fontsize=24)
            axLX.set_yticklabels([])

            for jj in range(ndim):
                kk = ndim - jj - 1 
                varjj = vardict[jj]
                if ii < jj:
                    ax2d = fig.add_axes([(8.8 + kk * 3) / width, ( 1 + ii * 3) / height, 3/width, 3/height])
                    ax2d.scatter(self.sdf[varjj['name']], self.sdf[varii['name']], s=1, marker='.', c=self.sdf['LogL'], cmap="plasma_r", alpha=0.7)
                    ax2d.set_xlim(varjj['lim'])
                    ax2d.set_ylim(varii['lim'])
                    ax2d.yaxis.set_minor_locator(AutoMinorLocator())
                    ax2d.xaxis.set_minor_locator(AutoMinorLocator())
                    ax2d.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
                    ax2d.tick_params(which='major', length=10)
                    ax2d.tick_params(which='minor', length=4)
                    ax2d.set_yticklabels([])
                    if ii != 0:
                        ax2d.set_xticklabels([])
                    else:
                        ax2d.set_xlabel(varjj['label'], fontsize=24)
                elif ii == jj:
                    ax1d = fig.add_axes([1 / width, ( 1 + ii * 3) / height, 1.95/width, 3/height])
                    wt_kde = gaussian_kde(np.array(self.sdf[varii['name']]), weights=np.array(self.sdf['weight']), bw_method="silverman")
                    yy = np.linspace(varii['lim'][0], varii['lim'][1], 500)
                    wtt = wt_kde(yy)
                    wtt = wtt / max(wtt)
                    ax1d.plot(wtt, yy, '-', linewidth=1.5, color="#283593")
                    ax1d.fill_betweenx(yy, wtt, 0., edgecolor=None, facecolor="#03A9F4", alpha=0.5)
                    ax1d.yaxis.set_minor_locator(AutoMinorLocator())
                    ax1d.xaxis.set_minor_locator(AutoMinorLocator())
                    ax1d.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
                    ax1d.tick_params(which='major', length=10)
                    ax1d.tick_params(which='minor', length=4)
                    ax1d.set_ylim(varii['lim'])
                    ax1d.set_xlim(1.2, 0)
                    if ii != 0:
                        ax1d.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
                        ax1d.set_xticklabels([])
                    else:  
                        ax1d.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
                    ax1d.grid(visible=True, which='major', axis="x")
                    # ax1d.set_yticklabels([])
                    ax1d.set_ylabel(varii['label'], fontsize=24)

        ax1d.xaxis.set_label_position("top")
        ax1d.set_xlabel(r"$\text{Post-PDF}$", fontsize=20)

        axc = fig.add_axes([3 / width, (1.05 + 3*ndim )/ height, 5.75 / width, 0.2 / height])
        plt.colorbar(a1, axc, orientation="horizontal")
        axc.tick_params(labelsize=11,  direction="in", bottom=False, left=False, top=True, right=False, which='both')
        axc.tick_params(which='major', length=7)
        axc.xaxis.set_minor_locator(AutoMinorLocator())
        axc.tick_params(which='minor', length=4)
        axc.xaxis.set_label_position('top')
        axc.xaxis.tick_top()
        axc.set_xlabel(r"$\log{\mathcal{L}}$", fontsize=24)
        plt.draw()
        
        image['name'] = f"{image['name']}_sample"
        image['fig'] = plt
        image['file'] = os.path.join(self.yaml['Plot_Config']['save_dir'], image['name'])
        self.savefig(image, plt)
        # plt.close()
        if self.yaml['Plot_Config']['screen_show']:
            plt.show(block=False)
            # plt.pause(1)
            input("\n\tPress 'Enter' to continue ...")
        plt.close()
            
    def plot_dynesty_results(self, name) -> None: 
        import matplotlib.pyplot as plt
        maxid = self.ddf.shape[0]
        nlive = self.ddf.iloc[0]['samples_nlive']
        
        fig = plt.figure(figsize=(10, 11))
        ax = fig.add_axes([0.01, 0.99-0.5/11., 0.05, 0.5/11.])
        draw_logo_in_square(ax)
        # from plot import draw_logo_in_square
        data = [
            {
                "y": self.ddf['samples_nlive'],
                "label":    "$\\text{Live points}$"
            },
            {   
                "y": np.exp(self.ddf['log_Like'] - max(self.ddf['log_Like'])),
                "label":    "Likelihood"
            },
            {
                "y": np.exp(self.ddf['log_weight']),
                "label":    "Importance\nweight PDF"
            },
            {
                "y": np.exp(self.ddf['log_Evidence']),
                "label":    "Evidence"
            },
            {
                "y": self.ddf['samples_it'],
                "label":    "Iters" 
            }
        ]
        from copy import deepcopy
        dtt = np.array(deepcopy(data[2]['y']))
        # endid = np.where(dtt > dtt[-1])[-1][-1]
        # lnx_end = self.ddf.iloc[endid]['log_PriorVolume']
        # print()

        from matplotlib.ticker import AutoMinorLocator
        ax1 = fig.add_axes([0.15, 0.8/11., 0.83, 2/11.])
        ax1.scatter( - self.ddf['log_PriorVolume'], data[0]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax1.set_xlim(0., max(- self.ddf['log_PriorVolume']))
        ax1.set_ylim(0., max(data[0]['y']) * 1.1)
        ax1.set_ylabel(data[0]['label'], fontsize=18)
        ax1.yaxis.set_label_coords(-0.08, 0.5)
        ax1.yaxis.set_minor_locator(AutoMinorLocator())
        ax1.xaxis.set_minor_locator(AutoMinorLocator())
        ax1.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax1.tick_params(which='major', length=7)
        ax1.tick_params(which='minor', length=4)
        # ax1.set_xticklabels([])
        plt.draw()


        ax2 = fig.add_axes([0.15, 2.8/11., 0.83, 2/11.])
        ax2.scatter( - self.ddf['log_PriorVolume'], data[1]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax2.plot( - self.ddf['log_PriorVolume'], data[1]['y'], '-', linewidth=0.8, alpha=0.4, color="#3f51b5" )
        ax2.set_xlim(0., max(- self.ddf['log_PriorVolume']))
        ax2.set_ylim(0., max(data[1]['y']) * 1.1)
        ax2.set_ylabel(data[1]['label'], fontsize=18)
        ax2.yaxis.set_label_coords(-0.08, 0.5)
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax2.tick_params(which='major', length=7)
        ax2.tick_params(which='minor', length=4)
        ax2.set_xticklabels([])


        plt.draw()

        from scipy.stats import gaussian_kde
        from scipy.interpolate import interp1d
        ax3 = fig.add_axes([0.15, 4.8/11., 0.83, 2/11.])
        dtt = dtt / max(dtt)
        ax3.scatter( - self.ddf['log_PriorVolume'], dtt, marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        wt_kde = gaussian_kde(np.array(- self.ddf['log_PriorVolume']), weights=np.array(data[2]['y']), bw_method="silverman")
        logvol = np.linspace(min(- self.ddf['log_PriorVolume']), max(- self.ddf['log_PriorVolume']), 1000)
        wt = wt_kde(logvol)
        wt = wt / wt.max()
        ax3.plot(logvol, wt, '-', linewidth=1.8, color='#8bc34a', alpha=0.8)
        ax3.set_xlim(0., max(- self.ddf['log_PriorVolume']))
        ax3.set_ylim(0., 1.1)
        ax3.set_ylabel(data[2]['label'], fontsize=18)
        ax3.yaxis.set_label_coords(-0.08, 0.5)
        ax3.yaxis.set_minor_locator(AutoMinorLocator())
        ax3.xaxis.set_minor_locator(AutoMinorLocator())
        ax3.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax3.tick_params(which='major', length=7)
        ax3.tick_params(which='minor', length=4)
        ax3.set_xticklabels([])
        ax3.text(0.02, 0.96, "Normalized to the maximum value", ha='left', va='top', transform=ax3.transAxes)

        plt.draw()

        ax4 = fig.add_axes([0.15, 6.8/11., 0.83, 2/11.])
        ax4.scatter( - self.ddf['log_PriorVolume'], data[3]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax4.set_xlim(0., max(- self.ddf['log_PriorVolume']))
        ax4.set_ylim(0., max(data[3]['y']) * 1.4)
        ax4.set_ylabel(data[3]['label'], fontsize=18)
        ax4.yaxis.set_label_coords(-0.08, 0.5)
        ax4.yaxis.set_minor_locator(AutoMinorLocator())
        ax4.xaxis.set_minor_locator(AutoMinorLocator())
        ax4.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax4.tick_params(which='major', length=7)
        ax4.tick_params(which='minor', length=4)
        ax4.set_xticklabels([])
        plt.draw()
        offset_text = ax4.yaxis.get_offset_text().get_text()
        maxx = str(max(data[3]['y'])).split("e")[0][0:5]
        upp = offset_text.split("e")[-1]
        txt = r"$\times 10^{" + upp + r"}$"
        ax4.text(0.4, max(data[3]['y']) * 0.95, maxx+txt, ha='left', va='top')
        ax4.plot([0, max(-self.ddf['log_PriorVolume'])], [max(data[3]['y']), max(data[3]['y'])], '-', linewidth=0.8, color="grey", alpha=0.4)

        logzerr = np.array(self.ddf['log_Evidence_err'])  
        logzerr[~np.isfinite(logzerr)] = 0.
        ax4.fill_between(- self.ddf['log_PriorVolume'], np.exp(self.ddf['log_Evidence'] - logzerr  ), np.exp(self.ddf['log_Evidence'] + logzerr), color='#8bc34a', alpha=0.2)

        ax5 = fig.add_axes([0.15, 8.8/11., 0.83, 2/11.])
        ax5.scatter( - self.ddf['log_PriorVolume'], data[4]['y'], marker='.', s=0.5, alpha=0.4, color="#3f51b5" )
        ax5.set_xlim(0., max(- self.ddf['log_PriorVolume']))
        ax5.set_ylim(0., max(data[4]['y']) * 1.1)
        ax5.set_ylabel(data[4]['label'], fontsize=18)
        ax5.yaxis.set_label_coords(-0.08, 0.5)
        ax5.yaxis.set_minor_locator(AutoMinorLocator())
        ax5.xaxis.set_minor_locator(AutoMinorLocator())
        ax5.tick_params(labelsize=11,  direction="in", bottom=True, left=True, top=True, right=True, which='both')
        ax5.tick_params(which='major', length=7)
        ax5.tick_params(which='minor', length=4)
        ax5.set_xticklabels([])

        # ax1.axvline( -lnx_end, color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        # ax2.axvline( -lnx_end, color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        # ax3.axvline( -lnx_end, color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        # ax4.axvline( -lnx_end, color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)
        # ax5.axvline( -lnx_end, color="#f44336", linestyle=":", linewidth=2.0, alpha=0.6)

        plt.draw()

        ax1.set_xlabel(r"$-\ln(X)$", fontsize=24)
        
        os.makedirs(self.yaml['Plot_Config']['save_dir'], exist_ok=True)
        
        img = {}
        img['fig'] = plt
        img['file'] = os.path.join(self.yaml['Plot_Config']['save_dir'], name)
        img['name'] = name
        self.savefig(img, plt)
        if self.yaml['Plot_Config']['screen_show']:
            plt.show(block=False)
        #     plt.pause(1)
            input("\n\tPress 'Enter' to continue ...")
        plt.close()
        
        # plt.savefig(savepath, dpi=300)

    def savefig(self, fig, plt):
        from matplotlib.backends.backend_pdf import PdfPages
        support_fmt_list = ['ps', 'eps', 'pdf', 'pgf', 'png', 'raw',
                            'rgba', 'svg', 'svgz', 'jpg', 'jpeg', 'tif', 'tiff']
        if "save_format" not in self.yaml["Plot_Config"]:
            self.yaml["Plot_Config"]['save_format'] = ['pdf']

        unsupport = []
        support = []
        file_list = []
        for fmt in self.yaml["Plot_Config"]['save_format']:
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
            elif fmt == "pdf" and ('ps' not in self.yaml["Plot_Config"]['save_format']):
                pp = plt
                pp.savefig("{}.{}".format(fig['file'], fmt))
            else:
                pp = plt
                pp.savefig("{}.{}".format(fig['file'], fmt), dpi=300)

        if ('pdf' not in support) and ('ps' in support):
            os.remove("{}.pdf".format(fig['file']))

        self.logger.warning("Figure {} saved in the path\n\t\t-> {} \n\t\t>> {}.".format(
            fig['name'], os.path.dirname(fig['file']), ", >> ".join(file_list)))
        if unsupport:
            self.logger.warning(emoji.emojize('\t:ghost::ghost::ghost: Figure format unsupport -> {}. '.format(
                ", ".join(unsupport)), use_aliases=True))

    def compress_figure_to_PS(self, figpath):
        os.system('pdf2ps {}.pdf {}.ps'.format(figpath, figpath))