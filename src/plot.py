#!/usr/bin/env python3

import os, io 
import sys
import pandas as pd
import time
import configparser
import math
import sympy
import json
import emoji
import numpy as np 
import matplotlib.pyplot as plt
import shapely as sp 
import contextlib

pwd = os.path.abspath(os.path.dirname(__file__))
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.abspath(os.path.join(pwd, "BudingPLOT")))


def draw_logo_in_square(ax):
    with io.StringIO() as buf, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):

        ax.axis("off")
        ax.text(0.5, -0.06, "Jarvis-HEP", fontsize=12, ha="center", va='top', transform=ax.transAxes)
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