#!/usr/bin/env python3 

import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import os, sys 
from _pytest import config
import json 
import emoji
import math 
import matplotlib.colors as mcolor 
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm

pwd = os.path.abspath(os.path.dirname(__file__))
jpath = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Update the colors 
color_card = os.path.join(pwd, "cards/named_color.json")
with open(color_card, 'r') as f1: 
    setting = json.loads(f1.read())
    matplotlib.colors.CSS4_COLORS.update(setting)
    named_colors = mcolor.get_named_colors_mapping()
    named_colors.update(setting)
    
# Modify the default color cycle 
from cycler import cycler
matplotlib.rcParams['axes.prop_cycle'] = cycler(color=['green', 'blue', 'red', 'yellow', 'indigo', "purple", "cyan", "amber", 'orange'])


plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'

import warnings
warnings.filterwarnings("ignore", message="There are no gridspecs with layoutgrids.")

with open(os.path.join(pwd, "cards/colormaps.json")) as f:
    cmap_data = json.load(f)
    
for cmap_name, color_list in cmap_data.items():
    new_cmap = LinearSegmentedColormap.from_list(cmap_name, color_list)
    new_cmap_r = new_cmap.reversed()

    matplotlib.colormaps.register(cmap=new_cmap)
    matplotlib.colormaps.register(cmap=new_cmap_r)
    