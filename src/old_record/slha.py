#!/usr/bin/env python3 

import os, sys 
import pyslha, xslha
from sympy.matrices.expressions import special
from sympy.functions.special.gamma_functions import uppergamma 

def get_variable_info(inf):
    res = {}
    if inf[3].strip().upper() == "BLOCK":
        if len(inf) == 6:
            res = {
                "patt": inf[4].strip().upper(),
                "code": int(inf[5].strip())
            }
        elif len(inf[4:]) > 2:
            res = {
                "patt": inf[4].strip().upper(),
                "code": tuple(list(map(int, inf[5:])))
            } 
    elif inf[3].strip().upper() == "DECAY":
        if len(inf) == 5:
            res = {
                "patt": "WIDTH",
                "code": (int(inf[4].strip()))
            }
        elif len(inf[6:]) >= 2:
            res = {
                "patt": "DECAY",
                "pdg":  int(inf[4].strip()),
                "finalstates":  [tuple(list(map(int, inf[6:])))]
            }
    return res 

def read_variable_value(var):
    readlib = "xSLHA"
    try:
        spectr = xslha.read(var['file'])
    except:
        spectr = pyslha.read(var['file'], ignorenomass=True)
        readlib = "PySLHA"
    if readlib == "xSLHA":
        value = read_variable_value_by_xSLHA(spectr, var)
    elif readlib == "PySLHA":
        value = read_variable_value_by_PySLHA(spectr, var)
    return value 
        
def read_variable_value_by_xSLHA(spectr, var):
    if var['patt'] == "WIDTH":
        try:
            width = spectr.Value("WIDTH", var['code'])
        except:
            width = 0. 
        return width 
    elif var['patt'] == "DECAY":
        br = 0. 
        from itertools import permutations
        try:
            for fstate in var['finalstates']:
                br += spectr.Value("BR", [var['pdg'], fstate])
        except:
            br += 0. 
        return br  
    else:
        if type(var['code']) == int:
            value = spectr.Value(var['patt'], [var['code']])
        else:
            value = spectr.Value(var['patt'], var['code'])
        return value

def read_variable_value_by_PySLHA(var):
    print(var)
    return 0.