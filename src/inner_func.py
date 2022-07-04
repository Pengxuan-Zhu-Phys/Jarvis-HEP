#!/usr/bin/env python3
import os, sys

def Gauss(xx, mean, err):
    from math import sqrt, pi, exp
    prob = 1./ (err * sqrt(2 * pi)) * exp(-0.5*((xx - mean)/err)**2)
    return prob

def updata_funcs(funcs):
    funcs['Gauss'] = Gauss
    return funcs