#! /usr/bin/env python3
import os, sys 
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append("{}/".format(pwd))
from SLHA_line import LoopLines

class special_blocks():
    @LoopLines
    def HIGGSCOUPLINGSBOSONS(line):
        s=line.split()
        v=float(s[0])
        np=int(s[1])
        tuple_p=tuple([int(i) for i in s[2:2+np]])
        return {tuple_p:v}
    @LoopLines
    def HIGGSCOUPLINGSFERMIONS(line):
        s=line.split()
        v=max(float(s[0]),float(s[1]))
        np=int(s[2])
        tuple_p=tuple([int(i) for i in s[3:3+np]])
        return {tuple_p:v}
    @LoopLines
    def HIGGSBOUNDSRESULTS(line):
        semanteme=line.split()
        if semanteme[1]=='1':
            return {'channel'+semanteme[0]:int(semanteme[2])}
        elif semanteme[1]=='2':
            return {'HBresult'+semanteme[0]:int(semanteme[2])}
        elif semanteme[1]=='3':
            return {'obsratio'+semanteme[0]:float(semanteme[2])}
        else: return {}
    @LoopLines
    def HIGGSSIGNALSRESULTS(line):
        semanteme=line.split()
        if semanteme[0]=='13':
            return {'Probability':float(semanteme[1])}
        elif semanteme[0]=='12':
            return {'X2_total':float(semanteme[1])}
        else: return {}

special_blocks.HIGGSBOUNDSINPUTHIGGSCOUPLINGSBOSONS=special_blocks.HIGGSCOUPLINGSBOSONS
special_blocks.HIGGSBOUNDSINPUTHIGGSCOUPLINGSFERMIONS=special_blocks.HIGGSCOUPLINGSFERMIONS
