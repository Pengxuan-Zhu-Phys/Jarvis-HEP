#!/usr/bin/env python3

import os, sys 
pwd = os.path.abspath(os.path.dirname(__file__))
sys.path.append("{}/".format(pwd))

def Dfloat(string):
    try:
        return float(string)
    except ValueError:
        return float(string.upper().replace('D','E',1))

from functools import wraps
# from ...operators.string import Dfloat
from SLHA_line import LoopLines

scalar_groups={
    'SUSY_input':   ['MINPAR','EXTPAR'],
    'additional':   ['NMSSMRUN','MSOFT'],
    'output'    :   ['MASS','SPhenoLowEnergy','FlavorKitQFV','LHCFIT','FINETUNING'],
    'omega'     :   ['ABUNDANCE','LSP','NDMCROSSSECT','INDIRECT_CHISQUARES']
}
matrix_groups={
    'Mass'      :   ['MSD2','MSE2','MSL2','MSQ2','MSU2'],
    'Mix'       :   ['NMHMIX','NMAMIX','STOPMIX','NMNMIX','UMIX','VMIX', 'SBOTMIX', 'STAUMIX'],
    'Triliner'  :   ['TD','TE','TU'],
    'SeeSaw'    :   ['MUX','MV2','MX2','YV','TV','BMUX','LAMN','TLAMN'],
    'output'    :   ['YE','YU','YD',
                     'HiggsLHC13','HiggsLHC14','REDCOUP']
}
scalar_list=[ i.upper() for j in scalar_groups.values()  for i in j ]
matrix_list=[ i.upper() for j in matrix_groups.values()  for i in j ]


@staticmethod
@LoopLines
def ReadScalar(line):
    s=line.split()
    code=int(s[0])
    value=Dfloat(s[1])
    return {code:value}
@staticmethod
@LoopLines
def ReadMatrix(line):
    s=line.split()
    code=(int(s[0]),int(s[1]))
    value=Dfloat(s[2])
    return {code:value}

from special_blocks import special_blocks
class ReadBlock(special_blocks): # Container of methods to read SLHA data
    pass

for name in scalar_list:
    setattr(ReadBlock,name,ReadScalar)
for name in matrix_list:
    setattr(ReadBlock,name,ReadMatrix)

class SLHA_block:
    '''block'''
    def __init__(self,block_text,block_format=ReadBlock):
        self.text_dict=block_text
        self.block_format=block_format
    def __getattr__(self,block_name):
        # print(f'find {block_name}')
        try: text=self.text_dict[block_name]
        except KeyError: 
            print(f'{block_name} not found in text')
            raise
        try: data=getattr(self.block_format,block_name)(text)
        except AttributeError:
            print(f'Load method for {block_name} not found in block_format')
            print(*self.text_dict[block_name])
            raise
        setattr(self,block_name,data)
        return data
