#!/usr/bin/env python3

from SLHA_block import  SLHA_block, ReadBlock, ReadScalar
from SLHA_decay import SLHA_decay
from object import lazyproperty
from SLHA_document import SLHA_text

class NTools_block_format(ReadBlock):
    LHCCROSSSECTIONS=ReadScalar
    DELTAMH=ReadScalar
    LOWEN=ReadScalar
    
    @staticmethod
    def SPINFO(lines):
        data={}
        for line in lines:
            code=line.strip()[0]
            if code in ['3','4']:
                data.setdefault(int(code),[]).append(line.split('#')[1])
        return data


class NToolsOutput(SLHA_text):
    def __init__(self,spectr_dir,*omega_dir,ignore=[]):
        self.spectr_dir=spectr_dir
        self.ignore=ignore
        output_text=[] 
        with open(spectr_dir,'r') as spectr:
            for line in spectr:
                if line=='# BLOCK FINETUNING\n':
                    output_text.append(line[2:])
                else:
                    output_text.append(line)
        try:
            with open(omega_dir[0],'r') as omega:
                output_text.extend(omega.readlines())
        except FileNotFoundError:
            pass
        except IndexError:
            pass
        else:
            self.omega_dir=omega_dir[0]
        super().__init__(output_text,block_format=NTools_block_format)
    @property
    def constraints(self):
        all_consts=[]
        for constraint in self.BLOCK.SPINFO[3]:
            for const in self.ignore:
                if const in constraint:
                    break
            else:
                all_consts.append(constraint)
        return all_consts
    @property
    def error(self):
        return not bool(self.constraints)

class SLHA_document(SLHA_text):
    def __init__(self,SLHA_document,block_format=ReadBlock):
        self.path=SLHA_document
        self.block_format=ReadBlock
    @lazyproperty
    def text(self):
        # print(f'Getting text from {self.path}')
        with open(self.path,'r') as SLHA:
            return SLHA.readlines()