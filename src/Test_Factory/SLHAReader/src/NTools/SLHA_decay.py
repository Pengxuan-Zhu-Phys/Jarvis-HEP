#!/usr/bin/env python3

from SLHA_block import Dfloat

def FormatTo_float_tuple(line):
    splited=line.split()
    value=float(splited[0])
    member_number=int(splited[1])
    member_code=tuple([int(i) for i in splited[2:2+member_number]])
    return {member_code: value}

def ReadDecay(text):
    data={}
    splited=text[0].split()
    if splited[0].upper()=='DECAY':
        data['WIDTH']=Dfloat(splited[2])
    for line in text:
        try:
            data.update(FormatTo_float_tuple(line))
        except ValueError: continue
        except IndexError: continue
    return data

class SLHA_decay:
    '''particals' decay information'''
    def __init__(self,decay_text):
        self.text_dict=decay_text
        self.data={}

    def ReadDecay(self,p_code=None):
        if p_code is None:
            for p_code in self.text_dict.keys():
                self.ReadDecay(p_code)
        else:
            try: text=self.text_dict[p_code]
            except KeyError:
                print(f'DECAY for particle-{code} not found in text')
                raise
            self.data[p_code]=ReadDecay(text)
    def __getitem__(self,p_code):
        try: return self.data[p_code]
        except KeyError:
            self.ReadDecay(p_code)
            return self[p_code]
    def __getattr__(self,p_str):
        if p_str[0]=='p':
            return self[int(p_str[1:])]
        else:
            print(f'Wrong Key: {p_str}')
            print(r'argument "p_str" should be p{PDG} like p24 for Z decay')