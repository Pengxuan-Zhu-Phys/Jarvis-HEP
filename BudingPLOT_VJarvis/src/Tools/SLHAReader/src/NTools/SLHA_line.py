#!/usr/bin/env python3

def GetBlockName(line):
    return line.split()[1].upper()

def GetDecayCode(line):
    return int(line.split()[1])

from functools import wraps
def LoopLines(func):
    '''
    Read SLHA informations line by line with wrapped functions.
    Return a dictionary
    '''
    @wraps(func)
    def wrapper(lines,*args,**kwargs):
        data=dict()
        for line in lines:
            try:
                data.update(func(line))
            except ValueError: pass
            except IndexError: pass
            except TypeError: print(line,func(line),sep='\n');raise
        return data
    return wrapper