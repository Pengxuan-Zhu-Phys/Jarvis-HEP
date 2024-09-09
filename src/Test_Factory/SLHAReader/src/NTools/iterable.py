#!/usr/bin/env python3

FTL=lambda x: sum(map(FTL,x),[]) if isinstance(x,(list,tuple)) else [x]
def FlatToList(*Lists):
    '''
    extract all elements in multilayer lists and dict.keys()s to one list.
    '''
    # DL=[list(i)for i in Lists]
    return FTL([list(i)for i in Lists])

def SortMixList(mix_list):
    try:
        return sorted(mix_list)
    except TypeError:
        type_list=set([type(i) for i in mix_list])
        deep_list=[]
        for i,type_i in enumerate(type_list):
            deep_list.extend(sorted([l_i for l_i in mix_list if type(l_i) is type_i]))
        return deep_list