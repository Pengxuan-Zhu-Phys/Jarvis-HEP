#!/usr/bin/env python3 

import xmltodict
import json 
import os, sys 

def read_variable_value(cont, var):
    flag = True
    value = None
    while flag:
        kk = var['code'].pop(0)
        if kk in cont.keys():
            cont = cont[kk]
        if not len(var['code']):
            flag = False
    if var['patt'] == "number":
        vls = []
        for vv in cont.split():
            try:
                float(vv)
                a = eval(vv)
                vls.append(a)
            except:
                pass 
    value = vls[var['posi']-1]
        
    return value
    
def load_yoda_data(yoda):
    import re
    plots = []
    from numpy import loadtxt
    from pandas import DataFrame
    with open(yoda, 'r') as f1:
        yd = f1.read()
        yd = yd.replace("\n", " endl " )
        p_rec = re.compile(r"BEGIN(.*?)END", re.M)
        a = re.findall(p_rec, yd)
        for item in a:
            item = item.replace(" endl ", "\n").split('\n')
            if ("Weight_MERGING" in item[0]) or ("RAW" in item[0]):
                continue
            elif "YODA_HISTO1D_V2" in item[0]:
                res = {
                    "name": item[0].split()[1].strip(),
                    "type": item[0].split()[0].strip(),
                    "sumw": float(item[item.index("---") + 4].split()[2]),
                    "data": DataFrame(loadtxt(item[item.index("---") + 8:]), columns=item[item.index("---") + 7][1:].split())
                }
                plots.append(res)
            elif "YODA_HISTO2D_V2" in item[0]:
                res = {
                    "name": item[0].split()[1].strip(),
                    "type": item[0].split()[0].strip(),
                    "sumw": float(item[item.index("---") + 4].split()[2]),
                    "data": DataFrame(loadtxt(item[item.index("---") + 7:]), columns=item[item.index("---") + 7][1:].split())
                }
                plots.append(res)
    return plots    

    
if __name__ == "__main__":
    var = {
        "expr": "XSect",
        "file": "run_01_tag_1_banner.txt",
        "meth": "xml",
        "patt": "number",
        "posi": 2,
        "code": ["LesHouchesEvents", "header", "MGGenerationInfo"]
    }
    res = read_variable_value(var)
    print(res)