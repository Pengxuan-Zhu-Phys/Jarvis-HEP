#!/usr/bin/env python3
import os, sys
import subprocess
from uuid import uuid4

from sympy.core.evalf import fastlog 
from scan import jpath
import json 
sys.path.append(os.path.join(os.path.dirname(__file__), "SLHAReader/src/"))

def set_inputvariable(var):
    if var['meth'] == "position":
        set_position_var(var)
    elif var['meth'] == "replace":
        set_replace_var(var)
        
def set_replace_var(var):
    with open(var['file'], 'r') as f1:
        cont = f1.read()
    cont = cont.replace(var['code'], str(var['expr']))
    with open(var['file'], 'w') as f1:
        f1.write(cont)

def set_position_var(var):
    from numpy import loadtxt, savetxt
    data = loadtxt(var['file'], ndmin=2)
    data[var['code'][0]-1][var['code'][1]-1] = var['expr']
    savetxt(var['file'], data)
    
def get_outputvariable(var):
    if var['meth'] == "json":
        res = get_Jsonout(var)
        return res
    elif var['meth'] == "file":
        res = get_Fileout(var)
        return res 
    elif var['meth'] == "slha":
        res = get_SLHAout(var)
        return res 
    elif var['meth'] == "xml":
        res = get_XMLout(var)
        return res 

def get_Fileout(var):
    res = {
        "error tag":    None,
        "error type":   None,
        "value":        None
    }
    if os.path.exists(var['file']):
        if os.path.exists(var['code']) and var['save']:
            savedir = os.path.join(var['code'], "{}_{}".format(var['tag'], os.path.basename(var['file'])))
            try:
                from shutil import copyfile
                copyfile(var['file'], savedir)
                res['error tag']    = False
                res['value']        = savedir
            except:
                res['error tag']    = True
                res['error type']   = "CopyError"
                res['value']        = savedir
        else:
            res['error tag']    = True
            res['error type']   = "CopyError"
            res['value']        = savedir
    else:
        res['error tag']    = True
        res['error type']   = "FileError"
        res['value']        = var['file']
    return res

def get_XMLout(var):
    res = {
        "error tag":    None,
        "error type":   None,
        "value":        None
    }
    if os.path.exists(var['file']):
        try:
            from xmlreader import read_variable_value
            res['value'] = read_variable_value(var)
            res['error tag'] = False
        except:
            res['error tag']    =   True
            res['error type']   =   "ValueError"
            res['value']        =   var['file']
    else:
        res['error tag']    =   True
        res['error type']   =   "FileError"
        res['value']        =   var['file']
    return res 

def get_SLHAout(var):
    res = {
        "error tag":    None,
        "error type":   None, 
        "value":        None
    }
    if os.path.exists(var['file']):
        try:
            from slha import read_variable_value
            res['value'] = read_variable_value(var)
            res['error tag'] = False
        except:
            res['error tag']    = True
            res['error type']   = "ValueError"
            res['value']        = var['file']
    else:
        res['error tag']    = True
        res['error type']   = "FileError"
        res['value']        = var['file']
    return res 

def get_Jsonout(var):
    res = {
        "error tag":    None,
        "error type":   None,
        "value":        None
    }
    if os.path.exists(var['file']):
        with open(var['file'], 'r') as f1:
            data = json.loads(f1.read())
            try:
                res['value'] = data[var['expr']]
                res['error tag']    = False
            except:
                res['error tag']    = True
                res['error type']   = "ValueError"
                res['value']        = var['file']
    else:
        res['error tag']    = True
        res['error type']   = "FileError"
        res['value']        = var['file']
    return res 

def decode_path_from_file(pathdir):
    if "&pwd" in pathdir:
        pwd = os.path.abspath(os.path.dirname(pathdir))
        pathdir = pathdir.replace("&pwd", pwd)
    if "&J" in pathdir:
        pathdir = pathdir.replace("&J", jpath)
    return pathdir

def ck_sect(cf, sect):
    if cf.has_section(sect):
        return True
    else:
        print("\tNo Section Found: -> {},  Terminating !!!".format(sect))
        sys.exit(0)
        return False

def checking_os():
    return sys.platform

def count_cores(ips):
    if ips == "cpu_count":
        return os.cpu_count()
    elif ips.isdigit():
        return int(ips)

def check_cern_root(ips):
    if ips[0:2] == "@C":
        import subprocess
        root_path = subprocess.getstatusoutput(ips[2:].strip())
        return not root_path[0], root_path[1]
    elif os.path.exists(ips):
        return True, ips
    else:
        return False, ""

def check_python2_version(ips):
    ips = ips.strip()
    if ips[0] == "2":
        cv = subprocess.getoutput("python2 -V").replace("Python ", "")
    ips = ips.split('.')
    pv  = cv.split('.')
    res = True
    if len(pv) == len(ips):
        for ii in range(len(pv)):
            if int(pv[ii]) < int(ips[ii]):
                res = False
    return res, cv

class InputTimeoutError(Exception):
    pass

def wait_for_input(info, timeout, default):
    from select import select
    print(info)
    ipt, o, e = select( [sys.stdin], [], [], timeout)
    if ipt:
        ipt = sys.stdin.readline().strip()
    else:
        ipt = default
    return ipt

def get_sample_id():
    return uuid4()

def Help():
    print("Help")





def decode_stdout(output):
    output = output.split("\n")
    output.remove("")
    output = "\n\t".join(output)
    output += "\n"
    output = "Screen Output -> \n\t" + output
    return output

def force_kill_process(proc_pid):
    import psutil
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        try:
            proc.kill()
        except:
            pass 
    try:
        process.kill()
    except:
        pass 

def get_time_lock(ts):
    from time import time 
    return int(time()) // ts 

if __name__ == "__main__":
    decode_path(sys.argv[1], "&pwd/makefile.json")
    decode_path(sys.argv[1], "&CL/makefile.json")

    from sympy import *
    expr = sympify("cos(x) + sin(y) + 1")
    expr.subs({symbols("x"): 0, symbols("y"): math.pi})


def format_PID(idx:int):
    return "{:0>7}".format(idx)