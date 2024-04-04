#!/usr/bin/env python3 
import os, sys 
import json

class Tree(object):
    def __init__(self, *args):
        self.layer = []
        self.nodes = {}
        self.info  = {}
        self.pars  = None
        self.plot_input = None
        self.BPlot_path = None
        
    def addNode(self, pkg, info):
        node = Node()
        node.name = pkg 
        node.config = info 
        node.parent = info['required package']
        self.nodes[pkg] = node
        
    def makeTree(self):
        tba = []
        rootlayer = []
        rba = []
        for pkg, node in self.nodes.items():
            if not node.config['multi']:
                if node.parent == None:
                    rootlayer.append(node.name)
                    rba.append(node.name)
                    node.depth = 0
                else:
                    tba.append(pkg)
            else:
                for pp in node.config['include']:
                    if node.parent == None:
                        rootlayer.append("{}\n@{}".format(pp, pkg))
                        rba.append("{}\n@{}".format(pp, pkg))
                        node.depth = 0
                    else:
                        tba.append(pkg)
        self.layer.append(rootlayer)
        nn = 1
        while tba:
            newlayer = []
            for pkg in tba:
                rbaflag = True
                rooflag = False
                for tt in self.nodes[pkg].parent:
                    if pkg not in self.nodes[tt].child:
                        self.nodes[tt].child.append(pkg)
                    if tt not in rba:
                        rbaflag = False
                        break
                    if tt in rootlayer:
                        rooflag = True
                if (rbaflag and rooflag):
                    newlayer.append(pkg)
            if newlayer:
                self.layer.append(newlayer)
                rootlayer = newlayer
                for pkg in newlayer:
                    tba.remove(pkg)
                    rba.append(pkg)
        for pkg, node in self.nodes.items():
            for ii in range(len(self.layer)):
                if pkg in self.layer[ii]:
                    node.depth = ii 
                    break
        
    def makeDrawInfo(self):
        self.info['plot'] = {
            "layer":    self.layer,
            "nodes":    {}
        }
        # print(self.layer)
        for pkg, node in self.nodes.items():
            # print(json.dumps(node.config, indent=4) )
            if not node.config['multi']:
                self.info['plot']['nodes'][pkg] = {
                    "input file":   [],
                    "output file":  []
                }
                for kk, ff in node.config['input file'].items():
                    filinfo = {
                        "name": os.path.basename(ff),
                        "vars": {}
                    }
                    for var in node.config['input variables']:
                        if var['file'] == ff:
                            filinfo['vars'][var['expr']] = var['meth']
                    self.info['plot']['nodes'][pkg]['input file'].append(filinfo)
                for kk, ff in node.config['output file'].items():
                    filinfo = {
                        "name": os.path.basename(ff),
                        "vars": {}
                    }
                    for var in node.config['output variables']:
                        if var['file'] == ff:
                            filinfo['vars'][var['expr']] = var['meth']
                    self.info['plot']['nodes'][pkg]['output file'].append(filinfo)           
            else:
                for pp in node.config['include']:
                    self.info['plot']['nodes']["{}\n@{}".format(pp, pkg)] = {
                        "input file":   [],
                        "output file":  []
                    }
                    for kk, ff in node.config['modes'][pp]['input file'].items():
                        filinfo = {
                            "name": os.path.basename(ff),
                            "vars": {}
                        }
                        for var in node.config['modes'][pp]['input variables']:
                            if var['file'] == ff:
                                filinfo['vars'][var['expr']] = var['meth']
                        self.info['plot']['nodes']["{}\n@{}".format(pp, pkg)]['input file'].append(filinfo)
                    for kk, ff in node.config['modes'][pp]['output file'].items():
                        filinfo = {
                            "name": os.path.basename(ff),
                            "vars": {}
                        }
                        for var in node.config['modes'][pp]['output variables']:
                            if var['file'] == ff:
                                filinfo['vars'][var['expr']] = var['meth']
                        self.info['plot']['nodes']["{}\n@{}".format(pp, pkg)]['output file'].append(filinfo)   



        self.info['scan variables'] = []
        for kk, vv in self.pars.items():
            self.info['scan variables'].append({
                "name":     kk,
                "prior":    vv['prior']
            })  


        

            
    
    
class Node(object):
    def __init__(self) -> None:
        self.name   = None
        self.depth  = 0
        self.parent = None
        self.child  = []
        self.config = None
    
    def print_info(self):
        print(self.name)
        
        