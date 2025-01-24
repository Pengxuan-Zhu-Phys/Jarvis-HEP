#!/usr/bin/env python3

import networkx as nx 
from base import Base
from Module.module import Module
import numpy as np 
import contextlib
import os, io 
import asyncio
from numpy.lib._type_check_impl import imag

class Workflow(Base):
    def __init__(self):
        super().__init__()
        self.modules = {}
        self.module_states = {}
        self.module_dependencies = {}
        self.parameter_module = None
        self.library_modules = {}
        self.graph = {}
        self.layer_info = {}
        self.layers = []
        self.workflow = {}

    def add_module(self, module):
        self.modules[module.name] = module
        self.module_states[module.name] = 'waiting'

    def add_library_module(self, module):
        self.library_modules[module.name] = module
        self.modules[module.name] = module

    def set_modules(self, modules):
        from Module.parameters import Parameters
        from Module.library import LibraryModule
        from Module.calculator import CalculatorModule
        
        parameter_module = Parameters("Parameters", modules['Parameter'])
        parameter_module.analyze_ios()
        self.add_module(parameter_module)
        self.parameter_module = parameter_module
        
        if hasattr(modules, "Library"):
            for lib in modules['Library']:
                module = LibraryModule(
                    name=lib['name'],
                    required_modules=lib.get("required_modules", []),
                    installed=lib.get("installed", False),
                    installation=lib.get("installation", {})
                )
                self.add_library_module(module)

        for calc in modules['Calculator']:
            module = CalculatorModule(
                name=calc['name'],
                config=calc
            )
            self.add_module(module)

    def get_workflow_dict(self):
        for kk, layer in self.layer_info.items():
            if kk > 1: 
                self.workflow[kk] = list(layer.keys())
    
    def resolve_layers(self):
        self.calc_layer = {}  
        self.libr_layer = {}

        unresolved_all = []
        unresolved_lib = []
        unresolved = []
        for name in self.modules.keys():
            if name in self.library_modules:
                unresolved_lib.append(name)
            else:
                unresolved.append(name)
            unresolved_all.append(name)

        # resolve the library layers 
        current_layer = 1
        resolved = []
        while unresolved_lib:
            new_resolved = []
            for ii in range(len(unresolved_lib)):
                if current_layer == 1:
                    if not self.modules[unresolved_lib[ii]].required_modules:
                        new_resolved.append(self.modules[unresolved_lib[ii]].name)
                else:
                    if set(self.modules[unresolved_lib[ii]].required_modules).issubset(set(resolved)):
                        new_resolved.append(self.modules[unresolved_lib[ii]].name)

            resolved.extend(new_resolved)
            unresolved_lib = [item for item in unresolved_lib if item not in resolved]
            self.libr_layer[current_layer] = new_resolved
            current_layer += 1         
        
        current_layer = 1
        while unresolved:
            # achieve the first layer module
            new_resolved = []
            for ii in range(len(unresolved)):
                if current_layer == 1:
                    if not self.modules[unresolved[ii]].required_modules:
                        new_resolved.append(self.modules[unresolved[ii]].name)
                else:
                    if set(self.modules[unresolved[ii]].required_modules).issubset(set(resolved)):
                        new_resolved.append(self.modules[unresolved[ii]].name) 
            
            resolved.extend(new_resolved)
            unresolved = [item for item in unresolved if item not in resolved]
            self.calc_layer[current_layer] = {"module": new_resolved, "ipfs": {}, "opfs": {}, "ipvs": {}, "opvs": {}, "mdls": {}}
            current_layer += 1 

    def resolve_dependencies(self):
        # Resolve dependencies based on the IOs 
        output_to_module = {} 
        for module in self.modules.values():
            for output in module.outputs:
                output_to_module[output] = module.name

        # Update required_modules via IOs 
        for module in self.modules.values():
            module.required_modules = module.required_modules or []
            for input_var in module.inputs:
                if input_var in output_to_module:
                    # if the output is the input of other module, then generate teh dependencies 
                    dep_module = output_to_module[input_var]
                    if dep_module not in module.required_modules:
                        module.required_modules.append(dep_module)

        # Resolve the layers using the current dependencies  
        self.resolve_layers()

    def analyze(self):
        # First, add edges based on module dependencies
        for module in self.modules:
            for output in module.outputs:
                for dependent_module in self.modules:
                    if output in dependent_module.inputs:
                        self.graph.add_edge(module.name, dependent_module.name)
        
        # Check for cycles which would indicate a problem in the dependency graph
        try:
            cycle = nx.find_cycle(self.graph)
            raise Exception(f"Cyclic dependency found: {cycle}")
        except nx.NetworkXNoCycle:
            pass  # All good if there are no cycles

    def run(self, sample_point):
        try:
            for module in self.modules:
                self.module_states[module.name] = 'running'
                result = module.run(sample_point)
                
                if not module.validate_result(result):
                    self.module_states[module.name] = 'failed'
                    print(f"Module {module.name} failed. Stopping workflow for this sample.")
                    # Stop the workflow for this sample point
                    return False  
                
                self.module_states[module.name] = 'completed'
            # Workflow completed successfully for this sample point    
            return True  
        except Exception as e:
            print(f"An error occurred: {e}")
            # An error occurred, stopping workflow for this sample point
            return False  

    def reset_states(self):
        for module in self.modules:
            self.module_states[module.name] = 'not started'
    
    def run_all_samples(self, sample_points):
        for sample in sample_points:
            self.reset_states()
            if not self.run(sample):
                print(f"Workflow failed for sample point {sample}. Moving to next sample.")
            else:
                print(f"Workflow completed successfully for sample point {sample}.")

    async def draw_flowchart(self, save_path="flowchart.png"):
        # with io.StringIO() as buf, contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
            import matplotlib.pyplot as plt 
            import logging
            logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
            logging.getLogger('PIL.PngImagePlugin').setLevel(logging.CRITICAL)
            from matplotlib.collections import LineCollection
            from plot import draw_logo_in_square
            from plot import create_round_square
            from PIL import Image

            layerW = 6.0 
            figL = layerW * len(self.calc_layer) - 3
            figH = 4

            def resolve_module(module, lid):
                if lid == 1:
                    res = {"name": module, "ipf": {}, "opf": {}, "ipv": {}, "opv": {}, 'mhh':0.9, 'ohh': 0., "ihh": 0., "bp": np.array([0.85, 0.]), "width": 0.9}
                else:     
                    res = {"name": module, "ipf": {}, "opf": {}, "ipv": {}, "opv": {}, 'mhh':0.9,'ohh': 0., "ihh": 0.0, "bp": np.array([0., 0.]), "width": 0.9}
                if self.modules[module].input:
                    ipfl = {}
                    nn = len(self.modules[module].input)
                    hh = 0.7 * nn - 0.2
                    for ii in range(nn):
                        ipf = self.modules[module].input[ii]
                        ipfl[ipf['name']] ={
                            "nam": ipf['name'],
                            "pos": np.array([-1.1, ((nn - 1)/2 - ii) * 0.7 + 0.2]) + res["bp"],
                            "inc": [],
                            "fil": os.path.basename(ipf['path'])
                        }
                        for var in ipf['variables'].values():
                            if "inc" in var.keys():
                                for kk in var['inc']:
                                    if kk not in ipfl[ipf['name']]['inc']:
                                        ipfl[ipf['name']]['inc'].append(
                                            {
                                                "name": str(kk)
                                            }
                                        )
                            else:
                                ipfl[ipf['name']]['inc'].append(var)
                    res["ipf"] = ipfl
                    if hh > 0.:
                        res['ihh'] = hh 

                if self.modules[module].output:
                    opfl = {}
                    nn = len(self.modules[module].output)
                    hh = 0.7 * nn - 0.2 
                    h0 = 0.
                    for ii in range(nn):
                        opf = self.modules[module].output[ii]
                        nv = len(opf['variables'])
                        if nv <= 2:
                            hl = 0.5 
                        else:
                            hl = 0.2 * nv
                        # print(opf['name'],)
                        op = {
                            "nam": opf['name'],
                            "pos": np.array([1.1, -h0 - 0.5*hl -0.2 ]) + res["bp"],
                            "inc": [var['name'] for var in opf['variables']],
                            'fil': os.path.basename(opf['path'])
                        }
                        h0 += hl 
                        h0 += 0.1
                        for jj in range(nv):
                            res['opv'][opf['variables'][jj]['name']] = {
                                "nam": opf['variables'][jj]['name'],
                                "pos": np.array([2.6, op["pos"][1] + ((nv-1)/2 - jj) * 0.2 ]) + res["bp"],
                                "wid": 0.
                            }
                        opfl[opf['name']] = (op)
                    h0 -= 0.1
                    res["opf"] = opfl 
                    if h0 > 0.:
                        res['ohh'] = h0
                    for op in res['opf']:
                        res['opf'][op]['pos'] += np.array([0., + 0.5 * h0])
                    for op in res["opv"]:
                        res['opv'][op]["pos"] += np.array([0., + 0.5 * h0])

                elif self.modules[module].outputs:
                    nv = (len(self.modules[module].outputs) - 1)/2
                    h0 = nv * 0.2
                    ii = 0
                    for op in self.modules[module].outputs:
                        res['opv'][op] = {
                            "nam": op, 
                            "pos": np.array([2.5, h0 - ii * 0.2]),
                            "wid": 0.
                        }
                        ii += 1
                
                
                if res['ipf']:
                    res['width'] += 0.65 
                if res['opf']:
                    res['width'] += 0.65
                if res['opv']:
                    res['width'] += 1.5 
                
                return res 

            def resolve_layer():
                for kk, item in self.calc_layer.items():
                    if kk == 1:
                        self.layer_info[1] = {
                            item["module"][0]: resolve_module(item["module"][0], kk)
                        }
                    else:
                        self.layer_info[kk] = {}
                        for module in item['module']:
                            self.layer_info[kk][module] = resolve_module(module, kk)

                layinfo = [None for ii in range(len(self.layer_info))]
                nlay = len(self.layer_info)
                for ii in range(nlay):
                    lid = len(self.layer_info) - ii
                    layer = self.layer_info[len(self.layer_info) - ii]
                    h0 = 0. 

                    y0 = {}
                    opvs = {}
                    for kk, item in layer.items():
                        mh = max([item['ihh'] + 0.4, item['ohh'] + 0.4, item['mhh']])
                        y00 = 0.5*mh + h0
                        item['bp'] += np.array([layerW * (lid - 1), y00])
                        h0 += mh 
                        y0[kk] = y00
                        for op, val in item['opv'].items():
                            opvs[op] = val
                        
                    layinfo[lid - 1] = {
                        "height": h0, 
                        "ipvars": {}, 
                        "y0":   y0, 
                        "opvars": opvs,
                        "bridge": {}
                    }
                for ii in range(nlay - 1):
                    clid = len(self.layer_info) - ii
                    # print(clid, "current layer")
                    layer = self.layer_info[clid]
                    layer_before = layinfo[clid - 2]
                    layer_current = layinfo[clid - 1]
                    # print(layer_before['opvars'])

                    bridge = {}
                    h0bf = layer_before['height']
                    for kk, item in layer.items():
                        # print(kk, item)
                        for iif, ipf in item['ipf'].items():
                            ipvs = {}
                            for ipvv in ipf['inc']:
                                # print(ipvv['name'], layer_before['opvars'], layer_before['bridge'])
                                # print("bg", ipvv, layer_before['opvars'].keys(), "end")
                                if ipvv["name"] not in layer_before['opvars'].keys() and ipvv['name'] not in layer_before['bridge'].keys():
                                    # print(ipvv, "Not Find in before in layer", clid - 1, len(layer_before["bridge"]))
                                    if not len(layer_before["bridge"]):
                                        h0bf += 1
                                    bridge[ipvv['name']] = {
                                        "nam": ipvv['name'],
                                        "pos": np.array([layerW*(clid - 2),  h0bf + 0.1])
                                    }
                                    layer_before['bridge'].update(bridge)
                                    ipvs[ipvv['name']] = {
                                        "nam": ipvv['name'],
                                        "pos": np.array([layerW*(clid - 2),  h0bf + 0.1])
                                    }
                                    h0bf += 0.2 
                                elif ipvv["name"] in layer_before['opvars'].keys(): 
                                    ipvs[ipvv['name']] = layer_before['opvars'][ipvv['name']]
                                elif ipvv['name'] in layer_before['bridge'].keys():
                                    # print(ipvv, "find in before", layer_before['bridge'])
                                    ipvs[ipvv['name']] = layer_before['bridge'][ipvv['name']]
                            
                            layer_current['ipvars'][iif] = [ipvs, ipf['pos']] 
                            layer_before['height'] = h0bf 
                
                for ii in range(nlay - 1):
                    clid = ii + 1 
                    layer = self.layer_info[clid]
                    layer_next = layinfo[clid]
                    if clid != 1:
                        for kk, mod in layer.items():
                            for off, opf in mod['opf'].items():
                                for op in opf['inc']:
                                    if op in layer_next['bridge'].keys():
                                        mod['opv'][op]['brg'] = layer_next['bridge'][op]['pos']
                    else: 
                        mod = layer['Parameters']
                        for op, opv in mod['opv'].items():
                            if op in layer_next['bridge'].keys():
                                opv['brg'] = layer_next['bridge'][op]['pos']
                
                maxh = max([layinfo[ii]['height'] for ii in range(nlay)]) + 0.4 
                figH = maxh
                for ii in range(nlay):
                    # print(layinfo[ii]['height'], layinfo[ii]['y0'], layinfo[ii]['bridge'])
                    lid =  ii + 1
                    layer = self.layer_info[lid]
                    for mm, mod in layer.items():
                        if ii != 0:
                            arrbp = mod['bp'] + np.array([0., 0.5 * (maxh - layinfo[ii]['height'])])
                        else: 
                            arrbp = np.array([0., 0.5 * (maxh - layinfo[ii]['height']) + mod['bp'][1]])
                        # update_module_pos(mod, arrbp)
                        mod['bp'] = mod['bp'] + np.array([0., 0.5 * (maxh - layinfo[ii]['height'])])
                        for oof, opf in mod['opf'].items():
                            opf['pos'] += arrbp
                        for oov, opv in mod['opv'].items():
                            opv['pos'] += arrbp 
                        for iif, ipf in mod['ipf'].items():
                            ipf['pos'] += arrbp
                for ii in range(nlay):
                    layer = self.layer_info[ii + 1]
                    layer_current = layinfo[ii]
                    for item in layer.values():
                        for kk, ipf in item['ipf'].items():
                            # pF = layer_current["ipvars"][kk][1]
                            # for op in layer_current["ipvars"][kk][0].values():
                            #     arraw_from_V2F(op['pos'], pF )
                            # print(kk, layer_current["ipvars"][kk])
                            ipf['link'] = layer_current["ipvars"][kk][0]
                            # break



                return layinfo, figH

            def update_module_pos(mod, arbp):
                from pprint import pprint
                pprint(mod)

            def draw_layer_module(klayer, mod):
                oh = np.array([0., + 0.5 * mod['ohh']])
                # ax.plot([(klayer-1)*layerW, (klayer-1)*layerW], [-100, 100], "-", c='grey', lw=0.6, alpha=0.2)
                if klayer != 1:
                    # print(mod)
                    module_at_pos(mod['bp'])
                    ax.text(mod['bp'][0], mod['bp'][1] - 0.4, mod['name'], ha="center", va='top', fontfamily="sans-serif", fontsize="medium", fontstyle="normal", fontweight="bold")
                    for kk, ipf in mod['ipf'].items():
                        # print(kk, ipf.keys())
                        input_file_at_pos(ipf['pos'] )
                        ax.text(
                            ipf['pos'][0] - 0.25, ipf['pos'][1] - 0.27, ipf['fil'], 
                            ha="center", va='top', 
                            fontfamily="sans-serif", fontsize="x-small", fontstyle="normal", fontweight="light"
                        )
                        line_from_F2M(ipf['pos'] , mod['bp'])
                        for link in ipf['link'].values():
                            # print(kk, link)
                            arraw_from_V2F(link['pos'], ipf['pos'])
                    for op in mod['opv']:
                        mod['opv'][op]['pos'] = mod['opv'][op]['pos'] 
                        mod['opv'][op] = outvar_at_OP(mod['opv'][op])
                    # ax.plot([bs[0] + 2.6, bs[0] + 2.6], [-oh[1] - 0.2, oh[1] -0.2], "-", c='grey', lw=0.6, alpha=0.2)
                    for kk, opf in mod['opf'].items():
                        # print(kk, opf.keys())
                        output_file_at_pos(opf['pos'])
                        ax.text(
                            opf['pos'][0] + 0.25, opf['pos'][1] - 0.27, opf['fil'], 
                            ha="center", va='top', 
                            fontfamily="sans-serif", fontsize="x-small", fontstyle="normal", fontweight="light"
                        )
                        
                        line_from_M2F(opf['pos'], mod['bp'])

                        for op in opf['inc']:
                            opp = mod['opv'][op]
                            line_from_F2V(opf['pos'], opp)
                else:
                    sampler_at_pos(mod['bp'])
                    ax.text(mod['bp'][0], mod['bp'][1] - 0.4, mod['name'], ha="center", va='top', fontfamily="sans-serif", fontsize="medium", fontstyle="normal", fontweight="bold")

                    for op in mod['opv']:
                        mod['opv'][op] = outvar_at_OP(mod['opv'][op])
                        line_from_F2V(mod['bp'], mod['opv'][op])
                        if "brg" in mod['opv'][op].keys():
                            line_from_V2B(mod['opv'][op]['pos'], mod['opv'][op]['brg'])

            def arraw_from_V2F(pA, pB):
                tt = np.linspace(-0.5 * np.pi, 0.4 * np.pi, 100)
                xx = np.linspace(pA[0] + 0.14, pB[0] - 0.4, 100)
                yy = pA[1] + (np.sin(np.sin(tt)) + 0.8414709848078965)/1.6555005931675704 * (pB[1]-pA[1])
                points = np.array([xx, yy]).T.reshape(-1, 1, 2)
                segments = np.concatenate([points[:-1], points[1:]], axis=1)
                cmap = plt.get_cmap('coolwarm')  
                norm = plt.Normalize(-1.5, 1.5)
                colors = cmap(norm(tt[:-1]))             
                lc = LineCollection(segments, colors=colors, linewidth=2)
                ax.add_collection(lc)
                ax.plot(pA[0] + 0.1, pA[1], 's', c="gray", markersize=4)
                ax.plot(pA[0] + 0.14, pA[1], 's', c="darkgray", markersize=6)

            def line_from_V2B(pV, pB):
                tt = np.linspace(-0.15 * np.pi, 0.5 * np.pi, 100)
                xx = np.linspace(pV[0] + 0.16, pB[0], 100)
                yy = pV[1] + (np.sin(np.sin(np.sin(np.sin(np.sin(tt))))) + 0.4004294116620898)/1.028001243711249 * (pB[1]-pV[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=2)

            def line_from_V2M(pA, pB):
                tt = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 100)
                xx = np.linspace(pA[0], pB[0], 100)
                yy = pA[1] + ((np.sin(np.sin(tt)) / 2 / 0.8414709848078965) + 0.5)  * (pB[1]-pA[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=2)

            def line_from_M2F(pF, pM):
                tt = np.linspace(0., 0.5*np.pi, 100)
                xx = np.linspace(0.3 + pM[0], 0.1 + pF[0], 100 )
                yy = pM[1] - 0.2 + (np.sin(np.sin(tt))) / 0.8414709848078965 * (pF[1] + 0.2 - pM[1])
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=3)

            def line_from_F2V(pF, op):
                tt = np.linspace(-0.5*np.pi, 0.5*np.pi, 100)
                xx = np.linspace(pF[0] + 0.3, op['pos'][0] - op['wid'] - 0.1 , 100)
                yy = pF[1] + (np.sin(np.sin(tt)) / 0.8414709848078965  + 1)/2 * (op['pos'][1] - pF[1])  
                ax.plot(xx, yy, "-", c="#3b4dc0", lw=0.8, alpha=0.7)
                ax.plot([op['pos'][0] - op['wid'] -0.1], [op['pos'][1]], "o", markersize=4, c="#3b4dc0", alpha=0.7)

            def line_from_F2M(pF, pM):
                tt = np.linspace(0., 0.5*np.pi, 100)
                xx = np.linspace(-0.3 + pM[0], -0.1 + pF[0] , 100 )
                yy = pM[1] + 0.2 + (np.sin(np.sin(tt))) / 0.8414709848078965 * (pF[1] - 0.2 - pM[1])
                ax.plot(xx, yy, "-", c="#d45040", lw=3)

            def input_file_at_pos(pos):
                image_path = "src/icons/inputfile.png" 
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.5, pos[0], pos[1]-0.25, pos[1]+0.25], zorder=100)

            def output_file_at_pos(pos):
                image_path = "src/icons/outputfile.png"  
                with Image.open(image_path) as image: 
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0], pos[0]+0.5, pos[1]-0.25, pos[1]+0.25], zorder=100)

            def outvar_at_OP(op):
                text = ax.text(op["pos"][0], op["pos"][1], op['nam'], ha="right", va="center", fontfamily="monospace", variant="small-caps", fontsize="x-small", fontstyle="normal", fontweight="bold")
                bbox = text.get_window_extent(renderer=fig.canvas.get_renderer())
                bbox_data = inv.transform([(bbox.x0, bbox.y0), (bbox.x1, bbox.y1)])
                op['wid'] = bbox_data[1][0] - bbox_data[0][0]
                return op 
                # ax.plot([op["pos"][0] - op['wid'], op["pos"][0]], [op["pos"][1], op["pos"][1]], "-")

            def module_at_pos(pos):
                image_path = "src/icons/calculator.png"  
                with Image.open(image_path) as image:
                    image = np.array(image.convert("RGBA"))  # Explicitly preserve transparency
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def lib_at_pos(pos):
                image_path = "src/icons/library.png"  
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def sampler_at_pos(pos):
                image_path = "src/icons/sampler.png"  
                with Image.open(image_path) as image:
                    image = np.array(image)
                    ax.imshow(image, extent=[pos[0]-0.45, pos[0]+0.45, pos[1]-0.45, pos[1]+0.45], zorder=100)

            def logo_at_pos(pos):
                image_path = "src/icons/JarvisHEP.png"
                with Image.open(image_path) as image:
                    image = np.array(image.convert("RGBA"))
                    ax.imshow(image, extent=[pos[0]-0.25, pos[0]+0.25, pos[1]-0.25, pos[1]+0.25], zorder=100)
                    ax.text(pos[0]+0.27, pos[1]+0.15, "Jarvis-HEP", ha="left", va='top', color="#0F66C3", fontfamily="sans-serif", fontsize="small", fontstyle="normal", fontweight="bold")


            layerInfo, figH = resolve_layer()
            if figH < 1.8: 
                figH = 1.8
            # figH += 
            fig = plt.figure(figsize=(figL, figH))

            ax  = fig.add_axes([0., 0., 1., 1.])
            ax.axis("off")
            logo_at_pos([0.7, figH - 0.35])
            ax.set_xlim([0, figL])
            ax.set_ylim([0, figH])
            inv = ax.transData.inverted()

            # ax.plot([-100, 100], [0., 0.], "-", c='grey', lw=0.4, alpha=0.2)
            # ax.plot([-100, 100], [-0.2, -0.2], "-", c='grey', lw=0.4, alpha=0.2)
            # ax.plot([-100, 100], [0.2, 0.2], "-", c='grey', lw=0.4, alpha=0.2)    
            outframe = create_round_square([2, 4,])
            x_coords = [point[0] for point in outframe.exterior.coords]
            y_coords = [point[1] for point in outframe.exterior.coords]  

            for kk, item in self.layer_info.items():
                for mod, info in item.items():
                    draw_layer_module(kk, info)
            
            plt.savefig(save_path, dpi=300)
            plt.close("all")
            await asyncio.sleep(2)