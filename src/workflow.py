#!/usr/bin/env python3

import networkx as nx 
from base import Base
from module import Module
import numpy as np 

class Workflow(Base):
    def __init__(self):
        self.modules = {}
        self.module_states = {}
        self.module_dependencies = {}
        self.parameter_module = None
        self.library_modules = {}
        self.graph = {}
        self.layer_info = {}

    def add_module(self, module):
        self.modules[module.name] = module
        self.module_states[module.name] = 'waiting'

    def add_library_module(self, module):
        self.library_modules[module.name] = module
        self.modules[module.name] = module

    def set_modules(self, modules):
        from module import Parameters, LibraryModule, CalculatorModule
        
        parameter_module = Parameters("Parameters", modules['Parameter'])
        parameter_module.analyze_ios()
        self.add_module(parameter_module)
        self.parameter_module = parameter_module
        
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

        print(unresolved_lib, unresolved)
        # resolve the library layers 
        current_layer = 1
        resolved = []
        while unresolved_lib:
            new_resolved = []
            for ii in range(len(unresolved_lib)):
                print(unresolved_lib[ii], self.modules[unresolved_lib[ii]].required_modules)
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
                print(unresolved[ii], self.modules[unresolved[ii]].required_modules)
                if current_layer == 1:
                    if not self.modules[unresolved[ii]].required_modules:
                        new_resolved.append(self.modules[unresolved[ii]].name)
                else:
                    if set(self.modules[unresolved[ii]].required_modules).issubset(set(resolved)):
                        new_resolved.append(self.modules[unresolved[ii]].name) 
            
            resolved.extend(new_resolved)
            unresolved = [item for item in unresolved if item not in resolved]
            self.calc_layer[current_layer] = new_resolved
            current_layer += 1 
        
        from pprint import pprint
        pprint(self.libr_layer)
        pprint(self.calc_layer)

    def resolve_dependencies(self):
        # Resolve dependencies based on the IOs 
        output_to_module = {} 
        for module in self.modules.values():
            print(module.outputs)
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

    def draw_flowchart(self):
        import matplotlib.pyplot as plt 
        import shapely as sp 
        import logging
        logging.getLogger('matplotlib').setLevel(logging.CRITICAL)
        from matplotlib.collections import LineCollection
        from plot import draw_logo_in_square
        from plot import create_round_square
        from PIL import Image
        
        fig = plt.figure(figsize=(10, 8))
        ax  = fig.add_axes([0., 0., 1., 1.])
        axlogo = fig.add_axes([0.02, 0.875, 0.08, 0.1])
        draw_logo_in_square(axlogo)
        ax.axis("off")
        
        def arraw_from_V2F(pA, pB, ax):
            tt = np.linspace(-0.5 * np.pi, 0.4 * np.pi, 100)
            xx = np.linspace(pA[0], pB[0], 100)
            yy = pA[1] + (np.sin(np.sin(tt)) + 0.8414709848078965)/1.6555005931675704 * (pB[1]-pA[1])
            points = np.array([xx, yy]).T.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            cmap = plt.get_cmap('coolwarm')  
            norm = plt.Normalize(1.5, 5)
            colors = cmap(norm(xx[:-1]))             
            lc = LineCollection(segments, colors=colors, linewidth=2)
            ax.add_collection(lc)
            ax.plot(pA[0]-0.04, pA[1], 's', c="gray", markersize=4)
            ax.plot(pA[0], pA[1], 's', c="darkgray", markersize=6)
        
        def line_from_V2M(pA, pB, ax):
            tt = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 100)
            xx = np.linspace(pA[0], pB[0], 100)
            yy = pA[1] + ((np.sin(np.sin(tt)) / 2 / 0.8414709848078965) + 0.5)  * (pB[1]-pA[1])
            
            ax.plot(xx, yy, "-", c="#3b4dc0", lw=2)
        
        def input_file_at_pos(pos, ax):
            image_path = "src/icons/inputfile.png"  
            image = Image.open(image_path)
            image = np.array(image)
            ax.imshow(image, extent=[pos[0]-0.5, pos[0], pos[1]-0.25, pos[1]+0.25])

        
        def module_at_pos(pos, ax):
            image_path = "src/icons/calculator.png"  
            image = Image.open(image_path)
            image = np.array(image)
            ax.imshow(image, extent=[pos[0]-0.65, pos[0]+0.65, pos[1]-0.65, pos[1]+0.65])
        
        outframe = create_round_square([2, 4,])
        x_coords = [point[0] for point in outframe.exterior.coords]
        y_coords = [point[1] for point in outframe.exterior.coords]  
        # ax.plot(x_coords, y_coords, '-', color='grey', lw=0.3)
        # ax.fill(x_coords, y_coords, fc='#072348', ec=None)
        ax.set_xlim([0, 10])
        ax.set_ylim([0, 8])
        arraw_from_V2F([1., 1.], [4., 3], ax)
        arraw_from_V2F([1., 1.3], [4., 3], ax)
        line_from_V2M([1., 1.3], [4.5, 7], ax)
        input_file_at_pos([4.5, 3], ax)
        module_at_pos([5.4, 2.7], ax)
        


        

        
              
        # plt.show()
        plt.savefig("flowchart.png", dpi=300)