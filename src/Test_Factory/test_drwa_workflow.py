#!/usr/bin/env python3 

import networkx as nx
import matplotlib.pyplot as plt

class Workflow:
    def __init__(self):
        self.graph = nx.DiGraph()

    def add_module(self, module_name, required_modules=None):
        self.graph.add_node(module_name)
        if required_modules:
            for rm in required_modules:
                self.graph.add_edge(rm, module_name)

    def draw_workflow(self):
        plt.figure(figsize=(10, 8))
        pos = nx.spring_layout(self.graph)
        pos = {
    'Delphes': (0, 0),
    'HepMC': (1, 0),
    'Pythia8': (2, 0),
    'SUSYHIT': (0.5, 1),
    'XSect_SmuonR': (1.5, 1)
        }


        nx.draw(self.graph, pos, with_labels=True, node_color='skyblue', node_size=2000, edge_color='k', linewidths=1, font_size=15, arrows=True)
        plt.show()

# 创建Workflow实例
workflow = Workflow()

# 添加模块和依赖关系
workflow.add_module('Delphes', required_modules=[])
workflow.add_module('HepMC', required_modules=[])
workflow.add_module('Pythia8', required_modules=['HepMC'])
workflow.add_module('SUSYHIT', required_modules=['Delphes'])
workflow.add_module('XSect_SmuonR', required_modules=['SUSYHIT'])

# 绘制工作流逻辑关系图
workflow.draw_workflow()
