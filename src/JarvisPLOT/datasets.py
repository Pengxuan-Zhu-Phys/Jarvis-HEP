# datasets.py
import pandas as pd, numpy as np, os

class DataEnv:
    def __init__(self, specs, yaml_dir):
        self.frames = {}
        self.yaml_dir = yaml_dir
        if specs:
            for ds in specs:
                path = ds.path if os.path.isabs(ds.path) else os.path.abspath(os.path.join(yaml_dir, ds.path))
                if ds.type == 'csv':
                    self.frames[ds.name] = pd.read_csv(path)
                # (extend: parquet/hdf5/etc.)

    def namespace(self):
        ns = {'np': np, 'pd': pd}
        ns.update(self.frames)
        return ns

    def eval_expr(self, expr: str):
        return eval(expr, self.namespace(), {})