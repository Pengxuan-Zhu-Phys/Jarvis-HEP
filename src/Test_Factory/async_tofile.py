#!/usr/bin/env python3 

import os, sys 

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from IOs.IOs import to_file_woblk

data = {'a': 5, "b": 3}
pt = "./test_woblk.json"

from pandas import Series
ds = Series(data)

# to_file_woblk(data=data, pth=pt, method="to_json")
to_file_woblk(data=ds, pth=pt, method="to_csv")
