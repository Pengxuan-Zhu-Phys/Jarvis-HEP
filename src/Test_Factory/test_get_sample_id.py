#!/usr/bin/env python3 

import numpy as np 
import time
import pandas as pd 
import os, sys

pwd = sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from Func_lib import get_sample_id

# t0 = time.time()
# dat = np.random.random_sample((1000000, 4))
# # print(dat)
# print(time.time() - t0)
# df = pd.DataFrame(dat, columns=["cube0", "cube1", "cube2", "cube3"])
# print(df)
# print(time.time() - t0)
# points = []
# for ids, item in df.iterrows():
#     raw = dict(item)
#     raw['ID'] = get_sample_id()
#     if ids % 10000 == 0:
#         print("Current time {:.3f} at {}".format(time.time() - t0, ids))
#     points.append(raw)


# points = pd.DataFrame(points)
# points.to_csv("test100w01.csv", index=False)

# print("Finishing at {:.3f}".format(time.time() - t0))


t0 = time.time()
dat = np.random.random_sample((1000000, 4))
# print(dat)
print("Sampling finished at {:.3f} second".format(time.time() - t0))
df = pd.DataFrame(dat, columns=["cube0", "cube1", "cube2", "cube3"])
from uuid import uuid4
df['ID'] = df.apply(lambda _: uuid4())
df['Status'] = "Free"
df['X0'] = pd.eval("3.0 + 100.0 * df['cube0']", target=df)
df['X1'] = pd.eval("10 ** (0.01 + 0.99 * df['cube1'])", target=df)
df['X2'] = pd.eval("10 ** (0.01 + 0.99 * df['cube2'])", target=df)
df['X3'] = pd.eval("10 ** (0.01 + 0.99 * df['cube3'])", target=df)
df = df.loc[pd.eval("df['X0'] > 2.0 * df['X3']")]
print("Evaluation at {:.3f}".format(time.time() - t0))
df.to_csv("test100w02.csv", index=False)

print("Finishing at {:.3f}".format(time.time() - t0))
