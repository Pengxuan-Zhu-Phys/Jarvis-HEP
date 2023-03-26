#!/usr/bin/env python3 

import time
import ipyparallel as ipp 

task_durations = [1] * 25
# print(task_durations)

cl = ipp.Cluster(n=5)
rc = cl.start_and_connect_sync()
nprocs = len(rc.ids)
print(rc.ids)
dview = rc[:]
dview.use_dill()



# print("{:.2f} -> starting tasks".format(time.time()))
# with ipp.Cluster(n=5) as rc:
#     view = rc.load_balanced_view()
#     print("{:.2f} -> get a view on the cluster".format(time.time()))
#     asyncresult = view.map_async(time.sleep, task_durations)
#     print("{:.2f} -> submit the tasks".format(time.time()))
#     asyncresult.wait_interactive()
#     print("{:.2f} -> wait interactively for results".format(time.time()))
#     result = asyncresult.get()
#     print("{:.2f} -> retrieve active results".format(time.time()))