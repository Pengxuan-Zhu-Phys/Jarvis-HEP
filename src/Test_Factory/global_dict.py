#!/usr/bin/env python3 

from multiprocessing import Manager
import multiprocessing
from random import random
import time

def func(ipdct, upt):
    print("appling func {} <- {} at {} ".format(ipdct, upt, time.time()))
    time.sleep(2.0 * random())
    
    ipdct.update(upt)
    print("closing {} at {}".format(ipdct, time.time()))


if __name__ == "__main__":
    glb_dict1 = Manager().dict()
    glb_dict1['a'] = 150
    print(glb_dict1)

    pool = multiprocessing.Pool(processes=3)
    pool.apply_async(func=func, args=(glb_dict1, {"b": 300}))
    pool.apply_async(func=func, args=(glb_dict1, {"a": 200}))
    pool.apply_async(func=func, args=(glb_dict1, {"b": 500}))
    print(glb_dict1)
    pool.close()
    pool.join()
    time.sleep(3)


    print(glb_dict1)