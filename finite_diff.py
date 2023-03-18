#!/usr/bin/env python3 

import numpy as np 
import time
import pandas as pd 
import matplotlib.pyplot as plt 

# Variables (x, y) is defined as (B, T) in your case
def Entropy(x, y):
    return np.sqrt(x**2 + y**2)

def GammaB(xx, yy, ff):
    dfx = np.full(ff.shape, np.NAN)
    dfy = np.full(ff.shape, np.NAN)
    dfxx = np.full(ff.shape, np.NAN)
    dfyy = np.full(ff.shape, np.NAN)
    dfxy = np.full(ff.shape, np.NAN)

    GammaB = np.full(ff.shape, np.NAN)
    for ii in range(ff.shape[0]):
        for jj in range(ff.shape[1]):
            if ii > 0 and ii < ff.shape[0] - 1:
                dfx[ii, jj] = (ff[ii+1, jj] - ff[ii-1, jj]) / (xx[ii+1] - xx[ii-1])
                dfxx[ii, jj] = 4.0 * (ff[ii+1, jj] - 2.0*ff[ii, jj] + ff[ii-1, jj]) / (xx[ii+1] - xx[ii-1])**2 
            if jj > 0 and jj < ff.shape[1] - 1:
                dfy[ii, jj] = (ff[ii, jj+1] - ff[ii, jj-1]) / (yy[jj+1] - yy[jj-1])
                dfyy[ii, jj] = 4.0 * (ff[ii, jj+1] - 2.0*ff[ii, jj] - ff[ii, jj-1]) / (yy[jj+1] - yy[jj-1])**2 
            if ii > 0 and ii < ff.shape[0] - 1 and jj > 0 and jj < ff.shape[1] - 1:
                dfxy[ii, jj] = (ff[ii+1, jj+1] - ff[ii-1, jj+1] - ff[ii+1, jj-1] + ff[ii-1, jj-1]) / (xx[ii+1] - xx[ii-1] + yy[jj+1] - yy[jj-1])
        
    for ii in range(ff.shape[0]):
        for jj in range(ff.shape[1]):
            GammaB[ii, jj] = - dfx[ii,jj] / (yy[jj] * dfy[ii, jj])
    return dfx, dfy, dfxx, dfyy, dfxy, GammaB


if __name__ == "__main__":
    # Here I chose the dx and the dy are constant 
    xx = np.linspace(2.0, 3.0, 1000)
    yy = np.linspace(3.0, 5.0, 6)
    
    # xx = 2.0 + np.random.random(1000)
    # yy = 3.0 + 2.0 * np.random.random(6)
    xv, yv = np.meshgrid(xx, yy)
    
    t0 = time.time()
    print(xx.shape, yy.shape, Entropy(xv.T, yv.T).shape)
    
    dfx, dfy, dfxx, dfyy, dfxy, gammb = GammaB(xx, yy, Entropy(xv.T, yv.T))
    # print(gammb[:,2])
    plt.plot(xx, dfx[:,2], "o")
    plt.show()
    
    
    # df = pd.DataFrame(gammb)
    # df.to_csv("GammaB.csv", index=False)
    print("dfx\n", dfx)
    # print(gammb)
    # print(time.time() - t0)
