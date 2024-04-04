#!/usr/bin/env python3 

from IOs.IOs import YodaFile
yd = YodaFile()
yd.file = "/home/buding/Workshop/CEPC/Signal/output/test03.dat"    
yd.vinf = [
    {"expr": "sr01", "type": "hist1d", "name": "/CEPC_SLEPTON/SR_low_DeltaM", "code": ["1.000000e+00", "2.000000e+00"]},
    {"expr": "sr02", "type": "hist2d", "name": "/CEPC_SLEPTON/hist2D_mll_cThetalp", "code": ["1.000000e+01", "2.000000e+01", "0.000000e+00", "1.000000e-01"]}
    ]     
yd.read()
print(yd.vars)