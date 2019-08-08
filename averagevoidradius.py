import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import h5py
import nbodykit.lab as nlab
import nbodykit.cosmology as ncosmo
from nbodykit import CurrentMPIComm, transform
import dask.array as da
from nbodykit.source.catalog import HaloCatalog, HDFCatalog

zlab = '1'
redshift = '1.00'
sim = 'Quijote'
#w = 'wiggles'
wiggles = ['wiggles','no-wiggles']
d = '_'
realization =1
scale = 'lin'
binlab = '_dr1'
#num_realizations = 1
c = '_halos'

Nmesh=256

BoxSize = 2000.0

plotarr={}
xi = {}
averagevoidradius = {}


for w in wiggles:
    averagevoidradius[w] = 0
    number = 0
    
    filenw = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_' +w + '_'+str(realization)+'_halos_z1.00/centers_central_Quijote_' + w+'_'+ str(realization) +'_halos_z1.00.out'
            
    gadget_format2 = pd.read_csv(filenw,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')

    gadget_format = gadget_format2[gadget_format2[4]>40]
    
    for i in range(0, len(gadget_format)-1):
        
        if i in gadget_format[4]:
            averagevoidradius[w] += gadget_format[4][i]
            number += 1

    print(averagevoidradius[w])

    averagevoidradius[w] = averagevoidradius[w]/(number)
        
    print(w +': ' + str(averagevoidradius[w]))
    
