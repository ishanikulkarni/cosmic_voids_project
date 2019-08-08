import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import h5py
import nbodykit.cosmology as ncosmo
from nbodykit import CurrentMPIComm, transform
import dask.array as da
from nbodykit.source.catalog import HaloCatalog, HDFCatalog

zlab = '1'
redshift = '1.00'
#zlabs = ['0','1']
#redshifts = ['0.00','1.00']
zlabs = ['1']
redshifts = ['1.00']
sim = 'Quijote'
w = 'wiggles'
wiggles = ['wiggles','no-wiggles']
d = '_'
#realization =0
scale = 'lin'
binlab = '_dr1'
num_realizations = 100
c = '_halos'# #'_cdm'
cola = '' #'COLA/'
crazy = ''

Nmesh=256

BoxSize = 2000.0


#----------------------------------------------CREATING VOID CUT CATALOGS-----------------------------------------------------------------------------------------------------
for (zlab,redshift) in zip(zlabs,redshifts):
    
    plotarr2 = {}

    for w in wiggles:
        for realization in range(0,num_realizations):
                        
            file2='/projects/SPERGEL/COLA_runs/voidcuts/z' + zlab + sim+ '/' +sim +'_'+ w + '_regular_large_rcut40_' + str(realization)+ '_voids_drdefault.dat'
            
            if not os.path.exists(file2):
                file = '/projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+zlab+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'/'+crazy+'centers_central_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'.out'

                gadget_format2 = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')

                heading = '# center x,y,z (Mpc/h), volume (normalized), radius (Mpc/h), redshift, volume (Mpc/h^3), void ID, density contrast, num part, parent ID, tree level, number of children, central density'

                gadget_format = gadget_format2[gadget_format2[4]>40]
 
                print "read in data"
        
                if not os.path.exists(os.path.dirname(file2)):
                    os.makedirs(os.path.dirname(file2))

                f = open(file2,'w+')
     
                np.savetxt(f, gadget_format, header = heading)
                f.close()
