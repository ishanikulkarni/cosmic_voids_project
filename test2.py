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
w = 'wiggles'
wiggles = ['wiggles','no-wiggles']
d = '_'
#realization =0
#scale = 'lin'
#binlab = '_dr1'
num_realizations = 100

Nmesh=256

BoxSize = 2000.0

plotarr={}
#fig, axs = plt.subplots(2, 1, sharex=True)

for w in wiggles:
    for realization in range(0,num_realizations):
        filewrite = '/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/' + w + '_Quijote_' + str(realization)+'_voids.dat'

        if not os.path.exists(filewrite):
            filenw = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_' +w + '_'+str(realization)+'_halos_z1.00/centers_central_Quijote_' + w+'_'+ str(realization) +'_halos_z1.00.out'
            #filenw = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_no-wiggles_0_halos_z1.00/centers_central_Quijote_no-wiggles_0_halos_z1.00.out'
            gadget_format = pd.read_csv(filenw,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
            
            print "read in data"

            data = {}
            
            data['Position'] = da.from_array(np.column_stack((gadget_format[0], gadget_format[1], gadget_format[2])), chunks=(100,3))
            data['Velocity'] = da.from_array(np.column_stack((np.zeros(len(gadget_format[0])), np.zeros(len(gadget_format[0])), np.zeros(len(gadget_format[0])))), chunks=(100,3)) #transform.StackColumns
            
            data = nlab.ArrayCatalog(data)
            
            print "stacked data"

            # Convert to MeshSource object, BoxSize in Mpc/h
            mesh = data.to_mesh(window='tsc', Nmesh=Nmesh, BoxSize = BoxSize, compensated=True, position='Position')
            print "converted data"

            # Get void correlation function from the simulation
            # Do the Fourier transform
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize, dr = 1.0)
            
            gg = t.corr
            Xigg = gg['corr']
            Xir = gg['r']

            # Voids
            # Do correlation function
            os.chdir('/')
	
            #print(type(realization))
            #if not (os.path.exists('/tigress/isk/Documents/research/projects/COLA_runs/plot_data/'+cola+'z'+redshift+sim+'/Xivv_'+crazy+w+d+'_nbodykit_'+sim+realization+c+'_'+scale+binlab+'.dat')):
       
            if not os.path.exists(os.path.dirname(filewrite)):
                os.makedirs(os.path.dirname(filewrite))

            f = open(filewrite,'w+')
            f.write('# r Xi(r) Xi(r)-ShotNoise\n')
            f.write('# z = '+zlab +'\n')

            datacc = np.array([gg['r'],gg['corr'].real,gg['corr'].real - gg.attrs['shotnoise']])
            datacc = datacc.T
            
            np.savetxt(f, datacc)
            f.close()

            # print out the meta-data
            for k in gg.attrs:
                print("%s = %s" %(k, str(gg.attrs[k])))

            print(filewrite)

        else:
            print(filewrite)
            data_array = []
            data_array =  pd.read_csv(filewrite,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
            gg={}
            #print(data_array)
            gg['r'] = data_array[0]
            gg['corr']=data_array[1]
            
            #gg[w+'avgcorr'] = np.zeros(len(gg[w+'corr']))
            #gg[w+'avgcorr'] += gg[w+'corr']
        
        if w+'avgcorr' not in plotarr.keys():
            plotarr[w+'avgcorr'] = np.zeros(len(gg['r']))
            plotarr[w+'r'] = gg['r']
            plotarr[w+'corr'] = gg['corr']
        #try:
        #   plotarr[w+'r'] == gg['r']
        #except ValueError:
        #   print('r values do not match previous r values')
        #print(gg.keys())

        print(plotarr.keys())
        plotarr[w+'avgcorr'] += gg['corr']

plt.figure()
plt.semilogx(plotarr['wigglesr'], plotarr['wigglescorr'], label = 'Wiggle, '+ redshift)
plt.legend()
plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')

plt.figure()
plt.semilogx(plotarr['no-wigglesr'], plotarr['no-wigglescorr'], label = 'no-wiggle, '+ redshift)
plt.legend()
plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')

plt.figure() 
plt.semilogx(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr'])/num_realizations)
plt.title('z = '+redshift)
plt.legend()
plt.hlines(y=0,xmin=0,xmax=1000, linestyles='solid')
plt.legend()
plt.show()

'''
ax = axs[1]
ax.plot(plotarr3['wigglesr'],plotarr3['wigglesr']**2.0*(plotarr3['wigglesavgcorr'] - plotarr3['no-wigglesavgcorr'])/num_realizations, label = 'wiggles - no-wiggles')
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.legend(loc = "lower left")
plt.show()
'''
