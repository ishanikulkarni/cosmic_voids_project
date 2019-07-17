import h5py
import nbodykit.lab as nlab
import nbodykit.cosmology as ncosmo
from nbodykit import CurrentMPIComm, transform
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dask.array as da
from nbodykit.source.catalog import HaloCatalog, HDFCatalog
import os

#Add two for loops to go through all the realizations and stuff

densities = ['']
wiggles = ['wiggles','no-wiggles']#['wiggles', 'no-wiggles']

sim = 'Quijote'
#redshift = '0'
# Should match redshift above; stands for z label in pathname
#zlab = '1.00'
zlabs = ['0','1']
redshifts = ['0.00','1.00']
scale = 'lin'
binlab = '_dr1'
c = '_halos'# #'_cdm'
cola = '' #'COLA/'
num_realizations = 1
crazy = '' #'untrimmed_dencut_'
#filewrites = ['/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat', '/tigress/isk/COLA_runs/trial_folder'+'z'+zlab+sim+'samplestrial'+w+d+'_nbodykit_' + sim + '_halos_' + str(realization)+'.dat', '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat' , '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat']
#filereads = ['/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/Quijote_'+w+d+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/random_sampling/z'+zlab+sim +'/'+sim+'_'+w+'_samples_'+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat']
Nmesh=256

BoxSize = 2000.0

d='_'
#realization=0
os.chdir('/')

### This section is for halos ###

# halos
plotarr2 = {}
for (redshift,zlab) in zip(redshifts,zlabs):
    for w in wiggles:
        for realization in range(0,num_realizations):
            
            #Original data:
            file4 = '/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat'
            file3 ='/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/Quijote_'+w+d+str(realization)+'_halos.dat'
            
            # Random sampling data:
            #file4 = '/tigress/isk/COLA_runs/trial_folder'+'z'+zlab+sim+'samplestrial'+w+d+'_nbodykit_' + sim + '_halos_' + str(realization)+'.dat'
            #file3 = '/projects/SPERGEL/COLA_runs/random_sampling/z'+zlab+sim +'/'+sim+'_'+w+'_samples_'+str(realization)+'_halos.dat'
            
            #Mass cuts larger than 10^13.5:
            #file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat'
            #file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat'
                
            #Mass cuts larger than 10^13:
            #file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_large_mcut1e13_'+str(realization)+'_halos.dat'
            #file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_large_mcut1e13_'+str(realization)+'_halos.dat'
                
            #Mass cuts smaller than 10^13.5:
            #file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat'
            #file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat'
                
            #Mass cuts smaller than 10^13:
            #file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_small_mcut1e13_'+str(realization)+'_halos.dat'
            #file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13_'+str(realization)+'_halos.dat'

            if not os.path.exists(file4):
                gadget_format = pd.read_csv(file3,skiprows=5, delim_whitespace=True,engine='python', header=None, index_col=False, comment='#')
                
                print "read in data"
                
                data = {}

                data['Position'] = da.from_array(np.column_stack((gadget_format[1], gadget_format[2], gadget_format[3])), chunks=(100,3))
                data['Velocity'] = da.from_array(np.column_stack((gadget_format[4], gadget_format[5], gadget_format[6])), chunks=(100,3)) #transform.StackColumns
                    
                data = nlab.ArrayCatalog(data)
                
                print "stacked data"

                # Convert to MeshSource object, BoxSize in Mpc/h
                mesh = data.to_mesh(window='tsc', Nmesh=Nmesh, BoxSize = BoxSize, compensated=True, position='Position')

                # Halos
                # Do correlation function
                #if not (os.path.exists(file4)):
                t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize, dr = 1.0)

                gg = t.corr
                Xigg = gg['corr']
                Xir = gg['r']

                if not os.path.exists(os.path.dirname(file4)):
                    os.makedirs(os.path.dirname(file4))
                f = open(file4, 'w+')
                f.write('# r Xi(r) Xi(r)-ShotNoise\n')
                f.write('# z = '+zlab +'\n')
                        
                datacc = np.array([gg['r'],gg['corr'].real,gg['corr'].real - gg.attrs['shotnoise']])
                datacc = datacc.T
                        
                np.savetxt(f, datacc)
                f.close()

                # print out the meta-data
                for k in gg.attrs:
                    print("%s = %s" %(k, str(gg.attrs[k])))
            else:
                data_array =  pd.read_csv(file4,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
                gg={}
                # print(data_array)
                gg['r'] = data_array[0]
                gg['corr']=data_array[1]
                
                #gg[w+'avgcorr'] = np.zeros(len(gg[w+'corr']))
                #gg[w+'avgcorr'] += gg[w+'corr']
                
            if w+'avgcorr' not in plotarr2.keys():
                plotarr2[w+'avgcorr'] = np.zeros(len(gg['r']))
                plotarr2[w+'r'] = gg['r']
                plotarr2[w+'corr'] = gg['corr']
            try:
                plotarr2[w+'r'] == gg['r']
            except ValueError:
                print('r values do not match previous r values')
            #print(gg.keys())
                            
            plotarr2[w+'avgcorr'] += gg['corr']
            #   plt.figure()
            #print(plotarr2['wigglesr'])
            #print(plotarr2['wigglesavgcorr'])
            #print(plotarr2['no-wigglesavgcorr'])
            
            #plt.figure()
            #plt.semilogx(plotarr2['wigglesr'], plotarr2['wigglescorr'], label = 'Wiggle, '+ redshift)
            #plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
            #plt.legend()
            #plt.show()
            
            #plt.figure()
            #plt.semilogx(plotarr2['no-wigglesr'],plotarr2['no-wigglescorr'], label = 'No Wiggles, '+redshift)
            #plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
            #plt.legend()
            #plt.show()
    plt.figure()
    plt.semilogx(plotarr2['wigglesr'],plotarr2['wigglesr']**2.0*(plotarr2['wigglesavgcorr'] - plotarr2['no-wigglesavgcorr'])/num_realizations, label = redshift)
    plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    plt.legend()
plt.show()
#plt.plot(datacc[0],datacc[1])
#plt.show()

