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
redshifts = ['1']
zlabs = ['1.00']
scale = 'lin'
binlab = '_dr1'
c = '_halos'# #'_cdm'
cola = '' #'COLA/'
num_realizations = 100
crazy = '' #'untrimmed_dencut_'

Nmesh=256

BoxSize = 2000.0

##### This section is for voids #####

# Voids
# Read in void catalog
#redshift='1'
#w='no-wiggles'
d='_'
#realization=0
os.chdir('/')
'''
print(type('projects/SPERGEL/COLA_runs/voidcatalogs/'))
print(type(cola))
print(type(redshift))
print(type(sim))
print(type(w))
print(type(d))
print(type(realization))
print(type(c))
print(type(crazy))
print(type(zlab))
'''

plotarr = {}

for (redshift,zlab) in zip(redshifts,zlabs):
    for w in wiggles:
        for realization in range(0,num_realizations): #range(0,5):

            file2='/tigress/isk/COLA_runs/plot_data/'+cola+'z'+redshift+sim+'/Xivv_'+crazy+w+d+'_nbodykit_'+sim+str(realization)+c+'_'+scale+binlab+'.dat'
       
            if not os.path.exists(file2):
                file = 'projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+redshift+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+str(realization)+c+'_z'+zlab+'/'+crazy+'centers_central_'+sim+'_'+w+d+str(realization)+c+'_z'+zlab+'.out'
                gadget_format = pd.read_csv(file,delim_whitespace=True, header=None,index_col=False,comment='#')
                #gadget_format = pd.read_csv('/projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+redshift+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+realization+c+'_z'+zlab+'/'+crazy+'centers_central_'+sim+'_'+w+d+realization+c+'_z'+zlab+'.out', delim_whitespace=True, header=None, index_col=False, comment='#')

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

                #print(type(cola))
                #print(type(redshift))
                #print(type(sim))
                #print(type(crazy))
                #print(type(w))
                #print(type(d))
                #print(type(realization))
                #print(type(c))
                #print(type(scale))
                #print(type(binlab))
                os.chdir('/')
	
                #print(type(realization))
                #if not (os.path.exists('/tigress/isk/Documents/research/projects/COLA_runs/plot_data/'+cola+'z'+redshift+sim+'/Xivv_'+crazy+w+d+'_nbodykit_'+sim+realization+c+'_'+scale+binlab+'.dat')):
       
                if not os.path.exists(os.path.dirname(file2)):
                    os.makedirs(os.path.dirname(file2))

                    #if not os.path.exists(file2):
                    #   pd.read_csv(file2,delim_whitespace=True, header=None,index_col=False,comment='#')
                    #else:
                f = open(file2,'w+')
                #f = open('/tigress/isk/Documents/research/projects/COLA_runs/plot_data/'+cola+'z'+redshift+sim+'/Xivv_'+crazy+w+d+'_nbodykit_'+sim+str(realization)+c+'_'+scale+binlab+'.dat', 'w+')
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
                data_array =  pd.read_csv(file2,delim_whitespace=True, header=None,index_col=False,comment='#')
                gg={}
                gg['r'] = data_array[0]
                gg['corr']=data_array[1]
            #if str(realization) == '0'
            if w+'avgcorr' not in plotarr.keys():
                plotarr[w+'avgcorr'] = np.zeros(len(gg['r']))
            plotarr[w+redshift+'r'] = gg['r']
            plotarr[w+redshift+'corr'] = gg['corr']
            #try:
            #    plotarr[w+redshift+'r'] == gg['r']
            #    plotarr[w+redshift+'corr'] = gg['corr']
            #except ValueError:
            #    print('r values do not match previous r values')

            plotarr[w+'avgcorr'] += gg['corr']
        #print(plotarr.keys())

    plt.figure()
    plt.semilogx(plotarr['wiggles'+redshift+'r'], plotarr['wiggles'+redshift+'r']**2*plotarr['wigglesavgcorr'], label = 'BAO')
    plt.legend()
    plt.title("Void Correlation Function")
    plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')

    #plt.figure()
    plt.semilogx(plotarr['no-wiggles'+redshift+'r'],plotarr['no-wiggles'+redshift+'r']**2*plotarr['no-wigglesavgcorr'], label = 'BAO-removed')
    plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    plt.title("Void Correlation Function")
    plt.legend()
    plt.xlabel(r"$r$ [$h^{-1} \ \mathrm{Mpc}$]", fontsize = 12)
    plt.ylabel(r"$r^2 \Delta \xi(r)$ [$h^{3}\mathrm{Mpc}^{-3}$]", fontsize = 12)
    plt.xlim(xmax = 250)

    plt.figure() 
    plt.semilogx(plotarr['wiggles'+redshift+'r'],plotarr['wiggles'+redshift+'r']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr'])/num_realizations, label = 'BAO- no-BAO')
    plt.title('Void Correlation Difference Function')
    plt.legend()
    plt.hlines(y=0,xmin=0,xmax=1000, linestyles='solid')
    plt.ylim(ymin = -150, ymax = 150)
    plt.xlim(xmin = 12, xmax = 250)
    plt.xlabel(r"$r$ [$h^{-1} \ \mathrm{Mpc}$]", fontsize = 12)
    plt.ylabel(r"$r^2 \Delta \xi(r)$ [$h^{3}\mathrm{Mpc}^{-3}$] BAO - no-BAO", fontsize = 12)
plt.legend()
plt.show()

#plt.savefig("test.png")
#fig1, ax1 = plt.subplots(figsize=(5,4))
#ax1 = plt.subplot(111)
#ax1.plot(gg['r'],gg['corr'])
#fig1.savefig('test.png')

'''
### This section is for halos ###

# halos
plotarr2 = {}
for (redshift,zlab) in zip(redshifts,zlabs):
    for w in wiggles:
        for realization in range(0,5):

            file4 = '/tigress/isk/COLA_runs/vide_files/plot_data/'+cola+'z'+redshift+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+scale+binlab+'.dat'
            
            if not os.path.exists(file4):
                file3 ='/projects/SPERGEL/COLA_runs/data/z'+redshift+sim+'/Quijote_'+w+d+str(realization)+'_halos.dat'
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
                data_array =  pd.read_csv(file2,delim_whitespace=True, header=None,index_col=False,comment='#')
                gg={}
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
            print(gg.keys())

            plotarr2[w+'avgcorr'] += gg['corr']
            #   plt.figure()
    #print(plotarr2['wigglesr'])
    #print(plotarr2['wigglesavgcorr'])
    #print(plotarr2['no-wigglesavgcorr'])
    
    plt.figure()
    plt.semilogx(plotarr2['wigglesr'], plotarr2['wigglescorr'], label = 'Wiggle, '+ redshift)
    plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    plt.legend()
    plt.show()

    plt.figure()
    plt.semilogx(plotarr2['no-wigglesr'],plotarr2['no-wigglescorr'], label = 'No Wiggles, '+redshift)
    plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    plt.legend()
plt.show()
#plt.figure()
#plt.semilogx(plotarr2['wigglesr'],plotarr2['wigglesr']**2.0*(plotarr2['wigglesavgcorr'] - plotarr2['no-wigglesavgcorr'])/num_realizations)
#plt.show()
#plt.plot(datacc[0],datacc[1])
#   plt.show()
'''
