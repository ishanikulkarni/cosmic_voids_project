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
realization =0
scale = 'lin'
binlab = '_dr1'
num_realizations = 100

Nmesh=256

BoxSize = 2000.0

#files for mass cuts at 10^13.5:
#files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out','/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_large_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_large_mcut1e13pt5_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_small_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_small_mcut1e13pt5_0_halos_z1.00.out']
#filecuts= ['ogwiggles','randomsampling', 'greaterthan10^13.5', 'lessthan10^13.5']

#files for mass cuts at 10^13:
#files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out','/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_large_mcut1e13_0_halos_z1.00/centers_central_Quijote_wiggles_large_mcut1e13_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_small_mcut1e13_0_halos_z1.00/centers_central_Quijote_wiggles_small_mcut1e13_0_halos_z1.00.out']

#filecuts= ['ogwiggles','randomsampling', 'greaterthan10^13', 'lessthan10^13']

#files for less than 10^12.55:
files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out','/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_small_mcut1e12.55_0_halos_z1.00/centers_central_Quijote_wiggles_small_mcut1e12.55_0_halos_z1.00.out']

filecuts= ['ogwiggles','randomsampling', 'lessthan10^12.55']

#files for mass cuts greater than 10^13.25:
#files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out','/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out', '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_large_mcut1e13pt25_0_halos_z1.00/centers_central_Quijote_wiggles_large_mcut1e13pt25_0_halos_z1.00.out']

#filecuts= ['ogwiggles','randomsampling', 'greaterthan10^13.25']

#ORIGINAL FOR TESTING
#files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out']

#filecuts= ['ogwiggles']

#CREATING PANEL OF CORRELATION FUNCTIONS. TOP IS CORR FUNC FOR WIGGLES FOR ALL MASS CUTS/SAMPLES AND BOTTOM IS ORIGINAL WIGGLES- NO-WIGGLES
fig, axs = plt.subplots(2, 1, sharex=True)
for (file, filecut) in zip(files, filecuts):
    plotarr2 = {}
    #for realization in range(0, 1): #num_realizations):

    file2='/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+filecut+'_Quijote_wiggles_'+str(realization)+'_voids.dat'

    if not os.path.exists(file2):
        #file = 'projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+redshift+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+str(realization)+c+'_z'+zlab+'/'+crazy+'centers_central_'+sim+'_'+w+d+str(realization)+c+'_z'+zlab+'.out'
        gadget_format = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
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
        os.chdir('/')
	
        #print(type(realization))
        #if not (os.path.exists('/tigress/isk/Documents/research/projects/COLA_runs/plot_data/'+cola+'z'+redshift+sim+'/Xivv_'+crazy+w+d+'_nbodykit_'+sim+realization+c+'_'+scale+binlab+'.dat')):
        
        if not os.path.exists(os.path.dirname(file2)):
            os.makedirs(os.path.dirname(file2))

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
        #file4 ='/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+filecut+'_Quijote_0_voids.dat'
        plotarr2= {}
        #print(file4)
        data_array = []
        data_array =  pd.read_csv(file2,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
        gg={}
        #print(data_array)
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
   
    #fig = plt.figure()
    #panel1 = plt.semilogx(plotarr2['wigglesr'], plotarr2['wigglescorr'], label = 'Wiggle, '+ redshift)
    #fig.add_subplot(111)
    #panel2 =  plt.semilogx(plotarr2['wigglesr'], plotarr2['wigglescorr'], label = 'Wiggle, '+ redshift)
    
    #print(plotarr2.keys())
    ax = axs[0]
    ax.plot(plotarr2['wigglesr'], plotarr2['wiggles'+'r']**2.0*plotarr2['wigglescorr'], label = filecut)
    ax.set_xscale('log')
    ax.legend(loc = "lower left")
ax.set_ylim(ymin = -5000, ymax =5000)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.set_title("Voids")
#plt.legend()

#---------------------------------CREATING WIGGLES AND NO-WIGGLES FILES AND PLOTTING-------------------------------------------------------------------------------------------------------------------

plotarr3={}

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

            #print(filewrite)

        else:
            #print(filewrite)
            data_array = []
            data_array =  pd.read_csv(filewrite,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
            gg={}
            #print(data_array)
            gg['r'] = data_array[0]
            gg['corr']=data_array[1]
            
            #gg[w+'avgcorr'] = np.zeros(len(gg[w+'corr']))
            #gg[w+'avgcorr'] += gg[w+'corr']
        
        if w+'avgcorr' not in plotarr3.keys():
            plotarr3[w+'avgcorr'] = np.zeros(len(gg['r']))
            plotarr3[w+'r'] = gg['r']
            plotarr3[w+'corr'] = gg['corr']
        #try:
        #   plotarr3[w+'r'] == gg['r']
        #except ValueError:
        #   print('r values do not match previous r values')
        #print(gg.keys())

        #print(plotarr3.keys())
        plotarr3[w+'avgcorr'] += gg['corr']

ax = axs[1]
ax.plot(plotarr3['wigglesr'],plotarr3['wigglesr']**2.0*(plotarr3['wigglesavgcorr'] - plotarr3['no-wigglesavgcorr'])/num_realizations, label = 'wiggles - no-wiggles')
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.set_ylim(ymin = -150, ymax =150)
ax.legend(loc = "lower left")
plt.show()

