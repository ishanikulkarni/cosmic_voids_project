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

fig, axs = plt.subplots(2, 1, sharex=True)
plotarr= {}

#CREATING PANEL OF CORRELATION FUNCTIONS. TOP IS CORR FUNC FOR WIGGLES FOR ALL MASS CUTS/SAMPLES AND BOTTOM IS ORIGINAL WIGGLES- NO-WIGGLES
for w in wiggles:
    for realization in range(0, num_realizations):
        file2='/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/Quijote_original_'+w+'_' + str(realization)+'_voids_drdefault.dat'

        if not os.path.exists(file2):
            file = '/projects/SPERGEL/COLA_runs/voidcatalogs/z'+zlab+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+'_'+str(realization)+'_halos_z'+redshift+'/'+ 'centers_central_'+sim+'_'+w+'_'+str(realization)+'_halos_z'+redshift+'.out'

            gadget_format = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
            
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
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize)#, dr = 5.0)
        
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
            #plotarr= {}
           
            data_array = []
            data_array =  pd.read_csv(file2,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
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
        try:
            plotarr[w+'r'] == gg['r']
        except ValueError:
            print('r values do not match previous r values')
        #print(gg.keys())
                            
        plotarr[w+'avgcorr'] += gg['corr']
   
ax = axs[0]
ax.plot(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglesavgcorr']/num_realizations, label = 'original wiggles')
ax.set_xscale('log')
ax.legend(loc = "lower left")
ax.set_ylim(ymin = -5000, ymax =5000)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')

ax = axs[0]
ax.plot(plotarr['no-wigglesr'], plotarr['no-wiggles'+'r']**2.0*plotarr['no-wigglesavgcorr']/num_realizations, label = 'original no-wiggles')
ax.set_xscale('log')
ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -5000, ymax =5000)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')

ax = axs[1]
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr'])/num_realizations, label = 'original: wiggles - no-wiggles')
ax.set_xscale('log')
ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -500, ymax =500)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')


#----------------------------------------------CREATING VOID CUT PLOT-----------------------------------------------------------------------------------------------------------
for (zlab,redshift) in zip(zlabs,redshifts):
    
    plotarr2 = {}

    for w in wiggles:
        for realization in range(0,num_realizations):
            
            #file2='/tigress/isk/COLA_runs/voidcuts/z' + zlab + sim+ '/' +sim +'_' + w + '_large_rcut40_' + str(realization)+ '_voids_drdefault.dat'
            
            file2='/projects/SPERGEL/COLA_runs/voidcuts/z' + zlab + sim+ '/' +sim +'_'+ w + '_large_rcut40_' + str(realization)+ '_voids_drdefault.dat'
           
            if not os.path.exists(file2):
                file = '/projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+zlab+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'/'+crazy+'centers_central_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'.out'

               # file = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_'+w+'_large_mcut1e13pt25_'+str(realization)+'_halos_z1.00/centers_central_Quijote_'+w+'_large_mcut1e13pt25_'+str(realization)+'_halos_z1.00.out'

                gadget_format2 = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
                
                #plt.hist(gadget_format2[4], bins = 50)
                #plt.show()
                gadget_format = gadget_format2[gadget_format2[4]>40]
 
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
                t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize)#, dr = 5.0)
                
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
                #plotarr2= {}
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
            
            
            plotarr2[w+'avgcorr'] += gg['corr']
    
    ax = axs[0]
    ax.plot(plotarr2['wigglesr'], plotarr2['wiggles'+'r']**2.0*plotarr2['wigglesavgcorr']/num_realizations, label = 'Radius cut Larger than 40 wiggles')
    ax.set_xscale('log')
    ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    ax.set_title("Voidcuts")
    ax.legend(loc = "lower left")

    ax = axs[0]
    ax.plot(plotarr2['no-wigglesr'], plotarr2['no-wiggles'+'r']**2.0*plotarr2['no-wigglesavgcorr']/num_realizations, label = 'Radius cut Larger than 40 no-wiggles')
    ax.set_xscale('log')
    ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    ax.legend(loc = "lower left")

    ax = axs[1]
    ax.plot(plotarr2['wigglesr'],plotarr2['wigglesr']**2.0*(plotarr2['wigglesavgcorr'] - plotarr2['no-wigglesavgcorr'])/num_realizations, label = 'Radius cut: wiggles - no-wiggles')
    ax.set_xscale('log')
    ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    ax.set_ylim(ymin = -500, ymax =500)
    ax.legend(loc = "lower left")

    
plt.show()


