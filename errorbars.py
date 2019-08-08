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
c = '_halos'

Nmesh=256

BoxSize = 2000.0

fig, axs = plt.subplots(2, 1, sharex=True)

#-----------------------------------------------TESTING ERROR BARS ---------------------------------------------------------------------------------

plotarr={}
xi = {}

for w in wiggles:
    for realization in range(0,num_realizations):
        filewrite = '/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+'sample_Quijote_'+w+'_'+str(realization)+'_voids_drdefault.dat'
        
        if not os.path.exists(filewrite):
            filenw = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_'+w+'_samples_'+str(realization)+'_halos_z1.00/centers_central_Quijote_'+w+'_samples_'+str(realization)+'_halos_z1.00.out'
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
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize)#, dr = 1.0)
            
            gg = t.corr
            Xigg = gg['corr']
            Xir = gg['r']
            
            xi[w+str(realization)] = Xigg
            
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
            xi[w+str(realization)] = gg['corr']
            #gg[w+'avgcorr'] = np.zeros(len(gg[w+'corr']))
            #gg[w+'avgcorr'] += gg[w+'corr']
        
        if w+'avgcorr' not in plotarr.keys():
            plotarr[w+'avgcorr'] = np.zeros(len(gg['r']))
            plotarr[w+'r'] = gg['r']
            plotarr[w+'corr'] = gg['corr']
        #try:
        #   plotarr3[w+'r'] == gg['r']
        #except ValueError:
        #   print('r values do not match previous r values')
        #print(gg.keys())

        #print(plotarr3.keys())
        plotarr[w+'avgcorr'] += gg['corr']
        
plotarr['wigglesavgcorr']= (1.0/num_realizations) *plotarr['wigglesavgcorr']
plotarr['no-wigglesavgcorr']= (1.0/num_realizations) *plotarr['no-wigglesavgcorr']

sigmaarr = {}
for w in wiggles:
    sigmaarr[w] = np.zeros(len(gg['r']))
    for realization in range(0,num_realizations):
        #realization = '_'+str(int(i))
        sigmaarr[w] += (xi[w+str(realization)] - plotarr[w+'avgcorr'])**2.0

    sigmaarr[w] /= num_realizations
    sigmaarr[w] = np.sqrt(sigmaarr[w])

    file = '/tigress/isk/COLA_runs/plot_data/z'+zlab+sim+'/Xivv_100AVG_ERROR_'+w+'_nbodykit_samples_halos_'+scale+binlab+'_drdefault.dat'
    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))
    if not (os.path.exists(file)):
        f = open(file,'w+')
        f.write('# r avgXi(r) 1Sigma\n')
        f.write('# z = '+zlab +'\n')
       
        datacc = np.array([gg['r'],plotarr[w+'avgcorr'],sigmaarr[w]])
        datacc = datacc.T
        
        np.savetxt(f, datacc)
        f.close()
    else:
        errordata = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
        plotarr[w + 'avgcorr'] = errordata[1]
        sigmaarr[w] = errordata[2]

print("Wrote avg and error file")

ax = axs[1]
ax.plot(plotarr['wigglesr'],np.zeros(len(plotarr['wigglesr'])), color='k')
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr']), label = 'Random Sampling: wiggles - no-wiggles')
ax.fill_between(plotarr['wigglesr'], plotarr['wigglesr']**2.0*(((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr']))+np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), plotarr['wigglesr']**2.0*((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr'])-np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), alpha=0.2)
ax.set_xscale('log')
ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -150, ymax =150)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
plt.show()


#THIS IS STUFF I DON'T NEED SORRY
'''
plotarr = {}
xi = {}

for w in wiggles:
    for realization in range(0,num_realizations):

        #file2='/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+'greaterthan13pt25_Quijote_'+w+'_'+str(realization)+'_voids.dat'
        file2 = '/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+'sample_Quijote_'+w+'_'+str(realization)+'_voids_dr1.dat'

        if not os.path.exists(file2):
            #file = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_'+w+'_large_mcut1e13pt25_'+str(realization)+'_halos_z1.00/centers_central_Quijote_'+w+'_large_mcut1e13pt25_'+str(realization)+'_halos_z1.00.out'
            file = '/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_'+w+'_samples_'+str(realization)+'_halos_z1.00/centers_central_Quijote_'+w+'_samples_'+str(realization)+'_halos_z1.00.out'

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
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize, dr = 1.0)
        
            gg = t.corr
            Xigg = gg['corr']
            Xir = gg['r']

            xi[w+str(realization)] = Xigg

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
            data_array = []
            data_array =  pd.read_csv(file2,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
            gg={}
            #print(data_array)
            gg['r'] = data_array[0]
            gg['corr']=data_array[1]
            xi[w+str(realization)] = gg['corr']

            #gg[w+'avgcorr'] = np.zeros(len(gg[w+'corr']))
            #gg[w+'avgcorr'] += gg[w+'corr']
                
        if w+'avgcorr' not in plotarr.keys():
            plotarr[w+'avgcorr'] = np.zeros(len(gg['r']))
            plotarr[w+'r'] = gg['r']
            plotarr[w+'corr'] = gg['corr']
       # else:
        #    try:
         #       plotarr[w+'r'] == gg['r']
          #  except ValueError:
           #     print('r values do not match previous r values')
                #print(gg.keys())
                
        plotarr[w+'avgcorr'] += gg['corr']
        
       
plotarr['wigglesavgcorr']= (1.0/num_realizations) *plotarr['wigglesavgcorr']
plotarr['no-wigglesavgcorr']= (1.0/num_realizations) *plotarr['no-wigglesavgcorr']

ax = axs[0]
ax.plot(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'], label = 'Random Sampling')
ax.set_xscale('log')
ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -5000, ymax =5000)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.set_title("Voids")

#plt.show()
#plt.legend()

sigmaarr = {}
for w in wiggles:
    sigmaarr[w] = np.zeros(len(gg['r']))
    for realization in range(0,num_realizations):
        #realization = '_'+str(int(i))
        sigmaarr[w] += (xi[w+str(realization)] - plotarr[w+'avgcorr'])**2.0

    sigmaarr[w] /= num_realizations
    sigmaarr[w] = np.sqrt(sigmaarr[w])

    file = '/tigress/isk/COLA_runs/plot_data/z'+zlab+sim+'/Xivv_100AVG_ERROR_'+w+'_nbodykit_samples_halos_'+scale+binlab+'_dr1.dat'
    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))
    if True: # not (os.path.exists(file)):
        f = open(file,'w+')
        f.write('# r avgXi(r) 1Sigma\n')
        f.write('# z = '+zlab +'\n')
       
        datacc = np.array([gg['r'],plotarr[w+'avgcorr'],sigmaarr[w]])
        datacc = datacc.T
        
        np.savetxt(f, datacc)
        f.close()
    else:
        errordata = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
        plotarr[w + 'avgcorr'] = errordata[1]
        sigmaarr[w] = errordata[2]


print("Wrote avg and error file")

ax = axs[1]
ax.plot(plotarr['wigglesr'],np.zeros(len(plotarr['wigglesr'])), color='k')
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr']), label = 'Random Sampling: wiggles - no-wiggles')
ax.fill_between(plotarr['wigglesr'], plotarr['wigglesr']**2.0*(((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr']))+np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), plotarr['wigglesr']**2.0*((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr'])-np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), alpha=0.2)
ax.set_xscale('log')
ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -150, ymax =150)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
plt.show()
'''
