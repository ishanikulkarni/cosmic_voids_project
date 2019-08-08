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

#files for less than 10^12.55:
files = ['/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_0_halos_z1.00/centers_central_Quijote_wiggles_0_halos_z1.00.out']#,'/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out']

filecuts= ['ogwiggles']#,'randomsampling']

#CREATING PANEL OF CORRELATION FUNCTIONS. TOP IS CORR FUNC FOR WIGGLES FOR ALL MASS CUTS/SAMPLES AND BOTTOM IS ORIGINAL WIGGLES- NO-WIGGLES

fig, axs = plt.subplots(2, 1, sharex=True, figsize = (8,10))
fig.tight_layout
fig.subplots_adjust(hspace=0)
fig.subplots_adjust(wspace=0)
plt.rc('legend', fontsize = 'medium')

#---------------------------------CREATING WIGGLES AND NO-WIGGLES FILES AND PLOTTING-------------------------------------------------------------------------------------------------------------------

plotarr={}
xi = {}

for w in wiggles:
    for realization in range(0,num_realizations):
        filewrite = '/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/' + w + '_Quijote_' + str(realization)+'_voids_drdefault.dat'

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
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize)#, dr = 9.0)
            
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

    file = '/tigress/isk/COLA_runs/plot_data/z'+zlab+sim+'/Xivv_100AVG_ERROR_'+w+'_nbodykit_original_halos_'+scale+'_drdefault.dat'
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


ax = axs[0]
ax.plot(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'], label = 'BAO')
ax.fill_between(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] +plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] -plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], alpha=0.2)
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
#ax.set_title("Voids")
ax.set_ylabel(r"$r^2 \Delta \xi(r)$ [$h^{3}\mathrm{Mpc}^{-3}$]", fontsize = 12)
ax.set_xlim(xmin = 12, xmax =250)
ax.set_ylim(ymin = -700, ymax =700)

ax.plot(plotarr['no-wigglesr'], plotarr['no-wiggles'+'r']**2.0*plotarr['no-wigglescorr'], label = 'BAO-removed')
ax.fill_between(plotarr['no-wigglesr'], plotarr['no-wiggles'+'r']**2.0*plotarr['no-wigglescorr'] +plotarr['no-wiggles'+'r']**2.0*sigmaarr['no-wiggles'], plotarr['no-wiggles'+'r']**2.0*plotarr['no-wigglescorr'] -plotarr['no-wiggles'+'r']**2.0*sigmaarr['no-wiggles'], alpha=0.2)
ax.legend(loc = "lower left")

ax = axs[1]
ax.plot(plotarr['wigglesr'],np.zeros(len(plotarr['wigglesr'])), color='k')
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr']), label = 'BAO - no-BAO')
ax.fill_between(plotarr['wigglesr'], plotarr['wigglesr']**2.0*(((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr']))+np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), plotarr['wigglesr']**2.0*((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr'])-np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), alpha=0.2)
ax.set_xscale('log')
ax.set_xlim(xmin = 12, xmax =250)
ax.set_ylim(ymin = -100, ymax = 100)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.legend(loc = "lower left")
ax.set_xlabel(r"$r$ [$h^{-1} \ \mathrm{Mpc}$]")
ax.set_ylabel(r"$r^2 \Delta \xi(r)$ [$h^{3}\mathrm{Mpc}^{-3}$] BAO - no-BAO", fontsize = 12)
plt.show()
