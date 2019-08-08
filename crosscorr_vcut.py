import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
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
plt.rc('legend', fontsize = 'large')

#---------------------------------CREATING WIGGLES AND NO-WIGGLES FILES AND PLOTTING-------------------------------------------------------------------------------------------------------------------

plotarr={}
xi = {}

for w in wiggles:
    for realization in range(0,num_realizations):

        filewrite = '/projects/SPERGEL/COLA_runs/plot_data/z1Quijote/Xivv_' + w + '_' + str(realization)+'_Quijote_halos_log_largebins.dat'

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

    file = '/tigress/isk/COLA_runs/plot_data/z'+zlab+sim+'/Xivv_100AVG_ERROR_'+w+'_vide_original_halos_'+scale+'_drdefault.dat'
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
ax.plot(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'])
ax.fill_between(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] +plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] -plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], alpha=0.2)
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
#ax.set_title("Voids")
ax.set_ylim(ymin = -2500, ymax =2500)

ax = axs[1]
ax.plot(plotarr['wigglesr'],np.zeros(len(plotarr['wigglesr'])), color='k')
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr']), label = 'Void autocorrelation function')
ax.fill_between(plotarr['wigglesr'], plotarr['wigglesr']**2.0*(((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr']))+np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), plotarr['wigglesr']**2.0*((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr'])-np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), alpha=0.2)
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')



#-----------------------------------------------CREATING 100 REALIZATIONS OF SPECIFIC FILE----------------------------------------------------------------------------------

plotarr={}
xi = {}

for w in wiggles:
    for realization in range(0,num_realizations):
      
        #file = '/projects/SPERGEL/COLA_runs/plot_data/z1Quijote/Xivh_' + w + '_large_voidcut40_' + str(realization)+'_Quijote_halos_log_largebins.dat'
        file = '/projects/SPERGEL/COLA_runs/plot_data/z1Quijote/Xivh_' + w+'_large_voidrcut40_' + str(realization)+'_Quijote_halos_log_largebins.dat'

        data_array = []
        data_array =  pd.read_csv(file,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
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
    print(len(plotarr[w+'avgcorr']))   
#print('r values = ')
#print(gg['r'][1]-gg['r'][0])

plotarr['wigglesavgcorr']= (1.0/num_realizations) *plotarr['wigglesavgcorr']
plotarr['no-wigglesavgcorr']= (1.0/num_realizations) *plotarr['no-wigglesavgcorr']
print(len(plotarr['wigglesavgcorr']))
print(len(plotarr['no-wigglesavgcorr']))
print(len(plotarr['wigglesr']))
print(len(plotarr['no-wigglesr']))


sigmaarr = {}
for w in wiggles:
    sigmaarr[w] = np.zeros(len(gg['r']))
    for realization in range(0,num_realizations):
        #realization = '_'+str(int(i))
        sigmaarr[w] += (xi[w+str(realization)] - plotarr[w+'avgcorr'])**2.0

    sigmaarr[w] /= num_realizations
    sigmaarr[w] = np.sqrt(sigmaarr[w])

    file = '/tigress/isk/COLA_runs/plot_data/z'+zlab+sim+'/Xivh_100AVG_ERROR_'+w+'_large_voidcut40_'+scale+binlab+'_drdefault.dat'
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
ax.plot(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'])
ax.fill_between(plotarr['wigglesr'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] +plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], plotarr['wiggles'+'r']**2.0*plotarr['wigglescorr'] -plotarr['wiggles'+'r']**2.0*sigmaarr['wiggles'], alpha=0.2)
ax.set_xscale('log')
ax.set_ylim(ymin = -640, ymax =640)
ax.tick_params(axis = 'y', labelsize = 14)
ax.set_ylabel(r"$r^2 \xi(r)$ [$h\mathrm{Mpc}^{-1}$]", fontsize = 15)

ax = axs[1]
ax.plot(plotarr['wigglesr'],np.zeros(len(plotarr['wigglesr'])), color='k')
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr']), label = 'Void halo cross-correlation for Radius >40')
ax.fill_between(plotarr['wigglesr'], plotarr['wigglesr']**2.0*(((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr']))+np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), plotarr['wigglesr']**2.0*((plotarr['wigglesavgcorr']-plotarr['no-wigglesavgcorr'])-np.sqrt(sigmaarr['wiggles']**2.0 + sigmaarr['no-wiggles']**2.0)), alpha=0.2)
ax.set_xscale('log')
ax.legend(loc = "lower left")
ax.set_ylim(ymin = -100, ymax = 100)
ax.set_xlim(xmin = 16, xmax = 250)
ax.tick_params(axis = 'x', labelsize = 14)
ax.tick_params(axis = 'y', labelsize = 14)
ax.set_xlabel(r"$r$ [$h^{-1} \ \mathrm{Mpc}$]", fontsize = 15)
ax.set_ylabel(r"$r^2 \Delta \xi(r)$ [$h\mathrm{Mpc}^{-1}$]", fontsize = 15)
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
plt.show()