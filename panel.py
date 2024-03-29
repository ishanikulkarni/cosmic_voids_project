import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

#Original data:
#file4 = '/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat'
#file3 ='/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/Quijote_'+w+d+str(realization)+'_halos.dat'

# Random sampling data:
#file4 = '/tigress/isk/COLA_runs/trial_folder'+'z'+zlab+sim+'samplestrial'+w+d+'_nbodykit_' + sim + '_halos_' + str(realization)+'.dat'
#file3 = '/projects/SPERGEL/COLA_runs/random_sampling/z'+zlab+sim +'/'+sim+'_'+w+'_samples_'+str(realization)+'_halos.dat'
            
#Mass cuts larger than 10^13.5:
#file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat'
#file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat'
                
#Mass cuts smaller than 10^13.5:
#file4 = '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat'
#file3 = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat'
#filereads = ['/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/Quijote_'+w+d+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/random_sampling/z'+zlab+sim +'/'+sim+'_'+w+'_samples_'+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat', '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat']

zlab = '1'
redshift = '1.00'
sim = 'Quijote'
w = 'wiggles'
wiggles = ['wiggles','no-wiggles']
d = '_'
realization =0
scale = 'lin'
binlab = '_dr1'
num_realizations = 5
files = ['/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat','/tigress/isk/COLA_runs/random_sampling/z'+zlab+sim+'/' + sim + '_samplestrial_'+w+d+'_nbodykit_halos_' + str(realization)+'.dat', '/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_large_mcut1e13pt5_'+str(realization)+'_halos.dat','/tigress/isk/COLA_runs/trial_folder/'+'z'+zlab +sim +'/' +sim+'_mcuttrial_'+w+'_small_mcut1e13_'+str(realization)+'_halos.dat']
filecuts= ['original', 'random sampling', 'greater than 10^13.5', 'less than 10^13.5']

#files = ['/projects/SPERGEL/COLA_runs/mass_data/z'+zlab +sim +'/' +sim+'_'+w+'_small_mcut1e13pt5_'+str(realization)+'_halos.dat']

fig, axs = plt.subplots(2, 1, sharex=True)

for (file4, filecut) in zip(files, filecuts):
    plotarr2= {}
    print(file4)
    data_array = []
    data_array =  pd.read_csv(file4,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
    gg={}
    print(data_array)
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
    
    ax = axs[0]
    ax.plot(plotarr2['wigglesr'], plotarr2['wiggles'+'r']**2.0*plotarr2['wigglescorr'], label = filecut)
    ax.set_xscale('log')
    ax.set_ylim(ymin = -200, ymax =200)
    ax.legend(loc = "lower left")
    #ax1.plot(plotarr2['wigglesr'], plotarr2['wigglescorr'])
    #ax2.plot(plotarr2['wigglesr'], plotarr2['wigglescorr'])
    #ax.Axes.set_xscale(fig,'log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.set_title("Halos")
#plt.legend()

plotarr={}

for w in wiggles:
    for realization in range(0,num_realizations):

        file5 = '/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat'
        data_array = []
        data_array =  pd.read_csv(file5,delim_whitespace=True, header=None,engine='python',index_col=False,comment='#')
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
                            
        plotarr[w+'avgcorr'] += gg['corr']

ax = axs[1]
ax.plot(plotarr['wigglesr'],plotarr['wigglesr']**2.0*(plotarr['wigglesavgcorr'] - plotarr['no-wigglesavgcorr'])/num_realizations, label = 'wiggles - no-wiggles')
ax.set_xscale('log')
ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
ax.legend(loc = "lower left")
plt.show()

#plt.figure()
    #plt.semilogx(plotarr2['wigglesr'],plotarr2['wigglesr']**2.0*(plotarr2['wigglesavgcorr'] - plotarr2['no-wigglesavgcorr'])/num_realizations, label = redshift)
    #plt.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
    #plt.legend()
