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

#files for mass cuts at 10^13.5:
files = ['/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat',/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out, /projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_large_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_large_mcut1e13pt5_0_halos_z1.00.out, /projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_small_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_small_mcut1e13pt5_0_halos_z1.00.out]
filecuts= ['original','random sampling', 'greater than 10^13.5', 'less than 10^13.5']

#files for mass cuts at 10^13:
#files = ['/tigress/isk/COLA_runs/plot_data/'+'z'+zlab+sim+'/Xihh_'+w+d+'_nbodykit_'+sim+'_halos_'+str(realization)+scale+binlab+'.dat',/projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_samples_0_halos_z1.00/centers_central_Quijote_wiggles_samples_0_halos_z1.00.out, /projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_large_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_large_mcut1e13_0_halos_z1.00.out, /projects/SPERGEL/COLA_runs/voidcatalogs/z1Quijote/Quijote_ss1.0/sample_Quijote_wiggles_small_mcut1e13pt5_0_halos_z1.00/centers_central_Quijote_wiggles_small_mcut1e13_0_halos_z1.00.out]

#filecuts= ['original','randomsampling', 'greaterthan10^13.5', 'lessthan10^13.5']


for (file, filecut) in zip(files, filecuts):

    file2='/tigress/isk/COLA_runs/plot_data/sampling_'+filecut+'_Quijote_1_voids.dat'

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

#CREATING PANEL OF CORRELATION FUNCTIONS. TOP IS CORR FUNC FOR WIGGLES FOR ALL MASS CUTS/SAMPLES AND BOTTOM IS ORIGINAL WIGGLES- NO-WIGGLES
fig, axs = plt.subplots(2, 1, sharex=True)

for filecut in filecuts:
    files4 ='/tigress/isk/COLA_runs/plot_data/sampling_'+filecut+'_Quijote_1_voids.dat'
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
