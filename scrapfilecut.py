
#----------------------------CREATING ORIGINAL WIGGLES AND RANDOM SAMPLING PLOTS---------------------------------------------------------------------
'''

for (file, filecut) in zip(files, filecuts):
    plotarr2 = {}
    for realization in range(0, 1): #num_realizations):

        file2='/tigress/isk/COLA_runs/plot_data/sampling/z1Quijote/'+filecut+'_Quijote_wiggles_'+str(realization)+'_voids_drdefault.dat'

        if not os.path.exists(file2):
            
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
            t = nlab.FFTCorr(mesh, mode='1d', Nmesh=Nmesh, BoxSize = BoxSize)#, dr = 1.0)
        
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
            #plotarr2= {}
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
    #ax = axs[0]
    #ax.plot(plotarr2['wigglesr'], plotarr2['wiggles'+'r']**2.0*plotarr2['wigglescorr'], label = filecut)
    #ax.set_xscale('log')
    #ax.legend(loc = "lower left")
#ax.set_ylim(ymin = -5000, ymax =5000)
#ax.hlines(y=0,xmin=0, xmax=1000, linestyles='solid')
#ax.set_title("Voids")
#plt.legend()
'''
