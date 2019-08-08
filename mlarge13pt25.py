import pandas as pd
import numpy as np
import os

#zlab='0'
redshifts = [0.0,1.0]
sim = 'Quijote'
#wiggle = 'wiggles'
#realization = 0
zlabs = ['0','1']
wiggles = ['wiggles', 'no-wiggles']

num_realizations = 100

boxlength=2000.0
Omega_m=0.3175
h=.6711

for (zlab,redshift) in zip(zlabs,redshifts):
    for wiggle in wiggles:
        for realization in range(0,num_realizations):
            newfile = '/projects/SPERGEL/COLA_runs/mass_data/z'+zlab + sim + '/'+sim+'_'+wiggle+'_large_mcut1e13pt25_'+str(realization)+'_halos.dat'

            if not os.path.exists(newfile):
                halo_file = '/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/'+sim+'_'+wiggle+'_'+str(realization)+'_halos.dat'
                #boxlength = pd.read_csv(halo_file,nrows=1)
                data_array = pd.read_csv(halo_file, delim_whitespace = True,skipfooter=1, engine='python',skiprows = 5,header=None,index_col = False)

                #print(data_array)
                #mass_array = data_array[7]
                #select_array = []
                #print(mass_array)
                #length = len(mass_array)
                #print(len(mass_array))
                #print(select_array)

         
                #print(newfile)
                select_array=data_array[data_array[7]>10**13.25]
                numhalos = str(len(select_array[0]))
                #print(select_array)
                os.chdir('/')
                
                if not os.path.exists(os.path.dirname(newfile)):
                    os.makedirs(os.path.dirname(newfile))
                f = open(newfile,'w+')
                data = np.array([select_array[0], select_array[1], select_array[2], select_array[3], select_array[4], select_array[5], select_array[6], select_array[7]])
                #next(f)
                data = data.T
                f.write(str(boxlength)+"\n"+ str(Omega_m)+"\n"+ str(h)+"\n"+ str(redshift) + "\n" + str(numhalos) + "\n")
                np.savetxt(f, data, fmt=('%d' ,'%.8f',' %.8f', '%.8f', '%.8f', '%.8f', '%.8f','%d'))
                #f.write(str(select_array))
                f.write("-99 -99 -99 -99 -99 -99 -99 -99")
                f.close()

'''
#if not os.path.exists(newfile):
    #print("in the first if statement")
    if not os.path.exists(os.path.dirname(newfile)):
        os.makedirs(os.path.dirname(newfile))
    for i in range(1,length):
        #print("in the for loop")
        if mass_array[i]>10**13.5:
            os.chdir('/')
            f = open(newfile,'a+')
            f.write(str(data_array.iloc[i:i+1]))
            f.close()
        #select_array.append(data_array[i,:])
'''
#print mass_array
