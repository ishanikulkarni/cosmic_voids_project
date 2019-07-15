import numpy as np
import pandas as pd
import os
import matplotlib as plt
import matplotlib.pyplot as plt

import seaborn as sns

import sklearn
from sklearn.model_selection import train_test_split

import sys
#zlab = "0"
sim = "Quijote"
wiggle = "no-wiggles"
realization = 1

redshift = 0.0
sim = 'Quijote'
zlabs = ['0']
#zlabs = ['0','1']
wiggles = ['wiggles'] # 'no-wiggles']

num_realizations = 2

boxlength=2000.0
Omega_m=0.3175
h=.6711

for zlab in zlabs:
    for wiggle in wiggles:
        for realization in range(1,num_realizations):

            newfile = '/tigress/isk/COLA_runs/random_sampling/z'+zlab+sim+'_'+wiggle+'_'+'_samples_'+str(realization)+'_halos.dat'

            if not os.path.exists(newfile):
                filename = '/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/' +sim+'_'+wiggle+'_'+str(realization)+'_halos.dat'

                data = pd.read_csv(filename, skiprows = 5,skipfooter=1, delim_whitespace=True, header= None)
                data.columns= ["ID", "x_pos", "y_pos", "z_pos", "x_vel", "y_vel", "z_vel", "mass"]
                #print(data.head(10))
                #data.info()
                #print(min(data["mass"]))
                #print(max(data["mass"]))
    
                #correlation_matrix = data.corr()
                #plt.subplots(figsize=(8,6))
                #sns.heatmap(correlation_matrix, center =0, annot=True,linewidths=.3)
                #plt.show()

                #corr = data.corr()
                #corr["mass"].sort_values(ascending=True)
                #print(corr)
    
                logdata = np.round(np.log10(data["mass"]),decimals=4)
                #logdata = np.log10(data["mass"])
    
                '''
                print(len(histogram))
                indices =np.where(bins<13.94)
                print(len(indices))
                print(indices)
                print(histogram[indices])
                '''
                logdata.where(logdata <14.06, 14.06,inplace=True)
                print(logdata.value_counts().sort_index())
                
                (histogram,bins,patches)= plt.hist(logdata, bins=25,range = [12.5,15.5],edgecolor='black', density=None, histtype='bar', align = 'mid')
    
    
                plt.show()
                print(histogram)
                print(bins)
    

                #sns.distplot(data.mass)
                #plt.show()


                from sklearn.model_selection import StratifiedShuffleSplit
                
                split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
                for train_index, test_index in split.split(data,logdata):
                    strat_train_set = data.loc[train_index]
                    strat_test_set = data.loc[test_index]
                print(strat_test_set.head(25))
                numhalos= len(strat_test_set)

                os.chdir('/')
                
                if not os.path.exists(os.path.dirname(newfile)):
                    os.makedirs(os.path.dirname(newfile))

                f = open(newfile,'w+')

                f.write(str(boxlength)+"\n"+ str(Omega_m)+"\n"+ str(h)+"\n"+ str(redshift) + "\n" + str(numhalos) + "\n")
                np.savetxt(f, strat_test_set, fmt=('%d' ,'%.8f',' %.8f', '%.8f', '%.8f', '%.8f', '%.8f','%d'))
                f.write("-99 -99 -99 -99 -99 -99 -99 -99")
                f.close()
