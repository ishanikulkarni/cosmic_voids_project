import numpy as np
import pandas as pd
import os
import matplotlib as plt
import matplotlib.pyplot as plt

import seaborn as sns

import sklearn
from sklearn.model_selection import train_test_split

import sys

d = '_'
scale = 'lin'
binlab = '_dr1'

c = '_halos'# #'_cdm'
cola = '' #'COLA/'
crazy = ''

#redshifts = [0.0,1.0]
sim = 'Quijote'
w = 'wiggles'
#realization = '0'
#zlabs = ['0','1']
#wiggles = ['wiggles', 'no-wiggles']

boxlength=2000.0
Omega_m=0.3175
h=.6711


zlab = '1'
wiggle = 'wiggles'
realization = '0'
redshift = '1.00'

file = '/projects/SPERGEL/COLA_runs/voidcatalogs/'+cola+'z'+zlab+sim+'/'+sim+'_ss1.0/sample_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'/'+crazy+'centers_central_'+sim+'_'+w+d+str(realization)+c+'_z'+redshift+'.out'

gadget_format2 = pd.read_csv(file,delim_whitespace=True,engine = 'python', header=None,index_col=False,comment='#')
                
plt.hist(gadget_format2[4], bins = 50, edgecolor = 'black')
plt.title("Histogram of Void Radius")
plt.xlabel("Radius (Mpc/h)")
plt.ylabel("Number of Voids")
plt.show()
