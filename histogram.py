import numpy as np
import pandas as pd
#import os
import matplotlib as plt
import matplotlib.pyplot as plt


sim = "Quijote"
wiggle = "no-wiggles"
#realization = 1

redshift = 1.0
sim = 'Quijote'
zlab = '1'
wiggle = 'wiggles'
realization = 0
#num_realizations = 2

boxlength=2000.0
Omega_m=0.3175
h=.6711

filename = '/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/' +sim+'_'+wiggle+'_'+str(realization)+'_halos.dat'

data = pd.read_csv(filename, skiprows = 5,skipfooter=1,engine = 'python', delim_whitespace=True, header= None)
data.columns= ["ID", "x_pos", "y_pos", "z_pos", "x_vel", "y_vel", "z_vel", "mass"]

logdata = np.log10(data["mass"])
plt.hist(logdata,bins = 50, range = [12.5,15.5], edgecolor = 'black', density = None, histtype = 'bar', align= 'mid')
plt.title('Histogram of Halo Mass')
plt.xlabel("Mass (Solar Masses)")
plt.ylabel("Number of Halos")
plt.show()
  
