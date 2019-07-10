import numpy as np
import pandas as pd

import matplotlib as plt
import matplotlib.pyplot as plt

import seaborn as sns

import sklearn
from sklearn.model_selection import train_test_split

import sys
zlab = "0"
sim = "Quijote"
wiggle = "no-wiggles"
realization = 1
filename = '/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/' +sim+'_'+wiggle+'_'+str(realization)+'_halos.dat'

with open(filename, 'r') as fn:
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

logdata = np.log10(data["mass"])

(histogram,bins,patches)= plt.hist(logdata, bins=25,range = [12.5,15.5],edgecolor='black', density=None, histtype='bar', align = 'mid')

print(len(histogram))
indices =np.where(bins<13.94)
print(len(indices))
print(indices)
print(histogram[indices])

#histogram.where(bins <13.94, 13.94,inplace=True)
#print(logdata.value_counts().sort_index())


plt.show()
print(histogram)
print(bins)


#sns.distplot(data.mass)
#plt.show()


from sklearn.model_selection import StratifiedShuffleSplit

split = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
for train_index, test_index in split.split(data,data["mass"]):
    strat_train_set = histogram.loc[train_index]
    strat_test_set = histogram.loc[test_index]
strat_test_set.head(25)
