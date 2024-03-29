import sys
sys.path.append('/tigress/isk/Cosmic_voids_project/cosmic_voids_project/Pylians/library/')
import readfof
import sys
import os.path
import os
import math
import numpy as np

#input files
wiggle_set = ["wiggles", "no-wiggles"]
realization= "0"

snapnum = 1			#redshift 0
redshift = 0.0
for wiggle in wiggle_set:
    for realization in range(0,100):
    	#read the halo catalogue
    	snapdir = '/projects/SPERGEL/COLA_runs/'+wiggle+'/'+str(realization)
    	FoF = readfof.FoF_catalog(snapdir, snapnum, long_ids=False, swap=False, SFR=False, read_IDs=False)

    	#get the properties of the halos
    	pos_h = FoF.GroupPos/1e3		#Halo positions in Mpc/h
    	mass = FoF.GroupMass*1e10 		#Halo masses in Msun/h
    	vel_h = FoF.GroupVel*(1.0+redshift)	#Halo peculiar velocities in km/s

    	dict = {}
    	dict['ID'] = np.arange(len(pos_h))
    	dict['x_pos'] = pos_h[:,0]
    	dict['y_pos'] = pos_h[:,1]
    	dict['z_pos'] = pos_h[:,2]
    	dict['x_vel'] = vel_h[:,0]
    	dict['y_vel'] = vel_h[:,1]
    	dict['z_vel'] = vel_h[:,2]

    	boxlength=2000.0
    	Omega_m=0.3175
    	h=.6711

    	zlab = "0"
    	sim = "Quijote"
    	numhalos = str(len(pos_h))

    	fpath = '/projects/SPERGEL/COLA_runs/data/z'+zlab+sim+'/'+sim+'_'+wiggle+'_'+str(realization)+'_halos.dat'
    	#Check if file exists, and if not, creates file for output 
    	if os.path.exists(fpath)==False:
       	   file = open(fpath, 'w+')
       	   file.write(str(boxlength)+"\n"+ str(Omega_m)+"\n"+ str(h)+"\n"+ str(redshift) + "\n" + str(numhalos) + "\n")
       	   data = np.array([dict['ID'], dict['x_pos'] , dict['y_pos'] , dict['z_pos'], dict['x_vel'], dict['y_vel'], dict['z_vel'], mass])
       	   data = data.T
       	   np.savetxt(file, data, fmt=('%d' ,'%.8f',' %.8f', '%.8f', '%.8f', '%.8f', '%.8f','%d'))
       	   file.write("-99 -99 -99 -99 -99 -99 -99 -99")
       	   file.close()
