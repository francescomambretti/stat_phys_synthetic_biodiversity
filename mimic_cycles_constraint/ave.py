# print the average over nruns of the histograms obtained at a given cycle
# written by Francesco Mambretti

import numpy as np
import sys
import string

ave_array=np.zeros(21)
ave_array2=np.zeros(21)

nruns=50

for c in np.arange (3,12,3):
	for i in np.arange (0,nruns):
		array=np.loadtxt("run_{}/cycle.{}.txt".format(i,c))
		ave_array+=array
		ave_array2+=array*array

	ave_array/=nruns
	ave_array2/=nruns

	for a in range (0,21):
		output="cycle.{}.txt".format(c)
		with open(output, mode='a') as file_object:
			print (a,ave_array[a],np.sqrt(np.abs(ave_array2[a]-ave_array[a]*ave_array[a])/(nruns)),file=file_object)
